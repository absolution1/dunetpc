// -*- mode: c++; c-basic-offset: 4; -*-
////////////////////////////////////////////////////////////////////////
//
//  \file LBNEOpDetResponse_service.cc
//
////////////////////////////////////////////////////////////////////////


#include "lbne/OpticalDetector/LBNEOpDetResponse.h"
#include "TGeoNode.h"
#include "TGeoBBox.h"
#include "Geometry/OpDetGeo.h"
#include "Utilities/LArProperties.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "CLHEP/Random/RandFlat.h"


namespace opdet{


    //--------------------------------------------------------------------
    LBNEOpDetResponse::LBNEOpDetResponse(fhicl::ParameterSet const& pset, 
                                         art::ActivityRegistry &/*reg*/)
    {
        this->doReconfigure(pset);
    }
    
    //--------------------------------------------------------------------
    LBNEOpDetResponse::~LBNEOpDetResponse() throw()
    { }


    //--------------------------------------------------------------------
    void LBNEOpDetResponse::doReconfigure(fhicl::ParameterSet const& pset)
    {
        double tempfQE =         pset.get<double>("QuantumEfficiency");
        fWavelengthCutLow =      pset.get<double>("WavelengthCutLow");
        fWavelengthCutHigh =     pset.get<double>("WavelengthCutHigh");
        fLightGuideAttenuation = pset.get<bool>("LightGuideAttenuation");
        fChannelConversion =     pset.get<std::string>("ChannelConversion");
        std::string tmpAxis =    pset.get<std::string>("LongAxis"); 

        boost::algorithm::to_lower(tmpAxis);

        if (tmpAxis == "x") fLongAxis = 0;
        if (tmpAxis == "y") fLongAxis = 1;
        if (tmpAxis == "z") fLongAxis = 2;
        
        // Only allow channel conversion once - so it must be set to happen
        // either during full simulation (library generation) or during
        // fast simulation (library use).
        
        boost::algorithm::to_lower(fChannelConversion);
        
        fFullSimChannelConvert = false;
        fFastSimChannelConvert = false;
        
        if (fChannelConversion == "full") fFullSimChannelConvert = true;
        if (fChannelConversion == "fast") fFastSimChannelConvert = true;

        // Correct out the prescaling applied during simulation
        art::ServiceHandle<util::LArProperties>   LarProp;
        fQE = tempfQE / LarProp->ScintPreScale();
        
        if (fQE > 1.0001 ) {
            mf::LogError("LBNEOpDetResponse_service") << "Quantum efficiency set in OpDetResponse_service, " << tempfQE
                                                      << " is too large.  It is larger than the prescaling applied during simulation, "
                                                      << LarProp->ScintPreScale()
                                                      << ".  Final QE must be equalt to or smaller than the QE applied at simulation time.";
            assert(false);
        }

    }


    //--------------------------------------------------------------------
    int  LBNEOpDetResponse::doNOpChannels() const
    {
        art::ServiceHandle<geo::Geometry> geom;
        if (fFastSimChannelConvert || fFullSimChannelConvert)
            return geom->NOpChannels();
        else
            return geom->NOpDets();

    }


    //--------------------------------------------------------------------
    bool LBNEOpDetResponse::doDetected(int OpChannel, const sim::OnePhoton& Phot, int &newOpChannel) const
    {
        
        // Find the Optical Detector using the geometry service
        art::ServiceHandle<geo::Geometry> geom;
        const TGeoNode* node = geom->OpDetGeoFromOpChannel(OpChannel).Node();

        // Identify the photon detector type
        int pdtype;
        std::string detname = node->GetName();
        boost::to_lower(detname);
        if      (detname.find("bar") != std::string::npos )   pdtype = 0;
        else if (detname.find("fiber") != std::string::npos ) pdtype = 1;
        else if (detname.find("plank") != std::string::npos ) pdtype = 2;
        else                                                  pdtype = -1;


        if (fFullSimChannelConvert){
            // Override default number of channels for Fiber and Plank
            float NOpHardwareChannels = geom->NOpHardwareChannels(OpChannel);
            if (pdtype == 1) NOpHardwareChannels = 3;
            if (pdtype == 2) NOpHardwareChannels = 2;
            
            int hardwareChannel = (int) ( CLHEP::RandFlat::shoot(1.0) * NOpHardwareChannels );
            newOpChannel = geom->OpChannel(OpChannel, hardwareChannel);
        }
        else{
            newOpChannel = OpChannel;
        }
        
        // Check QE
        if ( CLHEP::RandFlat::shoot(1.0) > fQE ) return false;

        double wavel = wavelength(Phot.Energy);
        // Check wavelength acceptance
        if (wavel < fWavelengthCutLow) return false;
        if (wavel > fWavelengthCutHigh) return false;

        if (fLightGuideAttenuation) {
            // Get the length of the photon detector
            TGeoBBox *box = (TGeoBBox*)node->GetVolume()->GetShape();
            double opdetLength = 0;
            double sipmDistance = 0;

            if (fLongAxis == 0) {
                opdetLength = box->GetDX();
                sipmDistance = opdetLength - Phot.FinalLocalPosition.x();
            }
            else if (fLongAxis == 1) {
                opdetLength = box->GetDY();
                sipmDistance = opdetLength - Phot.FinalLocalPosition.y();
            }
            else if (fLongAxis == 2) {
                opdetLength = box->GetDZ();
                sipmDistance = opdetLength - Phot.FinalLocalPosition.z();
            }
            else {
                mf::LogError("LBNEOpDetResponse") << "Unknown axis, fLongAxis = " << fLongAxis;
                assert(false);
            }



            if (pdtype == 0) {
                // Now assume Bar and Fiber have the same behavior, modify later when more is known
                double normalize   = 0.6719; // Normalize mean performance to be the same in all PDs
                double lambdaShort =  5.56; // cm
                double fracShort   =  0.40 * normalize;
                double lambdaLong  = 44.13; // cm
                double fracLong    =  0.60 * normalize;

                // Throw away some photons based on attenuation
                double AttenuationProb = fracShort*exp(-sipmDistance/lambdaShort) + fracLong*exp(-sipmDistance/lambdaLong);
                
                //mf::LogVerbatim("LBNEOpDetResponse") << "OpChannel: " << OpChannel << " is a " << pdtype 
                //                                     << " with length " << opdetLength << " in detector "
                //                                     << box->GetDX() << " x " << box->GetDY()  << " x " << box->GetDZ()
                //                                     << " named " << detname;
                //mf::LogVerbatim("LBNEOpDetResponse") << "   Local Position = (" << Phot.FinalLocalPosition.x() 
                //                                     << ", " << Phot.FinalLocalPosition.y() << ", " << Phot.FinalLocalPosition.z() << ")";
                //mf::LogVerbatim("LBNEOpDetResponse") << "   Distance to SiPM = " << sipmDistance << " along axis " << fLongAxis;
                //mf::LogVerbatim("LBNEOpDetResponse") << "   Attenuation Probability = " << AttenuationProb;

                if ( CLHEP::RandFlat::shoot(1.0) > AttenuationProb ) return false;

          
                
               
            }
            else if (pdtype == 1) {
                double normalize   = 1.0; // Normalize mean performance to be the same in all PDs
                double lambda      = 14.6; // cm

                // Throw away some photons based on attenuation
                double AttenuationProb = normalize*exp(-sipmDistance/lambda);
                if ( CLHEP::RandFlat::shoot(1.0) > AttenuationProb ) return false;
            }
            else if (pdtype == 2) {
                double normalize   = 0.4305; // Normalize mean performance to be the same in all PDs
                double lambda      = 48.4; // cm
                double altDistance = opdetLength - sipmDistance;
                double frac        = 0.5 * normalize;

                // Throw away some photons based on attenuation
                double AttenuationProb = frac*exp(-sipmDistance/lambda) + frac*exp(-altDistance/lambda);
                if ( CLHEP::RandFlat::shoot(1.0) > AttenuationProb ) return false;
            }
            else {
                mf::LogWarning("LBNEOpDetResponse") << "OpDet: " << OpChannel << " is an unknown PD type named: " << detname 
                                                    << ". Assuming no attenuation.";
            }

        }

        return true;
    }

    //--------------------------------------------------------------------
    bool LBNEOpDetResponse::doDetectedLite(int OpChannel, int &newOpChannel) const
    {
        if (fFastSimChannelConvert){

            // Find the Optical Detector using the geometry service
            art::ServiceHandle<geo::Geometry> geom;
            // Here OpChannel must be opdet since we are introducing
            // channel mapping here.
            const TGeoNode* node = geom->OpDetGeoFromOpDet(OpChannel).Node();

            
            // Identify the photon detector type
            int pdtype;
            std::string detname = node->GetName();
            boost::to_lower(detname);
            if      (detname.find("bar") != std::string::npos )   pdtype = 0;
            else if (detname.find("fiber") != std::string::npos ) pdtype = 1;
            else if (detname.find("plank") != std::string::npos ) pdtype = 2;
            else                                                  pdtype = -1;

            // Override default number of channels for Fiber and Plank
            float NOpHardwareChannels = geom->NOpHardwareChannels(OpChannel);
            if (pdtype == 1) NOpHardwareChannels = 3;
            if (pdtype == 2) NOpHardwareChannels = 2;

            int hardwareChannel = (int) ( CLHEP::RandFlat::shoot(1.0) * NOpHardwareChannels );
            newOpChannel = geom->OpChannel(OpChannel, hardwareChannel);
        }
        else{
            newOpChannel = OpChannel;
        }
        
        // Check QE
        if ( CLHEP::RandFlat::shoot(1.0) > fQE ) return false;

        
        return true;
    }

} // namespace

DEFINE_ART_SERVICE_INTERFACE_IMPL(opdet::LBNEOpDetResponse, opdet::OpDetResponseInterface)

