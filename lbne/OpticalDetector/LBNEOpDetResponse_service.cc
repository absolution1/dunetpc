////////////////////////////////////////////////////////////////////////
//
//  \file LBNEOpDetResponse_service.cc
//
////////////////////////////////////////////////////////////////////////


#include "lbne/OpticalDetector/LBNEOpDetResponse.h"
#include "TGeoNode.h"
#include "TGeoBBox.h"
#include "Geometry/OpDetGeo.h"
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
        fQE=                       pset.get<double>("QuantumEfficiency");
        fWavelengthCutLow=         pset.get<double>("WavelengthCutLow");
        fWavelengthCutHigh=        pset.get<double>("WavelengthCutHigh");
        fLightGuideAttenuation=    pset.get<bool>("LightGuideAttenuation");
    }


    //--------------------------------------------------------------------
    bool LBNEOpDetResponse::doDetected(int OpChannel, const sim::OnePhoton& Phot) const
    {
        // Check QE
        if ( CLHEP::RandFlat::shoot(1.0) > fQE ) return false;

        double wavel = wavelength(Phot.Energy);
        // Check wavelength acceptance
        if (wavel < fWavelengthCutLow) return false;
        if (wavel > fWavelengthCutHigh) return false;

        if (fLightGuideAttenuation) {
            // Find the Optical Detector using geom
            unsigned int cryostatID=0, opdetID=0;
            art::ServiceHandle<geo::Geometry> geom;
            geom->OpChannelToCryoOpDet(OpChannel, opdetID, cryostatID);
            const TGeoNode* node = geom->Cryostat(cryostatID).OpDet(opdetID).Node();

            // Get the length of the photon detector
            TGeoBBox *box = (TGeoBBox*)node->GetVolume()->GetShape();
            double opdetLength = box->GetDY();    // X is depth and Z is width

            // Identify the photon detector type
            int pdtype;
            std::string detname = node->GetName();
            if      (detname.find("Bar") != std::string::npos )   pdtype = 0;
            else if (detname.find("Fiber") != std::string::npos ) pdtype = 1;
            else if (detname.find("Plank") != std::string::npos ) pdtype = 2;
            else {
                pdtype = -1;
                mf::LogWarning("LBNEOpDetResponse") << "OpChannel: " << OpChannel << " is an unknown PD type named: " << detname 
                                                   << ". Assuming no attenuation.";
            }


            // Assume sipm is at the +y end of the sipm
            double sipmDistance = opdetLength - Phot.FinalLocalPosition.y();

            if (pdtype == 0) {
                // Now assume Bar and Fiber have the same behavior, modify later when more is known
                double normalize   = 0.6719; // Normalize mean performance to be the same in all PDs
                double lambdaShort =  5.56; // cm
                double fracShort   =  0.40 * normalize;
                double lambdaLong  = 44.13; // cm
                double fracLong    =  0.60 * normalize;

                // Throw away some photons based on attenuation
                double AttenuationProb = fracShort*exp(-sipmDistance/lambdaShort) + fracLong*exp(-sipmDistance/lambdaLong);
                if ( CLHEP::RandFlat::shoot(1.0) > AttenuationProb ) return false;

          
                
                  //mf::LogVerbatim("LBNEOpDetResponse") << "OpChannel: " << OpChannel << " is a " << typenames[pdtype] 
                  //<< " with length " << opdetLength << ", width " << opdetWidth << ", and depth " << opdetDepth;
                  //mf::LogVerbatim("LBNEOpDetResponse") << "   Local Position = (" << Phot.FinalLocalPosition.x() 
                  //<< ", " << Phot.FinalLocalPosition.y() << ", " << Phot.FinalLocalPosition.z() << ")";
                  //mf::LogVerbatim("LBNEOpDetResponse") << "   Distance to SiPM = " << sipmDistance;
                  //mf::LogVerbatim("LBNEOpDetResponse") << "   Attenuation Probability = " << AttenuationProb;
               
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

        }

        return true;
    }

    //--------------------------------------------------------------------
    bool LBNEOpDetResponse::doDetectedLite(int OpChannel) const
    {
        // Check QE
        if ( CLHEP::RandFlat::shoot(1.0) > fQE ) return false;

        return true;
    }



} // namespace

DEFINE_ART_SERVICE_INTERFACE_IMPL(opdet::LBNEOpDetResponse, opdet::OpDetResponseInterface)

