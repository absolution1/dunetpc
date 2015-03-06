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
        
        art::ServiceHandle<geo::Geometry> geom;
        Nchannels = geom->NOpChannels();


        if (fFullSimChannelConvert || fFastSimChannelConvert) {
            int nReadout = 0;
            
            for (int gChannel = 0; gChannel < Nchannels; gChannel++) {
                if (gChannel == 2) opChannelMap.push_back( {nReadout++, nReadout++} );  // Plank (channel 2) has only 2 SiPM's
                else               opChannelMap.push_back( {nReadout++, nReadout++, nReadout++} );
            }
            Nchannels = nReadout;

            // Print out the channel map.  Only once per job.
            PrintChannelMap();
        }
    }
    
    //--------------------------------------------------------------------
    LBNEOpDetResponse::~LBNEOpDetResponse() throw()
    { }


    //--------------------------------------------------------------------
    void LBNEOpDetResponse::doReconfigure(fhicl::ParameterSet const& pset)
    {
        fQE =                    pset.get<double>("QuantumEfficiency");
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
    }


    //--------------------------------------------------------------------
    int  LBNEOpDetResponse::doReadoutToGeoChannel(int readoutChannel) const
    {
        if (!fFullSimChannelConvert && !fFastSimChannelConvert) {
            mf::LogError("LBNEOpDetResponse") << "Trying to convert a readout channel to a geometry channel, but no conversion is turned on.";
            exit(1);
        }

        for (unsigned int g = 0; g < opChannelMap.size(); g++) {
            // Search for readoutChannel in this channel map vector
            auto itr = std::find(opChannelMap[g].begin(), opChannelMap[g].end(), readoutChannel);
            if (itr != opChannelMap[g].end())
                return g;
        }

        PrintChannelMap();
        mf::LogError("LBNEOpDetResponse") << "Readout channel " << readoutChannel << " was not found in the above channel map.";
        exit(2);
        return -1;
    }


    //--------------------------------------------------------------------
    bool LBNEOpDetResponse::doDetected(int OpChannel, const sim::OnePhoton& Phot, int &newOpChannel) const
    {
        if (fFullSimChannelConvert){
            int i = (int) ( CLHEP::RandFlat::shoot(1.0) * (float)opChannelMap[OpChannel].size() );
            newOpChannel = opChannelMap[OpChannel][i];
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
            // Find the Optical Detector using geom
            unsigned int cryostatID=0, opdetID=0;
            art::ServiceHandle<geo::Geometry> geom;
            geom->OpChannelToCryoOpDet(OpChannel, opdetID, cryostatID);
            const TGeoNode* node = geom->Cryostat(cryostatID).OpDet(opdetID).Node();

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

            // Identify the photon detector type
            int pdtype;
            std::string detname = node->GetName();
            boost::to_lower(detname);
            if      (detname.find("bar") != std::string::npos )   pdtype = 0;
            else if (detname.find("fiber") != std::string::npos ) pdtype = 1;
            else if (detname.find("plank") != std::string::npos ) pdtype = 2;
            else {
                pdtype = -1;
                mf::LogWarning("LBNEOpDetResponse") << "OpChannel: " << OpChannel << " is an unknown PD type named: " << detname 
                                                   << ". Assuming no attenuation.";
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

        }

        return true;
    }

    //--------------------------------------------------------------------
    bool LBNEOpDetResponse::doDetectedLite(int OpChannel, int &newOpChannel) const
    {
        if (fFastSimChannelConvert){
            int i = (int) ( CLHEP::RandFlat::shoot(1.0) * (float)opChannelMap[OpChannel].size() );
            newOpChannel = opChannelMap[OpChannel][i];
        }
        else{
            newOpChannel = OpChannel;
        }
        
        // Check QE
        if ( CLHEP::RandFlat::shoot(1.0) > fQE ) return false;

        return true;
    }


    //--------------------------------------------------------------------
    void LBNEOpDetResponse::PrintChannelMap() const
    {
        //mf::LogInfo("LBNEOpDetResponse") << "Converting optical channel numbers in " << fChannelConversion << " simulation" << std::endl;

        std::cout << "LBNE OpDetResponse channel map:" << std::endl;
        for (unsigned int g = 0; g < opChannelMap.size(); g++) {
            std::cout << "  " <<  g << " -> ";
            for (unsigned int d = 0; d < opChannelMap[g].size(); d++) {
                std::cout <<  opChannelMap[g][d] << " ";
            }
            std::cout << std::endl;
        }
    }



} // namespace

DEFINE_ART_SERVICE_INTERFACE_IMPL(opdet::LBNEOpDetResponse, opdet::OpDetResponseInterface)

