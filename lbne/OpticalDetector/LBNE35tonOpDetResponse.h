////////////////////////////////////////////////////////////////////////
// \file LBNE35tonOpDetResponse.h
//
// \brief service containing information about the response of optical detectors in LBNE35ton
//
// \author ahimmel@phy.duke.edu
//
////////////////////////////////////////////////////////////////////////

#ifndef LBNE35ton_OPDET_RESPONSE_H
#define LBNE35ton_OPDET_RESPONSE_H

// LArSoft includes
#include "Simulation/SimPhotons.h"
#include "OpticalDetector/OpDetResponseInterface.h"



namespace opdet
{
    class LBNE35tonOpDetResponse : public opdet::OpDetResponseInterface {
    public:

        LBNE35tonOpDetResponse(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
        ~LBNE35tonOpDetResponse() throw();



    private:

        virtual void doReconfigure(fhicl::ParameterSet const& p);

        virtual int  doNOpChannels() const;
        virtual bool doDetected(int OpChannel, const sim::OnePhoton& Phot, int &newOpChannel) const;
        virtual bool doDetectedLite(int OpChannel, int &newOpChannel) const;

        float fQE;                     // Quantum efficiency of tube
        
        float fWavelengthCutLow;       // Sensitive wavelength range 
        float fWavelengthCutHigh;      // 
        
        bool fLightGuideAttenuation;   // Flag to turn on position-dependent sensitivity

        std::string fChannelConversion;
        bool fFullSimChannelConvert;   // Flag to conver detector->electronics channels in full optical sim
        bool fFastSimChannelConvert;   // Flag to conver detector->electronics channels in fast optical sim

        int fLongAxis;                 // 0 = x, 1 = y, 2 = z

    }; // class LBNE35tonOpDetResponse

    
} //namespace opdet


DECLARE_ART_SERVICE_INTERFACE_IMPL(opdet::LBNE35tonOpDetResponse, opdet::OpDetResponseInterface, LEGACY)

#endif //OPDET_RESPONSE_LBNE35ton_H
