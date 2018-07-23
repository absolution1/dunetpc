///////////////////////////////////////////////////////////////////////
///
/// \file   SignalShapingServiceDUNEDPhase.h
///
/// \brief  Service to provide DUNE dual-phase signal shaping for
///         simulation (convolution) and reconstruction (deconvolution).
///
/// \author V. Galymov vgalymov@ipnl.in2p3.fr
///
/// This service inherits from SignalShaping and supplies
/// DUNE dual-phase specific configuration.
/// It is intended that SimWire and CalWire modules will access this service.
///
///
////////////////////////////////////////////////////////////////////////

#ifndef __SIGNALSHAPINGSERVICEDUNEDPHASE_H__
#define __SIGNALSHAPINGSERVICEDUNEDPHASE_H__

#include <vector>
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "lardata/Utilities/SignalShaping.h"
#include "TF1.h"

namespace util {

  class SignalShapingServiceDUNEDPhase {
  public:

    // Constructor, destructor.
    SignalShapingServiceDUNEDPhase(const fhicl::ParameterSet& pset,
				   art::ActivityRegistry& reg);
    ~SignalShapingServiceDUNEDPhase();

    // Update configuration parameters.
    void reconfigure(const fhicl::ParameterSet& pset);

    double GetASICGain(unsigned int const channel) const;
    double GetShapingTime(unsigned int const channel) const;

    double GetRawNoise(unsigned int const channel) const ;
    double GetDeconNoise(unsigned int const channel) const;

    double GetDeconNorm() const {return fDeconNorm;};

    // Accessors.
    //int FieldResponseTOffset(unsigned int const channel) const;

    const util::SignalShaping& SignalShaping(unsigned int channel) const;

    // Do convolution calcution (for simulation).
    template <class T> void Convolute(unsigned int channel, std::vector<T>& func) const;

    // Do deconvolution calcution (for reconstruction).
    template <class T> void Deconvolute(unsigned int channel, std::vector<T>& func) const;

  private:

    // Private configuration methods.

    // Post-constructor initialization.
    void init() const{const_cast<SignalShapingServiceDUNEDPhase*>(this)->init();}
    void init();


    // Calculate response functions.
    void SetElectResponse(std::vector<double> &eresp1meter, std::vector<double> &eresp3meter);
    void SetFieldResponse(std::vector<double> &fresp); //<- could be added later
    double FieldResponse(double tval_us); //harcoded function to calcuate the field response

    // Calculate filter functions.
    void SetFilters();

    // Sample the response function, including a configurable drift velocity of electrons
    void SetResponseSampling();

    // Attributes.
    bool fInit;               ///< Initialization flag.

    // Fcl parameters.
    double fDeconNorm;
    double fASICmVperfC1Meter;              ///< Amplifier gain per 1 fC input for 1 meter long readout channels
    double fASICmVperfC3Meter;              ///< Amplifier gain per 1 fC input for 3 meter long readout channels
    double fADCpermV;                       ///< Digitizer scale
    double fAmpENC;                         ///< Amplifier noise
    double fAmpENCADC;                      ///  Amplifier noise in ADC
    double fRespSamplingPeriod;             ///< Sampling period for response in ns

    std::string fCollection1MeterShape;
    std::string fCollection3MeterShape;

    std::vector<double> fCollection1MeterParams;
    std::vector<double> fCollection3MeterParams;

    TF1* fCollection1MeterFunction;             // Parameterized collection shaping function for 1 meter long readout channels
    TF1* fCollection3MeterFunction;             // Parameterized collection shaping function for 3 meter long readout channels

    TF1* fColFieldFunc;             ///<Parameterized collection field function.
    double fColFieldRespAmp;  			///< amplitude of response to field
    bool fAddFieldFunction;         ///< Enable the filed function

    TF1* fColFilterFunc;      			///< Parameterized collection filter function.
    // Following attributes hold the convolution and deconvolution kernels
    util::SignalShaping fColSignalShaping1Meter;
    util::SignalShaping fColSignalShaping3Meter;

    // Field response.

    // Filters.

    std::vector<TComplex> fColFilter;

    std::vector<double>fParArray;
  };
}

//----------------------------------------------------------------------
// Do convolution.
template <class T> inline void util::SignalShapingServiceDUNEDPhase::Convolute(unsigned int channel,     std::vector<T>& func) const
{
  SignalShaping(channel).Convolute(func);
}



//----------------------------------------------------------------------
// Do deconvolution.
template <class T> inline void util::SignalShapingServiceDUNEDPhase::Deconvolute(unsigned int channel,
std::vector<T>& func) const
{
  SignalShaping(channel).Deconvolute(func);
}


DECLARE_ART_SERVICE(util::SignalShapingServiceDUNEDPhase, LEGACY)
#endif
