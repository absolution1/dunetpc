///////////////////////////////////////////////////////////////////////
///
/// \file   SignalShapingServiceLBNE35t.h
///
/// \brief  Service to provide microboone-specific signal shaping for
///         simulation (convolution) and reconstruction (deconvolution).
///
/// \author H. Greenlee 
///
/// This service inherits from SignalShaping and supplies
/// microboone-specific configuration.  It is intended that SimWire and
/// CalWire modules will access this service.
///
/// FCL parameters:
///
/// FieldBins       - Number of bins of field response.
/// Col3DCorrection - 3D path length correction for collection plane.
/// Ind3DCorrection - 3D path length correction for induction plane.
/// ColFieldRespAmp - Collection field response amplitude.
/// IndFieldRespAmp - Induction field response amplitude.
/// ShapeTimeConst  - Time constants for exponential shaping.
/// ColFilter       - Root parameterized collection plane filter function.
/// ColFilterParams - Collection filter function parameters.
/// IndFilter       - Root parameterized induction plane filter function.
/// IndFilterParams - Induction filter function parameters.
///
////////////////////////////////////////////////////////////////////////

#ifndef SIGNALSHAPINGSERVICELBNE35T_H
#define SIGNALSHAPINGSERVICELBNE35T_H

#include <vector>
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "Utilities/SignalShaping.h"
#include "TF1.h"
#include "TH1D.h"

using DoubleVec = std::vector<double>;

namespace util {
  class SignalShapingServiceLBNE35t {
  public:

    // Constructor, destructor.

    SignalShapingServiceLBNE35t(const fhicl::ParameterSet& pset,
				   art::ActivityRegistry& reg);
    ~SignalShapingServiceLBNE35t();

    // Update configuration parameters.

    void reconfigure(const fhicl::ParameterSet& pset);

    std::vector<DoubleVec> GetNoiseFactVec()                { return fNoiseFactVec; }
    double GetASICGain(unsigned int const channel) const;
    double GetShapingTime(unsigned int const channel) const; 

    double GetRawNoise(unsigned int const channel) const ;
    double GetDeconNoise(unsigned int const channel) const;   
    
    // Accessors.

    const util::SignalShaping& SignalShaping(unsigned int channel) const;

    int FieldResponseTOffset(unsigned int const channel) const;

    // Do convolution calcution (for simulation).

    template <class T> void Convolute(unsigned int channel, std::vector<T>& func) const;

    // Do deconvolution calcution (for reconstruction).

    template <class T> void Deconvolute(unsigned int channel, std::vector<T>& func) const;
    
    double GetDeconNorm(){return fDeconNorm;};

  private:

    // Private configuration methods.

    // Post-constructor initialization.

    void init() const{const_cast<SignalShapingServiceLBNE35t*>(this)->init();}
    void init();

    // Calculate response functions.
    // Copied from SimWireLBNE35t.

    void SetFieldResponse();
    void SetElectResponse(double shapingtime, double gain);

    // Calculate filter functions.

    void SetFilters();

    // Attributes.

    bool fInit;               ///< Initialization flag.

    // Sample the response function, including a configurable
    // drift velocity of electrons
    
    void SetResponseSampling();

    // Fcl parameters.

    int fNFieldBins;         			///< number of bins for field response
    double fCol3DCorrection; 			///< correction factor to account for 3D path of 
						///< electrons thru wires
    double fInd3DCorrection;  			///< correction factor to account for 3D path of 
						///< electrons thru wires
    double fColFieldRespAmp;  			///< amplitude of response to field 
    double fIndFieldRespAmp;  			///< amplitude of response to field 
    
    std::vector<double> fFieldResponseTOffset;  ///< Time offset for field response in ns
    double fInputFieldRespSamplingPeriod;       ///< Sampling period in the input field response. 

    double fDeconNorm;
    double fADCPerPCAtLowestASICGain; ///< Pulse amplitude gain for a 1 pc charge impulse after convoluting it the with field and electronics response with the lowest ASIC gain setting of 4.7 mV/fC
    std::vector<DoubleVec> fNoiseFactVec; 
    
    std::vector<double> fASICGainInMVPerFC;    

    std::vector<double> fShapeTimeConst;  	///< time constants for exponential shaping
    TF1* fColFilterFunc;      			///< Parameterized collection filter function.
    TF1* fIndFilterFunc;      			///< Parameterized induction filter function.

    
    bool fUseFunctionFieldShape;   		///< Flag that allows to use a parameterized field response instead of the hardcoded version
    bool fUseHistogramFieldShape;               ///< Flag that turns on field response shapes from histograms
    bool fGetFilterFromHisto;   		///< Flag that allows to use a filter function from a histogram instead of the functional dependency
    TF1* fColFieldFunc;      			///< Parameterized collection field shape function.
    TF1* fIndFieldFunc;      			///< Parameterized induction field shape function.
    
    TH1F *fFieldResponseHist[3];                ///< Histogram used to hold the field response, hardcoded for the time being
    TH1D *fFilterHist[3];    			///< Histogram used to hold the collection filter, hardcoded for the time being
    
    // Following attributes hold the convolution and deconvolution kernels

    util::SignalShaping fColSignalShaping;
    util::SignalShaping fIndSignalShaping;

    // Field response.

    std::vector<double> fColFieldResponse;
    std::vector<double> fIndFieldResponse;

    // Electronics response.

    std::vector<double> fElectResponse;

    // Filters.

    std::vector<TComplex> fColFilter;
    std::vector<TComplex> fIndFilter;
  };
}
//----------------------------------------------------------------------
// Do convolution.
template <class T> inline void util::SignalShapingServiceLBNE35t::Convolute(unsigned int channel, std::vector<T>& func) const
{
  SignalShaping(channel).Convolute(func);
}


//----------------------------------------------------------------------
// Do deconvolution.
template <class T> inline void util::SignalShapingServiceLBNE35t::Deconvolute(unsigned int channel, std::vector<T>& func) const
{
  SignalShaping(channel).Deconvolute(func);
}

DECLARE_ART_SERVICE(util::SignalShapingServiceLBNE35t, LEGACY)
#endif
