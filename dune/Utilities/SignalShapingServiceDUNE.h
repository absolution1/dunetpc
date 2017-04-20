///////////////////////////////////////////////////////////////////////
///
/// \file   SignalShapingServiceDUNE.h
///
/// \brief  Service to provide microboone-specific signal shaping for
///         simulation (convolution) and reconstruction (deconvolution).
///
/// \author H. Greenlee (original class)
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

#ifndef SIGNALSHAPINGSERVICEDUNE_H
#define SIGNALSHAPINGSERVICEDUNE_H

#include "dune/DuneInterface/SignalShapingService.h"

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "lardata/Utilities/SignalShaping.h"
#include "TF1.h"
#include "TH1D.h"

using DoubleVec = std::vector<double>;
namespace util {

class SignalShapingServiceDUNE : public SignalShapingService {
public:

  // Constructor, destructor.

  SignalShapingServiceDUNE(const fhicl::ParameterSet& pset,
				   art::ActivityRegistry& reg);
  ~SignalShapingServiceDUNE();

  // Update configuration parameters.

  void reconfigure(const fhicl::ParameterSet& pset);

  std::vector<DoubleVec> GetNoiseFactVec() const override;
  double GetASICGain(Channel channel) const override;
  double GetShapingTime(Channel channel) const override; 
  double GetRawNoise(Channel channel) const override;
  double GetDeconNoise(Channel channel) const override;
  double GetDeconNorm() const override;
  unsigned int GetSignalSize() const override;

  // Accessors.
  int FieldResponseTOffset(Channel const channel) const override;
  const util::SignalShaping& SignalShaping(Channel channel) const override;

  // Do convolution calcution (for simulation).
  template <class T> void Convolute(Channel channel, std::vector<T>& func) const;
  void Convolute(Channel channel, FloatVector& func) const override;
  void Convolute(Channel channel, DoubleVector& func) const override;

  // Do deconvolution calcution (for reconstruction).
  template <class T> void Deconvolute(Channel channel, std::vector<T>& func) const;
  void Deconvolute(Channel channel, DoubleVector& func) const override;
  void Deconvolute(Channel channel, FloatVector& func) const override;

private:

  // Private configuration methods.

  // Post-constructor initialization.
  void init() const{const_cast<SignalShapingServiceDUNE*>(this)->init();}
  void init();

  // Calculate response functions.
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
  double fInd3DCorrection;  			///< correction factor to account for 3D path of 
  double fColFieldRespAmp;  			///< amplitude of response to field 
  double fIndUFieldRespAmp;  			///< amplitude of response to field  
  double fIndVFieldRespAmp;  			///< amplitude of response to field
  std::vector<double> fFieldResponseTOffset;  ///< Time offset for field response in ns
  std::vector<double> fCalibResponseTOffset;
  double fInputFieldRespSamplingPeriod;       ///< Sampling period in the input field response.
  double fDeconNorm;
  double fADCPerPCAtLowestASICGain; ///< Pulse amplitude gain for a 1 pc charge impulse after convoluting it the with field and electronics response with the lowest ASIC gain setting of 4.7 mV/fC
  std::vector<DoubleVec> fNoiseFactVec; 
  std::vector<double> fASICGainInMVPerFC;    
  std::vector<double> fShapeTimeConst;  	///< time constants for exponential shaping
  TF1* fColFilterFunc;      			///< Parameterized collection filter function.
  TF1* fIndUFilterFunc;      			///< Parameterized induction filter function.
  TF1* fIndVFilterFunc;      			///< Parameterized induction filter function.
  bool fUseFunctionFieldShape;   		///< Flag that allows to use a parameterized field response instead of the hardcoded version
  bool fUseHistogramFieldShape;               ///< Flag that turns on field response shapes from histograms
  bool fGetFilterFromHisto;   		///< Flag that allows to use a filter function from a histogram instead of the functional dependency
  TF1* fColFieldFunc;      			///< Parameterized collection field shape function.
  TF1* fIndUFieldFunc;      			///< Parameterized induction field shape function.
  TF1* fIndVFieldFunc;      			///< Parameterized induction field shape function.
  TH1F *fFieldResponseHist[3];                ///< Histogram used to hold the field response, hardcoded for the time being
  TH1D *fFilterHist[3];    			///< Histogram used to hold the collection filter, hardcoded for the time being
    
  // Following attributes hold the convolution and deconvolution kernels
  util::SignalShaping fColSignalShaping;
  util::SignalShaping fIndUSignalShaping;
  util::SignalShaping fIndVSignalShaping;

  // Field response.
  std::vector<double> fColFieldResponse;
  std::vector<double> fIndUFieldResponse;
  std::vector<double> fIndVFieldResponse;

  // Electronics response.
  std::vector<double> fElectResponse;

  // Filters.
  std::vector<TComplex> fColFilter;
  std::vector<TComplex> fIndUFilter;
  std::vector<TComplex> fIndVFilter;

};

}  // end namespace util

//----------------------------------------------------------------------

// Do convolution.
template <class T>
inline void util::SignalShapingServiceDUNE::
Convolute(unsigned int channel, std::vector<T>& func) const {
  SignalShaping(channel).Convolute(func);

  //negative number
  int time_offset = FieldResponseTOffset(channel);
  
  std::vector<T> temp;
  if (time_offset <=0){
    temp.assign(func.begin(),func.begin()-time_offset);
    func.erase(func.begin(),func.begin()-time_offset);
    func.insert(func.end(),temp.begin(),temp.end());
  }else{
    temp.assign(func.end()-time_offset,func.end());
    func.erase(func.end()-time_offset,func.end());
    func.insert(func.begin(),temp.begin(),temp.end());
  }
}

//----------------------------------------------------------------------

// Do deconvolution.
template <class T>
inline void util::SignalShapingServiceDUNE::
Deconvolute(unsigned int channel, std::vector<T>& func) const {
  SignalShaping(channel).Deconvolute(func);

  int time_offset = FieldResponseTOffset(channel);
  
  std::vector<T> temp;
  if (time_offset <=0){
    temp.assign(func.end()+time_offset,func.end());
    func.erase(func.end()+time_offset,func.end());
     func.insert(func.begin(),temp.begin(),temp.end());
  }else{
    temp.assign(func.begin(),func.begin()+time_offset);
    func.erase(func.begin(),func.begin()+time_offset);
    func.insert(func.end(),temp.begin(),temp.end());    
  }
}

//----------------------------------------------------------------------

DECLARE_ART_SERVICE(util::SignalShapingServiceDUNE, LEGACY)

#endif
