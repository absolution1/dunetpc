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
/// The general philosophy to follow is to split fC -> ADC conversion 
/// in two parts
///  - area normalization factor (typicaly in units ADC x us / fC)
///  - pulse shaping
/// Since the hit charge should be calculated from the integral, the norm
/// factor should be an area normalization constant. The shape function 
/// should then be area-normalized in order to eventually achive a proper 
/// normalization of the hit integral. Note that in this case even if the 
/// simulated pulse shape may not be exact, if correctly normalized the 
/// resulting simulated hit charge should be correct. 
/// The charge loss due to, for example, the detector capacitance should 
/// not be included in the calibration constant here but treated elsewhere 
/// to avoid lumping many effiiciency factors into one ambiguous number
///
/// changes:
///   11.20.2019, vgalymov
///     Changed pulse normalization to integral instead of amplitude. 
///     For 3x1x1 the parameters specified in fcl were computed based 
///     on what was originally put there for 3x1x1. Namely the pulse 
///     function was calculated with original parameters and its integral 
///     was then used to specify the area normalization in the service fcl.
///     Make the parameters for 1 m view optional, since this corresponds to
///     to a rather particular CRP configuration
///
///    02.17.2020, vgalymov
///     Updated filter function to sample up-to correct 1.25 MHz
///     This service now inherits from SignalShapingService
///
///        
////////////////////////////////////////////////////////////////////////

#ifndef __SIGNALSHAPINGSERVICEDUNEDPHASE_H__
#define __SIGNALSHAPINGSERVICEDUNEDPHASE_H__

#include "dune/DuneInterface/Service/SignalShapingService.h"

#include <vector>
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "lardata/Utilities/SignalShaping.h"
#include "TF1.h"

namespace detinfo {
  class DetectorClocksData;
}

using DoubleVec = std::vector<double>;

namespace util {

  class SignalShapingServiceDUNEDPhase : public SignalShapingService {
  public:

    // Constructor, destructor.
    SignalShapingServiceDUNEDPhase(const fhicl::ParameterSet& pset,
				   art::ActivityRegistry& reg);
    ~SignalShapingServiceDUNEDPhase();

    // Update configuration parameters.
    void reconfigure(const fhicl::ParameterSet& pset);

    double GetAreaNorm(unsigned int const channel) const;

    std::vector<DoubleVec> GetNoiseFactVec() const override;
    double GetASICGain(Channel channel) const override;
    double GetShapingTime(Channel channel) const override; 
    double GetRawNoise(Channel channel) const override;
    double GetDeconNoise(Channel channel) const override;
    double GetDeconNorm() const override;
    unsigned int GetSignalSize() const override;

    //double GetDeconNorm() const {return fDeconNorm;};

    // Accessors.
    int FieldResponseTOffset(detinfo::DetectorClocksData const& clockData,
                             unsigned int const channel) const override;
    const util::SignalShaping& SignalShaping(unsigned int channel) const override;

    // Do convolution calcution (for simulation).
    template <class T> void Convolute(unsigned int channel, std::vector<T>& func) const;
    void Convolute(detinfo::DetectorClocksData const& clockData, Channel channel, FloatVector& func) const override;
    void Convolute(detinfo::DetectorClocksData const& clockData, Channel channel, DoubleVector& func) const override;
    
    // Do deconvolution calcution (for reconstruction).
    template <class T> void Deconvolute(unsigned int channel, std::vector<T>& func) const;
    void Deconvolute(detinfo::DetectorClocksData const& clockData, Channel channel, DoubleVector& func) const override;
    void Deconvolute(detinfo::DetectorClocksData const& clockData, Channel channel, FloatVector& func) const override;

  private:

    // Private configuration methods.

    // Post-constructor initialization.
    void init() const{const_cast<SignalShapingServiceDUNEDPhase*>(this)->init();}
    void init();


    // Calculate response functions.
    double SetElectResponse( const TF1* fshape, util::SignalShaping &sig, double areanorm );
    void   SetFieldResponse(std::vector<double> &fresp); //<- could be added later
    double FieldResponse(double tval_us); //harcoded function to calcuate the field response

    // Calculate filter functions.
    void SetFilters(detinfo::DetectorClocksData const& clockData);

    // Sample the response function, including a configurable drift velocity of electrons
    void SetResponseSampling(detinfo::DetectorClocksData const& clockData,
                             util::SignalShaping &sig );

    // Attributes.
    bool fInit;               ///< Initialization flag.

    // Fcl parameters.
    double      fDeconNorm;
    double      fRespSamplingPeriod;            ///< Sampling period for response in ns
    double      fAreaNorm;                      ///< in units (ADC x us) / fC
    double      fPulseHeight;                   ///< pulse height in mV per fC
    std::string fShapeFunc;
    
    double fADCpermV;                           ///< Digitizer scale
    double fAmpENC;                             ///< noise in num of electrons used where???
    double fAmpENCADC;                          ///< noise rescaled to ADC

    
    // the 3x1x1 stuff
    bool        fHave1mView; 
    double      fAreaNorm1m;                ///< in units (ADC x us) / fC for 1m strips
    double      fPulseHeight1m;             ///< pulse amtplitide in mV per fC
    std::string fShapeFunc1m;

    
    /// 
    TF1* fColShapeFunc;      // Parameterized collection shaping function 
    TF1* fColShapeFunc1m;    // Parameterized collection shaping function for 1m 3x1x1 strips

    // Field response
    bool   fAddFieldFunction;         ///< Enable the filed function
    TF1*   fColFieldFunc;             ///< Parameterized collection field function.
    double fColFieldRespAmp;  	      ///< amplitude of response to field
    
    //
    // Following attributes hold the convolution and deconvolution kernels
    util::SignalShaping fColSignalShaping;
    
    //
    util::SignalShaping fColSignalShaping1m; /// for 3x1x1 1m strips
    
    // Filters.
    TF1* fColFilterFunc;      	      ///< Parameterized collection filter function.
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
