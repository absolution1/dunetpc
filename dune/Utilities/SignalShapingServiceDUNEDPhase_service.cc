////////////////////////////////////////////////////////////////////////
/// \file   SignalShapingServiceDUNEDPhase_service.cc
/// \author vgalymov@ipnl.in2p3.fr
////////////////////////////////////////////////////////////////////////

#include "dune/Utilities/SignalShapingServiceDUNEDPhase.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/TPCGeo.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/LArFFT.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "TSpline.h"

//#include "TFile.h"

//----------------------------------------------------------------------
// Constructor.
util::SignalShapingServiceDUNEDPhase::SignalShapingServiceDUNEDPhase(const fhicl::ParameterSet& pset,
								    art::ActivityRegistry& /* reg */) 
  : fInit(false)
{
  reconfigure(pset);
}


//----------------------------------------------------------------------
// Destructor.
util::SignalShapingServiceDUNEDPhase::~SignalShapingServiceDUNEDPhase()
{}


//----------------------------------------------------------------------
// Reconfigure method.
void util::SignalShapingServiceDUNEDPhase::reconfigure(const fhicl::ParameterSet& pset)
{
  // Reset initialization flag.
  fInit = false;

  // Reset kernels.
  fColSignalShaping.Reset();

  // Fetch fcl parameters.
  fDeconNorm            = pset.get<double>("DeconNorm");
  fASICmVperfC          = pset.get<double>("ASICmVperfC");
  fADCpermV             = pset.get<double>("ADCpermV");
  fAmpENC               = pset.get<double>("AmpENC");
  
  fRespSamplingPeriod   = pset.get<double>("RespSamplingPeriod");

  // amplifier noise in ADC
  fAmpENCADC            = fAmpENC * 1.60217657e-4 * fASICmVperfC * fADCpermV;

  // Construct parameterized collection filter function.
  
  //LOG_DEBUG("SignalShapingServiceDUNEDPhase") <<"ASIC Gain in mV per fC = "<<fASICmVperfC;
  mf::LogInfo("SignalShapingServiceDUNEDPhase") <<"ASIC Gain in mV per fC = "<<fASICmVperfC
						<<";  ADC conversion = "<<fADCpermV
						<<";  Amplifier ENC = "<<fAmpENC
						<<";  Amplifier ENC in ADC = "<<fAmpENCADC;
  
   mf::LogInfo("SignalShapingServiceDUNE") << "Getting Filter from .fcl file";
  std::string colFilt = pset.get<std::string>("ColFilter");
  std::vector<double> colFiltParams =
  pset.get<std::vector<double> >("ColFilterParams");
  fColFilterFunc = new TF1("colFilter", colFilt.c_str());
  for(unsigned int i=0; i<colFiltParams.size(); ++i)
    fColFilterFunc->SetParameter(i, colFiltParams[i]);
}

//----------------------------------------------------------------------
// Accessor for single-plane signal shaper.
const util::SignalShaping&
util::SignalShapingServiceDUNEDPhase::SignalShaping(unsigned int channel) const
{
  if(!fInit)
    init();

  auto const* geom = lar::providerFrom<geo::Geometry>();
  switch (geom->SignalType(channel)) {
    case geo::kCollection:
      return fColSignalShaping; //always collections in DP detector
    default:
      throw cet::exception("SignalShapingServiceDUNEDPhase")
        << "unexpected signal type " << geom->SignalType(channel)
        << " for channel #" << channel << "\n";
  } // switch
  
} // util::SignalShapingServiceDUNEDPhase::SignalShaping()

//----------------------------------------------------------------------
// Initialization method.
// Here we do initialization that can't be done in the constructor.
// All public methods should ensure that this method is called as necessary.
void util::SignalShapingServiceDUNEDPhase::init()
{
  if(!fInit) {
    fInit = true;
    
    // Calculate field and electronics response functions.
    std::vector<double> eresp;
    SetElectResponse(eresp);

    // Configure convolution kernels.
    fColSignalShaping.AddResponseFunction(eresp);
    fColSignalShaping.save_response();
    fColSignalShaping.set_normflag(false);
    
    
    // rebin to appropriate sampling rate of readout
    // NOTE: could have done it from the start in the eresp calculation
    //       but implemented it like this for future flexibility
    SetResponseSampling();

    // Calculate filter functions.

    SetFilters();

    // Configure deconvolution kernels.
    fColSignalShaping.AddFilterFunction(fColFilter);
    fColSignalShaping.CalculateDeconvKernel();
    
    // std::cout << " fColSignalShaping: " << fColSignalShaping
    //mf::LogInfo("SignalShapingServiceDUNEDPhase")<<"Done Init";
  }
}

//-----Give Gain Settings to SimWire-----//
double util::SignalShapingServiceDUNEDPhase::GetASICGain(unsigned int const channel) const
{
  art::ServiceHandle<geo::Geometry> geom;
  //geo::SigType_t sigtype = geom->SignalType(channel);

   // we need to distinguis between the U and V planes
  geo::View_t view = geom->View(channel); 
  
  double gain = 0;
  if(view == geo::kU || view == geo::kV || view == geo::kZ )
    gain = fASICmVperfC;
  else
    throw cet::exception("SignalShapingServiceDUNEDPhase")<< "can't determine"
							  << " View\n";
  return gain;
}


//-----Give Shaping time to SimWire-----//
double util::SignalShapingServiceDUNEDPhase::GetShapingTime(unsigned int const channel) const
{
  //art::ServiceHandle<geo::Geometry> geom;
  //geo::View_t view = geom->View(channel); 

  return 0.0; //not sure what this does in simwire
}

double util::SignalShapingServiceDUNEDPhase::GetRawNoise(unsigned int const channel) const
{
  return fAmpENCADC; 
}

double util::SignalShapingServiceDUNEDPhase::GetDeconNoise(unsigned int const channel) const
{
  return fAmpENC / fDeconNorm; //expected ENC
}



//----------------------------------------------------------------------
// Calculate field response. Maybe implemented at some point
//void util::SignalShapingServiceDUNEDPhase::SetFieldResponse()
//{}

//----------------------------------------------------------------------
// Calculate electronics response
void util::SignalShapingServiceDUNEDPhase::SetElectResponse(std::vector<double> &ElecResp)
{
  // Get services.

  art::ServiceHandle<util::LArFFT> fft;

  LOG_DEBUG("SignalShapingDUNEDPhase") << "Setting DUNEDPhase electronics response function...";
  
  int nticks = fft->FFTSize();
  std::vector<double> time(nticks,0.);
  ElecResp.resize(nticks, 0.);

  double max = 0;
  for(size_t i = 0; i < ElecResp.size(); ++i)
    {
      //convert time to microseconds, to match response function definition
      time[i]     = (1.*i)*fRespSamplingPeriod*1e-3; 
      ElecResp[i] = PreampETHZ(time[i]);
 
      if(ElecResp[i] > max) max = ElecResp[i];
    }// end loop over time buckets
    

  LOG_DEBUG("SignalShapingDUNEDPhase") << " Done.";
  
  //normalize to 1e charge before the convolution   
  for(auto& element : ElecResp)
    {
      element /= max;
      element *= fASICmVperfC * 1.60217657e-4; //mV
      element *= fADCpermV;                    //ADC
    }
  
  return;

}




//----------------------------------------------------------------------
// Calculate filter functions.
void util::SignalShapingServiceDUNEDPhase::SetFilters()
{ 
  // Get services.
  
  auto const *detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
  art::ServiceHandle<util::LArFFT> fft;

  double ts = detprop->SamplingRate();
  int n = fft->FFTSize() / 2;

  // Calculate collection filter.

  fColFilter.resize(n+1);

  fColFilterFunc->SetRange(0, double(n));

  for(int i=0; i<=n; ++i) 
  {
    	double freq = 400. * i / (ts * n);      // Cycles / microsecond. 
    	double f = fColFilterFunc->Eval(freq);
    	fColFilter[i] = TComplex(f, 0.);
  }

  return;
}


//----------------------------------------------------------------------
// (Re)sample electronics (+field?) response 
void util::SignalShapingServiceDUNEDPhase::SetResponseSampling()
{
  // Get services
  auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
  art::ServiceHandle<util::LArFFT> fft;

  // Operation permitted only if output of rebinning has a larger bin size
  if( fRespSamplingPeriod > detprop->SamplingRate() )
    throw cet::exception(__FUNCTION__) << "\033[93m"
				       << "Invalid operation: cannot rebin to a more finely binned vector!"
				       << "\033[00m" << std::endl;


  if( fRespSamplingPeriod == detprop->SamplingRate() ) //nothing to do
    {
      return;
    }
  
  ///
  // make a new response vector with different time sampling
  
  // default number of time ticks
  int nticks = fft->FFTSize();
  
  // resampling vectors
  std::vector<double> SamplingTime( nticks, 0. );
  std::vector<double> SamplingResp(nticks , 0. );
  for ( int itime = 0; itime < nticks; itime++ ) 
    SamplingTime[itime] = (1.*itime) * detprop->SamplingRate();
  
  // get response our old response vector
  const std::vector<double>* pResp = &(fColSignalShaping.Response_save());

  //
  int nticks_input = pResp->size();

  //old sampling time vector
  std::vector<double> InputTime(nticks_input, 0. );
  std::vector<double> InputResp(nticks_input, 0. );
  for ( int itime = 0; itime < nticks_input; itime++ )
    {
      InputResp[itime] = (*pResp)[itime];
      InputTime[itime] = (1.*itime) * fRespSamplingPeriod;
    }
  
  
  // build a spline for interpolation
  TSpline3 ispl("ispl", &InputTime[0], &InputResp[0], nticks_input);

  // interpolate
  int SamplingCount = 0;
  double maxtime = InputTime.back();
  for(int itime = 0; itime < nticks; itime++)
    {
      // don't go past max t value from old response
      if(SamplingTime[itime] > maxtime) break;
      
      SamplingResp[itime] = ispl.Eval(SamplingTime[itime]);
      SamplingCount++;
    }
  SamplingResp.resize( SamplingCount, 0.);    
  
  // set new response
  fColSignalShaping.AddResponseFunction( SamplingResp, true );

  return;
}

//----------------------------------------------------------------------
// response of ETHZ pre-amplifier
double util::SignalShapingServiceDUNEDPhase::PreampETHZ(double tval_us) 
{
  // parameters
  double dt = 3.5;  //us
  double T1 = 2.83; //us
  double T2 = 0.47; //us
  
  double fval = 1/(dt*(T1-T2)*(T1-T2)) * 
    (exp(-(tval_us)/T2)*( (tval_us)*(T1-T2)+(2*T1-T2)*T2-
		       ((tval_us-dt>=0)?1.:0.) * exp(dt/T2)*((tval_us)*(T1-T2)+2*T1*T2-T2*T2+dt*(T2-T1)) )+
     exp(-(tval_us)/T1)*( ((dt-tval_us> 0)?1.:0.) * exp((tval_us)/T1)*(T1-T2)*(T1-T2)+
		       T1*T1*(exp(dt/T1) * ((tval_us-dt>=0)?1:0) -1)));

  if(fval != fval) fval = 0;
  return fval;  
}



namespace util {

  DEFINE_ART_SERVICE(SignalShapingServiceDUNEDPhase)

}
