////////////////////////////////////////////////////////////////////////
/// \file   SignalShapingServiceDUNEDPhase_service.cc
/// \author vgalymov@ipnl.in2p3.fr
////////////////////////////////////////////////////////////////////////

#include "dune/Utilities/SignalShapingServiceDUNEDPhase.h"

#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
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

  // Fetch fcl parameters for normal 3m view
  fAreaNorm             = pset.get<double>("AreaNorm");
  fADCpermV             = pset.get<double>("ADCpermV");
  fAmpENC               = pset.get<double>("AmpENC");
  fDeconNorm            = pset.get<double>("DeconNorm");
  fRespSamplingPeriod   = pset.get<double>("RespSamplingPeriod");
  fShapeFunc            = pset.get<std::string>("ShapeFunc");
  
  fColShapeFunc = new TF1("fColShapeFunc", fShapeFunc.c_str());
  std::vector<double> params = pset.get<std::vector<double>>("ShapeFuncParams");
  for( unsigned i=0;i<params.size();++i ){
    //std::cout<<i<<" "<<params[i]<<std::endl;
    fColShapeFunc->SetParameter(i, params[i]);
  }
  
  // now check if there is anything specific for 3x1x1 1m view
  fColSignalShaping1m.Reset();
  fHave1mView          = (pset.has_key("AreaNorm1m") || pset.has_key( "ShapeFunc1m" ));
  fAreaNorm1m          = pset.get<double>("AreaNorm1m", fAreaNorm);
  fShapeFunc1m         = fShapeFunc;
  std::vector<double> params1 = params;
  if( pset.has_key( "ShapeFunc1m" ) )
    {
      fShapeFunc1m    = pset.get<std::string>("ShapeFunc1m");
      params1         = pset.get<std::vector<double>>("ShapeFuncParams1m");
    }
  fColShapeFunc1m = new TF1("fColShapeFunc1m", fShapeFunc1m.c_str());
  for( unsigned i=0;i<params1.size();++i ){
    //std::cout<<i<<" "<<params1[i]<<std::endl;
    fColShapeFunc1m->SetParameter(i, params1[i]);
  }
  
  //
  mf::LogInfo("SignalShapingServiceDUNEDPhase") << "  Area normalization : "<<fAreaNorm << " ADC x us\n"					
						<< "  Shape function : "<<fShapeFunc << "\n"
						<< "  Decon normalization : "<<fDeconNorm <<"\n"
						<< "  Response sampling period : "<< fRespSamplingPeriod<<" us \n"
						<< "  Noise : "<<fAmpENC<<" electrons\n"
						<< "  ADC per mV : "<<fADCpermV;

  if( fHave1mView )
    mf::LogInfo("SignalShapingServiceDUNEDPhase") << "  Parameters for 1m view :\n"
						  << "   Area normalization : "<<fAreaNorm1m << " ADC x us\n"					
						  << "   Shape function : "<<fShapeFunc1m;

  
    
  // Construct parameterized collection filter function.
  mf::LogInfo("SignalShapingServiceDUNE") << "Getting Field response from .fcl file"; //<<< to be changed if we keep the harcoded function
  fAddFieldFunction    = pset.get<bool>("AddFieldFunction");
  std::string colField = pset.get<std::string>("ColFieldFunction");
  std::vector<double> colFieldParams =
    pset.get<std::vector<double> >("ColFieldFunctionParams");
  fColFieldFunc = new TF1("colFieldFunction", colField.c_str());
  for(unsigned int i=0; i<colFieldParams.size(); ++i)
    fColFieldFunc->SetParameter(i, colFieldParams[i]);
  
  mf::LogInfo("SignalShapingServiceDUNE") << "Getting Filter from .fcl file";
  std::string colFilt = pset.get<std::string>("ColFilter");
  std::vector<double> colFiltParams =
    pset.get<std::vector<double> >("ColFilterParams");
  fColFilterFunc = new TF1("colFilter", colFilt.c_str());
  for(unsigned int i=0; i<colFiltParams.size(); ++i)
    fColFilterFunc->SetParameter(i, colFiltParams[i]);
  fColFieldRespAmp = pset.get<double>("ColFieldRespAmp");

  init();
}

//----------------------------------------------------------------------
// Accessor for single-plane signal shaper.
const util::SignalShaping&
util::SignalShapingServiceDUNEDPhase::SignalShaping(unsigned int channel) const
{
  if(!fInit)
    init();
  
  if( !fHave1mView ) return fColSignalShaping;
  
  //
  auto const* geom = lar::providerFrom<geo::Geometry>();
  std::vector<geo::WireID> Wires =  geom->ChannelToWire(channel);

  double wirestartpoint[3];
  double wireendpoint[3];
  geom->WireEndPoints(Wires[0],wirestartpoint,wireendpoint);
  double wirelength = sqrt(pow(wirestartpoint[0]-wireendpoint[0],2) + pow(wirestartpoint[1]-wireendpoint[1],2) + pow(wirestartpoint[2]-wireendpoint[2],2));

  if((int)wirelength == 300)
    return fColSignalShaping;
  else if((int)wirelength == 100)
    return fColSignalShaping1m;
  else
    throw cet::exception("SignalShapingServiceDUNEDPhase")
        << "unexpected signal type " << geom->SignalType(channel)
        << " for channel #" << channel << "\n";
/*
  switch (geom->SignalType(channel)) {
    case geo::kCollection:
      return fColSignalShaping; //always collections in DP detector
    default:
      throw cet::exception("SignalShapingServiceDUNEDPhase")
        << "unexpected signal type " << geom->SignalType(channel)
        << " for channel #" << channel << "\n";
  } // switch
*/
} // util::SignalShapingServiceDUNEDPhase::SignalShaping()

//----------------------------------------------------------------------
// Initialization method.
// Here we do initialization that can't be done in the constructor.
// All public methods should ensure that this method is called as necessary.
void util::SignalShapingServiceDUNEDPhase::init()
{    
  if( fInit ) return;
  
  fInit = true;
  
  // Lazy initialization is problematic wrt multi-threading.  For this
  // case, the best we can do is use the global data as provided by
  // the detector-clocks service.
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();

  //
  fPulseHeight = SetElectResponse( fColShapeFunc, fColSignalShaping, fAreaNorm );
  // noise in ADC
  fAmpENCADC = fAmpENC * 1.60217657e-4 * fPulseHeight * fADCpermV;
  MF_LOG_DEBUG("SignalShapingDUNEDPhase") << "Pulse height in mV 3m view : "<<fPulseHeight;
  MF_LOG_DEBUG("SignalShapingDUNEDPhase") << "Noise in ADC : "<<fAmpENCADC;
  // rebin to appropriate sampling rate of readout
  SetResponseSampling(clockData, fColSignalShaping );

  // 3x1x1 1m view
  fPulseHeight1m = SetElectResponse( fColShapeFunc1m, fColSignalShaping1m, fAreaNorm1m );
  MF_LOG_DEBUG("SignalShapingDUNEDPhase") << "Pulse height in mV 1m view : "<<fPulseHeight1m;
  // rebin to appropriate sampling rate of readout
  SetResponseSampling(clockData, fColSignalShaping1m );


  // Calculate field and electronics response functions.
  if(fAddFieldFunction){
    std::vector<double> fresp;
    SetFieldResponse( fresp );
    fColSignalShaping.AddResponseFunction(fresp);
    fColSignalShaping1m.AddResponseFunction(fresp);
  }

  //
  // Calculate filter functions.
  SetFilters(clockData);

  // Configure deconvolution kernels.
  fColSignalShaping.AddFilterFunction(fColFilter);
  fColSignalShaping.CalculateDeconvKernel();
  fColSignalShaping1m.AddFilterFunction(fColFilter);
  fColSignalShaping1m.CalculateDeconvKernel();
}

//-----Give Gain Settings to SimWire-----//
double util::SignalShapingServiceDUNEDPhase::GetASICGain(unsigned int const channel) const
{
  //
  if( !fHave1mView ) return fPulseHeight;

  //
  double gain = 0;
  art::ServiceHandle<geo::Geometry> geom;
  std::vector<geo::WireID> Wires =  geom->ChannelToWire(channel);
  double wirestartpoint[3];
  double wireendpoint[3];
  geom->WireEndPoints(Wires[0],wirestartpoint,wireendpoint);
  double wirelength = sqrt(pow(wirestartpoint[0]-wireendpoint[0],2) + pow(wirestartpoint[1]-wireendpoint[1],2) + pow(wirestartpoint[2]-wireendpoint[2],2));

  if((int)wirelength == 300)
    gain = fPulseHeight;
  else if((int)wirelength == 100)
    gain = fPulseHeight1m;
  else
    throw cet::exception("SignalShapingServiceDUNEDPhase")<< "can't determine"
							  << " View\n";
  return gain;
}

// area normalization factor should be in (ADC x us) / fC
double util::SignalShapingServiceDUNEDPhase::GetAreaNorm(unsigned int const channel) const
{
  //
  if( !fHave1mView ) return fAreaNorm;

  //
  double gain = 0;
  art::ServiceHandle<geo::Geometry> geom;
  std::vector<geo::WireID> Wires =  geom->ChannelToWire(channel);
  double wirestartpoint[3];
  double wireendpoint[3];
  geom->WireEndPoints(Wires[0],wirestartpoint,wireendpoint);
  double wirelength = sqrt(pow(wirestartpoint[0]-wireendpoint[0],2) + pow(wirestartpoint[1]-wireendpoint[1],2) + pow(wirestartpoint[2]-wireendpoint[2],2));

  if((int)wirelength == 300)
    gain = fAreaNorm;
  else if((int)wirelength == 100)
    gain = fAreaNorm1m;
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
void util::SignalShapingServiceDUNEDPhase::SetFieldResponse(std::vector<double> &fresp)
{
  // set the response for the collection plane first
  // the first entry is 0
  
  //double integral = 0;
  art::ServiceHandle<util::LArFFT> fft;
  MF_LOG_DEBUG("SignalShapingDUNEDPhase") << "Setting DUNEDPhase field response function...";
  
  int nticks = fft->FFTSize();
  std::vector<double> time(nticks,0.);
  fresp.resize(nticks, 0.);
  
  double integral =0;
  for(size_t i = 0; i < fresp.size(); ++i)
    {
      //convert time to microseconds, to match response function definition
      time[i]     = (1.*i)*fRespSamplingPeriod*1e-3;
      fresp[i] = FieldResponse(time[i]);
      integral +=fresp[i];
    }// end loop over time buckets
  
  //insert normalization
  for(size_t i = 0; i < fresp.size(); ++i)
    fresp[i] *= fColFieldRespAmp/integral;
  
  MF_LOG_DEBUG("SignalShapingDUNEDPhase") << " Done.";
  
  return;
}

//----------------------------------------------------------------------
// Calculate electronics response
double util::SignalShapingServiceDUNEDPhase::SetElectResponse( const TF1* fshape, 
							       util::SignalShaping &sig,
							       double areanorm )
{
  MF_LOG_DEBUG("SignalShapingDUNEDPhase") << "Setting DUNEDPhase electronics response function...";
  art::ServiceHandle<util::LArFFT> fft;
  int nticks = fft->FFTSize();
  std::vector<double> resp(nticks, 0.);
  
  // convert to us
  double dt  = fRespSamplingPeriod * 1.0E-3; 
  double max = 0.;
  double sum = 0.;
  for( auto i=0;i<nticks; ++i )
    {
      double t = i*dt;
      double v = fshape->Eval( t );
      if( t == 0 || v != v || v < 0 ) v = 0;
      
      if( v > max ) max = v;
      sum    += v;
      resp[i] = v;
    }
  
  //
  double norm    = areanorm / (sum * dt);
  double pheight = norm * max / fADCpermV;
  //normalize to 1e charge before the convolution
  for(auto& element : resp)
    {
      element *= norm * 1.60217657e-4; // per fC -> per electron
    }
  
  MF_LOG_DEBUG("SignalShapingDUNEDPhase") << "Area norm "<<areanorm<<" -> norm "<<norm <<"\n"
					  << "Pulse height "<<pheight;
  
  //
  sig.AddResponseFunction(resp);

  // Configure convolution kernels.
  sig.save_response();
  sig.set_normflag(false);

  MF_LOG_DEBUG("SignalShapingDUNEDPhase") << " Done.";

  return pheight;
}

//----------------------------------------------------------------------
// (Re)sample electronics (+field?) response
void util::SignalShapingServiceDUNEDPhase::SetResponseSampling(detinfo::DetectorClocksData const& clockData,
                                                               util::SignalShaping &sig )
{
  // Get services
  art::ServiceHandle<util::LArFFT> fft;

  // Operation permitted only if output of rebinning has a larger bin size
  if( fRespSamplingPeriod > sampling_rate(clockData) )
    throw cet::exception(__FUNCTION__) << "\033[93m"
				       << "Invalid operation: cannot rebin to a more finely binned vector!"
				       << "\033[00m" << std::endl;
  
  
  if( fRespSamplingPeriod == sampling_rate(clockData) ) //nothing to do
    {
      return;
    }
  
  //
  // make a new response vector with different time sampling
    
  const std::vector<double>* pResp = &(sig.Response_save());
  int nticks_input = pResp->size();
  
  //old sampling time vector
  std::vector<double> InputTime(nticks_input, 0. );
  std::vector<double> InputResp(nticks_input, 0. );
  for ( int itime = 0; itime < nticks_input; itime++ )
    {
      InputResp[itime] = (*pResp)[itime];
      InputTime[itime] = (1.*itime) * fRespSamplingPeriod; // in ns
    }

  // build a spline for interpolation
  TSpline3 ispl("ispl", &InputTime[0], &InputResp[0], nticks_input);
  
  
  //
  // default number of time ticks
  int nticks = fft->FFTSize();
  
  // resampling vectors
  std::vector<double> SamplingResp( nticks , 0. );
  double dt      = sampling_rate(clockData);
  double maxtime = InputTime.back();
  //
  for ( int itime = 0; itime < nticks; itime++ )
    {
      double t = itime * dt; // in ns
      // don't go past max t value from old response
      if(t > maxtime) break;
      SamplingResp[itime] = ispl.Eval(t);
    }
  
  sig.AddResponseFunction( SamplingResp, true );
  return;
}
  

//----------------------------------------------------------------------
// Calculate filter functions.
void util::SignalShapingServiceDUNEDPhase::SetFilters(detinfo::DetectorClocksData const& clockData)
{
  // Get services.
  art::ServiceHandle<util::LArFFT> fft;
  
  double ts = sampling_rate(clockData); // should return in ns
  int nsize = fft->FFTSize();
  int nhalf = nsize / 2;
  fColFilterFunc->SetRange(0, double(nhalf));
  
  //
  // Calculate collection filter.
  
  // the frequency step is given by fMhz/nsamples
  double df = 1.0/(ts*1.0E-3 * nsize);
  fColFilter.resize(nhalf+1);
  for(int i=0; i<=nhalf; ++i)
    {
      //double freq = 400. * i / (ts * nhalf);      // Cycles / microsecond.
      double freq   = i * df; 
      double f      = fColFilterFunc->Eval(freq);
      fColFilter[i] = TComplex(f, 0.);
    }
  
  return;
}

//----------------------------------------------------------------------
// Field response function (Harcoded, but it shouldn't)
double util::SignalShapingServiceDUNEDPhase::FieldResponse(double tval_us)
{
  /*
  //Studies for the Slow component function
  //reaad parameters from configuration
  double fWidth  = fParArray[0]; //Fast component time width
  double fHeight = fParArray[1]; //Qfast/Qextract
  double fNorm   = fParArray[2]; //Qextract/Qall
  double fTau 	 = fParArray[3]; //Trapping time slow component
  double f1=0;
  double f2=0;
  if(tval_us >= 0 && tval_us < fWidth) {f1 = 1/fWidth;}
  
  //f1 = (1/fWidth)*exp(-tval_us/fWidth);
  //f1 = TMath::Gaus(0, fWidth, true);
  //f1 = 1;
  f2 = (1/fTau)*exp(-tval_us/fTau);
  return fNorm*(fHeight*f1+(1-fHeight)*f2); //normalized field response
  */
  return fColFieldFunc->Eval(tval_us);
}

//
double util::SignalShapingServiceDUNEDPhase::GetDeconNorm() const {
  return fDeconNorm;
}

//
unsigned int util::SignalShapingServiceDUNEDPhase::GetSignalSize() const {
  return SignalShaping(0).Response().size();
}

// n/a
std::vector<DoubleVec> util::SignalShapingServiceDUNEDPhase::GetNoiseFactVec() const {
  return std::vector<DoubleVec>();
}


// n/a
int util::SignalShapingServiceDUNEDPhase::
FieldResponseTOffset(detinfo::DetectorClocksData const& clockData,
                     unsigned int const channel) const {
  return 0;
}

void util::SignalShapingServiceDUNEDPhase::
Convolute(detinfo::DetectorClocksData const&,
          unsigned int channel, FloatVector& func) const {
  return Convolute<float>(channel, func);
}

void util::SignalShapingServiceDUNEDPhase::
Convolute(detinfo::DetectorClocksData const&,
          unsigned int channel, DoubleVector& func) const {
  return Convolute<double>(channel, func);
}

void util::SignalShapingServiceDUNEDPhase::
Deconvolute(detinfo::DetectorClocksData const&,
            unsigned int channel, FloatVector& func) const {
  return Deconvolute<float>(channel, func);
}

void util::SignalShapingServiceDUNEDPhase::
Deconvolute(detinfo::DetectorClocksData const&,
            unsigned int channel, DoubleVector& func) const {
  return Deconvolute<double>(channel, func);
}


namespace util {

  DEFINE_ART_SERVICE(SignalShapingServiceDUNEDPhase)

}
