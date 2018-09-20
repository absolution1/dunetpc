////////////////////////////////////////////////////////////////////////
/// \file   SignalShapingServiceDUNEDPhase_service.cc
/// \author vgalymov@ipnl.in2p3.fr
////////////////////////////////////////////////////////////////////////

#include "dune/Utilities/SignalShapingServiceDUNEDPhase.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
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
  fColSignalShaping1Meter.Reset();
  fColSignalShaping3Meter.Reset();

  // Fetch fcl parameters.
  fDeconNorm            = pset.get<double>("DeconNorm");
  fASICmVperfC1Meter    = pset.get<double>("ASICmVperfC1Meter");
  fASICmVperfC3Meter    = pset.get<double>("ASICmVperfC3Meter");
  fADCpermV             = pset.get<double>("ADCpermV");
  fAmpENC               = pset.get<double>("AmpENC");
  fRespSamplingPeriod   = pset.get<double>("RespSamplingPeriod");

  fCollection1MeterShape = pset.get<std::string>("Collection1MeterShape");
  fCollection3MeterShape = pset.get<std::string>("Collection3MeterShape");

  fCollection1MeterParams = pset.get<std::vector<double>>("Collection1MeterParams");
  fCollection3MeterParams = pset.get<std::vector<double>>("Collection3MeterParams");

  fCollection1MeterFunction = new TF1("fCollection1MeterFunction", fCollection1MeterShape.c_str());
  for(unsigned int i=0; i<fCollection1MeterParams.size(); ++i)
  {
    fCollection1MeterFunction->SetParameter(i, fCollection1MeterParams[i]);
  }
  fCollection1MeterFunction->SetRange(0,fRespSamplingPeriod*1e-3);

  fCollection3MeterFunction = new TF1("fCollection3MeterFunction", fCollection3MeterShape.c_str());
  for(unsigned int i=0; i<fCollection3MeterParams.size(); ++i)
  {
    fCollection3MeterFunction->SetParameter(i, fCollection3MeterParams[i]);
  }
  fCollection3MeterFunction->SetRange(0,fRespSamplingPeriod*1e-3);

  // amplifier noise in ADC
  fAmpENCADC            = fAmpENC * 1.60217657e-4 * fASICmVperfC1Meter * fADCpermV;

  // Construct parameterized collection filter function.

  //LOG_DEBUG("SignalShapingServiceDUNEDPhase") <<"ASIC Gain in mV per fC = "<<fASICmVperfC1Meter;
  mf::LogInfo("SignalShapingServiceDUNEDPhase") <<"ASIC Gain in mV per fC for 1 meter long readout channeks = "<<fASICmVperfC1Meter
						<<"ASIC Gain in mV per fC for 3 meter long readout channeks = "<<fASICmVperfC3Meter
						<<";  ADC conversion = "<<fADCpermV
						<<";  Amplifier ENC = "<<fAmpENC
						<<";  Amplifier ENC in ADC = "<<fAmpENCADC;


	mf::LogInfo("SignalShapingServiceDUNE") << "Getting Field response from .fcl file"; //<<< to be changed if we keep the harcoded function
	std::string colField = pset.get<std::string>("ColFieldFunction");
	std::vector<double> colFieldParams =
  pset.get<std::vector<double> >("ColFieldFunctionParams");
	fColFieldFunc = new TF1("colFieldFunction", colField.c_str());
	for(unsigned int i=0; i<colFieldParams.size(); ++i)
		fColFieldFunc->SetParameter(i, colFieldParams[i]);

  mf::LogInfo("SignalShapingServiceDUNE") << "Getting Filter from .fcl file";
	fAddFieldFunction = pset.get<bool>("AddFieldFunction");
  std::string colFilt = pset.get<std::string>("ColFilter");
  std::vector<double> colFiltParams =
  	pset.get<std::vector<double> >("ColFilterParams");
  fColFilterFunc = new TF1("colFilter", colFilt.c_str());
  for(unsigned int i=0; i<colFiltParams.size(); ++i)
    fColFilterFunc->SetParameter(i, colFiltParams[i]);
	fColFieldRespAmp = pset.get<double>("ColFieldRespAmp");

	//fParArray = pset.get< std::vector<double> >("FieldFunctionParameters");
}

//----------------------------------------------------------------------
// Accessor for single-plane signal shaper.
const util::SignalShaping&
util::SignalShapingServiceDUNEDPhase::SignalShaping(unsigned int channel) const
{
  if(!fInit)
    init();
  auto const* geom = lar::providerFrom<geo::Geometry>();
  std::vector<geo::WireID> Wires =  geom->ChannelToWire(channel);

  double wirestartpoint[3];
  double wireendpoint[3];
  geom->WireEndPoints(Wires[0],wirestartpoint,wireendpoint);
  double wirelength = sqrt(pow(wirestartpoint[0]-wireendpoint[0],2) + pow(wirestartpoint[1]-wireendpoint[1],2) + pow(wirestartpoint[2]-wireendpoint[2],2));

  if((int)wirelength == 300)
    return fColSignalShaping3Meter;
  else if((int)wirelength == 100)
    return fColSignalShaping1Meter;
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
  if(!fInit) {
    fInit = true;

    // Calculate field and electronics response functions.
		if(fAddFieldFunction){
			std::vector<double> fresp;
			SetFieldResponse(fresp);
			fColSignalShaping1Meter.AddResponseFunction(fresp);
			fColSignalShaping3Meter.AddResponseFunction(fresp);
		}

		std::vector<double> eresp1meter, eresp3meter;
    SetElectResponse(eresp1meter,eresp3meter); //Step 5 and 6
		fColSignalShaping1Meter.AddResponseFunction(eresp1meter);
		fColSignalShaping3Meter.AddResponseFunction(eresp3meter);

    // Configure convolution kernels.
    fColSignalShaping1Meter.save_response();
    fColSignalShaping1Meter.set_normflag(false);
    fColSignalShaping3Meter.save_response();
    fColSignalShaping3Meter.set_normflag(false);


    // rebin to appropriate sampling rate of readout
    // NOTE: could have done it from the start in the eresp calculation
    //       but implemented it like this for future flexibility
    SetResponseSampling(); //Step 7

    // Calculate filter functions.

    SetFilters(); //Step 8

    // Configure deconvolution kernels.
    fColSignalShaping1Meter.AddFilterFunction(fColFilter);
    fColSignalShaping1Meter.CalculateDeconvKernel();
    fColSignalShaping3Meter.AddFilterFunction(fColFilter);
    fColSignalShaping3Meter.CalculateDeconvKernel();

    // std::cout << " fColSignalShaping1Meter: " << fColSignalShaping1Meter
    //mf::LogInfo("SignalShapingServiceDUNEDPhase")<<"Done Init";
  }
}

//-----Give Gain Settings to SimWire-----//
double util::SignalShapingServiceDUNEDPhase::GetASICGain(unsigned int const channel) const
{
  art::ServiceHandle<geo::Geometry> geom;
  //geo::SigType_t sigtype = geom->SignalType(channel);

  double gain = 0;

   // we need to distinguis between the U and V planes
//  geo::View_t view = geom->View(channel);
//  if(view == geo::kU || view == geo::kV || view == geo::kZ )
//    gain = fASICmVperfC;


  std::vector<geo::WireID> Wires =  geom->ChannelToWire(channel);
  double wirestartpoint[3];
  double wireendpoint[3];
  geom->WireEndPoints(Wires[0],wirestartpoint,wireendpoint);
  double wirelength = sqrt(pow(wirestartpoint[0]-wireendpoint[0],2) + pow(wirestartpoint[1]-wireendpoint[1],2) + pow(wirestartpoint[2]-wireendpoint[2],2));

  if((int)wirelength == 300)
    gain = fASICmVperfC3Meter;
  else if((int)wirelength == 100)
    gain = fASICmVperfC1Meter;
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
	LOG_DEBUG("SignalShapingDUNEDPhase") << "Setting DUNEDPhase field response function...";

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

	LOG_DEBUG("SignalShapingDUNEDPhase") << " Done.";

   return;
}

//----------------------------------------------------------------------
// Calculate electronics response
void util::SignalShapingServiceDUNEDPhase::SetElectResponse(std::vector<double> &ElecResp1Meter, std::vector<double> &ElecResp3Meter)
{
  // Get services.

  art::ServiceHandle<util::LArFFT> fft;

  LOG_DEBUG("SignalShapingDUNEDPhase") << "Setting DUNEDPhase electronics response function...";

  int nticks = fft->FFTSize();
  std::vector<double> time(nticks,0.);
  ElecResp1Meter.resize(nticks, 0.);
  double max = 0.;

  for(size_t i = 0; i < ElecResp1Meter.size(); ++i)
    {
      time[i] = (1.*i)*fRespSamplingPeriod*1e-3; //convert time to microseconds, to match response function definition

      ElecResp1Meter[i] = fCollection1MeterFunction->Eval(time[i]); 

      if(time[i] == 0 || ElecResp1Meter[i] != ElecResp1Meter[i]) ElecResp1Meter[i] = 0.; //force function value to be 0 for t = 0 or if function value = nan

      if(ElecResp1Meter[i] > max) max = ElecResp1Meter[i];
    }// end loop over time buckets


  LOG_DEBUG("SignalShapingDUNEDPhase") << " Done.";

  //normalize to 1e charge before the convolution
  for(auto& element : ElecResp1Meter)
    {
      element /= max;
      element *= fASICmVperfC1Meter * 1.60217657e-4; //mV
      element *= fADCpermV;                    //ADC
    }
//  std::cout << "ElecResp1Meter after normalization" << std::endl;
//  std::cout << "tval_us" << "\t" << "fval" << std::endl;
//  for(size_t i = 0; i <= ElecResp1Meter.size(); ++i)
//  {
//    std::cout << time[i] << "\t" << ElecResp1Meter[i] << std::endl;
//  }



  ElecResp3Meter.resize(nticks, 0.);
  time.resize(nticks, 0.);
  max = 0.;

  for(size_t i = 0; i < ElecResp3Meter.size(); ++i)
    {
      time[i] = (1.*i)*fRespSamplingPeriod*1e-3; //convert time to microseconds, to match response function definition

      ElecResp3Meter[i] = fCollection3MeterFunction->Eval(time[i]);

      if(time[i] == 0 || ElecResp3Meter[i] != ElecResp3Meter[i]) ElecResp3Meter[i] = 0.; //force function value to be 0 for t=0 or if function value = nan

      if(ElecResp3Meter[i] > max) max = ElecResp3Meter[i];

    }// end loop over time buckets


  LOG_DEBUG("SignalShapingDUNEDPhase") << " Done.";

  //normalize to 1e charge before the convolution
  for(auto& element : ElecResp3Meter)
    {
      element /= max;
      element *= fASICmVperfC3Meter * 1.60217657e-4; //mV
      element *= fADCpermV;                    //ADC
    }
//  std::cout << "ElecResp3Meter after normalization" << std::endl;
//  std::cout << "tval_us" << "\t" << "fval" << std::endl;
//  for(size_t i = 0; i <= ElecResp3Meter.size(); ++i)
//  {
//    std::cout << time[i] << "\t" << ElecResp3Meter[i] << std::endl;
//  }

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

  for ( int iplane = 0; iplane < 2; iplane++ )
  {
    const std::vector<double>* pResp;
    switch ( iplane ) 
    {
    case 0: pResp = &(fColSignalShaping1Meter.Response_save()); break;
    case 1: pResp = &(fColSignalShaping3Meter.Response_save()); break;
    }
  // get response our old response vector
  //  const std::vector<double>* pResp = &(fColSignalShaping1Meter.Response_save());

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
    switch ( iplane )
    {
    case 0: fColSignalShaping1Meter.AddResponseFunction( SamplingResp, true ); break;
    case 1: fColSignalShaping3Meter.AddResponseFunction( SamplingResp, true ); break;
    }
  } // for ( int iplane = 0; iplane < 2; iplane++ )


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

namespace util {

  DEFINE_ART_SERVICE(SignalShapingServiceDUNEDPhase)

}
