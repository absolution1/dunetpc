// AdcMultiThreshSignalFinder_tool.cc
//
// This is a more basic version of the algorithm by 
// Christoph Alt for the 3x1x1 dual-phase TPC detector
// 
//  
//********************************
//
// A ROI is defined when these criteria are met:
// 1. a minimum number of consecutive bins above a relatively low threhsold 
// 2. inside temporary ROI: a minimum number of consecutive bins above a medium threshold
// 3. at least one sample above a maximum threshold
//
// Putting settings for condition 2. to those of 1. effectively removes condition 2.
//
// All thresholds are relative to 0. Threy can be expressed either
// in absolute ADC counts or as multiples of the standard deviation. 
// The standard deviation is obtained from pedestal RMS so the pedestal 
// fit tool should be run first in this case
//
// The ROI is build by expanding the window of found ROI until one goes
// below some minimum value specified via ThresholdMin. After this it is padded 
// with additional samples on either side of the widnow (BinsBefore, BinsAfter)
//
//********************************

#include "AdcMultiThreshSignalFinder.h"
#include <iostream>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include "lardata/Utilities/LArFFT.h"

using std::vector;
using std::string;
using std::ostream;
using std::cout;
using std::endl;
using std::setw;

//**********************************************************************
// Class methods.
//**********************************************************************

AdcMultiThreshSignalFinder::AdcMultiThreshSignalFinder(fhicl::ParameterSet const& ps): 
  m_LogLevel(ps.get<int>("LogLevel")), 
  m_UseStd(ps.get<bool>("UseStandardDeviation")),
  m_Thresh1(ps.get<AdcSignal>("Threshold1")),
  m_NsaAbove1(ps.get<AdcIndex>("SamplesAboveThreshold1")),
  m_Thresh2(ps.get<AdcSignal>("Threshold2")),
  m_NsaAbove2(ps.get<AdcIndex>("SamplesAboveThreshold2")),
  m_ThreshMax(ps.get<AdcSignal>("ThresholdMax")),
  m_ThreshMin(ps.get<AdcSignal>("ThresholdMin")),
  m_BinsBefore(ps.get<AdcIndex>("BinsBefore")),
  m_BinsAfter(ps.get<AdcIndex>("BinsAfter")) {
  const string myname = "AdcMultiThreshSignalFinder::ctor: ";
  if ( m_LogLevel > 0 ) {
    cout << myname << "          Configuration: " << endl;
    cout << myname << "               LogLevel: " << m_LogLevel   << endl;
    cout << myname << "   UseStandardDeviation: " << m_UseStd     << endl;
    cout << myname << "             Threshold1: " << m_Thresh1    << endl;
    cout << myname << " SamplesAboveThreshold1: " << m_NsaAbove1  << endl;
    cout << myname << "             Threshold2: " << m_Thresh2    << endl;
    cout << myname << " SamplesAboveThreshold2: " << m_NsaAbove2  << endl;
    cout << myname << "           ThresholdMax: " << m_ThreshMax  << endl;
    cout << myname << "           ThresholdMin: " << m_ThreshMin  << endl;
    cout << myname << "             BinsBefore: " << m_BinsBefore << endl;
    cout << myname << "              BinsAfter: " << m_BinsAfter  << endl;
  }
}

//**********************************************************************


AdcMultiThreshSignalFinder::~AdcMultiThreshSignalFinder() {
}

//**********************************************************************

DataMap AdcMultiThreshSignalFinder::update(AdcChannelData& data) const {
  const string myname = "AdcMultiThreshSignalFinder:build: ";
  if ( m_LogLevel >= 2 ) cout << myname << "Building ROIs for channel "
                              << data.channel() << "." << endl;
  data.rois.clear();
  
  //Create dummy DataMap to return
  DataMap res(0);
  res.setInt("Test", 0);

  //Raw ADC.
  AdcSignalVector sigs;
  sigs = data.samples;

  // Build ROIs before padding and merging.
  AdcFilterVector& signal = data.signal;
  AdcRoiVector& rois      = data.rois;
  signal.clear();
  signal.resize(sigs.size(), false);

  AdcSignal sigth1   = m_Thresh1;
  AdcSignal sigth2   = m_Thresh2;
  AdcSignal sigthmax = m_ThreshMax;
  AdcSignal sigthmin = m_ThreshMin;
  AdcIndex nsig      = sigs.size();

  if ( nsig < 1 ) {
    if ( m_LogLevel >= 2 ) cout << myname << "Channel " << data.channel()
                                << " has no samples." << endl;
    return res;
  }

  AdcSignal ped    = data.pedestal;
  AdcSignal pedrms = data.pedestalRms;

  if( m_LogLevel >= 2 )
    cout << myname << "Channel "<< data.channel()
	 << " pedestal "<<ped<<" "<<pedrms << endl;

  if(m_UseStd)
  {
    if( ped == AdcChannelData::badSignal() )
      {
	cout << myname << "  Channel "<< data.channel()
	     <<" pedestal is not valid" << endl;
	return res;

      }

    // sligthly increase if 0
    if( pedrms == 0 ) pedrms = 0.1;
    
    sigth1   = m_Thresh1 * pedrms;
    sigth2   = m_Thresh2 * pedrms;
    sigthmax = m_ThreshMax * pedrms;
    sigthmin = m_ThreshMin * pedrms;
  }

  if( m_LogLevel >= 3 )
    {
      cout << myname << "  sigth1:   " << sigth1   << endl;
      cout << myname << "  sigth2:   " << sigth2   << endl;
      cout << myname << "  sigthmax: " << sigthmax << endl;
      cout << myname << "  sigthmin: " << sigthmin << endl;
    }
  
  ROICandidate_t roiCandidate;
  roiCandidate.init();
  AdcIndex nsaTh2 = 0;
  AdcIndex nsaTh3 = 0;
  
  // loop over samples
  for( AdcIndex isig=0; isig<nsig; ++isig ) 
    {
      AdcSignal sig = sigs[isig];

      if( sig >= sigth1 )
	{
	  if( !roiCandidate.isRoi )
	    {
	      roiCandidate.isRoi    = true;
	      roiCandidate.StartRoi = isig;
	    }
	  
	  if( sig > roiCandidate.MaxValue )
	    roiCandidate.MaxValue = sig;

	  roiCandidate.NsaTh1++;
	  roiCandidate.EndRoi = isig;
	  
	  //
	  if( sig >= sigth2 )
	    {
	      nsaTh2++;
	      if( nsaTh2 > roiCandidate.NsaTh2 )
		roiCandidate.NsaTh2 = nsaTh2;
	    }
	  else
	    {
	      if( nsaTh2 > roiCandidate.NsaTh2 )
		roiCandidate.NsaTh2 = nsaTh2;
	      nsaTh2 = 0;
	    }
	  
	  //
	  if( sig >= sigthmax )
	    {
	      nsaTh3++;
	      if( nsaTh3 > roiCandidate.NsaTh3 )
		roiCandidate.NsaTh3 = nsaTh3;
	    }
	  else
	    {
	      if( nsaTh3 > roiCandidate.NsaTh3 )
		roiCandidate.NsaTh3 = nsaTh3;
	      nsaTh3 = 0;
	    }
	  
	  continue; // go to a next ADC sample
	}
      
      // else ... check and finalize an ROI candidate
      if( roiCandidate.isRoi ) 
	{
	  // check that we meet all the cut requirements
	  bool ok = (( roiCandidate.NsaTh1 >= m_NsaAbove1 ) && 
		     ( roiCandidate.NsaTh2 >= m_NsaAbove2 ) && 
		     ( roiCandidate.NsaTh3 >= 1 ));
	  
	  
	  if( ok ) // finalize ROI
	    {
	      if( m_LogLevel >= 3 )
		{
		  cout<< myname <<" Candidate: "
		      <<roiCandidate.NsaTh1<<" "
		      <<roiCandidate.NsaTh2<<" "
		      <<roiCandidate.NsaTh3<<" "
		      <<roiCandidate.MaxValue<<endl;
		}
	  
	      // expand roi until we pass the min theshold
	      // use only signed int to move backwards to avoid infinit loop when at 0
	      int istart = (int)roiCandidate.StartRoi;
	      for( int ii=istart; ii>=0; --ii )
		{
		  if( sigs[ii] <= m_ThreshMin ) break;
		  if( roiCandidate.StartRoi>0 ) roiCandidate.StartRoi--;
		} 
	      for( AdcIndex ii = roiCandidate.EndRoi; ii < nsig; ++ii)
		{
		  if( sigs[ii] <= m_ThreshMin ) break;
		  roiCandidate.EndRoi++;
		}

	      
	      // pad ROI
	      roiCandidate.StartRoi = roiCandidate.StartRoi > m_BinsBefore ? roiCandidate.StartRoi-m_BinsBefore : 0;
	      
	      roiCandidate.EndRoi += m_BinsAfter;
	      if( roiCandidate.EndRoi >= nsig )
		roiCandidate.EndRoi = nsig - 1;
	      
	      // flag the ROIs
	      for(AdcIndex isigroi = roiCandidate.StartRoi; 
		  isigroi <= roiCandidate.EndRoi; ++isigroi)
		{
		  signal[isigroi] = true;
		}
	    }
	}
      
      // clear everything
      roiCandidate.init(); 
      nsaTh2 = 0; nsaTh3 = 0;
    } // end of sample loop
      
  
  // Fill from ROIs.
  data.roisFromSignal();
  
  // Display final ROIs.
  if ( m_LogLevel >= 3 ) 
    {
      cout << myname << "  ROIs (size = " << rois.size() << "):" << endl;
      for ( const AdcRoi& roi : rois ) {
	cout << myname << setw(8) << roi.first << " " << setw(8) << roi.second << endl;
      }
    } 
  else if ( m_LogLevel >= 2 ) 
    {
      cout << myname << "  ROI count: " << data.rois.size() << endl;
    }
  
  //
  return res;
}

DEFINE_ART_CLASS_TOOL(AdcMultiThreshSignalFinder)
