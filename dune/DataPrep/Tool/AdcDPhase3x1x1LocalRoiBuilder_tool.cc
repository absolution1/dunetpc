// AdcDPhase3x1x1LocalRoiBuilder_tool.cc
// christoph.alt@cern.ch

//********************************
// This code is intended to be used on raw data. It determines ROI that should be excluded from pedestal and noise pattern calculation.
//
// Thode calculates a local pedestal in a sliding window and looks for ROI in the bins that follow this window.
// A ROI is defined when three criteria are met:
// 1. a minimum number of consecutive bins above a relatively low threhsold -> temporary ROI
// 2. inside temporary ROI: a minimum number of consecutive bins above a medium threshold
// 3. inside temporary ROI: at least one bin above a high threshold
// All thresholds are relative to the local pedestal of the sliding window. Threshold 1 and 2 can be expressed
// both in absolute ADC counts or as multiples of the standard deviation. 
// The standard deviation is only calculated for the first window while the pedestal is re-calculated for each window.
//
//********************************

#include "AdcDPhase3x1x1LocalRoiBuilder.h"
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
using art::ServiceHandle;

//**********************************************************************
// Class methods.
//**********************************************************************

AdcDPhase3x1x1LocalRoiBuilder::AdcDPhase3x1x1LocalRoiBuilder(fhicl::ParameterSet const& ps): 
  m_LogLevel(ps.get<int>("LogLevel")), 
  m_BinsToAverageForPedestal(ps.get<AdcIndex>("BinsToAverageForPedestal")),
  m_BinsToSkip(ps.get<AdcIndex>("BinsToSkip")),
  m_UseStandardDeviation(ps.get<bool>("UseStandardDeviation")),
  m_NConsecBinsAboveThreshold1(ps.get<AdcIndex>("NConsecBinsAboveThreshold1")),
  m_NSigmaStart1(ps.get<AdcSignal>("NSigmaStart1")),
  m_NSigmaEnd1(ps.get<AdcSignal>("NSigmaEnd1")),
  m_NConsecBinsAboveThreshold2(ps.get<AdcIndex>("NConsecBinsAboveThreshold2")),
  m_NSigmaStart2(ps.get<AdcSignal>("NSigmaStart2")),
  m_NSigmaMax(ps.get<AdcSignal>("NSigmaMax")),
  m_PadLow(ps.get<AdcIndex>("PadLow")),
  m_PadHigh(ps.get<AdcIndex>("PadHigh")) {
  const string myname = "AdcDPhase3x1x1LocalRoiBuilder::ctor: ";
  if ( m_LogLevel > 0 ) {
    cout << myname << "              Configuration: " << endl;
    cout << myname << "                   LogLevel: " << m_LogLevel << endl;
    cout << myname << "   BinsToAverageForPedestal: " << m_BinsToAverageForPedestal << endl;
    cout << myname << "                 BinsToSkip: " << m_BinsToSkip << endl;
    cout << myname << "       UseStandardDeviation: " << m_UseStandardDeviation << endl;
    cout << myname << " NConsecBinsAboveThreshold1: " << m_NConsecBinsAboveThreshold1 << endl;
    cout << myname << "               NSigmaStart1: " << m_NSigmaStart1 << endl;
    cout << myname << "                 NSigmaEnd1: " << m_NSigmaEnd1 << endl;
    cout << myname << " NConsecBinsAboveThreshold2: " << m_NConsecBinsAboveThreshold2 << endl;
    cout << myname << "               NSigmaStart2: " << m_NSigmaStart2 << endl;
    cout << myname << "                  NSigmaMax: " << m_NSigmaMax << endl;
    cout << myname << "                     PadLow: " << m_PadLow << endl;
    cout << myname << "                    PadHigh: " << m_PadHigh << endl;
  }
}

//**********************************************************************


AdcDPhase3x1x1LocalRoiBuilder::~AdcDPhase3x1x1LocalRoiBuilder() {
}

//**********************************************************************

DataMap AdcDPhase3x1x1LocalRoiBuilder::update(AdcChannelData& data) const {
  const string myname = "AdcDPhase3x1x1LocalRoiBuilder:build: ";
  if ( m_LogLevel >= 2 ) cout << myname << "Building ROIs for channel "
                              << data.channel << "." << endl;
  data.rois.clear();

  //Create dummy DataMap to return
  DataMap res(0);
  res.setInt("Test", 0);

  //Raw ADC.
  AdcSignalVector sigs;
  sigs = data.samples;

  // Build ROIS before padding and merging.
  AdcFilterVector& signal = data.signal;
  AdcRoiVector& rois = data.rois;
  signal.clear();
  signal.resize(sigs.size(), false);
  bool inroi = false;
  AdcSignal siglow1 = m_NSigmaEnd1;
  AdcSignal sighigh1 = m_NSigmaStart1;
  AdcSignal sighigh2 = m_NSigmaStart2;
  AdcSignal sigmax = m_NSigmaMax;
  AdcIndex nsig = sigs.size();

  int ROICount=0;
  std::vector<AdcIndex> ROIStart;
  std::vector<AdcIndex> ROIEnd;
  ROIStart.clear();
  ROIEnd.clear();

  if ( nsig < 1 ) {
    if ( m_LogLevel >= 2 ) cout << myname << "Channel " << data.channel
                                << " has no samples." << endl;
    return res;
  }

  //*****************************************
  //  FIRST ITERATION: BIN 0 TO LAST BIN  ***
  //*****************************************
  double SumPedestal=0.;
  double SumADCMinusPedestalSquared=0.;
  double StandardDeviationPedestal=0.;
  AdcIndex FirstEntryInPedestalSum=0;

  //Calculate pedestal for first m_BinsToAverageForPedestal ticks
  for ( AdcIndex isig=m_BinsToSkip; isig<m_BinsToAverageForPedestal+m_BinsToSkip; ++isig ) 
  {
    SumPedestal+=sigs[isig];
  }
  FirstEntryInPedestalSum=m_BinsToSkip;

  //Calculate standard deviation for first m_BinsToAverageForPedestal ticks
  if(m_UseStandardDeviation)
  {
    for ( AdcIndex isig=m_BinsToSkip; isig<m_BinsToAverageForPedestal+m_BinsToSkip; ++isig ) 
    {
      SumADCMinusPedestalSquared+=pow(sigs[isig]-SumPedestal/m_BinsToAverageForPedestal,2);
    }
  StandardDeviationPedestal = sqrt(SumADCMinusPedestalSquared/(m_BinsToAverageForPedestal-1));
  siglow1 = m_NSigmaEnd1*StandardDeviationPedestal;
  sighigh1 = m_NSigmaStart1*StandardDeviationPedestal;
  sighigh2 = m_NSigmaStart2*StandardDeviationPedestal;
  }

  for ( AdcIndex isig = m_BinsToAverageForPedestal+m_BinsToSkip; isig<sigs.size(); ++isig ) 
  {
    AdcSignal sig = sigs[isig];
    if ( inroi ) 
    {
      if ( sig > siglow1 + SumPedestal/m_BinsToAverageForPedestal &&  isig < sigs.size()-1)
      {
        signal[isig] = true;
      } 
      else  
      {
	ROIEnd.push_back(isig-1);
        ROICount++;

	//check second criterion in this temporary ROI.
	bool KeepThisROI = false;
	AdcIndex m_NConsecBinsAboveThreshold2Temp = std::min(m_NConsecBinsAboveThreshold2,(AdcIndex)(ROIEnd[ROICount-1]-ROIStart[ROICount-1]+1));
	AdcIndex NConsecBinsAboveThreshold2Count=0;

	for(AdcIndex isigroi = ROIStart[ROICount-1]; isigroi <= ROIEnd[ROICount-1]; isigroi++)
	{
	  if( sigs[isigroi] >= sighigh2 + SumPedestal/m_BinsToAverageForPedestal ) 
	  {
	    NConsecBinsAboveThreshold2Count++;
	    if(NConsecBinsAboveThreshold2Count == m_NConsecBinsAboveThreshold2Temp)
	    {
	      KeepThisROI = true;
	      break;
	    }
	  }
	  else 
	  {
	    NConsecBinsAboveThreshold2Count = 0;
	  }
	}

	//if second criterion for this temporary ROI is met, check third criterion
	if(KeepThisROI)
	{
	  KeepThisROI = false;
	  for(AdcIndex isigroi = ROIStart[ROICount-1]; isigroi <= ROIEnd[ROICount-1]; isigroi++)
	  {
	    if( sigs[isigroi] >= sigmax + SumPedestal/m_BinsToAverageForPedestal) 
	    {
	      KeepThisROI = true;
	      break;
	    }
	  }
	}

	//if second or third criteria is not met, delete temporary ROI. Otherwise keep it.
	if(!KeepThisROI)
	{
	  for(AdcIndex isigroi = ROIStart[ROICount-1]; isigroi <= ROIEnd[ROICount-1]; isigroi++)
	  {
	    signal[isigroi] = false;
	    //recalculate pedesta, including the ones in the deleted ROI
	    SumPedestal -= sigs[FirstEntryInPedestalSum]; //remove last entry in pedestal sum
	    FirstEntryInPedestalSum++; //increase index for last entry in pedestal sum by 1
	    while(signal[FirstEntryInPedestalSum]) //check if the increased index was a signal. if yes, increase until index with no signal was found.
	    {
	      FirstEntryInPedestalSum++;
	    }
	    SumPedestal += sigs[isigroi]; //add current ADC count to pedestal sum
	  }

	  ROIStart.pop_back();
	  ROIEnd.pop_back();
	  ROICount--;
	}

        inroi = false;

	SumPedestal -= sigs[FirstEntryInPedestalSum];//remove last entry in pedestal sum

	FirstEntryInPedestalSum++; //increase index for last entry in pedestal sum by 1
	while(signal[FirstEntryInPedestalSum]) //check if the increased index was a signal. if yes, increase until index with no signal was found.
	{
	  FirstEntryInPedestalSum++;
	}

	SumPedestal += sigs[isig]; //add current ADC count to sum
      }
    } 
    else 
    {
      bool ROIStartIsHere = true;
      for(AdcIndex isignext = isig; isignext < std::min((AdcIndex)sigs.size(),isig+m_NConsecBinsAboveThreshold1); isignext++)
      {
	if(sigs[isignext] < sighigh1 + SumPedestal/m_BinsToAverageForPedestal)
	{
	  ROIStartIsHere = false;
	  break;
	}
      }

      if(ROIStartIsHere) 
      {
        ROIStart.push_back(isig);
        signal[isig] = true; 
        inroi = true;
      }
      else
      {
	SumPedestal -= sigs[FirstEntryInPedestalSum]; //remove last entry in pedestal sum

	FirstEntryInPedestalSum++; //increase index for last entry in pedestal sum by 1
	while(signal[FirstEntryInPedestalSum]) //check if the increased index was a signal. if yes, increase until index with no signal was found.
	{
	  FirstEntryInPedestalSum++;
	}

	SumPedestal += sigs[isig]; //add current ADC count to sum
      }
    }
  }

//removal isolated signal at the end of the waveform
if(signal[sigs.size()-1] && !signal[sigs.size()-2])
{
signal[sigs.size()-1] = false;
}

  //******************************************
  //  SECOND ITERATION: LAST BIN TO BIN 0  ***
  //******************************************
  AdcIndex PedestalIndex=sigs.size()-1;
  AdcIndex PedestalCounter=0;
  inroi = false;
  siglow1 = m_NSigmaEnd1;
  sighigh1 = m_NSigmaStart1;
  sighigh2 = m_NSigmaStart2;

  SumPedestal=0.;
  SumADCMinusPedestalSquared=0.;
  StandardDeviationPedestal=0.;

  FirstEntryInPedestalSum=0;

  ROIStart.clear();
  ROIEnd.clear();
  ROICount=0;

  //Calculate pedestal for last m_BinsToAverageForPedestal ticks
  while(PedestalCounter < m_BinsToAverageForPedestal)
  {
    if(!signal[PedestalIndex])
    {
       if(PedestalCounter == 0) FirstEntryInPedestalSum=PedestalIndex; //remember first enntry of pedestal sum
       SumPedestal+=sigs[PedestalIndex];
       PedestalCounter++;
    }
    PedestalIndex--;
  }

  //Calculate standard deviation for last m_BinsToAverageForPedestal ticks
  if(m_UseStandardDeviation)
  {
    PedestalCounter=0;
    PedestalIndex=sigs.size()-1;
    while(PedestalCounter < m_BinsToAverageForPedestal)
    {
      if(!signal[PedestalIndex])
      {
        SumADCMinusPedestalSquared+=pow(sigs[PedestalIndex]-SumPedestal/m_BinsToAverageForPedestal,2);
        PedestalCounter++;
      }
    PedestalIndex--;
    }
  StandardDeviationPedestal = sqrt(SumADCMinusPedestalSquared/(m_BinsToAverageForPedestal-1));
  siglow1 = m_NSigmaEnd1*StandardDeviationPedestal;
  sighigh1 = m_NSigmaStart1*StandardDeviationPedestal;
  sighigh2 = m_NSigmaStart2*StandardDeviationPedestal;
  }


  for ( AdcIndex isig = PedestalIndex; isig >= m_BinsToSkip && isig <= PedestalIndex; --isig ) 
  {
    if(signal[isig]) continue;
    AdcSignal sig = sigs[isig];
    if ( inroi ) 
    {
      if ( sig > siglow1 + SumPedestal/m_BinsToAverageForPedestal && isig > m_BinsToSkip)
      {
        signal[isig] = true;
      } 
      else  
      {
	ROIEnd.push_back(isig+1);
        ROICount++;

	//check second criterion in this temporary ROI.
	bool KeepThisROI = false;
	AdcIndex m_NConsecBinsAboveThreshold2Temp = std::min(m_NConsecBinsAboveThreshold2,(AdcIndex)(ROIEnd[ROICount-1]-ROIStart[ROICount-1]+1));
	AdcIndex NConsecBinsAboveThreshold2Count=0;

	for(AdcIndex isigroi = ROIStart[ROICount-1]; isigroi >= ROIEnd[ROICount-1]; isigroi--)
	{
	  if( sigs[isigroi] >= sighigh2 + SumPedestal/m_BinsToAverageForPedestal ) 
	  {
	    NConsecBinsAboveThreshold2Count++;
	    if(NConsecBinsAboveThreshold2Count == m_NConsecBinsAboveThreshold2Temp)
	    {
	      KeepThisROI = true;
	      break;
	    }
	  }
	  else 
	  {
	    NConsecBinsAboveThreshold2Count = 0;
	  }
	}

	//if second criterion for this temporary ROI is met, check third criterion
	if(KeepThisROI)
	{
	  KeepThisROI = false;
	  for(AdcIndex isigroi = ROIStart[ROICount-1]; isigroi >= ROIEnd[ROICount-1]; isigroi--)
	  {
	    if( sigs[isigroi] >= sigmax + SumPedestal/m_BinsToAverageForPedestal ) 
	    {
	      KeepThisROI = true;
	      break;
	    }
	  }
	}

	//if second or third criteria is not met, delete temporary ROI. Otherwise keep it.
	if(!KeepThisROI)
	{
	  for(AdcIndex isigroi = ROIStart[ROICount-1]; isigroi >= ROIEnd[ROICount-1]; isigroi--)
	  {
	    signal[isigroi] = false;
	    //recalculate pedestal, including the ones in the deleted ROI
	    SumPedestal -= sigs[FirstEntryInPedestalSum]; //remove last entry in pedestal sum
	    FirstEntryInPedestalSum--; //decreas index for last entry in pedestal sum by 1
	    while(signal[FirstEntryInPedestalSum]) //check if the decreased index was a signal. if yes, increase until index with no signal was found.
	    {
	      FirstEntryInPedestalSum--;
	    }
	    SumPedestal += sigs[isigroi]; //add current ADC count to sum
	  }

	  ROIStart.pop_back();
	  ROIEnd.pop_back();
	  ROICount--;
	}

        inroi = false;
	SumPedestal -= sigs[FirstEntryInPedestalSum];//remove last entry in pedestal sum

	FirstEntryInPedestalSum--; //increase index for last entry in pedestal sum by 1
	while(signal[FirstEntryInPedestalSum]) //check if the increased index was a signal. if yes, increase until index with no signal was found.
	{
	  FirstEntryInPedestalSum--;
	}

	SumPedestal += sigs[isig]; //add current ADC count to sum

      }
    } 
    else 
    {
      bool ROIStartIsHere = true;
      for(AdcIndex isignext = isig; isignext >= (AdcIndex)std::max((int)m_BinsToSkip,(int)isig-(int)m_NConsecBinsAboveThreshold1+1) && isignext <= PedestalIndex ; isignext--)
      {
	if(sigs[isignext] < sighigh1 + SumPedestal/m_BinsToAverageForPedestal)
	{
	  ROIStartIsHere = false;
	  break;
	}
      }

      if(ROIStartIsHere) 
      {
        ROIStart.push_back(isig);
        signal[isig] = true; 
        inroi = true;
      }
      else
      {
	SumPedestal -= sigs[FirstEntryInPedestalSum]; //remove first entry in pedestal sum
	FirstEntryInPedestalSum--; //decrease index for last entry in pedestal sum by 1

	while(signal[FirstEntryInPedestalSum]) //check if the increased index was a signal. if yes, increase until index with no signal was found.
	{
	  FirstEntryInPedestalSum--;
	}
	SumPedestal += sigs[isig]; //add current ADC count to sum
      }
    }
  }

//removal isolated signal at the beginning of the waveform
if(signal[m_BinsToSkip] && !signal[m_BinsToSkip+1])
{
signal[m_BinsToSkip] = false;
}

  // Fill the unpadded ROIs.
  data.roisFromSignal();
  // Display ROIs before padding and merging.
  if ( m_LogLevel >= 3 ) {
    cout << myname << "  ROIs before merge (size = " << rois.size() << "):" << endl;
    for ( const AdcRoi& roi : rois ) {
      cout << myname << setw(8) << roi.first << " " << setw(8) << roi.second << endl;
    }
  } else if ( m_LogLevel >= 2 ) {
    cout << myname << "  ROI count before merge: " << data.rois.size() << endl;
  }
  if ( rois.size() == 0 ) return res;
  // Pad ROIs.
  unsigned int isig1 = 0;
  unsigned int isig2 = 0;
  for ( AdcRoi roi : rois ) {
    isig2 = roi.first;
    isig1 = isig2 > m_PadLow ? isig2 - m_PadLow : 0;
    for ( unsigned int isig=isig1; isig<isig2; ++isig ) signal[isig] = true;
    isig1 = roi.second + 1;
    isig2 = isig1 + m_PadHigh;
    if ( isig2 > nsig ) isig2 = nsig;
    for ( unsigned int isig=isig1; isig<isig2; ++isig ) signal[isig] = true;
  }
  // Fill the final ROIs.
  data.roisFromSignal();
  // Display final ROIs.
  if ( m_LogLevel >= 3 ) {
    cout << myname << "  ROIs after merge (size = " << rois.size() << "):" << endl;
    for ( const AdcRoi& roi : rois ) {
      cout << myname << setw(8) << roi.first << " " << setw(8) << roi.second << endl;
    }
  } else if ( m_LogLevel >= 2 ) {
    cout << myname << "  ROI count after merge: " << data.rois.size() << endl;
  }
  return res;
}
//**********************************************************************
