////////////////////////////////////////////////////////////////////////
// \file SlowControlsProtoDUNE.cxx
//
// \brief implementation of class for slow controls data for ProtoDUNE
//
// \author jpaley@fnal.gov
// 
////////////////////////////////////////////////////////////////////////

// C++ language includes
#include <iostream>
#include <fstream>
#include "math.h"
#include "stdio.h"

// LArSoft includes
#include "dune/SlowControls/SlowControlsProtoDUNE.h"

// nutools includes
#include "nutools/IFDatabase/Table.h"

// Framework includes
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//-----------------------------------------------
slowctrls::SlowControlsProtoDUNE::SlowControlsProtoDUNE()
{
  fUseCondb = true;
  fIsLoaded = false;
  fTimeWindow = 3600;
  fCurrentTS = 0;
  fVerbosity = 0;
  fNameToIdMap.clear();
  fIdToNameMap.clear();
  fValues.clear();
  fTimes.clear();
}

//-----------------------------------------------
slowctrls::SlowControlsProtoDUNE::SlowControlsProtoDUNE(
  fhicl::ParameterSet const& pset)
{
  fUseCondb = true;
  fIsLoaded = false;
  fTimeWindow = 3600;
  fCurrentTS = 0;
  fVerbosity = 0;
  fNameToIdMap.clear();
  fIdToNameMap.clear();
  fValues.clear();
  fTimes.clear();

  Configure(pset);
}

//------------------------------------------------
bool slowctrls::SlowControlsProtoDUNE::Configure(fhicl::ParameterSet const& pset)
{  
  fUseCondb = pset.get<bool>("UseCondb");
  fSlowCtrlFileName = pset.get<std::string>("SlowCtrlFileName");
  fTimeWindow = pset.get<int>("TimeWindow");
  fVerbosity = pset.get<int>("Verbosity");
  fNameToIdMap = pset.get< std::map<std::string,int> >("ChannelMap");
  std::map<std::string,int>::iterator itr = fNameToIdMap.begin();
  for (; itr != fNameToIdMap.end(); ++itr) {
    fIdToNameMap[itr->second] = itr->first;
  }
  
  this->SetTimeWindow(fTimeWindow);
  return true;
}

//------------------------------------------------
// convert time window from minutes to seconds, 
// in five-minute steps:
//------------------------------------------------

void slowctrls::SlowControlsProtoDUNE::SetTimeWindow(int dt)
{
  fTimeWindow = int((float(dt) + 2.5)/5.)*5;
  fTimeWindow *= 60;
}

//------------------------------------------------
bool slowctrls::SlowControlsProtoDUNE::Update(float ts) 
{
  if (fCurrentTS == 0 || !fIsLoaded) {
    int ti = ts;
    ti -= 150.;
    ti /= fTimeWindow;
    ti *= fTimeWindow;
    fIsLoaded = false;
    fValues.clear();
    fTimes.clear();
    fCurrentTS = float(ti);
    return true;
  }
  
  if (ts > fCurrentTS && ts < fCurrentTS+fTimeWindow)
    return true;

  // else, outside of current time window, so get new window of data
  fCurrentTS = 0;
  return this->Update(ts);
}

//------------------------------------------------
double slowctrls::SlowControlsProtoDUNE::GetValue(std::string& chanId,
						  float ts) 
{
  if (fNameToIdMap.empty()) {
    mf::LogError("SlowControlsProtoDUNE") << "No slow controls channels defined!";
    return 0.;
  }

  if (fNameToIdMap.find(chanId) == fNameToIdMap.end()) {
    mf::LogError("SlowControlsProtoDUNE") << "Unknown slow controls channel!";
    return 0.;
  }

  this->Update(ts);

  if (! this->Load()) {
    mf::LogError("SlowControlsProtoDUNE") << "Failed to load table!";
    return 0;
  }

  std::vector<float>& tVec = fTimes[fNameToIdMap[chanId]];
  std::vector<double>& values = fValues[fNameToIdMap[chanId]];

  if (values.empty()) {
    mf::LogWarning("SlowControlsProtoDUNE") << "No data found for channel: "
					    << chanId 
					    << ". Returning a value of 0." 
					    << std::endl;
    return 0.;
  }
  
  unsigned int i = 0;
  
  for (; i< tVec.size() && ts >= tVec[i]; ++i);
  if (i==0) {
    mf::LogWarning("SlowControlsProtoDUNE") << "No data found for channel "
					    << chanId << " at time "
					    << ts << ".  Returning 0."
					    << std::endl;
    return 0.;
  }
  
  return values[i-1];

}

//------------------------------------------------
bool slowctrls::SlowControlsProtoDUNE::Load()
{
  if (fIsLoaded) return true;

  nutools::dbi::Table SlowCtrlTable;
  
  SlowCtrlTable.SetDetector("pdunesp");
  SlowCtrlTable.SetTableName("slowctrls");
  SlowCtrlTable.SetTableType(nutools::dbi::kConditionsTable);
  SlowCtrlTable.SetDataTypeMask(nutools::dbi::kDataOnly);
  SlowCtrlTable.SetMinTSVld(fCurrentTS);
  SlowCtrlTable.SetMaxTSVld(fCurrentTS+fTimeWindow);
  SlowCtrlTable.SetVerbosity(fVerbosity);

  int valueIdx = SlowCtrlTable.AddCol("value","double");
  
  bool readOk = false;
  if (!fSlowCtrlFileName.empty()) 
    readOk = SlowCtrlTable.LoadFromCSV(fSlowCtrlFileName);
  else
    readOk = SlowCtrlTable.Load();

  if (! readOk) {
    mf::LogError("SlowControlsProtoDUNE") << "Load from slow controls database table failed."; 
   
    return false; //std::abort();
  }
  
  if (SlowCtrlTable.NRow() == 0) {
    mf::LogError("SlowControlsProtoDUNE") << "Number of rows in slow controls table is 0.  This should never be the case!";
    return false;
  }
  
  nutools::dbi::Row* row;
  int chan;
  double value;
  float ts;
  for (int i=0; i<SlowCtrlTable.NRow(); ++i) {
    row = SlowCtrlTable.GetRow(i);      
    chan = row->Channel();
    ts = row->VldTime();
    row->Col(valueIdx).Get(value);
    fValues[chan].push_back(value);
    fTimes[chan].push_back(ts);
  }    

  fIsLoaded = true;
  return true;
}


