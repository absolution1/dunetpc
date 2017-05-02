#ifndef DETPEDESTALDUNE_CXX
#define DETPEDESTALDUNE_CXX

#include "dune/RunHistory/DetPedestalDUNE.h"
#include <string>
#include <iostream>
#include "nutools/IFDatabase/Table.h"

using std::string;
using std::cout;
using std::endl;

namespace dune {

  DetPedestalDUNE::DetPedestalDUNE(std::string detName) {
    fDetName = detName;
    fVldTime = fVldTimeUsed = 0;   
    fMeanMap.clear();
    fRmsMap.clear();
    fMeanErrMap.clear();
    fRmsErrMap.clear();
    fUseDB = false;
    fAbortIfNoPeds = false;
  }

  //------------------------------------------------------------

  DetPedestalDUNE::DetPedestalDUNE(fhicl::ParameterSet const& pset) {
    fVldTime = fVldTimeUsed = 0;   
    fMeanMap.clear();
    fRmsMap.clear();
    fMeanErrMap.clear();
    fRmsErrMap.clear();
    fUseDB = false;
    fAbortIfNoPeds = false;
    fDetName = "";
    Configure(pset);
  }
  
  //------------------------------------------------------------

  bool DetPedestalDUNE::Configure(fhicl::ParameterSet const& p) {
    const string myname = "DetPedestalDUNE::Configure: ";
    fUseDB = p.get<bool>("UseDB",false);
    if (!fUseDB) {
      fVldTime = p.get<int>("Run",0);
    }
    fAbortIfNoPeds = p.get<bool>("AbortIfNoPeds",false);
    fCSVFileName = p.get<std::string>("CSVFile","");    
    fUseDefaults = p.get<bool>("UseDefaults",false);
    fDefaultMean[geo::kU] = p.get<float>("DefaultMeanU",0.);
    fDefaultMeanErr[geo::kU] = p.get<float>("DefaultMeanErrU",0.);
    fDefaultRms[geo::kU] = p.get<float>("DefaultRmsU",0.);
    fDefaultRmsErr[geo::kU] = p.get<float>("DefaultRmsErrU",0.);    
    fDefaultMean[geo::kV] = p.get<float>("DefaultMeanV",0.);
    fDefaultMeanErr[geo::kV] = p.get<float>("DefaultMeanErrV",0.);
    fDefaultRms[geo::kV] = p.get<float>("DefaultRmsV",0.);
    fDefaultRmsErr[geo::kV] = p.get<float>("DefaultRmsErrV",0.);    
    fDefaultMean[geo::kZ] = p.get<float>("DefaultMeanZ",0.);
    fDefaultMeanErr[geo::kZ] = p.get<float>("DefaultMeanErrZ",0.);
    fDefaultRms[geo::kZ] = p.get<float>("DefaultRmsZ",0.);
    fDefaultRmsErr[geo::kZ] = p.get<float>("DefaultRmsErrZ",0.);    
    fDefaultMean[geo::kUnknown] = p.get<float>("DefaultMean",0.);
    fDefaultMeanErr[geo::kUnknown] = p.get<float>("DefaultMeanErr",0.);
    fDefaultRms[geo::kUnknown] = p.get<float>("DefaultRms",0.);
    fDefaultRmsErr[geo::kUnknown] = p.get<float>("DefaultRmsErr",0.);    
    
    fLogLevel = p.get<int>("LogLevel", 1);    
    if ( fLogLevel > 0 ) {
      std::cout << myname << "            UseDB: " << fUseDB << std::endl;
      std::cout << myname << "              Run: " << fVldTime << std::endl;
      std::cout << myname << "          CSVFile: " << fCSVFileName << std::endl;
      std::cout << myname << "      UseDefaults: " << fUseDefaults << std::endl;
      std::cout << myname << "   DefaultMean[U]: " << fDefaultMean[geo::kU] << std::endl;
      std::cout << myname << "DefaultMeanErr[U]: " << fDefaultMeanErr[geo::kU] << std::endl;
      std::cout << myname << "    DefaultRms[U]: " << fDefaultRms[geo::kU] << std::endl;
      std::cout << myname << " DefaultRmsErr[U]: " << fDefaultRmsErr[geo::kU] << std::endl;
      std::cout << myname << "   DefaultMean[V]: " << fDefaultMean[geo::kV] << std::endl;
      std::cout << myname << "DefaultMeanErr[V]: " << fDefaultMeanErr[geo::kV] << std::endl;
      std::cout << myname << "    DefaultRms[V]: " << fDefaultRms[geo::kV] << std::endl;
      std::cout << myname << " DefaultRmsErr[V]: " << fDefaultRmsErr[geo::kV] << std::endl;
      std::cout << myname << "   DefaultMean[Z]: " << fDefaultMean[geo::kZ] << std::endl;
      std::cout << myname << "DefaultMeanErr[Z]: " << fDefaultMeanErr[geo::kZ] << std::endl;
      std::cout << myname << "    DefaultRms[Z]: " << fDefaultRms[geo::kZ] << std::endl;
      std::cout << myname << " DefaultRmsErr[Z]: " << fDefaultRmsErr[geo::kZ] << std::endl;
      std::cout << myname << "         LogLevel: " << fLogLevel << std::endl;
    }

    return true;
  }

  //------------------------------------------------------------

  bool DetPedestalDUNE::Update(uint64_t ts)
  {
    if (fUseDefaults) return true;

    if ( fLogLevel > 0 ) 
      std::cout << __PRETTY_FUNCTION__ << " Called with run " << ts << std::endl;

    if (!fUseDB && fCSVFileName.empty()) {
      std::cout << __PRETTY_FUNCTION__ << " Method of determining pedestals is undefined! Either set UseDB or CSVFileName in fhicl. ";
      if (fAbortIfNoPeds) {
	std::cout << "Aborting per request." << std::endl;
	abort();
      }
      return false;
    }

    if (ts == fVldTime)
      return true;
    else
      fVldTime = ts;
    
    fMeanMap.clear();
    fMeanErrMap.clear();
    fRmsMap.clear();
    fRmsErrMap.clear();
    
    std::string tableName = "pedestals";
    nutools::dbi::Table t;

    if (fDetName.empty()) {
      std::cerr << "Detector name is undefined.  Aborting." << std::endl;
      std::abort();
    }
    
    t.SetDetector(fDetName);
    t.SetTableName(tableName);
    t.SetTableType(nutools::dbi::kConditionsTable);
    t.SetDataTypeMask(nutools::dbi::kDataOnly);

    int meanIdx = t.AddCol("mean","float");
    int rmsIdx = t.AddCol("rms","float");
    int meanErrIdx = t.AddCol("meanerr","float");
    int rmsErrIdx = t.AddCol("rmserr","float");
    
    t.SetMinTSVld(ts);
    t.SetMaxTSVld(ts);

    t.SetVerbosity(100);
    if (fUseDB)
      t.Load();
    else 
      t.LoadFromCSV(fCSVFileName);
    
    if (t.NRow() == 0) {
      std::cout << "Number of pedestals from database/CSV file is 0.  This should never be the case!  ";
      if (fAbortIfNoPeds) {
	std::cout << "Aborting, per request." << std::endl;
	abort();
      }
      return false;
    }

    nutools::dbi::Row* row;
    float mean, rms, meanerr, rmserr;
    uint64_t chan, offlineChan;
    for (int i=0; i<t.NRow(); ++i) {
      mean = rms = meanerr = rmserr = 0.;
      row = t.GetRow(i);
      chan = row->Channel();                     // MW: as of 3/1/16, this is now an online channel number
      offlineChan = fChannelMap->Offline(chan);  // Need to map to an offline channel at this point
      //offlineChan = chan;
      row->Col(meanIdx).Get(mean);
      row->Col(meanErrIdx).Get(meanerr);
      row->Col(rmsIdx).Get(rms);
      row->Col(rmsErrIdx).Get(rmserr);
      fMeanMap[offlineChan] = mean;
      fMeanErrMap[offlineChan] = meanerr;
      fRmsMap[offlineChan] = rms;
      fRmsErrMap[offlineChan] = rmserr;
      if (i==0) { // print out
	fVldTimeUsed = row->VldTime();
	if (fLogLevel > 0) {
	  std::cout << __PRETTY_FUNCTION__ << ": using run " << row->VldTime()
		    << " for pedestals." << std::endl;
	}
      }
    }

    return true;
  }

  
  //------------------------------------------------------------

  float DetPedestalDUNE::DefaultMean(geo::View_t v) const {
    float retVal=-1.;
    auto it = fDefaultMean.find(v);
    if (it != fDefaultMean.end())
      return it->second;
    else
      return retVal;

  }
  
  //------------------------------------------------------------

  float DetPedestalDUNE::DefaultMeanErr(geo::View_t v) const {
    float retVal=-1.;
    auto it = fDefaultMeanErr.find(v);
    if (it != fDefaultMeanErr.end())
      return it->second;
    else
      return retVal;

  }
  
  //------------------------------------------------------------

  float DetPedestalDUNE::DefaultRms(geo::View_t v) const {
    float retVal=-1.;
    auto it = fDefaultRms.find(v);
    if (it != fDefaultRms.end())
      return it->second;
    else
      return retVal;

  }
  
  //------------------------------------------------------------

  float DetPedestalDUNE::DefaultRmsErr(geo::View_t v) const {
    float retVal=-1.;
    auto it = fDefaultRmsErr.find(v);
    if (it != fDefaultRmsErr.end())
      return it->second;
    else
      return retVal;

  }
  
  //------------------------------------------------------------

  float DetPedestalDUNE::PedMean(raw::ChannelID_t ch) const { 
    float retVal=-1.;
    if (fUseDefaults) {
      auto it = fDefaultMean.find(geo::kUnknown);
      if (it != fDefaultMean.end())
	return it->second;
      else
	return retVal;
    }
    if (fVldTime == 0) {
      std::cerr << "DetPedestalDUNE: Validity time is not set!  Aborting." << std::endl;
      abort();
    }
    auto it = fMeanMap.find(ch);
    if (it != fMeanMap.end())
      retVal = it->second;
    
    return retVal;
  }

  //------------------------------------------------------------

  float DetPedestalDUNE::PedMeanErr(raw::ChannelID_t ch) const { 
    float retVal=0.;
    if (fUseDefaults) {
      auto it = fDefaultMeanErr.find(geo::kUnknown);
      if (it != fDefaultMeanErr.end())
	return it->second;
      else
	return retVal;
    }

    if (fVldTime == 0) {
      std::cerr << "DetPedestalDUNE: Validity time is not set!  Aborting." << std::endl;
      abort();
    }
    auto it = fMeanErrMap.find(ch);
    if (it != fMeanErrMap.end())
      retVal = it->second;
    
    return retVal;
  }

    //------------------------------------------------------------

  float DetPedestalDUNE::PedRms(raw::ChannelID_t ch) const { 
    float retVal=-1;
    if (fUseDefaults) {
      auto it = fDefaultRms.find(geo::kUnknown);
      if (it != fDefaultRms.end())
	return it->second;
      else
	return retVal;
    }
    if (fVldTime == 0) {
      std::cerr << "DetPedestalDUNE: Validity time is not set!  Aborting." << std::endl;
      abort();
    }
    auto it = fRmsMap.find(ch);
    if (it != fRmsMap.end())
      retVal = it->second;
    
    return retVal;
  }

  //------------------------------------------------------------

  float DetPedestalDUNE::PedRmsErr(raw::ChannelID_t ch) const { 
    float retVal=0.;
    if (fUseDefaults) {
      auto it = fDefaultRmsErr.find(geo::kUnknown);
      if (it != fDefaultRmsErr.end())
	return it->second;
      else
	return retVal;
    }
    if (fVldTime == 0) {
      std::cerr << "DetPedestalDUNE: Validity time is not set!  Aborting." << std::endl;
      abort();
    }
    auto it = fRmsErrMap.find(ch);
    if (it != fRmsErrMap.end())
      retVal = it->second;
    
    return retVal;
  }

  //------------------------------------------------------------

  void DetPedestalDUNE::PrintAllValues() {
    auto itMean = fMeanMap.begin();
    for (; itMean != fMeanMap.end(); ++itMean) {
      std::cout << "Channel: " << itMean->first << ", Mean = " << itMean->second
		<< ", RMS = " << fRmsMap[itMean->first] << std::endl;
    }	
    
  }
  
} // end namespace dune

#endif
