#ifndef DETPEDESTALDUNE_CXX
#define DETPEDESTALDUNE_CXX

#include "DetPedestalDUNE.h"
#include "IFDatabase/Table.h"

namespace dune {

  DetPedestalDUNE::DetPedestalDUNE(std::string detName)
  {
    fDetName = detName;
    fVldTime = 0;   
    fMeanMap.clear();
    fRmsMap.clear();
    fMeanErrMap.clear();
    fRmsErrMap.clear();
    fUseDB = false;
  }

  //------------------------------------------------------------

  DetPedestalDUNE::DetPedestalDUNE(fhicl::ParameterSet const& pset)
  {
    fVldTime = 0;   
    fMeanMap.clear();
    fRmsMap.clear();
    fMeanErrMap.clear();
    fRmsErrMap.clear();
    fUseDB = false;
    fDetName = "";
    Configure(pset);
  }
  
  //------------------------------------------------------------

  bool DetPedestalDUNE::Configure(fhicl::ParameterSet const& p)
  {
    fUseDB = p.get<bool>("UseDB",false);
    if (!fUseDB) {
      fVldTime = p.get<bool>("Run",0);
    }
    fCSVFileName = p.get<std::string>("CSVFile","");    
    fUseDefaults = p.get<bool>("UseDefaults",false);
    fDefaultMean = p.get<float>("DefaultMean",0.);
    fDefaultMeanErr = p.get<float>("DefaultMeanErr",0.);
    fDefaultRms = p.get<float>("DefaultRms",0.);
    fDefaultRmsErr = p.get<float>("DefaultRmsErr",0.);    

    return true;
  }

  //------------------------------------------------------------

  bool DetPedestalDUNE::Update(uint64_t ts)
  {
    if (!fUseDB) return true;

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
    if (fCSVFileName.empty())
      t.Load();
    else
      t.LoadFromCSV(fCSVFileName);

    if (t.NRow() == 0) {
      std::cerr << "Number of pedestals from database/CSV file is 0.  This should never be the case, aborting." << std::endl;
      abort();
    }

    nutools::dbi::Row* row;
    float mean, rms, meanerr, rmserr;
    uint64_t chan;
    for (int i=0; i<t.NRow(); ++i) {
      mean = rms = meanerr = rmserr = 0.;
      row = t.GetRow(i);      
      chan = row->Channel();
      row->Col(meanIdx).Get(mean);
      row->Col(meanErrIdx).Get(meanerr);
      row->Col(rmsIdx).Get(rms);
      row->Col(rmsErrIdx).Get(rmserr);
      fMeanMap[chan] = mean;
      fMeanErrMap[chan] = meanerr;
      fRmsMap[chan] = rms;
      fRmsErrMap[chan] = rmserr;      
    }
    
    return true;
  }

  
  //------------------------------------------------------------

  float DetPedestalDUNE::PedMean(raw::ChannelID_t ch) const { 
    if (fUseDefaults)
      return fDefaultMean;
    if (fVldTime == 0) {
      std::cerr << "DetPedestalDUNE: Validity time is not set!  Aborting." << std::endl;
      abort();
    }
    float retVal=-1.;
    auto it = fMeanMap.find(ch);
    if (it != fMeanMap.end())
      retVal = it->second;
    
    return retVal;
  }

  //------------------------------------------------------------

  float DetPedestalDUNE::PedMeanErr(raw::ChannelID_t ch) const { 
    if (fUseDefaults)
      return fDefaultMeanErr;
    if (fVldTime == 0) {
      std::cerr << "DetPedestalDUNE: Validity time is not set!  Aborting." << std::endl;
      abort();
    }
    float retVal=0.;
    auto it = fMeanErrMap.find(ch);
    if (it != fMeanErrMap.end())
      retVal = it->second;
    
    return retVal;
  }

    //------------------------------------------------------------

  float DetPedestalDUNE::PedRms(raw::ChannelID_t ch) const { 
    if (fUseDefaults)
      return fDefaultRms;
    if (fVldTime == 0) {
      std::cerr << "DetPedestalDUNE: Validity time is not set!  Aborting." << std::endl;
      abort();
    }
    float retVal=0.;
    auto it = fRmsMap.find(ch);
    if (it != fRmsMap.end())
      retVal = it->second;
    
    return retVal;
  }

  //------------------------------------------------------------

  float DetPedestalDUNE::PedRmsErr(raw::ChannelID_t ch) const { 
    if (fUseDefaults)
      return fDefaultRmsErr;
    if (fVldTime == 0) {
      std::cerr << "DetPedestalDUNE: Validity time is not set!  Aborting." << std::endl;
      abort();
    }
    float retVal=0.;
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
