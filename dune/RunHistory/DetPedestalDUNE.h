/**
 *
 *
 * \brief Pedestal provider class for DUNE
 *
 * @author jpaley@fnal.gov
 */

#ifndef DETPEDESTALDUNE_H
#define DETPEDESTALDUNE_H

#include "larevt/CalibrationDBI/IOVData/DetPedestal.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalProvider.h"
#include "fhiclcpp/ParameterSet.h"
#include <unordered_map>

namespace dune {

  class DetPedestalDUNE : public lariov::DetPedestalProvider
  {
  public:

    DetPedestalDUNE(std::string detName="");
    DetPedestalDUNE(fhicl::ParameterSet const& pset);
    
    virtual float PedMean(raw::ChannelID_t ch) const;    
    virtual float PedRms(raw::ChannelID_t ch) const;
    virtual float PedMeanErr(raw::ChannelID_t ch) const;
    virtual float PedRmsErr(raw::ChannelID_t ch) const;
    
    bool Configure(fhicl::ParameterSet const& pset);
    bool Update(uint64_t ts);

    void SetUseDefaults(bool f) { fUseDefaults = f; }
    void SetUseDB(bool f) { fUseDB = f;}
    void PrintAllValues();

    void SetDetName(std::string detName) { fDetName = detName;}
    std::string DetName() { return fDetName; }

    void SetDefaults(float mean, float meanerr, float rms, float rmserr) {
      fDefaultMean = mean; fDefaultMeanErr = meanerr;
      fDefaultRms  = rms ; fDefaultRmsErr  = rmserr;
    }

    void SetCSVFileName(std::string fname) { fCSVFileName = fname;}
    std::string CSVFileName() const {return fCSVFileName; }
    
  private:
    void LoadFromCSV();
    
    bool fUseDB;
    bool fUseDefaults;
    float fDefaultMean;
    float fDefaultMeanErr;
    float fDefaultRms;
    float fDefaultRmsErr;
    uint64_t fVldTime;   
    std::string fCSVFileName;
    std::string fDetName;
    std::unordered_map<raw::ChannelID_t,float> fMeanMap;
    std::unordered_map<raw::ChannelID_t,float> fRmsMap;
    std::unordered_map<raw::ChannelID_t,float> fMeanErrMap;
    std::unordered_map<raw::ChannelID_t,float> fRmsErrMap;
  };
  
} // end namespace dune

#endif
