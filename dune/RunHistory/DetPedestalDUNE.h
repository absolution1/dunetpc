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
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "lbne-raw-data/Services/ChannelMap/ChannelMapService.h"
#include "fhiclcpp/ParameterSet.h"
#include <unordered_map>
#include <map>

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

    void SetDefaults(geo::View_t v, float mean, float meanerr, float rms, float rmserr)
    {
      fDefaultMean[v] = mean; fDefaultMeanErr[v] = meanerr;
      fDefaultRms[v]  = rms ; fDefaultRmsErr[v]  = rmserr;
    }
    float DefaultMean(geo::View_t v) const; 
    float DefaultMeanErr(geo::View_t v) const; 
    float DefaultRms(geo::View_t v) const; 
    float DefaultRmsErr(geo::View_t v) const; 
    
    void SetCSVFileName(std::string fname) { fCSVFileName = fname;}
    std::string CSVFileName() const {return fCSVFileName; }

    uint64_t VldTimeUsed() const { return fVldTimeUsed; }
    
  private:
    void LoadFromCSV();
    
    bool fUseDB;
    bool fUseDefaults;
    bool fAbortIfNoPeds;
    std::map<geo::View_t,float> fDefaultMean;
    std::map<geo::View_t,float> fDefaultMeanErr;
    std::map<geo::View_t,float> fDefaultRms;
    std::map<geo::View_t,float> fDefaultRmsErr;
    uint64_t fVldTime;
    uint64_t fVldTimeUsed;
    std::string fCSVFileName;
    std::string fDetName;
    int fLogLevel;
    std::unordered_map<raw::ChannelID_t,float> fMeanMap;
    std::unordered_map<raw::ChannelID_t,float> fRmsMap;
    std::unordered_map<raw::ChannelID_t,float> fMeanErrMap;
    std::unordered_map<raw::ChannelID_t,float> fRmsErrMap;

    art::ServiceHandle<lbne::ChannelMapService> fChannelMap;
  };
  
} // end namespace dune

#endif
