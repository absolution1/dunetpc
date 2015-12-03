/**
 *
 *
 * \brief Pedestal provider class for DUNE
 *
 * @author jpaley@fnal.gov
 */

#ifndef DETPEDESTALDUNE_H
#define DETPEDESTALDUNE_H

#include "CalibrationDBI/IOVData/DetPedestal.h"
#include "CalibrationDBI/Interface/IDetPedestalProvider.h"
#include "fhiclcpp/ParameterSet.h"
#include <unordered_map>

namespace dune {

  class DetPedestalDUNE : public lariov::IDetPedestalProvider
  {
  public:

    DetPedestalDUNE(int detid);
    DetPedestalDUNE(fhicl::ParameterSet const& pset);
    
    virtual float PedMean(raw::ChannelID_t ch) const;    
    virtual float PedRms(raw::ChannelID_t ch) const;
    virtual float PedMeanErr(raw::ChannelID_t ch) const;
    virtual float PedRmsErr(raw::ChannelID_t ch) const;
    
    bool Configure(fhicl::ParameterSet const& pset);
    bool Update(uint64_t ts);

    void SetUseDB(bool f) { fUseDB = f;}
    void PrintAllValues();
    
  private:
    void LoadFromCSV();
    
    bool fUseDB;
    bool fUseDefaults;
    int  fDetId;
    float fDefaultMean;
    float fDefaultMeanErr;
    float fDefaultRms;
    float fDefaultRmsErr;
    uint64_t fVldTime;   
    std::string fCSVFileName;
    std::unordered_map<raw::ChannelID_t,float> fMeanMap;
    std::unordered_map<raw::ChannelID_t,float> fRmsMap;
    std::unordered_map<raw::ChannelID_t,float> fMeanErrMap;
    std::unordered_map<raw::ChannelID_t,float> fRmsErrMap;
  };
  
} // end namespace dune

#endif
