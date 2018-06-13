////////////////////////////////////////////////////////////////////////
// RunHistoryDUNE.h
//
//  Data provider class for run history
//
// jpaley@fnal.gov
//
////////////////////////////////////////////////////////////////////////
#ifndef RUNHISTORYDUNE_H
#define RUNHISTORYDUNE_H

#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <memory>

#include "fhiclcpp/ParameterSet.h"

#include "lardataalg/DetectorInfo/RunHistory.h"
#include "nutools/IFDatabase/Table.h"

///General LArSoft Utilities
namespace dune {

  class ASICSetting {
  public:
    ASICSetting() {};
    ASICSetting(float g, float s, int b) {gain=g; shape=s; base=b;};
    ~ASICSetting() {};
    float gain;
    float shape;
    int base;
  }; 
  
  class SubRunDUNE : public detinfo::SubRun {
  public:
    SubRunDUNE() : fTStart(0) {};
    SubRunDUNE(SubRunDUNE const&) = delete;
    virtual ~SubRunDUNE() {};
    
    virtual uint64_t TStart() const override { return fTStart; }
    void SetTStart(uint64_t t) { fTStart = t; }
    
  private:
    uint64_t fTStart;
  };
    
  class RunHistoryDUNE : public detinfo::RunHistory {
  public: 
    typedef enum {
      kUnknownDet = 0,
      k35t,
      kProtoDUNE,
      kFarDet,
      kNearDet,
      kNDUNEDetectors
    } DetId_t;

    static std::string DetIdToString(int detId) {
      if (detId == k35t) return std::string("dune35t");
      else return std::string("UNKNOWN_DET");
    }

    static int DetNameToId(std::string detStr) {
      if (detStr == "dune35t" || detStr == "DUNE35t" ||
	  detStr == "35t") return k35t;
      else return kUnknownDet;
    }


  public:
    RunHistoryDUNE(int detid, int runnum);
    RunHistoryDUNE(RunHistoryDUNE const&) = delete;
    virtual ~RunHistoryDUNE();
      
    virtual bool Update(uint64_t ts=0) override;
       
    virtual int RunNumber() const override{ return fRun; }
    int DetId() const { return fDetId; }
    virtual int NSubruns() const override{ return fNSubruns; }
    virtual int RunType() const override{ return fRunType; }
    virtual std::string RunTypeAsString() const override;
    virtual uint64_t TStart() const override { return fTStart; }
    virtual uint64_t TStop()  const override { return fTStop; }
    virtual uint64_t Duration() const override { return fTStop-fTStart; }

    std::vector<std::string> Shifters() { return fShifter; }     
    std::vector<std::string> Components() {return fComponents; }
    
    void SetNSubruns(int nsr) { fNSubruns = nsr;}
    void SetRunType(int rt) { fRunType = rt; }
    void SetTStart(uint64_t t) { fTStart = t; }
    void SetTStop(uint64_t t) { fTStop = t; }
    void AddShifter(std::string sh) { fShifter.push_back(sh); }
    void SetShifters(std::vector<std::string> sh) { fShifter = sh; }

    std::string CfgLabel() const { return fCfgLabel; }
    std::string TStartAsString() const { return fTStartStr; }
    std::string TStopAsString() const { return fTStopStr; }
    
    void DumpSCData();
    void DumpASICSettings();
    
  private:

    bool LoadSCChanMap();
    bool LoadSCData();
    bool LoadASICSettings();
    
  protected:
    int    fRun;
    int    fNSubruns;
    int    fRunType;
    int    fDetId;
    
    uint64_t  fTStart;
    uint64_t  fTStop;      
    
    std::vector<std::string> fShifter;
    std::vector<std::string> fComponents;
    std::string fDetName;
    std::string fCfgLabel;
    std::string fTStartStr;
    std::string fTStopStr;
    
    std::vector<SubRunDUNE> fSubrun;
    std::unordered_map<std::string,int> fSCChanMap;
    std::unordered_map<int,std::string> fSCInvChanMap;
    std::unordered_map<int,ASICSetting> fASICSettingsMap;

    std::unique_ptr<nutools::dbi::Table> fSCDataTable;
    std::unique_ptr<nutools::dbi::Table> fASICSettingsTable;
    
  }; // class RunHistoryDUNE
} //namespace dune
#endif // RUNHISTORYDUNE_H
