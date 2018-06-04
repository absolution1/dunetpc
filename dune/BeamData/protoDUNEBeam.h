////////////////////////////////////////////////////////////////////////
//  protoDUNEBeam.h
//
//  Data provider class for beam data
//
// jpaley@fnal.gov
//
////////////////////////////////////////////////////////////////////////
#ifndef PROTODUNEBEAM_H
#define PROTODUNEBEAM_H

#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <memory>

#include "fhiclcpp/ParameterSet.h"
#include "ifbeam.h"

//#include "lardataalg/DetectorInfo/RunHistory.h"
//#include "nutools/IFDatabase/Table.h"

///General LArSoft Utilities
namespace proto {

  class ProtoDUNEBeam  {
  public: 
    ProtoDUNEBeam(fhicl::ParameterSet const& pset, ifbeam_ns::BeamFolder* ptr);
    ~ProtoDUNEBeam() {};

    bool Update(uint64_t ts = 0);

    void SetLoadFromDB(bool v) {_loadFromDB = v;}
    void SetCSVFileName(std::string s) {_csvFileName=s;}
    void SetTimeWindow(int tw) {_timeWindow = tw;}
    void SetBundleName(std::string bn) {_bundleName = bn; }
    bool LoadBeamFolder();
    void LoadData();
    
  private:
    
    bool        _loadFromDB;
    double      _timeWindow; // time window in seconds to search for data
    std::vector<double> _xy1; 
    std::vector<double> _xy2; 
    std::vector<double> _xy3;
    std::vector<double> _xy4; 
    std::vector<double> _ckov1; 
    std::vector<double> _ckov2; 
 
    uint64_t    _currentTS;
    std::string _csvFileName;
    std::string _bundleName;

    std::unique_ptr<ifbeam_ns::BeamFolder> _bfp;
    
  }; // class ProtoDUNEBeam
} //namespace dune
#endif // RUNHISTORYDUNE_H

