////////////////////////////////////////////////////////////////////////
//
//  protoDUNEBeam
//
//  jpaley@fnal.gov
//
////////////////////////////////////////////////////////////////////////
// Framework includes

// C++ language includes
#include <iostream>

// LArSoft includes
#include "protoDUNEBeam.h"

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

//#include "nutools/IFDatabase/Util.h"
//#include <boost/tokenizer.hpp>
//#include <ctime>

namespace proto {
  //-----------------------------------------------
  ProtoDUNEBeam::ProtoDUNEBeam(fhicl::ParameterSet const& pset, ifbeam_ns::BeamFolder* ptr) :
    _loadFromDB   (pset.get<bool        >("LoadFromDB")),
    _timeWindow   (pset.get<double      >("TimeWindow")),
    _csvFileName  (pset.get<std::string >("CSVFileName")),
    _bundleName   (pset.get<std::string >("BundleName")),
    _bfp(ptr)
  {

  }
  
  //------------------------------------------------
  bool ProtoDUNEBeam::Update(uint64_t ts) 
  {
    if (_currentTS == ts) return true;

    //    double trtgtd=0.;

    _xy1.clear();
    _xy2.clear();
    _xy3.clear();
    _xy4.clear();
    _ckov1.clear();
    _ckov2.clear();
    _xy1 = _bfp->GetNamedVector(ts,"E:NT04XY1[]");
    _xy2 = _bfp->GetNamedVector(ts,"E:NT04XY2[]");
    _xy3 = _bfp->GetNamedVector(ts,"E:NT04XY3[]");
    _xy4 = _bfp->GetNamedVector(ts,"E:NT04XY4[]");
    _ckov1 = _bfp->GetNamedVector(ts,"E:NT04Ckov1[]");
    _ckov2 = _bfp->GetNamedVector(ts,"E:NT04Ckov2[]");

    return true;
  }

}
