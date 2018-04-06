////////////////////////////////////////////////////////////////////////
// Class:       BeamAna
// Plugin Type: analyzer (art v2_08_03)
// File:        BeamAna_module.cc
//
// Generated at Thu Nov  2 22:57:41 2017 by Jonathan Paley using cetskelgen
// from cetlib version v3_01_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "IFBeam_service.h"

namespace proto {
  class BeamAna;
}


class proto::BeamAna : public art::EDAnalyzer {
public:
  explicit BeamAna(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  BeamAna(BeamAna const &) = delete;
  BeamAna(BeamAna &&) = delete;
  BeamAna & operator = (BeamAna const &) = delete;
  BeamAna & operator = (BeamAna &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;
  void beginRun(art::Run const & r) override;
  void beginSubRun(art::SubRun const & sr) override;
  void endJob() override;
  void endRun(art::Run const & r) override;
  void endSubRun(art::SubRun const & sr) override;
  void reconfigure(fhicl::ParameterSet const & p) override;

private:

  // Declare member data here.
  bool fLoadFromDB;
  double  fTimeWindow;
  std::string fCSVFileName;
  std::string fBundleName;
  std::string fURLStr;
  uint64_t fFixedTime;
};


proto::BeamAna::BeamAna(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{
  this->reconfigure(p);
}

void proto::BeamAna::analyze(art::Event const & e)
{

  std::cerr << "%%%%%%%%%% Getting ifbeam service handle %%%%%%%%%%" << std::endl;
  art::ServiceHandle<ifbeam_ns::IFBeam> ifb;
  std::cerr << "%%%%%%%%%% Got ifbeam handle %%%%%%%%%%" << std::endl;

  std::unique_ptr<ifbeam_ns::BeamFolder> bfp = ifb->getBeamFolder(fBundleName,fURLStr,fTimeWindow);
  std::cerr << "%%%%%%%%%% Got beam folder %%%%%%%%%%" << std::endl;

  uint64_t ts = fFixedTime;
  if (!fFixedTime)
    ts = uint64_t(e.time().timeLow());

  std::vector<double> x1 = bfp->GetNamedVector(ts,"E:NP04BPROF1X[]");
  std::vector<double> y1 = bfp->GetNamedVector(ts,"E:NP04BPROF1Y[]");
  
  std::vector<double> x2 = bfp->GetNamedVector(ts,"E:NP04BPROF2X[]");  
  std::vector<double> y2 = bfp->GetNamedVector(ts,"E:NP04BPROF2Y[]");  

  std::vector<double> x3 = bfp->GetNamedVector(ts,"E:NP04BPROF3X[]");
  std::vector<double> y3 = bfp->GetNamedVector(ts,"E:NP04BPROF3Y[]");
  
  std::vector<double> x4 = bfp->GetNamedVector(ts,"E:NP04BPROF4X[]");
  std::vector<double> y4 = bfp->GetNamedVector(ts,"E:NP04BPROF4Y[]");

  //  std::vector<double> x5 = bfp->GetNamedVector(ts,"E:NP04BPROF5X[]");
  //  std::vector<double> y5 = bfp->GetNamedVector(ts,"E:NP04BPROF5Y[]");
  
  std::vector<double> ckov1 = bfp->GetNamedVector(ts,"E:NP04CKOV1[]");
  std::vector<double> ckov2 = bfp->GetNamedVector(ts,"E:NP04CKOV2[]");

  std::vector<double> timeList = bfp->GetTimeList();
  std::vector<std::string> deviceList = bfp->GetDeviceList();
  /*  
  for (size_t i=0; i<timeList.size(); ++i)
    std::cout << "time[" << i << "] = " << timeList[i] << std::endl;
  for (size_t i=0; i<deviceList.size(); ++i)
    std::cout << "device[" << i << "] = " << deviceList[i] << std::endl;
  */

  std::cout << "x1.size() = " << x1.size() 
	    << ", y1.size() = " << y1.size() << std::endl;
  std::cout << "x2.size() = " << x2.size() 
	    << ", y2.size() = " << y2.size() << std::endl;
  std::cout << "x3.size() = " << x3.size() 
	    << ", y3.size() = " << y3.size() << std::endl;
  std::cout << "x4.size() = " << x4.size() 
	    << ", y4.size() = " << y4.size() << std::endl;
}

void proto::BeamAna::beginJob()
{
  // Implementation of optional member function here.
}

void proto::BeamAna::beginRun(art::Run const & r)
{
  // Implementation of optional member function here.
}

void proto::BeamAna::beginSubRun(art::SubRun const & sr)
{
  // Implementation of optional member function here.
}

void proto::BeamAna::endJob()
{
  // Implementation of optional member function here.
}

void proto::BeamAna::endRun(art::Run const & r)
{
  // Implementation of optional member function here.
}

void proto::BeamAna::endSubRun(art::SubRun const & sr)
{
  // Implementation of optional member function here.
}

void proto::BeamAna::reconfigure(fhicl::ParameterSet const & p)
{
  // Implementation of optional member function here.
  fBundleName  = p.get<std::string>("BundleName");
  fURLStr      = "";
  fTimeWindow  = p.get<double>("TimeWindow");
  fFixedTime   = p.get<uint64_t>("FixedTime");
}

DEFINE_ART_MODULE(proto::BeamAna)
