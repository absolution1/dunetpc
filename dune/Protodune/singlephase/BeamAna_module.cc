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
//#include "dune/BeamData/ProtoDUNEBeamSpill/ProtoDUNEBeamSpill.h"
#include "dune/DuneObj/ProtoDUNEBeamSpill.h"
#include <bitset>
#include <iomanip>

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
  void reconfigure(fhicl::ParameterSet const & p);
  std::bitset<sizeof(double)*CHAR_BIT> toBinary(const long num);
  void toBits(double);

private:

  // Declare member data here.
  //bool fLoadFromDB; // unused
  double  fTimeWindow;
  std::string fCSVFileName;
  std::string fBundleName;
  std::string fURLStr;
  uint64_t fFixedTime;
  std::vector<uint64_t> fMultipleTimes;
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

  std::cout << "TESTING BINARY" << std::endl;
  std::cout << toBinary(long(5.)) << std::endl;

  std::unique_ptr<ifbeam_ns::BeamFolder> bfp = ifb->getBeamFolder(fBundleName,fURLStr,fTimeWindow);
  std::cerr << "%%%%%%%%%% Got beam folder %%%%%%%%%%" << std::endl;

  //Use multiple times provided to fcl
  if( fMultipleTimes.size() ){
    std::cout << "Using multiple times: " << std::endl;
    for(size_t it = 0; it < fMultipleTimes.size(); ++it){
      std::cout << fMultipleTimes[it] << std::endl;
    }
  }
  //Use singular time provided to fcl
  else if(fFixedTime){
    std::cout << "Using singular time: " << fFixedTime << std::endl;
    fMultipleTimes.push_back(fFixedTime);
  }
  //Use time of event
  else{
    std::cout <<" Using Event Time: " << uint64_t(e.time().timeLow()) << std::endl;
    fMultipleTimes.push_back(uint64_t(e.time().timeLow()));
  }
  std::string prefix = "dip/acc/NORTH/NP04/BI/XBPF/";
  std::string devices[6] = {"XBPF022697_V",
                            "XBPF022697_H",
                            "XBPF022707_V",
                            "XBPF022707_H",
                            "XBPF022716_V",
                            "XBPF022716_H"};
  size_t nDev = 6;
  beamspill::ProtoDUNEBeamSpill * spill = new beamspill::ProtoDUNEBeamSpill();
  spill->InitFBMs(nDev);


  for(size_t it = 0; it < fMultipleTimes.size(); ++it){
    std::cout << "Time: " << fMultipleTimes[it] << std::endl;
    for(int d = 0; d < 6; ++d){
      std::cout <<"Device: " << devices[d] << std::endl;
      std::vector<double> data = bfp->GetNamedVector(fMultipleTimes[it], prefix + devices[d] + ":eventsData[]");
      std::vector<double> counts = bfp->GetNamedVector(fMultipleTimes[it], prefix + devices[d] + ":countsRecords[]");
    
  
      std::cout << "Data: " << data.size() << std::endl;
      std::cout << "Counts: " << counts.size() << std::endl;
      for(size_t i = 0; i < counts.size(); ++i){
        std::cout << counts[i] << std::endl;
      }
      if(counts[1] > data.size())continue;

      beamspill::FBM fbm;
      fbm.ID = d;

      for(size_t i = 0; i < counts[1]; ++i){
        std::cout << "Count: " << i << std::endl;
        for(int j = 0; j < 10; ++j){
          double theData = data[20*i + (2*j + 1)];
          std::cout << std::setw(12) << theData ;
          if(j < 4){
            fbm.timeData[j] = theData;           
          }
          else{
            fbm.fiberData[j - 4] = theData;
          }
        }
        spill->AddFBMTrigger(d,fbm);
        std::cout << std::endl;
      }

      for(size_t i = 0; i < spill->GetNFBMTriggers(d); ++i){
        spill->DecodeFibers(d,i);
      }

      for(size_t i = 0; i < counts[1]; ++i){
        std::cout << "Count: " << i << std::endl;
        for(int j = 0; j < 10; ++j){
          double theData = data[20*i + (2*j + 1)];
          toBinary(long(theData));        

        }
        std::cout <<std::endl;
      }
    }
  }

  for(size_t id = 0; id < 6; ++id){
    std::cout << "FBM: " << id << std::endl;
    std::cout << "N Triggers: " << spill->GetNFBMTriggers(id) << std::endl;

  }


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
  fMultipleTimes = p.get< std::vector<uint64_t> >("MultipleTimes");
}

std::bitset<sizeof(long)*CHAR_BIT> proto::BeamAna::toBinary(long num){
/*  std::cout << "double Size: " << sizeof(double)*CHAR_BIT << std::endl;
  std::cout << "long size: " << sizeof(unsigned long)*CHAR_BIT << std::endl;*/

//  /*unsigned*/ long * p = reinterpret_cast</*unsigned*/ long*>(&num);
/*  std::cout << "P! " << p << std::endl;
  std::cout << "bits: " << std::bitset<32>(p[0]) << " " << std::bitset<32>(p[1]) << std::endl;
  std::cout << "MyBits: " << std::bitset<64>(*p) << std::endl;*/
//  std::bitset<sizeof(double)*CHAR_BIT> mybits(*p);
    
  std::bitset<sizeof(double)*CHAR_BIT> mybits(num);  
  std::bitset<32> upper, lower;
  for(int i = 0; i < 32; ++i){
    lower[i] = mybits[i];
    upper[i] = mybits[i + 32];   
  }
  if(upper.any()) std::cout << "WARNING: NONZERO HALF" << std::endl;
  /*char * p = reinterpret_cast<char*>(&num);
  for (int i = sizeof(double)*CHAR_BIT-1 ; i >= 0 ; --i){
    mybits.set(i, (*(p)&(1<<i) ));
  }*/
  return mybits;
}

void proto::BeamAna::toBits(double input){

  std::cout << "size: " << sizeof(double) << std::endl;
  unsigned char rawBytes[sizeof(double)];  
  memcpy(rawBytes,&input,sizeof(double));

  unsigned char startMask=1;
  while (0!=static_cast<unsigned char>(startMask<<1)) {
    startMask<<=1;
  }

  bool hasLeadBit=false;

  size_t byteIndex;
  for (byteIndex=0; byteIndex<sizeof(double); ++byteIndex){
    unsigned char bitMask=startMask;
    while (0!=bitMask) {
      if (0!=(bitMask&rawBytes[byteIndex])) {
        std::cout<<"1";
        hasLeadBit=true;
      } 
      else if (hasLeadBit) {
        std::cout<<"0";
      }
      bitMask>>=1;
    }
  }
  if(!hasLeadBit){
    std::cout << "0";
  }
}

DEFINE_ART_MODULE(proto::BeamAna)
