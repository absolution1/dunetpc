////////////////////////////////////////////////////////////////////////
//
// DisambigAlg35t.cxx
//
// trj@fnal.gov
// tjyang@fnal.gov
//
// description
//
// Based on Tom Junk's idea of triplets matching
//
//
//
////////////////////////////////////////////////////////////////////////


//Framework includes:
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/AssociationUtil.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Wire.h"
#include "Geometry/CryostatGeo.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/WireGeo.h"
#include "DisambigAlg35t.h"

#include <map>
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>

namespace lbne{

DisambigAlg35t::DisambigAlg35t(fhicl::ParameterSet const& pset) 
{
  this->reconfigure(pset); 
}

//----------------------------------------------------------
void DisambigAlg35t::reconfigure(fhicl::ParameterSet const& p)
{

}


//----------------------------------------------------------
//----------------------------------------------------------
void DisambigAlg35t::RunDisambig( const std::vector< art::Ptr<recob::Hit> > &OrigHits   )
{
  fDisambigHits.clear();

  art::ServiceHandle<util::DetectorProperties> detprop;
  art::ServiceHandle<geo::Geometry> geo;

  std::vector<art::Ptr<recob::Hit> > hitsU;
  std::vector<art::Ptr<recob::Hit> > hitsV;
  std::vector<art::Ptr<recob::Hit> > hitsZ;

  for (size_t i = 0; i<OrigHits.size(); ++i){
    switch (OrigHits[i]->View()){
    case geo::kU:
      hitsU.push_back(OrigHits[i]);
      break;
    case geo::kV:
      hitsV.push_back(OrigHits[i]);
      break;
    case geo::kZ:
      hitsZ.push_back(OrigHits[i]);
      break;
    default:
      throw cet::exception("DisambigAlg35t")
	<<": hit view unkonw. \n";
    }
  }
  std::map<unsigned int, unsigned int> fHasBeenDisambigedU;
  std::map<unsigned int, unsigned int> fHasBeenDisambigedV;

  //Look for triplets of U,V,Z hits that are common in time
  //Try all possible wire segments for U and V hits and 
  //see if U, V, Z wires cross
  for (size_t z = 0; z<hitsZ.size(); ++z){//loop over z hits
    //find triplets of U,V,Z hits that are common in time
    unsigned int apaz(0), cryoz(0);
    fAPAGeo.ChannelToAPA(hitsZ[z]->Channel(), apaz, cryoz);
    double tz = hitsZ[z]->PeakTime()
      - detprop->GetXTicksOffset(hitsZ[z]->WireID().Plane,
				 hitsZ[z]->WireID().TPC,
				 hitsZ[z]->WireID().Cryostat);
    double mintvz = 1e10;
    unsigned int tmatchv = 0;
    for (size_t v = 0; v<hitsV.size(); ++v){
      if (fHasBeenDisambigedV.find(v)!=fHasBeenDisambigedV.end()) continue;
      unsigned int apav(0), cryov(0);
      fAPAGeo.ChannelToAPA(hitsV[v]->Channel(), apav, cryov);
      if (apav!=apaz) continue;
      double tv = hitsV[v]->PeakTime()
	- detprop->GetXTicksOffset(1,
				   hitsZ[z]->WireID().TPC,
				   hitsZ[z]->WireID().Cryostat);
      if (std::abs(tv-tz)<mintvz){
	mintvz = std::abs(tv-tz);
	tmatchv = v;
      }
    }
    if (mintvz>0.5) continue;
    double mintuz = 1e10;
    unsigned int tmatchu = 0;
    for (size_t u = 0; u<hitsU.size(); ++u){
      if (fHasBeenDisambigedU.find(u)!=fHasBeenDisambigedU.end()) continue;
      unsigned int apau(0), cryou(0);
      fAPAGeo.ChannelToAPA(hitsU[u]->Channel(), apau, cryou);
      if (apau!=apaz) continue;
      double tu = hitsU[u]->PeakTime()
	- detprop->GetXTicksOffset(0,
				   hitsZ[z]->WireID().TPC,
				   hitsZ[z]->WireID().Cryostat);
      if (std::abs(tu-tz)<mintuz){
	mintuz = std::abs(tu-tz);
	tmatchu = u;
      }
    }
    if (mintuz>0.5) continue;
    //found out which 3 wire segments cross
    geo::WireID zwire = geo->ChannelToWire(hitsZ[z]->Channel())[0];
    std::vector<geo::WireID> uwires = geo->ChannelToWire(hitsU[tmatchu]->Channel());
    std::vector<geo::WireID> vwires = geo->ChannelToWire(hitsV[tmatchv]->Channel());

    unsigned int totalintersections = 0;
    unsigned int bestu = 0;
    unsigned int bestv = 0;
    for (size_t uw = 0; uw<uwires.size(); ++uw){
      for (size_t vw = 0; vw<vwires.size(); ++vw){
	geo::WireIDIntersection widiuv;
	geo::WireIDIntersection widiuz;
	geo::WireIDIntersection widivz;
	if (!geo->WireIDsIntersect(zwire,uwires[uw],widiuz)) continue;
	if (!geo->WireIDsIntersect(zwire,vwires[vw],widivz)) continue;
	if (!geo->WireIDsIntersect(uwires[uw],vwires[vw],widiuv)) continue;
	double dis1 = sqrt(pow(widiuz.y-widivz.y,2)+pow(widiuz.z-widivz.z,2));
	double dis2 = sqrt(pow(widiuz.y-widiuv.y,2)+pow(widiuz.z-widiuv.z,2));
	double dis3 = sqrt(pow(widiuv.y-widivz.y,2)+pow(widiuv.z-widivz.z,2));
	double maxdis = std::max(dis1,dis2);
	maxdis = std::max(maxdis,dis3);
	if (maxdis<1){
	  ++totalintersections;
	  bestu = uw;
	  bestv = vw;
	}
      }
    }
    if (totalintersections==1){
      fDisambigHits.push_back(std::pair<art::Ptr<recob::Hit>, geo::WireID>(hitsU[tmatchu],uwires[bestu]));
      fDisambigHits.push_back(std::pair<art::Ptr<recob::Hit>, geo::WireID>(hitsV[tmatchv],vwires[bestv]));
      fHasBeenDisambigedU[tmatchu] = 1+bestu;

      fHasBeenDisambigedV[tmatchv] = 1+bestv;
      //std::cout<<"*** "<<hitsU[tmatchu]->Wire()<<" "<<hitsV[tmatchv]<<" "<<fDisambigHits.size()<<std::endl;
    }
  }
  //Done finding trivial disambiguated hits
  //loop over undisambiguated hits, find the nearest channel of disambiguated hits and determine the correct wire segment.
  for (size_t u1 = 0; u1<hitsU.size(); ++u1){
    if (fHasBeenDisambigedU.find(u1)!=fHasBeenDisambigedU.end()) continue;
    unsigned int apau1(0), cryou1(0);
    fAPAGeo.ChannelToAPA(hitsU[u1]->Channel(), apau1, cryou1);
    unsigned int channdiff = 100000;
    geo::WireID nearestwire;
    for (auto& u2 : fHasBeenDisambigedU){
      unsigned int apau2(0), cryou2(0);
      fAPAGeo.ChannelToAPA(hitsU[u2.first]->Channel(), apau2, cryou2);
      if (apau1!=apau2) continue;
      if (std::max(hitsU[u2.first]->Channel(),hitsU[u1]->Channel())
	  -std::min(hitsU[u2.first]->Channel(),hitsU[u1]->Channel())<channdiff){
	channdiff = std::max(hitsU[u2.first]->Channel(),hitsU[u1]->Channel())
	  -std::min(hitsU[u2.first]->Channel(),hitsU[u1]->Channel());
	nearestwire = geo->ChannelToWire(hitsU[u2.first]->Channel())[u2.second-1];
      }
    }
    if (nearestwire.isValid){
      unsigned wirediff = 1000000;
      unsigned bestwire = 0;
      std::vector<geo::WireID> wires = geo->ChannelToWire(hitsU[u1]->Channel());	
      for (size_t w = 0; w<wires.size(); ++w){
	if (wires[w].TPC!=nearestwire.TPC) continue;
	if (std::max(nearestwire.Wire,wires[w].Wire)-
	    std::min(nearestwire.Wire,wires[w].Wire)<wirediff){
	  wirediff = std::max(nearestwire.Wire,wires[w].Wire)-
	    std::min(nearestwire.Wire,wires[w].Wire);
	  bestwire = w;
	}
      }
      if (wirediff<1000000){
	fDisambigHits.push_back(std::pair<art::Ptr<recob::Hit>, geo::WireID>(hitsU[u1],wires[bestwire]));
	fHasBeenDisambigedU[u1] = 1+bestwire;
      }
    }
  }

  for (size_t v1 = 0; v1<hitsV.size(); ++v1){
    if (fHasBeenDisambigedV.find(v1)!=fHasBeenDisambigedV.end()) continue;
    unsigned int apav1(0), cryov1(0);
    fAPAGeo.ChannelToAPA(hitsV[v1]->Channel(), apav1, cryov1);
    unsigned int channdiff = 100000;
    geo::WireID nearestwire;
    for (auto& v2 : fHasBeenDisambigedV){
      unsigned int apav2(0), cryov2(0);
      fAPAGeo.ChannelToAPA(hitsV[v2.first]->Channel(), apav2, cryov2);
      if (apav1!=apav2) continue;
      if (std::max(hitsV[v2.first]->Channel(),hitsV[v1]->Channel())
	  -std::min(hitsV[v2.first]->Channel(),hitsV[v1]->Channel())<channdiff){
	channdiff = std::max(hitsV[v2.first]->Channel(),hitsV[v1]->Channel())
	  -std::min(hitsV[v2.first]->Channel(),hitsV[v1]->Channel());
	nearestwire = geo->ChannelToWire(hitsV[v2.first]->Channel())[v2.second-1];
      }
    }
    if (nearestwire.isValid){
      unsigned wirediff = 1000000;
      unsigned bestwire = 0;
      std::vector<geo::WireID> wires = geo->ChannelToWire(hitsV[v1]->Channel());	
      for (size_t w = 0; w<wires.size(); ++w){
	if (wires[w].TPC!=nearestwire.TPC) continue;
	if (std::max(nearestwire.Wire,wires[w].Wire)-
	    std::min(nearestwire.Wire,wires[w].Wire)<wirediff){
	  wirediff = std::max(nearestwire.Wire,wires[w].Wire)-
	    std::min(nearestwire.Wire,wires[w].Wire);
	}
      }
      if (wirediff<1000000){
	fDisambigHits.push_back(std::pair<art::Ptr<recob::Hit>, geo::WireID>(hitsV[v1],wires[bestwire]));
	fHasBeenDisambigedV[v1] = 1+bestwire;
      }
    }
  }

}

} //end namespace apa
