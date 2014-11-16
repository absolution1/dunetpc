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

  fTimeCut = p.get<double>("TimeCut");
  fDistanceCut = p.get<double>("DistanceCut");

}


//----------------------------------------------------------
//----------------------------------------------------------
void DisambigAlg35t::RunDisambig( const std::vector< art::Ptr<recob::Hit> > &OrigHits   )
{
  fDisambigHits.clear();

  art::ServiceHandle<util::DetectorProperties> detprop;
  art::ServiceHandle<geo::Geometry> geo;

  std::vector<std::vector<art::Ptr<recob::Hit> > > hitsUV(2);
  std::vector<art::Ptr<recob::Hit> > hitsZ;

  for (size_t i = 0; i<OrigHits.size(); ++i){
    switch (OrigHits[i]->View()){
    case geo::kU:
      hitsUV[0].push_back(OrigHits[i]);
      break;
    case geo::kV:
      hitsUV[1].push_back(OrigHits[i]);
      break;
    case geo::kZ:
      hitsZ.push_back(OrigHits[i]);
      break;
    default:
      throw cet::exception("DisambigAlg35t")
	<<": hit view unkonw. \n";
    }
  }
  std::vector<std::map<unsigned int, unsigned int> >fHasBeenDisambigedUV(2);

//  for (int i = 0; i<8; ++i){
//    std::cout<<detprop->GetXTicksOffset(0,i,0)<<" "
//	     <<detprop->GetXTicksOffset(1,i,0)<<" "
//	     <<detprop->GetXTicksOffset(2,i,0)<<std::endl;
//  }

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
    //if (hitsZ[z]->WireID().TPC!=0) continue;
    //if (geo->ChannelToWire(hitsZ[z]->Channel())[0].Wire!=60) continue;
    //if (z!=0) continue;
    //loop over u hits
    bool findmatch = false;
    for (size_t u = 0; u<hitsUV[0].size() && !findmatch; ++u){
      if (fHasBeenDisambigedUV[0].find(u)!=fHasBeenDisambigedUV[0].end()) continue;
      unsigned int apau(0), cryou(0);
      fAPAGeo.ChannelToAPA(hitsUV[0][u]->Channel(), apau, cryou);
      if (apau!=apaz) continue;
      double tu = hitsUV[0][u]->PeakTime()
	- detprop->GetXTicksOffset(0,
				   hitsZ[z]->WireID().TPC,
				   hitsZ[z]->WireID().Cryostat);
      //std::cout<<"u "<<tz<<" "<<tu<<" "<<geo->ChannelToWire(hitsUV[0][u]->Channel())[0].Wire<<" "<<geo->ChannelToWire(hitsUV[0][u]->Channel())[0].TPC<<std::endl;
      //std::cout<<z<<" "<<u<<" "<<std::abs(tu-tz)<<std::endl;
      if (std::abs(tu-tz)<fTimeCut){
	//if (hitsZ[z]->WireID().TPC == 0) std::cout<<z<<" "<<u<<" "<<std::abs(tu-tz)<<std::endl;
	//find a matched u hit, loop over v hits
	//std::cout<<hitsUV[1].size()<<std::endl;
	for (size_t v = 0; v<hitsUV[1].size(); ++v){
	  if (fHasBeenDisambigedUV[1].find(v)!=fHasBeenDisambigedUV[1].end()) continue;
	  unsigned int apav(0), cryov(0);
	  fAPAGeo.ChannelToAPA(hitsUV[1][v]->Channel(), apav, cryov);
	  if (apav!=apaz) continue;
	  double tv = hitsUV[1][v]->PeakTime()
	    - detprop->GetXTicksOffset(1,
				       hitsZ[z]->WireID().TPC,
				       hitsZ[z]->WireID().Cryostat);
	  //std::cout<<std::abs(tv-tz)<<std::endl;
	  //std::cout<<z<<" "<<u<<" "<<v<<" "<<std::abs(tu-tz)<<" "<<std::abs(tv-tz)<<std::endl;
	  //std::cout<<"v "<<tz<<" "<<tv<<" "<<geo->ChannelToWire(hitsUV[1][v]->Channel())[0].Wire<<" "<<geo->ChannelToWire(hitsUV[1][v]->Channel())[0].TPC<<std::endl;
	  if (std::abs(tv-tz)<fTimeCut){
	    //std::cout<<"triplets "<<z<<" "<<u<<" "<<v<<std::endl;
	    //find a matched v hit, see if the 3 wire segments cross
	    geo::WireID zwire = geo->ChannelToWire(hitsZ[z]->Channel())[0];
	    std::vector<geo::WireID>  uwires = geo->ChannelToWire(hitsUV[0][u]->Channel());
	    std::vector<geo::WireID>  vwires = geo->ChannelToWire(hitsUV[1][v]->Channel());
	    
	    unsigned int totalintersections = 0;
	    unsigned int bestu = 0;
	    unsigned int bestv = 0;
	    for (size_t uw = 0; uw<uwires.size(); ++uw){
	      for (size_t vw = 0; vw<vwires.size(); ++vw){
		geo::WireIDIntersection widiuv;
		geo::WireIDIntersection widiuz;
		geo::WireIDIntersection widivz;
		if (uwires[uw].TPC!=vwires[vw].TPC) continue;
		if (uwires[uw].TPC!=zwire.TPC) continue;
		if (vwires[vw].TPC!=zwire.TPC) continue;
		//std::cout<<"! "<<uwires[uw].Wire<<" "<<vwires[vw].Wire<<" "<<zwire.Wire<<" "<<tu<<" "<<tv<<" "<<tz<<" "<<uwires[uw].TPC<<" "<<vwires[vw].TPC<<" "<<zwire.TPC<<std::endl;
		if (!geo->WireIDsIntersect(zwire,uwires[uw],widiuz)) continue;
		if (!geo->WireIDsIntersect(zwire,vwires[vw],widivz)) continue;
		if (!geo->WireIDsIntersect(uwires[uw],vwires[vw],widiuv)) continue;
		double dis1 = sqrt(pow(widiuz.y-widivz.y,2)+pow(widiuz.z-widivz.z,2));
		double dis2 = sqrt(pow(widiuz.y-widiuv.y,2)+pow(widiuz.z-widiuv.z,2));
		double dis3 = sqrt(pow(widiuv.y-widivz.y,2)+pow(widiuv.z-widivz.z,2));
		double maxdis = std::max(dis1,dis2);
		maxdis = std::max(maxdis,dis3);
		//std::cout<<z<<" "<<uw<<" "<<vw<<" "<<dis1<<" "<<dis2<<" "<<dis3<<" "<<zwire.Wire<<" "<<uwires[uw].Wire<<" "<<vwires[vw].Wire<<std::endl;
		if (maxdis<fDistanceCut){
		  //std::cout<<"*** "<<totalintersections<<std::endl;
		  ++totalintersections;
		  bestu = uw;
		  bestv = vw;
		}
	      }
	    }
	    if (totalintersections==1){
	      fDisambigHits.push_back(std::pair<art::Ptr<recob::Hit>, geo::WireID>(hitsUV[0][u],uwires[bestu]));
	      fDisambigHits.push_back(std::pair<art::Ptr<recob::Hit>, geo::WireID>(hitsUV[1][v],vwires[bestv]));
	      fHasBeenDisambigedUV[0][u] = 1+bestu;
	      fHasBeenDisambigedUV[1][v] = 1+bestv;
	      findmatch = true;
	      break;
	    }
	  }
	}
      }
    }
  }
  //Done finding trivial disambiguated hits

  
  //loop over undisambiguated hits, find the nearest channel of disambiguated hits and determine the correct wire segment.
  for (size_t i = 0; i<2; ++i){//loop over U and V hits
    for (size_t hit = 0; hit<hitsUV[i].size(); ++hit){
      if (fHasBeenDisambigedUV[i].find(hit)!=fHasBeenDisambigedUV[i].end()) continue;
      unsigned int apa1(0), cryo1(0);
      fAPAGeo.ChannelToAPA(hitsUV[i][hit]->Channel(), apa1, cryo1);
      unsigned int channdiff = 100000;
      geo::WireID nearestwire;
      for (auto& u2 : fHasBeenDisambigedUV[i]){
	unsigned int apa2(0), cryo2(0);
	fAPAGeo.ChannelToAPA(hitsUV[i][u2.first]->Channel(), apa2, cryo2);
	if (apa1!=apa2) continue;
	if (std::max(hitsUV[i][u2.first]->Channel(),hitsUV[i][hit]->Channel())
	    -std::min(hitsUV[i][u2.first]->Channel(),hitsUV[i][hit]->Channel())<channdiff){
	  channdiff = std::max(hitsUV[i][u2.first]->Channel(),hitsUV[i][hit]->Channel())
	    -std::min(hitsUV[i][u2.first]->Channel(),hitsUV[i][hit]->Channel());
	  nearestwire = geo->ChannelToWire(hitsUV[i][u2.first]->Channel())[u2.second-1];
	}
      }
      if (nearestwire.isValid){
	unsigned wirediff = 1000000;
	unsigned bestwire = 0;
	std::vector<geo::WireID> wires = geo->ChannelToWire(hitsUV[i][hit]->Channel());	
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
	  fDisambigHits.push_back(std::pair<art::Ptr<recob::Hit>, geo::WireID>(hitsUV[i][hit],wires[bestwire]));
	  fHasBeenDisambigedUV[i][hit] = 1+bestwire;
	}
      }
    }
  }
 
}

} //end namespace apa
