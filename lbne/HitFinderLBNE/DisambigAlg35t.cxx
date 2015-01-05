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
#include "Filters/ChannelFilter.h"

#include <map>
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>

#include "TH1D.h"

namespace lbne{

DisambigAlg35t::DisambigAlg35t(fhicl::ParameterSet const& pset)
  : fDBScan(pset.get< fhicl::ParameterSet >("DBScanAlg"))
{
  this->reconfigure(pset); 
}

//----------------------------------------------------------
void DisambigAlg35t::reconfigure(fhicl::ParameterSet const& p)
{

  fTimeCut = p.get<double>("TimeCut");
  fDistanceCut = p.get<double>("DistanceCut");
  fDistanceCutClu = p.get<double>("DistanceCutClu");
  fDBScan.reconfigure(p.get< fhicl::ParameterSet >("DBScanAlg"));
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

//  TH1D *histu = new TH1D("histu","histu",4000,0,4000);
//  TH1D *histv = new TH1D("histv","histv",4000,0,4000);
//  TH1D *histz = new TH1D("histz","histz",4000,0,4000);

  for (size_t i = 0; i<OrigHits.size(); ++i){
    switch (OrigHits[i]->View()){
    case geo::kU:
      hitsUV[0].push_back(OrigHits[i]);
//      if (OrigHits[i]->WireID().TPC==1) 
//	histu->Fill(OrigHits[i]->PeakTime()
//		    - detprop->GetXTicksOffset(0,
//					       OrigHits[i]->WireID().TPC,
//					       OrigHits[i]->WireID().Cryostat)
//		    ,OrigHits[i]->Charge());
      break;
    case geo::kV:
      hitsUV[1].push_back(OrigHits[i]);
//      if (OrigHits[i]->WireID().TPC==1) 
//	histv->Fill(OrigHits[i]->PeakTime()
//		    - detprop->GetXTicksOffset(1,
//					       OrigHits[i]->WireID().TPC,
//					       OrigHits[i]->WireID().Cryostat)
//		    ,OrigHits[i]->Charge());
      break;
    case geo::kZ:
      hitsZ.push_back(OrigHits[i]);
//      if (OrigHits[i]->WireID().TPC==1) 
//	histz->Fill(OrigHits[i]->PeakTime()
//		    - detprop->GetXTicksOffset(2,
//					       OrigHits[i]->WireID().TPC,
//					       OrigHits[i]->WireID().Cryostat)
//		    ,OrigHits[i]->Charge());
      break;
    default:
      throw cet::exception("DisambigAlg35t")
	<<": hit view unkonw. \n";
    }
  }
//  std::cout<<histu->GetMean()<<" "<<histv->GetMean()<<" "<<histz->GetMean()<<std::endl;
//  delete histu;
//  delete histv;
//  delete histz;
  std::vector<std::map<unsigned int, unsigned int> >fHasBeenDisambigedUV(2);
  const unsigned int ntpc = geo->NTPC();
  //hits and wireids for DBScan
  std::vector<std::vector<art::Ptr<recob::Hit> > > allhitsu(ntpc);
  std::vector<std::vector<art::Ptr<recob::Hit> > > allhitsv(ntpc);
  std::vector<std::vector<geo::WireID> > wireidsu(ntpc);
  std::vector<std::vector<geo::WireID> > wireidsv(ntpc);
  
  std::vector<std::vector<unsigned int> > cluidu(ntpc);
  std::vector<std::vector<unsigned int> > cluidv(ntpc);
  std::vector<std::vector<unsigned int> > bestwireidu(ntpc);
  std::vector<std::vector<unsigned int> > bestwireidv(ntpc);
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
	    unsigned int bestu = 0;  //index to wires associated with channel
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
		//if (hitsZ[z]->WireID().TPC==1) std::cout<<z<<" "<<uw<<" "<<vw<<" "<<dis1<<" "<<dis2<<" "<<dis3<<" "<<zwire.Wire<<" "<<uwires[uw].Wire<<" "<<vwires[vw].Wire<<std::endl;
		if (maxdis<fDistanceCut){
		  ++totalintersections;
		  bestu = uw;
		  bestv = vw;
		}
	      }
	    }
	    if (totalintersections==1){
//	      fDisambigHits.push_back(std::pair<art::Ptr<recob::Hit>, geo::WireID>(hitsUV[0][u],uwires[bestu]));
//	      fDisambigHits.push_back(std::pair<art::Ptr<recob::Hit>, geo::WireID>(hitsUV[1][v],vwires[bestv]));
//	      fHasBeenDisambigedUV[0][u] = 1+bestu;
//	      fHasBeenDisambigedUV[1][v] = 1+bestv;
	      allhitsu[uwires[bestu].TPC].push_back(hitsUV[0][u]);
	      allhitsv[vwires[bestv].TPC].push_back(hitsUV[1][v]);
	      wireidsu[uwires[bestu].TPC].push_back(uwires[bestu]);
	      wireidsv[vwires[bestv].TPC].push_back(vwires[bestv]);
	      cluidu[uwires[bestu].TPC].push_back(u);
	      cluidv[vwires[bestv].TPC].push_back(v);
	      bestwireidu[uwires[bestu].TPC].push_back(1+bestu);
	      bestwireidv[vwires[bestv].TPC].push_back(1+bestv);	      
	      findmatch = true;
	      break;
	    }
	  }//find v hit consistent in time
	}//loop over all v hits
      }//find u hit consistent in time
    }//loop over all u hits
  }//loop over all z hits
  //Done finding trivial disambiguated hits

  //running DB scan to identify and remove outlier hits
  // get the ChannelFilter
  filter::ChannelFilter chanFilt;
  for (size_t i = 0; i<ntpc; ++i){//loop over all TPCs
    if (!allhitsu[i].size()) continue;
    fDBScan.InitScan(allhitsu[i], chanFilt.SetOfBadChannels(), wireidsu[i]);
    fDBScan.run_cluster();
    if (allhitsu[i].size()!=fDBScan.fpointId_to_clusterId.size())
      throw cet::exception("DisambigAlg35t") <<"DBScan hits do not match input hits"<<allhitsu[i].size()<<" "<<fDBScan.fpointId_to_clusterId.size()<<"\n";
    //find hits associated with the biggest cluster
    int dbclu = -1;
    unsigned int maxhit = 0;
    std::vector<unsigned int> hitstoaddu;
    std::vector<unsigned int> dbcluhits(fDBScan.fclusters.size(),0);
    for(size_t j = 0; j < fDBScan.fpointId_to_clusterId.size(); ++j){
      if (fDBScan.fpointId_to_clusterId[j]>=0&&fDBScan.fpointId_to_clusterId[j]<fDBScan.fclusters.size())
      ++dbcluhits[fDBScan.fpointId_to_clusterId[j]];
    }
    for(size_t j = 0; j<dbcluhits.size(); ++j){
      if (dbcluhits[j]>maxhit){
	dbclu = j;
	maxhit = dbcluhits[j];
      }
    }
    for(size_t j = 0; j < fDBScan.fpointId_to_clusterId.size(); ++j){
      if (int(fDBScan.fpointId_to_clusterId[j]) == dbclu) hitstoaddu.push_back(j);
    }
    fDBScan.InitScan(allhitsv[i], chanFilt.SetOfBadChannels(), wireidsv[i]);
    fDBScan.run_cluster();
    if (allhitsv[i].size()!=fDBScan.fpointId_to_clusterId.size())
      throw cet::exception("DisambigAlg35t") <<"DBScan hits do not match input hits"<<allhitsv[i].size()<<" "<<fDBScan.fpointId_to_clusterId.size()<<"\n";
    dbclu = -1;
    maxhit = 0;
    std::vector<unsigned int> hitstoaddv;
    dbcluhits.clear();
    dbcluhits.resize(fDBScan.fclusters.size(),0);
    for(size_t j = 0; j < fDBScan.fpointId_to_clusterId.size(); ++j){
      if (fDBScan.fpointId_to_clusterId[j]>=0&&fDBScan.fpointId_to_clusterId[j]<fDBScan.fclusters.size())
	++dbcluhits[fDBScan.fpointId_to_clusterId[j]];
    }
    for(size_t j = 0; j<dbcluhits.size(); ++j){
      if (dbcluhits[j]>maxhit){
	dbclu = j;
	maxhit = dbcluhits[j];
      }
    }
    for(size_t j = 0; j < fDBScan.fpointId_to_clusterId.size(); ++j){
      if (int(fDBScan.fpointId_to_clusterId[j]) == dbclu) hitstoaddv.push_back(j);
    }
    for (size_t j = 0; j < allhitsu[i].size(); ++j){
      bool foundu = false;
      bool foundv = false;
      for (size_t k = 0; k<hitstoaddu.size(); ++k){
	if (j==hitstoaddu[k]) foundu = true;
      }
      for (size_t k = 0; k<hitstoaddv.size(); ++k){
	if (j==hitstoaddv[k]) foundv = true;
      }
      if (foundu&&foundv){
	fDisambigHits.push_back(std::pair<art::Ptr<recob::Hit>, geo::WireID>(allhitsu[i][j],wireidsu[i][j]));
	fDisambigHits.push_back(std::pair<art::Ptr<recob::Hit>, geo::WireID>(allhitsv[i][j],wireidsv[i][j]));
	fHasBeenDisambigedUV[0][cluidu[i][j]] = bestwireidu[i][j];
	fHasBeenDisambigedUV[1][cluidv[i][j]] = bestwireidv[i][j];
      }
    }
  }

  //loop over undisambiguated hits, find the nearest channel of disambiguated hits and determine the correct wire segment.
  for (size_t i = 0; i<2; ++i){//loop over U and V hits
    for (size_t hit = 0; hit<hitsUV[i].size(); ++hit){
      if (fHasBeenDisambigedUV[i].find(hit)!=fHasBeenDisambigedUV[i].end()) continue;
      unsigned int apa1(0), cryo1(0);
      fAPAGeo.ChannelToAPA(hitsUV[i][hit]->Channel(), apa1, cryo1);
      //unsigned int channdiff = 100000;
      //geo::WireID nearestwire;
      std::vector<geo::WireID> wires = geo->ChannelToWire(hitsUV[i][hit]->Channel());
      //std::vector<double> disttoallhits(wires.size(),-1);
      std::vector<int> nearbyhits(wires.size(),-1);
      for (auto& u2 : fHasBeenDisambigedUV[i]){
	unsigned int apa2(0), cryo2(0);
	fAPAGeo.ChannelToAPA(hitsUV[i][u2.first]->Channel(), apa2, cryo2);
	if (apa1!=apa2) continue;
	geo::WireID hitwire = geo->ChannelToWire(hitsUV[i][u2.first]->Channel())[u2.second-1];
	double wire_pitch = geo->WirePitch(0,1,hitwire.Plane,hitwire.TPC,hitwire.Cryostat);    //wire pitch in cm
	for (size_t w = 0; w<wires.size(); ++w){
	  if (wires[w].TPC!= hitwire.TPC) continue;
	  if (nearbyhits[w]<0) nearbyhits[w] = 0;
//	  if (disttoallhits[w]<0) disttoallhits[w] = 0;
//	  //std::cout<<hitwire.Wire<<" "<<wires[w].Wire<<" "<<double(hitwire.Wire)-double(wires[w].Wire)<<std::endl;
	  double distance = sqrt(pow((double(hitwire.Wire)-double(wires[w].Wire))*wire_pitch,2)+pow(detprop->ConvertTicksToX(hitsUV[i][u2.first]->PeakTime(),hitwire.Plane,hitwire.TPC,hitwire.Cryostat)-detprop->ConvertTicksToX(hitsUV[i][hit]->PeakTime(),hitwire.Plane,hitwire.TPC,hitwire.Cryostat),2));
//	  std::cout<<distance<<std::endl;
//	  disttoallhits[w] += distance;
	  if (distance<fDistanceCutClu) ++nearbyhits[w];
	}
//	if (std::max(hitsUV[i][u2.first]->Channel(),hitsUV[i][hit]->Channel())
//	    -std::min(hitsUV[i][u2.first]->Channel(),hitsUV[i][hit]->Channel())<channdiff){
//	  channdiff = std::max(hitsUV[i][u2.first]->Channel(),hitsUV[i][hit]->Channel())
//	    -std::min(hitsUV[i][u2.first]->Channel(),hitsUV[i][hit]->Channel());
//	  nearestwire = geo->ChannelToWire(hitsUV[i][u2.first]->Channel())[u2.second-1];
//	}
      }
      unsigned bestwire = 0;
      int maxnumhits = 0;
      for (size_t w = 0; w<wires.size(); ++w){
	if (nearbyhits[w]<0) continue;
	if (nearbyhits[w]>maxnumhits){
	  bestwire = w;
	  maxnumhits = nearbyhits[w];
	}
      }
      if (maxnumhits){
	fDisambigHits.push_back(std::pair<art::Ptr<recob::Hit>, geo::WireID>(hitsUV[i][hit],wires[bestwire]));
	fHasBeenDisambigedUV[i][hit] = 1+bestwire;
      }
      else{
	mf::LogWarning("DisambigAlg35t")<<"Could not find disambiguated hit for "<<hitsUV[i][hit]<<"/n";
      }
//      if (nearestwire.isValid){
//	unsigned wirediff = 1000000;
//	unsigned bestwire = 0;
//	std::vector<geo::WireID> wires = geo->ChannelToWire(hitsUV[i][hit]->Channel());	
//	for (size_t w = 0; w<wires.size(); ++w){
//	  if (wires[w].TPC!=nearestwire.TPC) continue;
//	  if (std::max(nearestwire.Wire,wires[w].Wire)-
//	      std::min(nearestwire.Wire,wires[w].Wire)<wirediff){
//	    wirediff = std::max(nearestwire.Wire,wires[w].Wire)-
//	      std::min(nearestwire.Wire,wires[w].Wire);
//	    bestwire = w;
//	  }
//	}
//	if (wirediff<1000000){
//	  fDisambigHits.push_back(std::pair<art::Ptr<recob::Hit>, geo::WireID>(hitsUV[i][hit],wires[bestwire]));
//	  fHasBeenDisambigedUV[i][hit] = 1+bestwire;
//	}
//      }

    }
  }
  
}

} //end namespace apa
