////////////////////////////////////////////////////////////////////////
/// \file    RegPixelMapProducer.h
/// \brief   RegPixelMapProducer for RegCVN modified from PixelMapProducer.h
/// \author  Ilsoo Seong - iseong@uci.edu
//
//  Modifications to allow unwrapped collection view
//   - Leigh Whitehead - leigh.howard.whitehead@cern.ch
////////////////////////////////////////////////////////////////////////

#include  <iostream>
#include  <ostream>
#include  <list>
#include  <algorithm>

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "dune/RegCVN/art/RegPixelMapProducer.h"
#include  "TVector2.h"

namespace cvn
{

  RegPixelMapProducer::RegPixelMapProducer(unsigned int nWire, unsigned int nTdc, double tRes):
  fNWire(nWire),
  fNTdc(nTdc),
  fTRes(tRes),
  fOffset{0,0}
  {}


  RegPixelMap RegPixelMapProducer::CreateMap(std::vector< art::Ptr< recob::Hit > >& cluster, art::FindManyP<recob::Wire> fmwire)
  {

    hitwireidx.clear();
    tmin_each_wire.clear();
    tmax_each_wire.clear();
    trms_max_each_wire.clear();
    fOffset[0] = 0; fOffset[1] = 0;

    RegCVNBoundary bound = DefineBoundary(cluster);

    return CreateMapGivenBoundary(cluster, bound, fmwire);


  }

  RegPixelMap RegPixelMapProducer::CreateMapGivenBoundary(std::vector< art::Ptr< recob::Hit > >& cluster,
                                                    const RegCVNBoundary& bound,
                                                    art::FindManyP<recob::Wire> fmwire)
  {

    RegPixelMap pm(fNWire, fNTdc, fTRes, bound);

    if (!fmwire.isValid()) return pm;

    auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();

    // get all raw adc of every hit wire
    for (size_t iwire = 0; iwire < hitwireidx.size(); ++iwire)
    {
      unsigned int iHit = hitwireidx[iwire];
      std::vector< art::Ptr<recob::Wire> > wireptr = fmwire.at(iHit);
      geo::WireID wireid = cluster[iHit]->WireID();

      for (size_t iwireptr = 0; iwireptr < wireptr.size(); ++iwireptr){
        std::vector<geo::WireID> wireids = geom->ChannelToWire(wireptr[iwireptr]->Channel());
        bool goodWID = false; 
        for (auto const & wid:wireids){ 
          if (wid.Plane == wireid.Plane &&   
              wid.Wire  == wireid.Wire &&
              wid.TPC   == wireid.TPC &&
              wid.Cryostat == wireid.Cryostat) goodWID = true;
        }
        if (!goodWID) continue;


        int t0_hit = (int)( tmin_each_wire[iwire] - 3 * (trms_max_each_wire[iwire]) );
        int t1_hit = (int)( tmax_each_wire[iwire] + 3 * (trms_max_each_wire[iwire]) );
	t0_hit = (t0_hit < 0) ? 0 : t0_hit;
	t1_hit = (t1_hit > 4491) ? 4491 : t1_hit;

	const std::vector<float>& signal = wireptr[0]->Signal();
        double globalWire = GetGlobalWire(wireid);
	if (wireid.TPC%2 == 1) {
	  if (wireid.Plane == 0) globalWire += fOffset[0];
	  if (wireid.Plane == 1) globalWire += fOffset[1];
	}
	for (int tt = t0_hit; tt <= t1_hit; ++tt)
        {
	  //double correctedadc = (double) signal[tt];
	  double correctedadc = ( (double)signal[tt] * TMath::Exp( (detprop->SamplingRate() * tt) / (detprop->ElectronLifetime()*1.e3) ) );
	  int tdc = tt;
  	  if (wireid.TPC%2 == 0) tdc = -tdc;
	  pm.Add((int)globalWire, tdc, wireid.Plane, correctedadc);
    
        } // end of tt
      } // end of iwireptr
    } // end of iwire

    //std::cout<< "===============> Offsets: " << fOffset[0] << " " << fOffset[1] << std::endl;
    return pm;

  }



  std::ostream& operator<<(std::ostream& os, const RegPixelMapProducer& p)
  {
    os << "RegPixelMapProducer: "
       << p.NTdc()  <<" tdcs X  " <<  p.NWire() << " wires";
    return os;
  }



  RegCVNBoundary RegPixelMapProducer::DefineBoundary(std::vector< art::Ptr< recob::Hit > >& cluster)
  {
    ShiftGlobalWire(cluster);

    std::vector<int> time_0;
    std::vector<int> time_1;
    std::vector<int> time_2;

    std::vector<int> wire_0;
    std::vector<int> wire_1;
    std::vector<int> wire_2;

    unsigned int temp_wire = cluster[0]->WireID().Wire;
    int temp_time_min = 99999;
    int temp_time_max = -99999;
    float temp_trms_max = -99999;
    for(size_t iHit = 0; iHit < cluster.size(); ++iHit)
    {
        geo::WireID wireid = cluster[iHit]->WireID();
	if (temp_wire != wireid.Wire){
		temp_wire = wireid.Wire;
		hitwireidx.push_back(iHit-1);
		tmin_each_wire.push_back(temp_time_min);
		tmax_each_wire.push_back(temp_time_max);
		trms_max_each_wire.push_back(temp_trms_max);
		temp_time_min = 99999; temp_time_max = -99999, temp_trms_max = -99999;
	}
	if (temp_time_min > cluster[iHit]->PeakTime()) temp_time_min = (int)cluster[iHit]->PeakTime();
	if (temp_time_max < cluster[iHit]->PeakTime()) temp_time_max = (int)cluster[iHit]->PeakTime();
	if (temp_trms_max < cluster[iHit]->RMS()) temp_trms_max = (float)cluster[iHit]->RMS();

	unsigned int planeid = wireid.Plane;
	int peaktime = (int)cluster[iHit]->PeakTime() ;
	if (wireid.TPC%2 == 0) peaktime = -peaktime;

        double globalWire = GetGlobalWire(wireid);
	if (wireid.TPC%2 == 1) {
	  if (wireid.Plane == 0) globalWire += fOffset[0];
	  if (wireid.Plane == 1) globalWire += fOffset[1];
	}

        if(planeid==0){
          time_0.push_back(peaktime);
          wire_0.push_back((int)globalWire);
        }
        if(planeid==1){
          time_1.push_back(peaktime);
          wire_1.push_back((int)globalWire);
        }
        if(planeid==2){
          time_2.push_back(peaktime);
          wire_2.push_back((int)globalWire);
        }
    }
  

    double tsum_0 = std::accumulate(time_0.begin(), time_0.end(), 0.0);
    double tmean_0 = tsum_0 / time_0.size();

    double tsum_1 = std::accumulate(time_1.begin(), time_1.end(), 0.0);
    double tmean_1 = tsum_1 / time_1.size();

    double tsum_2 = std::accumulate(time_2.begin(), time_2.end(), 0.0);
    double tmean_2 = tsum_2 / time_2.size();


    double wiresum_0 = std::accumulate(wire_0.begin(), wire_0.end(), 0.0);
    double wiremean_0 = wiresum_0 / wire_0.size();

    double wiresum_1 = std::accumulate(wire_1.begin(), wire_1.end(), 0.0);
    double wiremean_1 = wiresum_1 / wire_1.size();

    double wiresum_2 = std::accumulate(wire_2.begin(), wire_2.end(), 0.0);
    double wiremean_2 = wiresum_2 / wire_2.size();

    //std::cout << "TDC ===> " << (int)tmean_0 << " " << (int)tmean_1 << " " << (int)tmean_2 << std::endl;
    //std::cout << "Wire ==> " << (int)wiremean_0 << " " << (int)wiremean_1 << " " << (int)wiremean_2 << std::endl;

    //auto minwireelement_0= std::min_element(wire_0.begin(), wire_0.end());
    //std::cout<<"minwire 0: "<<*minwireelement_0<<std::endl;
    //auto minwireelement_1= std::min_element(wire_1.begin(), wire_1.end());
    //std::cout<<"minwire 1: "<<*minwireelement_1<<std::endl;
    //auto minwireelement_2= std::min_element(wire_2.begin(), wire_2.end());
    //std::cout<<"minwire 2: "<<*minwireelement_2<<std::endl;

    //auto maxwireelement_0= std::max_element(wire_0.begin(), wire_0.end());
    //std::cout<<"maxwire 0: "<<*maxwireelement_0<<std::endl;
    //auto maxwireelement_1= std::max_element(wire_1.begin(), wire_1.end());
    //std::cout<<"maxwire 1: "<<*maxwireelement_1<<std::endl;
    //auto maxwireelement_2= std::max_element(wire_2.begin(), wire_2.end());
    //std::cout<<"maxwire 2: "<<*maxwireelement_2<<std::endl;


    //int minwire_0 = *minwireelement_0-1;
    //int minwire_1 = *minwireelement_1-1;
    //int minwire_2 = *minwireelement_2-1;

    RegCVNBoundary bound(fNWire,fNTdc,fTRes,wiremean_0,wiremean_1,wiremean_2,tmean_0,tmean_1,tmean_2);

    return bound;
  }

  double RegPixelMapProducer::GetGlobalWire(const geo::WireID& wireID){
    // Get Global Wire Coordinate for RegCVN
    double globalWire = -9999;
    unsigned int nwires = geom->Nwires(wireID.Plane, 0, wireID.Cryostat); 
    // Induction
    if (geom->SignalType(wireID) == geo::kInduction) {
      double WireCentre[3] = {0};
      geom->WireIDToWireGeo(wireID).GetCenter(WireCentre);
      geo::PlaneID p1;
      int temp_tpc = 0;
      if (wireID.TPC % 2 == 0) { temp_tpc = 0; }
      else { temp_tpc = 1;  }
      p1 = geo::PlaneID(wireID.Cryostat, temp_tpc, wireID.Plane);
      globalWire = geom->WireCoordinate(WireCentre[1], WireCentre[2], p1);
    }
    // Collection
    else {
       int block = wireID.TPC / 4;
       globalWire = (double)( ((int)nwires*block) + wireID.Wire );
    }
    return round(globalWire);

  }

  std::vector<unsigned int> getUniques(std::vector<unsigned int> coll)
  {
    std::vector<unsigned int> uniques;
    for (unsigned int tpc : coll) 
    {
      if (std::find(uniques.begin(), uniques.end(), tpc) == uniques.end())
          uniques.push_back(tpc);
    }
   return uniques;
  }

  void RegPixelMapProducer::ShiftGlobalWire(std::vector< art::Ptr< recob::Hit > >& cluster)
  {
    // find hits passing through an APA
    std::vector<unsigned int> list_TPCs;
    std::map<unsigned int, double> map_mean_U, map_mean_V, map_Twei_U, map_Twei_V;
    for (unsigned int iHit = 0; iHit < cluster.size(); ++iHit)
    {
	geo::WireID wireid = cluster[iHit]->WireID();
        if (wireid.Plane == 2) continue; // skip collection plane
	// select 30 ticks
	if (cluster[iHit]->PeakTime() < 30){
		double Twei = cluster[iHit]->PeakTime()+1.; //
		list_TPCs.push_back(wireid.TPC);
		if (wireid.Plane == 0){
		    map_mean_U[wireid.TPC] += GetGlobalWire(wireid)/Twei;
		    map_Twei_U[wireid.TPC] += 1.0/Twei;
		}
		if (wireid.Plane == 1){
		    map_mean_V[wireid.TPC] += GetGlobalWire(wireid)/Twei;
		    map_Twei_V[wireid.TPC] += 1.0/Twei;
		}
	}
    }
    std::vector<unsigned int> unique_TPCs = getUniques(list_TPCs);
    for (unsigned int tmp_idx : unique_TPCs)
    {
      if (map_Twei_U[tmp_idx]>0) map_mean_U[tmp_idx] /= map_Twei_U[tmp_idx];
      if (map_Twei_V[tmp_idx]>0) map_mean_V[tmp_idx] /= map_Twei_V[tmp_idx];
    }
    // obtain offsets for U and V planes
    if (unique_TPCs.size() > 1){
      // shift global wire for an event passing through an APA
      double offsetU = 0, offsetV = 0;
      for (unsigned int ii = 0; ii < unique_TPCs.size()-1; ++ii)
      {
        unsigned int tpc1 = unique_TPCs[ii];
        unsigned int tpc2 = unique_TPCs[ii+1];
        if (tpc1/4 != tpc2/4) continue;
        if (tpc1%2 == tpc2%2) continue;
        if (map_Twei_U[tpc1] > 0 && map_Twei_U[tpc2] > 0) offsetU = map_mean_U[tpc1] - map_mean_U[tpc2];
        if (map_Twei_V[tpc1] > 0 && map_Twei_V[tpc2] > 0) offsetV = map_mean_V[tpc1] - map_mean_V[tpc2];
      }
      if (abs(offsetU) > 0 && abs(offsetV) > 0){ //FIXME
        fOffset[0] = round(offsetU);
        fOffset[1] = round(offsetV);
      }
    }

  }

}
