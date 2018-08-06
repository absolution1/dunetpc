////////////////////////////////////////////////////////////////////////
/// \file  ChannelMapCRMAlg.cxx
/// \brief Interface to algorithm class for a dual-phase detector channel mapping
///
/// \version $Id:  $
/// \author  brebel@fnal.gov vgalymov@ipnl.in2p3.fr
////////////////////////////////////////////////////////////////////////

#include "dune/Geometry/ChannelMapCRMAlg.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/AuxDetGeo.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"

#include "messagefacility/MessageLogger/MessageLogger.h" 

namespace geo{

  //----------------------------------------------------------------------------
  ChannelMapCRMAlg::ChannelMapCRMAlg(fhicl::ParameterSet const& p)
    : fSorter(geo::GeoObjectSorterCRM(p))
  {
  }

  //----------------------------------------------------------------------------
  void ChannelMapCRMAlg::Initialize( GeometryData_t const& geodata )
  {
    // start over:
    Uninitialize();
    
    std::vector<geo::CryostatGeo> const& cgeo = geodata.cryostats;
    
    fNcryostat = cgeo.size();
    
    mf::LogInfo("ChannelMapCRMAlg") << "Initializing CRM ChannelMap...";

    fNTPC.resize(fNcryostat);
    fWireCounts.resize(fNcryostat);
    fNPlanes.resize(fNcryostat);
    fFirstWireProj.resize(fNcryostat);
    fOrthVectorsY.resize(fNcryostat);
    fOrthVectorsZ.resize(fNcryostat);
    fPlaneBaselines.resize(fNcryostat);
    fWiresPerPlane.resize(fNcryostat);
    fFirstChannelInNextPlane.resize(fNcryostat);
    fFirstChannelInThisPlane.resize(fNcryostat);
    fPlaneIDs.clear();
    fTopChannel = 0;

    int RunningTotal = 0;

    for(unsigned int cs = 0; cs != fNcryostat; ++cs){
      
      fNTPC[cs] = cgeo[cs].NTPC();
      
      // Size up all the vectors 
      fWireCounts[cs]             .resize(fNTPC[cs]);
      fFirstWireProj[cs]          .resize(fNTPC[cs]);
      fOrthVectorsY[cs]           .resize(fNTPC[cs]);
      fOrthVectorsZ[cs]           .resize(fNTPC[cs]);
      fPlaneBaselines[cs]         .resize(fNTPC[cs]);
      fWiresPerPlane[cs]          .resize(fNTPC[cs]);
      fNPlanes[cs]                .resize(fNTPC[cs]);
      fFirstChannelInThisPlane[cs].resize(fNTPC[cs]);
      fFirstChannelInNextPlane[cs].resize(fNTPC[cs]);

      for(unsigned int TPCCount = 0; TPCCount != fNTPC[cs]; ++TPCCount)
	{
	  unsigned int PlanesThisTPC = cgeo[cs].TPC(TPCCount).Nplanes();
	  fWireCounts[cs][TPCCount]   .resize(PlanesThisTPC);
	  fFirstWireProj[cs][TPCCount].resize(PlanesThisTPC);
	  fOrthVectorsY[cs][TPCCount] .resize(PlanesThisTPC);
	  fOrthVectorsZ[cs][TPCCount] .resize(PlanesThisTPC);
	  fNPlanes[cs][TPCCount]=PlanesThisTPC;
	  
	  for(unsigned int PlaneCount = 0; PlaneCount != PlanesThisTPC; ++PlaneCount){
	    
	    if(PlaneCount >= 2) //should only find two views for dual-phase
	      throw cet::exception("Geometry") <<"CANNOT HAVE more than two readout planes for dual-phase" 
					       << "\n";
	    /*
	    if(cgeo[cs].TPC(TPCCount).Plane(PlaneCount).View() == geo::kU)
	      std::cout<<" View is kU "<<std::endl;
	    else if(cgeo[cs].TPC(TPCCount).Plane(PlaneCount).View() == geo::kV)
	      std::cout<<" View is kV "<<std::endl;
	    else if(cgeo[cs].TPC(TPCCount).Plane(PlaneCount).View() == geo::kZ)
	      std::cout<<" View is kZ "<<std::endl;

	    if(cgeo[cs].TPC(TPCCount).Plane(PlaneCount).SignalType() == geo::kInduction)
	      std::cout<<" View is kInduction "<<std::endl;
	    if(cgeo[cs].TPC(TPCCount).Plane(PlaneCount).SignalType() == geo::kCollection)
	      std::cout<<" View is kCollection "<<std::endl;
	    */

	    fPlaneIDs.emplace(PlaneID(cs, TPCCount, PlaneCount));
	    double ThisWirePitch = cgeo[cs].TPC(TPCCount).WirePitch(PlaneCount);
	    fWireCounts[cs][TPCCount][PlaneCount] = cgeo[cs].TPC(TPCCount).Plane(PlaneCount).Nwires();
	    
	    double  WireCentre1[3] = {0.,0.,0.};
	    double  WireCentre2[3] = {0.,0.,0.};
	  
	    const geo::WireGeo& firstWire = cgeo[cs].TPC(TPCCount).Plane(PlaneCount).Wire(0);
	    const double sth = firstWire.SinThetaZ(), cth = firstWire.CosThetaZ();

	    firstWire.GetCenter(WireCentre1,0);
	    cgeo[cs].TPC(TPCCount).Plane(PlaneCount).Wire(1).GetCenter(WireCentre2,0);
	    
	    // figure out if we need to flip the orthogonal vector 
	    // (should point from wire n -> n+1)
	    double OrthY = cth, OrthZ = -sth;
	    if(((WireCentre2[1] - WireCentre1[1])*OrthY 
		+ (WireCentre2[2] - WireCentre1[2])*OrthZ) < 0){
	      OrthZ *= -1;
	      OrthY *= -1;
	    }
          
          // Overall we are trying to build an expression that looks like
          //   int NearestWireNumber = round((worldPos.OrthVector - FirstWire.OrthVector)/WirePitch);     
          // That runs as fast as humanly possible.
          // We predivide everything by the wire pitch so we don't do this in the loop.
          //
          // Putting this together into the useful constants we will use later per plane and tpc:
          fOrthVectorsY[cs][TPCCount][PlaneCount] = OrthY / ThisWirePitch;
          fOrthVectorsZ[cs][TPCCount][PlaneCount] = OrthZ / ThisWirePitch;
          
          fFirstWireProj[cs][TPCCount][PlaneCount]  = WireCentre1[1]*OrthY + WireCentre1[2]*OrthZ;
          fFirstWireProj[cs][TPCCount][PlaneCount] /= ThisWirePitch;
          
          // now to count up wires in each plane and get first channel in each plane
          int WiresThisPlane = cgeo[cs].TPC(TPCCount).Plane(PlaneCount).Nwires();        
          fWiresPerPlane[cs] .at(TPCCount).push_back(WiresThisPlane);
          fPlaneBaselines[cs].at(TPCCount).push_back(RunningTotal);
          
          RunningTotal += WiresThisPlane;

          fFirstChannelInThisPlane[cs].at(TPCCount).push_back(fTopChannel);
          fTopChannel += WiresThisPlane;
          fFirstChannelInNextPlane[cs].at(TPCCount).push_back(fTopChannel);

        }// end loop over planes
      }// end loop over TPCs
    }// end loop over cryostats

    // calculate the total number of channels in the detector
    fNchannels = fTopChannel;

    LOG_DEBUG("ChannelMapCRM") << "# of channels is " << fNchannels;


    return;
  }
   
  //----------------------------------------------------------------------------
  void ChannelMapCRMAlg::Uninitialize()
  {
  }

  //----------------------------------------------------------------------------
  std::vector<geo::WireID> ChannelMapCRMAlg::ChannelToWire(raw::ChannelID_t channel)  const
  {
    std::vector< geo::WireID > AllSegments; 
    unsigned int cstat = 0;
    unsigned int tpc   = 0;
    unsigned int plane = 0;
    unsigned int wire  = 0;
    
    // first check if this channel ID is legal
    if(channel > fTopChannel)
      throw cet::exception("Geometry") << "ILLEGAL CHANNEL ID for channel " << channel << "\n";
    
    // then go find which plane, tpc and cryostat it is in from the information we stored earlier
    bool foundWid(false);
    for(unsigned int csloop = 0; csloop != fNcryostat; ++csloop){
      for(unsigned int tpcloop = 0; tpcloop != fNTPC[csloop]; ++tpcloop){
        for(unsigned int planeloop = 0;
	    planeloop!=fFirstChannelInNextPlane[csloop][tpcloop].size();
	    ++planeloop)
	  {
	    if(channel < fFirstChannelInNextPlane[csloop][tpcloop][planeloop]){
	      cstat = csloop;
	      tpc   = tpcloop;
	      plane = planeloop;
	      wire  = channel - fFirstChannelInThisPlane[cstat][tpcloop][planeloop];
	      foundWid = true;
	      break;
	    }     
          if(foundWid) break;
        }// end plane loop
        if(foundWid) break;
      }// end tpc loop
      if(foundWid) break;
    }// end cryostat loop

    geo::WireID CodeWire(cstat, tpc, plane, wire);
    
    AllSegments.push_back(CodeWire);

    return AllSegments;
  }

  //----------------------------------------------------------------------------
  raw::ChannelID_t ChannelMapCRMAlg::Nchannels() const
  {
    return fNchannels;
  }

  //----------------------------------------------------------------------------
  unsigned int ChannelMapCRMAlg::Nchannels(readout::ROPID const& ropid) const {
    if (!HasROP(ropid)) return 0;
    // The number of channels matches the number of wires. Life is easy.
    return WireCount(FirstWirePlaneInROP(ropid));
  } // ChannelMapCRMAlg::Nchannels()
  
  
  //----------------------------------------------------------------------------
  double ChannelMapCRMAlg::WireCoordinate
    (double YPos, double ZPos, geo::PlaneID const& planeID) const
  {
    // Returns the wire number corresponding to a (Y,Z) position in PlaneNo 
    // with float precision.
    // B. Baller August 2014
    return YPos*AccessElement(fOrthVectorsY, planeID) 
         + ZPos*AccessElement(fOrthVectorsZ, planeID)
         - AccessElement(fFirstWireProj, planeID);
  }

  
  //----------------------------------------------------------------------------
  WireID ChannelMapCRMAlg::NearestWireID
    (const TVector3& worldPos, geo::PlaneID const& planeID) const
  {

    // This part is the actual calculation of the nearest wire number, where we assume
    //  uniform wire pitch and angle within a wireplane
    
    // add 0.5 to have the correct rounding
    int NearestWireNumber = int
      (0.5 + WireCoordinate(worldPos.Y(), worldPos.Z(), planeID));
    
    // If we are outside of the wireplane range, throw an exception
    // (this response maintains consistency with the previous
    // implementation based on geometry lookup)
    if(NearestWireNumber < 0 ||
       (unsigned int) NearestWireNumber >= WireCount(planeID))
    {
      int wireNumber = NearestWireNumber; // save for the output
      
      if(NearestWireNumber < 0 ) NearestWireNumber = 0;
      else                       NearestWireNumber = WireCount(planeID) - 1;
      
      throw InvalidWireIDError("Geometry", wireNumber, NearestWireNumber)
        << "Can't Find Nearest Wire for position (" 
        << worldPos[0] << "," << worldPos[1] << "," << worldPos[2] << ")"
        << " in plane " << std::string(planeID) << " approx wire number # "
        << wireNumber << " (capped from " << NearestWireNumber << ")\n";
    }

    return geo::WireID(planeID, (geo::WireID::WireID_t) NearestWireNumber);
  }
  
  //----------------------------------------------------------------------------
  // This method returns the channel number, assuming the numbering scheme
  // is heirachical - that is, channel numbers run in order, for example:
  //                                             (Ben J Oct 2011)                   
  //                    Wire1     | 0
  //           Plane1 { Wire2     | 1
  //    TPC1 {          Wire3     | 2
  //           Plane2 { Wire1     | 3   increasing channel number
  //                    Wire2     | 4     (with no gaps)
  //    TPC2 { Plane1 { Wire1     | 5
  //           Plane2 { Wire1     | 6
  //                    Wire2     v 7
  //
  raw::ChannelID_t ChannelMapCRMAlg::PlaneWireToChannel
    (geo::WireID const& wireID) const
  {
    unsigned int const* pBaseLine = GetElementPtr(fPlaneBaselines, wireID);
    // This is the actual lookup part - first make sure coordinates are legal
    if (pBaseLine) {
      // if the channel has legal coordinates, its ID is given by the wire
      // number above the number of wires in lower planes, tpcs and cryostats
      return *pBaseLine + wireID.Wire;
    }
    else{  
      // if the coordinates were bad, throw an exception
      throw cet::exception("ChannelMapCRMAlg")
        << "NO CHANNEL FOUND for " << std::string(wireID);
    }
    
    // made it here, that shouldn't happen, return raw::InvalidChannelID
    mf::LogWarning("ChannelMapCRMAlg") << "should not be at the point in the function, returning "
                                            << "invalid channel";
    return raw::InvalidChannelID;

  }


  //----------------------------------------------------------------------------
  SigType_t ChannelMapCRMAlg::SignalTypeForChannelImpl(raw::ChannelID_t const channel) const
  {

    // still assume one cryostat for now -- faster
    unsigned int nChanPerTPC = fNchannels/fNTPC[0];
    // casting wil trunc towards 0 -- faster than floor
    unsigned int tpc = channel / nChanPerTPC;  
    //need number of planes to know Collection 
    unsigned int PlanesThisTPC = fNPlanes[0][tpc];
    
      
    //for dual-phase force a return of two collection views
    SigType_t sigt = geo::kMysteryType;
    if(      (channel >= fFirstChannelInThisPlane[0][tpc][0]) &&
             (channel <  fFirstChannelInNextPlane[0][tpc][PlanesThisTPC-1])    )
      { sigt = geo::kCollection; }
    else
      mf::LogWarning("BadChannelSignalType") << "Channel " << channel
                                             << " not given signal type." << std::endl;
            

    return sigt;
  }


  //----------------------------------------------------------------------------
  /* // this code should be equivalent to the logic implemented in geo::PlaneGeo::UpdateView()
  View_t ChannelMapCRMAlg::View(raw::ChannelID_t const channel) const
  {

    // still assume one cryostat for now -- faster
    unsigned int nChanPerTPC = fNchannels/fNTPC[0];
    // casting wil trunc towards 0 -- faster than floor
    unsigned int tpc = channel / nChanPerTPC; 
    //unsigned int PlanesThisTPC = fNPlanes[0][tpc];

    View_t view = geo::kUnknown; 
    
    if(      (channel >= fFirstChannelInThisPlane[0][tpc][0]) &&
             (channel <  fFirstChannelInNextPlane[0][tpc][0])    ){ view = geo::kY; } 
    else if( (channel >= fFirstChannelInThisPlane[0][tpc][1]) &&
             (channel <  fFirstChannelInNextPlane[0][tpc][1])    ){ view = geo::kZ; } 
    //else if( (channel >= fFirstChannelInThisPlane[0][tpc][2]) &&
    //(channel <  fFirstChannelInNextPlane[0][tpc][2])    ){ view = geo::kZ; }
    else
      mf::LogWarning("BadChannelSignalType") << "Channel " << channel
                                             << " not given view type.";

    return view;
  }  
  */

  //----------------------------------------------------------------------------
  std::set<PlaneID> const& ChannelMapCRMAlg::PlaneIDs() const
  {
    return fPlaneIDs;
  }

  
  //----------------------------------------------------------------------------
  /**
   * @brief Returns the total number of TPC sets in the specified cryostat
   * @param cryoid cryostat ID
   * @return number of TPC sets in the cryostat, or 0 if no cryostat found
   */
  unsigned int ChannelMapCRMAlg::NTPCsets
    (readout::CryostatID const& cryoid) const
  {
    // return the same number as the number of TPCs
    return (cryoid.isValid && cryoid.Cryostat < fNTPC.size())?
      fNTPC[cryoid.Cryostat]: 0;
  } // ChannelMapCRMAlg::NTPCsets()
  
  
  //----------------------------------------------------------------------------
  /// Returns the largest number of TPC sets any cryostat in the detector has
  unsigned int ChannelMapCRMAlg::MaxTPCsets() const {
    return MaxTPCs();
  } // ChannelMapCRMAlg::MaxTPCsets()
  
  
  //----------------------------------------------------------------------------
  /// Returns whether we have the specified TPC set
  /// @return whether the TPC set is valid and exists
  bool ChannelMapCRMAlg::HasTPCset(readout::TPCsetID const& tpcsetid) const {
    return tpcsetid.TPCset < NTPCsets(tpcsetid);
  } // ChannelMapCRMAlg::HasTPCset()
  
  
  //----------------------------------------------------------------------------
  /// Returns the ID of the TPC set tpcid belongs to
  readout::TPCsetID ChannelMapCRMAlg::TPCtoTPCset
    (geo::TPCID const& tpcid) const
  {
    return ConvertTPCtoTPCset(tpcid);
  } // ChannelMapCRMAlg::TPCtoTPCset()
  
  
  //----------------------------------------------------------------------------
  /**
   * @brief Returns a list of ID of TPCs belonging to the specified TPC set
   * @param tpcsetid ID of the TPC set to convert into TPC IDs
   * @return the list of TPCs, empty if TPC set is invalid
   *
   * Note that the check is performed on the validity of the TPC set ID, that
   * does not necessarily imply that the TPC set specified by the ID actually
   * exists. Check the existence of the TPC set first (HasTPCset()).
   * Behaviour on valid, non-existent TPC set IDs is undefined.
   */
  std::vector<geo::TPCID> ChannelMapCRMAlg::TPCsetToTPCs
    (readout::TPCsetID const& tpcsetid) const
  {
    std::vector<geo::TPCID> IDs;
    if (tpcsetid.isValid) IDs.emplace_back(ConvertTPCsetToTPC(tpcsetid));
    return IDs;
  } // ChannelMapCRMAlg::TPCsetToTPCs()
  
  
  //----------------------------------------------------------------------------
  /// Returns the ID of the first TPC belonging to the specified TPC set
  geo::TPCID ChannelMapCRMAlg::FirstTPCinTPCset
    (readout::TPCsetID const& tpcsetid) const
  {
    return ConvertTPCsetToTPC(tpcsetid);
  } // ChannelMapCRMAlg::FirstTPCinTPCset()
  
  
  //----------------------------------------------------------------------------
  unsigned int ChannelMapCRMAlg::MaxTPCs() const
  {
    unsigned int max = 0;
    for (unsigned int nTPCs: fNTPC) if (nTPCs > max) max = nTPCs;
    return max;
  } // ChannelMapCRMAlg::MaxTPCs()
  
  
  //----------------------------------------------------------------------------
  /**
   * @brief Returns the total number of ROP in the specified TPC set
   * @param tpcsetid TPC set ID
   * @return number of readout planes in the TPC set, or 0 if no TPC set found
   * 
   * Note that this methods explicitly check the existence of the TPC set.
   */
  unsigned int ChannelMapCRMAlg::NROPs(readout::TPCsetID const& tpcsetid) const
  {
    if (!HasTPCset(tpcsetid)) return 0;
    return AccessElement(fNPlanes, FirstTPCinTPCset(tpcsetid));
  } // ChannelMapCRMAlg::NROPs()
  
  
  //----------------------------------------------------------------------------
  /// Returns the largest number of ROPs a TPC set in the detector has
  unsigned int ChannelMapCRMAlg::MaxROPs() const {
    unsigned int max = 0;
    for (auto const& cryo_tpc: fNPlanes)
      for (unsigned int nPlanes: cryo_tpc)
        if (nPlanes > max) max = nPlanes;
    return max;
  } // ChannelMapCRMAlg::MaxROPs()
  
  
  //----------------------------------------------------------------------------
  /// Returns whether we have the specified ROP
  /// @return whether the readout plane is valid and exists
  bool ChannelMapCRMAlg::HasROP(readout::ROPID const& ropid) const {
    return ropid.ROP < NROPs(ropid);
  } // ChannelMapCRMAlg::HasROP()
  
  
  //----------------------------------------------------------------------------
  /// Returns the ID of the ROP planeid belongs to
  readout::ROPID ChannelMapCRMAlg::WirePlaneToROP
    (geo::PlaneID const& planeid) const
  {
    return ConvertWirePlaneToROP(planeid);
  } // ChannelMapCRMAlg::WirePlaneToROP()
  
  
  //----------------------------------------------------------------------------
  /// Returns a list of ID of planes belonging to the specified ROP
  std::vector<geo::PlaneID> ChannelMapCRMAlg::ROPtoWirePlanes
    (readout::ROPID const& ropid) const
  {
    std::vector<geo::PlaneID> IDs;
    if (ropid.isValid) IDs.emplace_back(FirstWirePlaneInROP(ropid));
    return IDs;
  } // ChannelMapCRMAlg::ROPtoWirePlanes()
  
  
  //----------------------------------------------------------------------------
  /// Returns the ID of the first plane belonging to the specified ROP
  geo::PlaneID ChannelMapCRMAlg::FirstWirePlaneInROP
    (readout::ROPID const& ropid) const
  {
    return ConvertROPtoWirePlane(ropid);
  } // ChannelMapCRMAlg::FirstWirePlaneInROP()
  
  
  //----------------------------------------------------------------------------
  /**
   * @brief Returns a list of ID of TPCs the specified ROP spans
   * @param ropid ID of the readout plane
   * @return the list of TPC IDs, empty if readout plane ID is invalid
   *
   * Note that this check is performed on the validity of the readout plane
   * ID, that does not necessarily imply that the readout plane specified by
   * the ID actually exists. Check if the ROP exists with HasROP().
   * The behaviour on non-existing readout planes is undefined.
   */
  std::vector<geo::TPCID> ChannelMapCRMAlg::ROPtoTPCs
    (readout::ROPID const& ropid) const
  {
    std::vector<geo::TPCID> IDs;
    // we take the TPC set of the ROP and convert it straight into a TPC ID
    if (ropid.isValid) IDs.emplace_back(ConvertTPCsetToTPC(ropid.asTPCsetID()));
    return IDs;
  } // ChannelMapCRMAlg::ROPtoTPCs()
  
  
  //----------------------------------------------------------------------------
  /// Returns the ID of the ROP the channel belongs to
  /// @throws cet::exception (category: "Geometry") if non-existent channel
  readout::ROPID ChannelMapCRMAlg::ChannelToROP(raw::ChannelID_t channel) const
  {
    if (!raw::isValidChannelID(channel)) return {}; // invalid ROP returned
    
    // which wires does the channel cover?
    std::vector<geo::WireID> wires = ChannelToWire(channel);
    
    // - none:
    if (wires.empty()) return {}; // default-constructed ID, invalid
    
    // - one: maps its plane ID into a ROP ID
    return WirePlaneToROP(wires[0]);
  } // ChannelMapCRMAlg::ChannelToROP()
  
  
  //----------------------------------------------------------------------------
  /**
   * @brief Returns the ID of the first channel in the specified readout plane
   * @param ropid ID of the readout plane
   * @return ID of first channel, or raw::InvalidChannelID if ID is invalid
   * 
   * Note that this check is performed on the validity of the readout plane
   * ID, that does not necessarily imply that the readout plane specified by
   * the ID actually exists. Check if the ROP exists with HasROP().
   * The behaviour for non-existing readout planes is undefined.
   */
  raw::ChannelID_t ChannelMapCRMAlg::FirstChannelInROP
    (readout::ROPID const& ropid) const
  {
    if (!ropid.isValid) return raw::InvalidChannelID;
    return (raw::ChannelID_t)
      AccessElement(fPlaneBaselines, ConvertROPtoWirePlane(ropid));
  } // ChannelMapCRMAlg::FirstChannelInROP()
  
  
  //----------------------------------------------------------------------------
  readout::TPCsetID ChannelMapCRMAlg::ConvertTPCtoTPCset
    (geo::TPCID const& tpcid)
  {
    if (!tpcid.isValid) return {}; // invalid ID, default-constructed
    return {
      (readout::CryostatID::CryostatID_t) tpcid.Cryostat,
      (readout::TPCsetID::TPCsetID_t) tpcid.TPC
      };
  } // ChannelMapCRMAlg::ConvertTPCtoTPCset()
  
  
  //----------------------------------------------------------------------------
  geo::TPCID ChannelMapCRMAlg::ConvertTPCsetToTPC
    (readout::TPCsetID const& tpcsetid)
  {
    if (!tpcsetid.isValid) return {};
    return {
      (geo::CryostatID::CryostatID_t) tpcsetid.Cryostat,
      (geo::TPCID::TPCID_t) tpcsetid.TPCset
      };
  } // ChannelMapCRMAlg::ConvertTPCsetToTPC()
  
  
  //----------------------------------------------------------------------------
  readout::ROPID ChannelMapCRMAlg::ConvertWirePlaneToROP
    (geo::PlaneID const& planeid)
  {
    if (!planeid.isValid) return {}; // invalid ID, default-constructed
    return {
      (readout::CryostatID::CryostatID_t) planeid.Cryostat,
      (readout::TPCsetID::TPCsetID_t) planeid.TPC,
      (readout::ROPID::ROPID_t) planeid.Plane
      };
    
  } // ChannelMapCRMAlg::ConvertWirePlaneToROP()
  
  
  //----------------------------------------------------------------------------
  geo::PlaneID ChannelMapCRMAlg::ConvertROPtoWirePlane
    (readout::ROPID const& ropid)
  {
    if (!ropid.isValid) return {};
    return {
      (geo::CryostatID::CryostatID_t) ropid.Cryostat,
      (geo::TPCID::TPCID_t) ropid.TPCset,
      (geo::PlaneID::PlaneID_t) ropid.ROP
      };
  } // ChannelMapCRMAlg::ConvertROPtoWirePlane()
  
  
  //----------------------------------------------------------------------------
  
} // namespace
