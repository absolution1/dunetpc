////////////////////////////////////////////////////////////////////////
/// \file  ProtoDUNEChannelMapAlg.cxx
/// \brief Interface to algorithm class for a specific detector channel mapping
///
/// \version $Id:  $
/// \author  tylerdalion@gmail.com
////////////////////////////////////////////////////////////////////////

#include "dune/Geometry/ProtoDUNEChannelMapAlg.h"
#include "dune/Geometry/DuneApaChannelMapAlg.h"
#include "larcorealg/Geometry/GeoObjectSorter.h"
#include "dune/Geometry/GeoObjectSorterAPA.h"

#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/AuxDetGeo.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h" 

#include <iomanip>

using std::map;
using geo::ProtoDUNEChannelMapAlg;
using std::cout;
using std::endl;
using std::setw;

typedef unsigned int Index;

namespace geo {

//----------------------------------------------------------------------------

  ProtoDUNEChannelMapAlg::ProtoDUNEChannelMapAlg(const fhicl::ParameterSet& p)
    : DuneApaChannelMapAlg(p), 
      fMaxOpChannel(0), fNOpChannels(0)
  {
    fSorter = new GeoObjectSorterAPA(p);
  }


//----------------------------------------------------------------------------

  void ProtoDUNEChannelMapAlg::Initialize(GeometryData_t const& geodata) {
    Uninitialize();
  
    // Call the base-class version to perform TPC-related actions
    DuneApaChannelMapAlg::Initialize(geodata);

    // Manually entered based on maps from Chris Macias
    // http://indico.fnal.gov/event/14582/session/30/contribution/30/material/slides/0.pdf
    //
    // That gives SSP and channel ranges from OpDet. These offline channels correspond
    // to the APAs listed. Channel numbers increase by 2 going from top to bottom
    //
    //   USDaS   MSDaS  DSDaS
    //   41-59   21-39  1-19
    //
    //   USRaS   MSRaS  DSRaS
    //   40-58   20-38  0-18


    // DSDaS -- "normal" pattern
    fSSP[ 1] = 41; fSSPChOne[ 1] = 0;
    fSSP[ 3] = 41; fSSPChOne[ 3] = 4;
    fSSP[ 5] = 41; fSSPChOne[ 5] = 8;
    fSSP[ 7] = 42; fSSPChOne[ 7] = 0;
    fSSP[ 9] = 42; fSSPChOne[ 9] = 4;
    fSSP[11] = 42; fSSPChOne[11] = 8;
    fSSP[13] = 43; fSSPChOne[13] = 0;
    fSSP[15] = 43; fSSPChOne[15] = 4;
    fSSP[17] = 43; fSSPChOne[17] = 8;
    fSSP[19] = 44; fSSPChOne[19] = 0;

    // MSDaS -- unusual pattern for ARAPUCA
    fSSP[21] = 61; fSSPChOne[21] = 0;
    fSSP[23] = 61; fSSPChOne[23] = 4;
    fSSP[25] = 61; fSSPChOne[25] = 8;
    fSSP[27] = 62; fSSPChOne[27] = 0;
    fSSP[29] = 62; fSSPChOne[29] = 4;
    fSSP[31] = 63; fSSPChOne[31] = 0; // ARAPUCA
    fSSP[33] = 62; fSSPChOne[33] = 8;
    fSSP[35] = 64; fSSPChOne[35] = 0;
    fSSP[37] = 64; fSSPChOne[37] = 4;
    fSSP[39] = 64; fSSPChOne[39] = 8;

    // USDaS -- unusual pattern
    fSSP[41] = 53; fSSPChOne[41] = 0;
    fSSP[43] = 53; fSSPChOne[43] = 4;
    fSSP[45] = 53; fSSPChOne[45] = 8;
    fSSP[47] = 54; fSSPChOne[47] = 0;
    fSSP[49] = 51; fSSPChOne[49] = 0;
    fSSP[51] = 54; fSSPChOne[51] = 4;
    fSSP[53] = 51; fSSPChOne[53] = 4;
    fSSP[55] = 51; fSSPChOne[55] = 8;
    fSSP[57] = 54; fSSPChOne[57] = 8;
    fSSP[59] = 52; fSSPChOne[59] = 0;
  
    // DSRaS -- "normal" pattern
    fSSP[ 0] = 11; fSSPChOne[ 0] = 0;
    fSSP[ 2] = 11; fSSPChOne[ 2] = 4;
    fSSP[ 4] = 11; fSSPChOne[ 4] = 8;
    fSSP[ 6] = 12; fSSPChOne[ 6] = 0;
    fSSP[ 8] = 12; fSSPChOne[ 8] = 4;
    fSSP[10] = 12; fSSPChOne[10] = 8;
    fSSP[12] = 13; fSSPChOne[12] = 0;
    fSSP[14] = 13; fSSPChOne[14] = 4;
    fSSP[16] = 13; fSSPChOne[16] = 8;
    fSSP[18] = 14; fSSPChOne[18] = 0;

    // MSRaS -- "normal" pattern
    fSSP[20] = 21; fSSPChOne[20] = 0;
    fSSP[22] = 21; fSSPChOne[22] = 4;
    fSSP[24] = 21; fSSPChOne[24] = 8;
    fSSP[26] = 22; fSSPChOne[26] = 0;
    fSSP[28] = 22; fSSPChOne[28] = 4;
    fSSP[30] = 22; fSSPChOne[30] = 8;
    fSSP[32] = 23; fSSPChOne[32] = 0;
    fSSP[34] = 23; fSSPChOne[34] = 4;
    fSSP[36] = 23; fSSPChOne[36] = 8;
    fSSP[38] = 24; fSSPChOne[38] = 0;
  
    // USRaS -- unusual pattern for ARAPUCA
    fSSP[40] = 31; fSSPChOne[40] = 0;
    fSSP[42] = 31; fSSPChOne[42] = 4;
    fSSP[44] = 31; fSSPChOne[44] = 8;
    fSSP[46] = 34; fSSPChOne[46] = 0; // ARAPUCA
    fSSP[48] = 32; fSSPChOne[48] = 0;
    fSSP[50] = 32; fSSPChOne[50] = 4;
    fSSP[52] = 32; fSSPChOne[52] = 8;
    fSSP[54] = 33; fSSPChOne[54] = 0;
    fSSP[56] = 33; fSSPChOne[56] = 4;
    fSSP[58] = 33; fSSPChOne[58] = 8;


    // The above enables OpDet + HW channel -> OpChannel
    //
    // Fill the maps below to do the reverse by looping through
    // all possible OpDet and HW Channel combinations

    fMaxOpChannel = 0;
    fNOpChannels = 0;
  
    for (Index opDet = 0; opDet < 60; opDet++) {
      for (Index hwCh = 0; hwCh < NOpHardwareChannels(opDet); hwCh++) {

        // Find the channel number for this opDet and hw channel
        Index opChannel = OpChannel(opDet, hwCh);

        // Count channels and record the maximum possible channel
        if (opChannel > fMaxOpChannel) fMaxOpChannel = opChannel;
        fNOpChannels++;

        // Fill maps for opChannel -> opDet and hwChannel
        fOpDet[opChannel] = opDet;
        fHWChannel[opChannel] = hwCh;
      }
    }

  }


//----------------------------------------------------------------------------

  void ProtoDUNEChannelMapAlg::Uninitialize() {
    // Call the base-class version to perform TPC-related actions
    DuneApaChannelMapAlg::Uninitialize();
  }

//----------------------------------------------------------------------------

  Index ProtoDUNEChannelMapAlg::NOpChannels(Index /*NOpDets*/) const {
    return fNOpChannels;
  }

//----------------------------------------------------------------------------

  Index ProtoDUNEChannelMapAlg::MaxOpChannel(Index /*NOpDets*/) const {
    return fMaxOpChannel;
  }

//----------------------------------------------------------------------------

  Index ProtoDUNEChannelMapAlg::NOpHardwareChannels(Index opDet) const {

    // ARAPUCAs
    if (opDet == 31 or opDet == 46)
      return 12;
    else
      return 4;  
  }

//----------------------------------------------------------------------------

  Index ProtoDUNEChannelMapAlg::OpChannel(Index detNum, Index channel) const {

    Index sspch = fSSPChOne.at(detNum) + channel;
    Index ssp   = fSSP.at(detNum);
    
    // Special handling of ARAPUCA in MSDaS which cross between SSP 62 and 63
    if (sspch > 12) {
        mf::LogError("ProtoDUNEChannelMapAlg") << "Invalid address: SSP #" << ssp << ", SSP channel" << sspch;
    }

    return OpChannelFromSSP(ssp, sspch);
  }

//----------------------------------------------------------------------------

  Index ProtoDUNEChannelMapAlg::OpDetFromOpChannel(Index opChannel) const {
    if (!IsValidOpChannel(opChannel, 60)) {
      mf::LogWarning("ProtoDUNEChannelMapAlg") << "Requesting an OpDet number for an uninstrumented channel, " << opChannel;
      return 99999;
    }
    return fOpDet.at(opChannel);
  }

//----------------------------------------------------------------------------

  Index ProtoDUNEChannelMapAlg::HardwareChannelFromOpChannel(Index opChannel) const {
    if (!IsValidOpChannel(opChannel, 60)) {
      mf::LogWarning("ProtoDUNEChannelMapAlg") << "Requesting an OpDet number for an uninstrumented channel, " << opChannel;
      return 99999;
    }
    return fHWChannel.at(opChannel);
  }


//----------------------------------------------------------------------------

  Index ProtoDUNEChannelMapAlg::OpChannelFromSSP(Index ssp, Index sspch) const {
    // Expects SSP #'s of the from NM where N is APA number and M is SSP within the APA
    // So, IP 504 -> AP # 54

    //         ( (  APA # - 1 )    )*4 + SSP per APA)*12 + SSP channel
    Index ch = ( (trunc(ssp/10) - 1)*4 + ssp%10 - 1 )*12 + sspch;
    return ch;

  }


//----------------------------------------------------------------------------


  void ProtoDUNEChannelMapAlg::PrintChannelMaps() const {


    cout << "---------------------------------------------------------------" << endl;
    cout << "---------------------------------------------------------------" << endl;
    cout << "---------------------------------------------------------------" << endl;
    cout << "---------------------------------------------------------------" << endl;
    cout << "---------------------------------------------------------------" << endl;
    cout << "---------------------------------------------------------------" << endl;
    cout << "---------------------------------------------------------------" << endl;
    cout << endl << endl;

    std::vector<Index> ssps = { 11, 12, 13, 14,
                                21, 22, 23, 24,
                                31, 32, 33, 34,
                                41, 42, 43, 44,
                                51, 52, 53, 54,
                                61, 62, 63, 64 };
    
    cout << endl << "By SSP" << endl;
    for (Index ssp : ssps) {
      for (Index sspch = 0; sspch < 12; sspch++) {
        cout << setw(2) << ssp << " " << setw(2) << sspch << ": " << setw(3) << OpChannelFromSSP(ssp, sspch) << endl;
      }
    }

    cout << endl << "Beam side" << endl;
    for (Index opDet = 1; opDet < 60; opDet += 2) {
      cout << setw(2) << opDet << ":";
      for (Index hwCh = 0; hwCh < NOpHardwareChannels(opDet); hwCh++) {
        cout << " " << setw(2) << OpChannel(opDet, hwCh);
      }
      cout << endl;
    }


    cout << endl << "Non-Beam side" << endl;
    for (Index opDet = 0; opDet < 60; opDet += 2) {
      cout << setw(2) << opDet << ":";
      for (Index hwCh = 0; hwCh < NOpHardwareChannels(opDet); hwCh++) {
        cout << " " << setw(2) << OpChannel(opDet, hwCh);
      }
      cout << endl;
    }


    cout << endl << "Online -> offline" << endl;
    for (Index opCh = 0; opCh < MaxOpChannel(60); opCh++) {
      cout << setw(3) << opCh << ": ";
      if ( IsValidOpChannel(opCh, 60) ) {
        cout << setw(2) << OpDetFromOpChannel(opCh) << ", "
             << setw(2) << HardwareChannelFromOpChannel(opCh) << endl;
      }
      else {
        cout << "empty channel" << endl;
      }
    }

    cout << endl << endl;
    cout << "---------------------------------------------------------------" << endl;
    cout << "---------------------------------------------------------------" << endl;
    cout << "---------------------------------------------------------------" << endl;
    cout << "---------------------------------------------------------------" << endl;
    cout << "---------------------------------------------------------------" << endl;
    cout << "---------------------------------------------------------------" << endl;
    cout << "---------------------------------------------------------------" << endl;

  }
}
