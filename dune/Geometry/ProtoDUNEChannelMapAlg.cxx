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
    fSSP[ 0] = 41; fSSPChOne[ 0] = 0;
    fSSP[ 2] = 41; fSSPChOne[ 2] = 4;
    fSSP[ 4] = 41; fSSPChOne[ 4] = 8;
    fSSP[ 6] = 42; fSSPChOne[ 6] = 0;
    fSSP[ 8] = 42; fSSPChOne[ 8] = 4;
    fSSP[10] = 42; fSSPChOne[10] = 8;
    fSSP[12] = 43; fSSPChOne[12] = 0;
    fSSP[14] = 43; fSSPChOne[14] = 4;
    fSSP[16] = 43; fSSPChOne[16] = 8;
    fSSP[18] = 44; fSSPChOne[18] = 0;

    // MSDaS -- unusual patternfor ARAPUCA
    fSSP[20] = 61; fSSPChOne[20] = 0;
    fSSP[22] = 61; fSSPChOne[22] = 4;
    fSSP[24] = 61; fSSPChOne[24] = 8;
    fSSP[26] = 62; fSSPChOne[26] = 0;
    fSSP[28] = 62; fSSPChOne[28] = 4;
    fSSP[30] = 63; fSSPChOne[30] = 0; // ARAPUCA
    fSSP[32] = 62; fSSPChOne[32] = 8;
    fSSP[34] = 64; fSSPChOne[34] = 0;
    fSSP[36] = 64; fSSPChOne[36] = 4;
    fSSP[38] = 64; fSSPChOne[38] = 8;

    // USDaS -- unusual pattern
    fSSP[40] = 53; fSSPChOne[40] = 0;
    fSSP[42] = 53; fSSPChOne[42] = 4;
    fSSP[44] = 53; fSSPChOne[44] = 8;
    fSSP[46] = 54; fSSPChOne[46] = 0;
    fSSP[48] = 51; fSSPChOne[48] = 0;
    fSSP[50] = 54; fSSPChOne[50] = 4;
    fSSP[52] = 51; fSSPChOne[52] = 4;
    fSSP[54] = 51; fSSPChOne[54] = 8;
    fSSP[56] = 54; fSSPChOne[56] = 8;
    fSSP[58] = 52; fSSPChOne[58] = 0;

    // DSRaS -- "normal" pattern
    fSSP[ 1] = 11; fSSPChOne[ 1] = 0;
    fSSP[ 3] = 11; fSSPChOne[ 3] = 4;
    fSSP[ 5] = 11; fSSPChOne[ 5] = 8;
    fSSP[ 7] = 12; fSSPChOne[ 7] = 0;
    fSSP[ 9] = 12; fSSPChOne[ 9] = 4;
    fSSP[11] = 12; fSSPChOne[11] = 8;
    fSSP[13] = 13; fSSPChOne[13] = 0;
    fSSP[15] = 13; fSSPChOne[15] = 4;
    fSSP[17] = 13; fSSPChOne[17] = 8;
    fSSP[19] = 14; fSSPChOne[19] = 0;

    // MSRaS -- "normal" pattern
    fSSP[21] = 21; fSSPChOne[21] = 0;
    fSSP[23] = 21; fSSPChOne[23] = 4;
    fSSP[25] = 21; fSSPChOne[25] = 8;
    fSSP[27] = 22; fSSPChOne[27] = 0;
    fSSP[29] = 22; fSSPChOne[29] = 4;
    fSSP[31] = 22; fSSPChOne[31] = 8;
    fSSP[33] = 23; fSSPChOne[33] = 0;
    fSSP[35] = 23; fSSPChOne[35] = 4;
    fSSP[37] = 23; fSSPChOne[37] = 8;
    fSSP[39] = 24; fSSPChOne[39] = 0;

    // USRaS -- unusual patternfor ARAPUCA
    fSSP[41] = 31; fSSPChOne[41] = 0;
    fSSP[43] = 31; fSSPChOne[43] = 4;
    fSSP[45] = 31; fSSPChOne[45] = 8;
    fSSP[47] = 34; fSSPChOne[47] = 0; // ARAPUCA
    fSSP[49] = 32; fSSPChOne[49] = 0;
    fSSP[51] = 32; fSSPChOne[51] = 4;
    fSSP[53] = 32; fSSPChOne[53] = 8;
    fSSP[55] = 33; fSSPChOne[55] = 0;
    fSSP[57] = 33; fSSPChOne[57] = 4;
    fSSP[59] = 33; fSSPChOne[59] = 8;

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
    if (opDet == 30 or opDet == 47)
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
