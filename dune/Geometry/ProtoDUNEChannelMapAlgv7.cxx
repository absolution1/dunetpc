////////////////////////////////////////////////////////////////////////
/// \file  ProtoDUNEChannelMapAlgv7.cxx
/// \brief Interface to algorithm class for a specific detector channel mapping
///
/// \version $Id:  $
/// \author  tylerdalion@gmail.com
////////////////////////////////////////////////////////////////////////

#include "dune/Geometry/ProtoDUNEChannelMapAlgv7.h"
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
using geo::ProtoDUNEChannelMapAlgv7;
using std::cout;
using std::endl;
using std::setw;

typedef unsigned int Index;

namespace geo {

//----------------------------------------------------------------------------

  ProtoDUNEChannelMapAlgv7::ProtoDUNEChannelMapAlgv7(const fhicl::ParameterSet& p)
    : DuneApaChannelMapAlg(p), 
      fMaxOpChannel(0), fNOpChannels(0)
  {
    fSorter = new GeoObjectSorterAPA(p);
  }


//----------------------------------------------------------------------------

  void ProtoDUNEChannelMapAlgv7::Initialize(GeometryData_t const& geodata) {
    Uninitialize();
  
    // Call the base-class version to perform TPC-related actions
    DuneApaChannelMapAlg::Initialize(geodata);

    // Manually entered based on maps from Chris Macias
    // http://indico.fnal.gov/event/14582/session/30/contribution/30/material/slides/0.pdf
    //
    // That gives SSP and channel ranges from OpDet. These offline channels correspond
    // to the APAs listed. This is for the geometry that includes the arapuca bar with its
    // 16 optical windows (protodune geometry v7).
    //
    //   USDaS(APA5)   MSDaS(APA6)  DSDaS(APA4)
    //   19-28         10-18+29-44   0-9
    //
    //   USRaS(APA3)   MSRaS(APA2)  DSRaS(APA1)
    //   45-60+81-89     71-80       61-70


    // DSDaS -- "normal" pattern
    fSSP[ 0] = 41; fSSPChOne[ 0] = 0;
    fSSP[ 1] = 41; fSSPChOne[ 1] = 4;
    fSSP[ 2] = 41; fSSPChOne[ 2] = 8;
    fSSP[ 3] = 42; fSSPChOne[ 3] = 0;
    fSSP[ 4] = 42; fSSPChOne[ 4] = 4;
    fSSP[ 5] = 42; fSSPChOne[ 5] = 8;
    fSSP[ 6] = 43; fSSPChOne[ 6] = 0;
    fSSP[ 7] = 43; fSSPChOne[ 7] = 4;
    fSSP[ 8] = 43; fSSPChOne[ 8] = 8;
    fSSP[ 9] = 44; fSSPChOne[ 9] = 0;

    // MSDaS -- unusual patternfor ARAPUCA
    fSSP[10] = 61; fSSPChOne[10] = 0;
    fSSP[11] = 61; fSSPChOne[11] = 4;
    fSSP[12] = 61; fSSPChOne[12] = 8;
    fSSP[13] = 62; fSSPChOne[13] = 0;
    fSSP[14] = 62; fSSPChOne[14] = 4;
    fSSP[29] = 63; fSSPChOne[29] = 0; // ARAPUCA
    fSSP[30] = 63; fSSPChOne[30] = 0; // ARAPUCA
    fSSP[31] = 63; fSSPChOne[31] = 1; // ARAPUCA
    fSSP[32] = 63; fSSPChOne[32] = 1; // ARAPUCA
    fSSP[33] = 63; fSSPChOne[33] = 2; // ARAPUCA
    fSSP[34] = 63; fSSPChOne[34] = 2; // ARAPUCA
    fSSP[35] = 63; fSSPChOne[35] = 3; // ARAPUCA
    fSSP[36] = 63; fSSPChOne[36] = 3; // ARAPUCA
    fSSP[37] = 63; fSSPChOne[37] = 4; // ARAPUCA
    fSSP[38] = 63; fSSPChOne[38] = 5; // ARAPUCA
    fSSP[39] = 63; fSSPChOne[39] = 6; // ARAPUCA
    fSSP[40] = 63; fSSPChOne[40] = 7; // ARAPUCA
    fSSP[41] = 63; fSSPChOne[41] = 8; // ARAPUCA
    fSSP[42] = 63; fSSPChOne[42] = 9; // ARAPUCA
    fSSP[43] = 63; fSSPChOne[43] = 10; // ARAPUCA
    fSSP[44] = 63; fSSPChOne[44] = 11; // ARAPUCA
    fSSP[15] = 62; fSSPChOne[15] = 8;
    fSSP[16] = 64; fSSPChOne[16] = 0;
    fSSP[17] = 64; fSSPChOne[17] = 4;
    fSSP[18] = 64; fSSPChOne[18] = 8;

    // USDaS -- unusual pattern
    fSSP[19] = 53; fSSPChOne[19] = 0;
    fSSP[20] = 53; fSSPChOne[20] = 4;
    fSSP[21] = 53; fSSPChOne[21] = 8;
    fSSP[22] = 54; fSSPChOne[22] = 0;
    fSSP[23] = 51; fSSPChOne[23] = 0;
    fSSP[24] = 54; fSSPChOne[24] = 4;
    fSSP[25] = 51; fSSPChOne[25] = 4;
    fSSP[26] = 51; fSSPChOne[26] = 8;
    fSSP[27] = 54; fSSPChOne[27] = 8;
    fSSP[28] = 52; fSSPChOne[28] = 0;

    // DSRaS -- "normal" pattern
    fSSP[61] = 11; fSSPChOne[61] = 0;
    fSSP[62] = 11; fSSPChOne[62] = 4;
    fSSP[63] = 11; fSSPChOne[63] = 8;
    fSSP[64] = 12; fSSPChOne[64] = 0;
    fSSP[65] = 12; fSSPChOne[65] = 4;
    fSSP[66] = 12; fSSPChOne[66] = 8;
    fSSP[67] = 13; fSSPChOne[67] = 0;
    fSSP[68] = 13; fSSPChOne[68] = 4;
    fSSP[69] = 13; fSSPChOne[69] = 8;
    fSSP[70] = 14; fSSPChOne[70] = 0;

    // MSRaS -- "normal" pattern
    fSSP[71] = 21; fSSPChOne[71] = 0;
    fSSP[72] = 21; fSSPChOne[72] = 4;
    fSSP[73] = 21; fSSPChOne[73] = 8;
    fSSP[74] = 22; fSSPChOne[74] = 0;
    fSSP[75] = 22; fSSPChOne[75] = 4;
    fSSP[76] = 22; fSSPChOne[76] = 8;
    fSSP[77] = 23; fSSPChOne[77] = 0;
    fSSP[78] = 23; fSSPChOne[78] = 4;
    fSSP[79] = 23; fSSPChOne[79] = 8;
    fSSP[80] = 24; fSSPChOne[80] = 0;

    // USRaS -- unusual patternfor ARAPUCA
    fSSP[81] = 31; fSSPChOne[81] = 0;
    fSSP[82] = 31; fSSPChOne[82] = 4;
    fSSP[83] = 31; fSSPChOne[83] = 8;
    fSSP[45] = 34; fSSPChOne[45] = 0; // ARAPUCA
    fSSP[46] = 34; fSSPChOne[46] = 0; // ARAPUCA
    fSSP[47] = 34; fSSPChOne[47] = 1; // ARAPUCA
    fSSP[48] = 34; fSSPChOne[48] = 1; // ARAPUCA
    fSSP[49] = 34; fSSPChOne[49] = 2; // ARAPUCA
    fSSP[50] = 34; fSSPChOne[50] = 2; // ARAPUCA
    fSSP[51] = 34; fSSPChOne[51] = 3; // ARAPUCA
    fSSP[52] = 34; fSSPChOne[52] = 3; // ARAPUCA
    fSSP[53] = 34; fSSPChOne[53] = 4; // ARAPUCA
    fSSP[54] = 34; fSSPChOne[54] = 5; // ARAPUCA
    fSSP[55] = 34; fSSPChOne[55] = 6; // ARAPUCA
    fSSP[56] = 34; fSSPChOne[56] = 7; // ARAPUCA
    fSSP[57] = 34; fSSPChOne[57] = 8; // ARAPUCA
    fSSP[58] = 34; fSSPChOne[58] = 9; // ARAPUCA
    fSSP[59] = 34; fSSPChOne[59] = 10; // ARAPUCA
    fSSP[60] = 34; fSSPChOne[60] = 11; // ARAPUCA
    fSSP[84] = 32; fSSPChOne[84] = 0;
    fSSP[85] = 32; fSSPChOne[85] = 4;
    fSSP[86] = 32; fSSPChOne[86] = 8;
    fSSP[87] = 33; fSSPChOne[87] = 0;
    fSSP[88] = 33; fSSPChOne[88] = 4;
    fSSP[89] = 33; fSSPChOne[89] = 8;

    // The above enables OpDet + HW channel -> OpChannel
    //
    // Fill the maps below to do the reverse by looping through
    // all possible OpDet and HW Channel combinations

    fMaxOpChannel = 0;
    fNOpChannels = 0;
  
    for (Index opDet = 0; opDet < 90; opDet++) {
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

//    PrintChannelMaps();
  }


//----------------------------------------------------------------------------

  void ProtoDUNEChannelMapAlgv7::Uninitialize() {
    // Call the base-class version to perform TPC-related actions
    DuneApaChannelMapAlg::Uninitialize();
  }

//----------------------------------------------------------------------------

  Index ProtoDUNEChannelMapAlgv7::NOpChannels(Index /*NOpDets*/) const {
    return fNOpChannels;
  }

//----------------------------------------------------------------------------

  Index ProtoDUNEChannelMapAlgv7::MaxOpChannel(Index /*NOpDets*/) const {
    return fMaxOpChannel;
  }

//----------------------------------------------------------------------------

  Index ProtoDUNEChannelMapAlgv7::NOpHardwareChannels(Index opDet) const {

    // ARAPUCAs
    if ((opDet > 28 && opDet < 61))
      return 1;
    else
      return 4;  
  }

//----------------------------------------------------------------------------

  Index ProtoDUNEChannelMapAlgv7::OpChannel(Index detNum, Index channel) const {

    Index sspch = fSSPChOne.at(detNum) + channel;
    Index ssp   = fSSP.at(detNum);
    
    // Special handling of ARAPUCA in MSDaS which cross between SSP 62 and 63
    if (sspch > 12) {
        mf::LogError("ProtoDUNEChannelMapAlgv7") << "Invalid address: SSP #" << ssp << ", SSP channel" << sspch;
    }

    return OpChannelFromSSP(ssp, sspch);
  }

//----------------------------------------------------------------------------

  Index ProtoDUNEChannelMapAlgv7::OpDetFromOpChannel(Index opChannel) const {
    if (!IsValidOpChannel(opChannel, 90)) {
      mf::LogWarning("ProtoDUNEChannelMapAlgv7") << "Requesting an OpDet number for an uninstrumented channel, " << opChannel;
      return 99999;
    }
    return fOpDet.at(opChannel);
  }

//----------------------------------------------------------------------------

  Index ProtoDUNEChannelMapAlgv7::HardwareChannelFromOpChannel(Index opChannel) const {
    if (!IsValidOpChannel(opChannel, 90)) {
      mf::LogWarning("ProtoDUNEChannelMapAlgv7") << "Requesting an OpDet number for an uninstrumented channel, " << opChannel;
      return 99999;
    }
    return fHWChannel.at(opChannel);
  }


//----------------------------------------------------------------------------

  Index ProtoDUNEChannelMapAlgv7::OpChannelFromSSP(Index ssp, Index sspch) const {
    // Expects SSP #'s of the from NM where N is APA number and M is SSP within the APA
    // So, IP 504 -> AP # 54

    //         ( (  APA # - 1 )    )*4 + SSP per APA)*12 + SSP channel
    Index ch = ( (trunc(ssp/10) - 1)*4 + ssp%10 - 1 )*12 + sspch;
    return ch;

  }


//----------------------------------------------------------------------------


  void ProtoDUNEChannelMapAlgv7::PrintChannelMaps() const {


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

    cout << endl << "Non-Beam side" << endl;
    for (Index opDet = 0; opDet < 45; opDet ++) {
      cout << setw(2) << opDet << ":";
      for (Index hwCh = 0; hwCh < NOpHardwareChannels(opDet); hwCh++) {
        cout << " " << setw(2) << OpChannel(opDet, hwCh);
      }
      cout << endl;
    }


    cout << endl << "Beam side" << endl;
    for (Index opDet = 45; opDet < 90; opDet ++) {
      cout << setw(2) << opDet << ":";
      for (Index hwCh = 0; hwCh < NOpHardwareChannels(opDet); hwCh++) {
        cout << " " << setw(2) << OpChannel(opDet, hwCh);
      }
      cout << endl;
    }


    cout << endl << "Online -> offline" << endl;
    for (Index opCh = 0; opCh < MaxOpChannel(90); opCh++) {
      cout << setw(3) << opCh << ": ";
      if ( IsValidOpChannel(opCh, 90) ) {
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
