// ProtoduneOnlineChannel.cxx

#include "ProtoduneOnlineChannel.h"
#include <iostream>

using std::string;
using std::cout;
using std::endl;
using Index = ProtoduneOnlineChannel::Index;

//**********************************************************************

ProtoduneOnlineChannel::ProtoduneOnlineChannel(const fhicl::ParameterSet&) 
: m_LogLevel(2) {
    Index val = 0;
    // Create the maps from the to CCW wire numbers for each plane to
    // the FEMB channel number.
    // We follow DUNE DocDB 4064 except numberings start from 0 instead of 1.
    // ASIC 0
    uch[18] = val++;
    uch[16] = val++;
    uch[14] = val++;
    uch[12] = val++;
    uch[10] = val++;
    vch[18] = val++;
    vch[16] = val++;
    vch[14] = val++;
    vch[12] = val++;
    vch[10] = val++;
    zch[22] = val++;
    zch[20] = val++;
    zch[18] = val++;
    zch[16] = val++;
    zch[14] = val++;
    zch[12] = val++;
    // ASIC 1
    uch[ 8] = val++;
    uch[ 6] = val++;
    uch[ 4] = val++;
    uch[ 2] = val++;
    uch[ 0] = val++;
    vch[ 8] = val++;
    vch[ 6] = val++;
    vch[ 4] = val++;
    vch[ 2] = val++;
    vch[ 0] = val++;
    zch[10] = val++;
    zch[ 8] = val++;
    zch[ 6] = val++;
    zch[ 4] = val++;
    zch[ 2] = val++;
    zch[ 0] = val++;
    // ASIC 2
    zch[13] = val++;
    zch[15] = val++;
    zch[17] = val++;
    zch[19] = val++;
    zch[21] = val++;
    zch[23] = val++;
    vch[11] = val++;
    vch[13] = val++;
    vch[15] = val++;
    vch[17] = val++;
    vch[19] = val++;
    uch[11] = val++;
    uch[13] = val++;
    uch[15] = val++;
    uch[17] = val++;
    uch[19] = val++;
    // ASIC 3
    zch[ 1] = val++;
    zch[ 3] = val++;
    zch[ 5] = val++;
    zch[ 7] = val++;
    zch[ 9] = val++;
    zch[11] = val++;
    vch[ 1] = val++;
    vch[ 3] = val++;
    vch[ 5] = val++;
    vch[ 7] = val++;
    vch[ 9] = val++;
    uch[ 1] = val++;
    uch[ 3] = val++;
    uch[ 5] = val++;
    uch[ 7] = val++;
    uch[ 9] = val++;
    // ASIC 4
    uch[28] = val++;
    uch[26] = val++;
    uch[24] = val++;
    uch[22] = val++;
    uch[20] = val++;
    vch[28] = val++;
    vch[26] = val++;
    vch[24] = val++;
    vch[22] = val++;
    vch[20] = val++;
    zch[34] = val++;
    zch[32] = val++;
    zch[30] = val++;
    zch[28] = val++;
    zch[26] = val++;
    zch[24] = val++;
    // ASIC 5
    uch[38] = val++;
    uch[36] = val++;
    uch[34] = val++;
    uch[32] = val++;
    uch[30] = val++;
    vch[38] = val++;
    vch[36] = val++;
    vch[34] = val++;
    vch[32] = val++;
    vch[30] = val++;
    zch[46] = val++;
    zch[44] = val++;
    zch[42] = val++;
    zch[40] = val++;
    zch[38] = val++;
    zch[36] = val++;
    // ASIC 6
    zch[25] = val++;
    zch[27] = val++;
    zch[29] = val++;
    zch[31] = val++;
    zch[33] = val++;
    zch[35] = val++;
    vch[21] = val++;
    vch[23] = val++;
    vch[25] = val++;
    vch[27] = val++;
    vch[29] = val++;
    uch[21] = val++;
    uch[23] = val++;
    uch[25] = val++;
    uch[27] = val++;
    uch[29] = val++;
    // ASIC 7
    zch[37] = val++;
    zch[39] = val++;
    zch[41] = val++;
    zch[43] = val++;
    zch[45] = val++;
    zch[47] = val++;
    vch[31] = val++;
    vch[33] = val++;
    vch[35] = val++;
    vch[37] = val++;
    vch[39] = val++;
    uch[31] = val++;
    uch[33] = val++;
    uch[35] = val++;
    uch[37] = val++;
    uch[39] = val++;
    if ( val != 128 ) abort();
}

//**********************************************************************

Index ProtoduneOnlineChannel::get(Index chanOff) const {
  const string myname = "ProtoduneOnlineChannel::get: ";
  if ( chanOff >= 15360 ) {
    if ( m_LogLevel > 1 ) cout << myname << "Invalid offline channel: " << chanOff << endl;
    return badIndex();
  }
  // Get the APA.
  Index iapa = chanOff/2560;
  bool beamLeft = iapa & 1;
  bool beamRight = ! beamLeft;
  // Get the channel in the apa.
  Index ichApa = chanOff%2560;
  // Get the plane (ipla) and wire number in the plane (ichPla).
  Index ipla = 0;
  Index ichPla = ichApa;
  while ( ichPla >= nwirPlane[ipla] ) {
    ichPla -= nwirPlane[ipla];
    ++ipla;
  }
  // Get the FEMB # in the detector ifmbDet.
  Index ifmbApa = ichPla/nwirFemb[ipla];
  if ( ipla == 0 ) {
    if ( beamRight ) ifmbApa = (ifmbApa + 10) % 20;
  } else if ( ipla == 1 ) {
    ifmbApa = 19 - ifmbApa;
    if ( beamLeft ) ifmbApa = (ifmbApa + 10) % 20;
  } else if ( ipla == 2 ) {
    if ( beamLeft ) ifmbApa = 9 - ifmbApa;
    else ifmbApa = 19 - ifmbApa;
  } else {
    if ( beamLeft ) ifmbApa += 10;
  }
  Index ifmbDet = 20*iapa + ifmbApa;
  // Get the wire number in the FEMB.
  Index iwchFemb = ichPla % nwirFemb[ipla];  // Wire number in the plane and FEMB.
  if ( iwchFemb > nwirFemb[ipla] ) {
    if ( m_LogLevel > 1 ) cout << myname << "ERROR: Invalid FEMB channel: " << chanOff
                               << " --> " << iwchFemb << endl;
    return badIndex();
  }
  // Note the wire numbers in the channel-to-wire tables increase CCW while
  // offline is CW for u and z2 and CCW for v and z1.
  // Flip the wire numbers for the former.
  if ( ipla == 0 || ipla == 3 ) iwchFemb = nwirFemb[ipla] - 1 - iwchFemb;
  // Find the FEMB channel for the wire.
  Index ichFemb = ipla == 0 ? uch[iwchFemb] :
                  ipla == 1 ? vch[iwchFemb] :
                              zch[iwchFemb];
  // Build offline index.
  Index ichOn = 128*ifmbDet + ichFemb;
  return ichOn;
}

//**********************************************************************
