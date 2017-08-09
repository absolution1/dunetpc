// AcdWireReader_tool.cc

#include "AcdWireReader.h"
#include <iostream>
#include "lardataobj/RecoBase/Wire.h"

using std::string;
using std::cout;
using std::endl;
using raw::RawDigit;
using recob::Wire;

//**********************************************************************

AcdWireReader::AcdWireReader(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")) {
  const string myname = "AcdWireReader::ctor: ";
  if ( m_LogLevel > 0 ) {
    cout << myname << "    LogLevel: " << m_LogLevel << endl;
  }
}

//**********************************************************************

int AcdWireReader::update(AdcChannelData& acd) const {
  const string myname = "AcdWireReader::update: ";
  // Take the wire from the channel data.
  const Wire* pwir = acd.wire;
  if ( pwir == nullptr ) {
    cout << myname << "ERROR: Wire is null." << endl;
    return 1;
  }
  const Wire& wir = *pwir;
  // Check the input data is empty.
  if ( acd.samples.size() ) {
    cout << myname << "ERROR: ADC channel has prepared data." << endl;
    return 2;
  }
  // Check the signal flags record is empty.
  if ( acd.signal.size() ) {
    cout << myname << "ERROR: ADC channel has signal flags." << endl;
    return 3;
  }
  // Check the ROI data is empty.
  if ( acd.rois.size() ) {
    cout << myname << "ERROR: ADC channel has ROIs." << endl;
    return 4;
  }
  // Set or check the wire index.
  if ( acd.wireIndex != AdcChannelData::badIndex ) {
    cout << myname << "ERROR: ADC channel has a wire index." << endl;
    return 5;
  }
  // Set or check the channel.
  if ( acd.channel == AdcChannelData::badChannel ) {
    acd.channel = wir.Channel();
  } else {
    if ( acd.channel != wir.Channel() ) {
      cout << myname << "ERROR: Wire has inconsistent channel number." << endl;
      return 6;
    }
  }
  // Add code to copy wire ROIs to the channel data.
  return 99;
  return 0;
}

//**********************************************************************
