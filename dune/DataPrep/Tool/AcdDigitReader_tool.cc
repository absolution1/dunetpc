// AcdDigitReader_tool.cc

#include "AcdDigitReader.h"
#include <iostream>
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"   // compression fns are here

using std::string;
using std::cout;
using std::endl;
using raw::RawDigit;

//**********************************************************************

AcdDigitReader::AcdDigitReader(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")) {
  const string myname = "AcdDigitReader::ctor: ";
  if ( m_LogLevel > 0 ) {
    cout << myname << "    LogLevel: " << m_LogLevel << endl;
  }
}

//**********************************************************************

DataMap AcdDigitReader::update(AdcChannelData& acd) const {
  const string myname = "AcdDigitReader::update: ";
  // Take the digit from the channel data.
  const RawDigit* pdig = acd.digit;
  if ( pdig == nullptr ) {
    cout << myname << "ERROR: Digit is null." << endl;
    return DataMap(1);
  }
  const RawDigit& dig = *pdig;
  // Check the input data is empty.
  if ( acd.raw.size() ) {
    cout << myname << "ERROR: ADC channel has raw data." << endl;
    return DataMap(2);
  }
  if ( acd.flags.size() ) {
    cout << myname << "ERROR: ADC channel has flag data." << endl;
    return DataMap(3);
  }
  if ( acd.pedestal != AdcChannelData::badSignal ) {
    cout << myname << "ERROR: ADC channel has a pedestal." << endl;
    return DataMap(4);
  }
  // Set or check the channel number.
  if ( acd.channel == AdcChannelData::badChannel ) {
    acd.channel = dig.Channel();
  } else {
    if ( acd.channel != dig.Channel() ) {
      cout << myname << "ERROR: Raw digit has inconsistent channel number." << endl;
      return DataMap(5);
    }
  }
  // Copy pedestal.
  acd.pedestal = dig.GetPedestal();
  // Copy raw data.
  unsigned int nsig = dig.Samples();
  acd.raw.resize(nsig, -999);  // See https://cdcvs.fnal.gov/redmine/issues/11572.
  raw::Uncompress(dig.ADCs(), acd.raw, dig.GetPedestal(), dig.Compression());
  if ( m_LogLevel >= 3 ) cout << myname << "Channel " << acd.channel << " raw count: " << acd.raw.size() << endl;
  return DataMap(0);
}

//**********************************************************************
