// ConfigurableAdcChannelDataCopyService_service.cc

#include "ConfigurableAdcChannelDataCopyService.h"
#include <iostream>
#include "art/Framework/Services/Registry/ServiceHandle.h"

using std::string;
using std::cout;
using std::endl;

//**********************************************************************

ConfigurableAdcChannelDataCopyService::
ConfigurableAdcChannelDataCopyService(fhicl::ParameterSet const& pset, art::ActivityRegistry&) {
  const string myname = "ConfigurableAdcChannelDataCopyService::ctor: ";
  pset.get_if_present<int>("LogLevel", m_LogLevel);
  pset.get_if_present<bool>("CopyChannel", m_CopyChannel);
  pset.get_if_present<bool>("CopyPedestal", m_CopyPedestal);
  pset.get_if_present<bool>("CopyRaw", m_CopyRaw);
  pset.get_if_present<bool>("CopySamples", m_CopySamples);
  pset.get_if_present<bool>("CopyFlags", m_CopyFlags);
  pset.get_if_present<bool>("CopySignal", m_CopySignal);
  pset.get_if_present<bool>("CopyRois", m_CopyRois);
  pset.get_if_present<bool>("CopyDigit", m_CopyDigit);
  pset.get_if_present<bool>("CopyWire", m_CopyWire);
  pset.get_if_present<bool>("CopyDigitIndex", m_CopyDigitIndex);
  pset.get_if_present<bool>("CopyWireIndex", m_CopyWireIndex);
  print(cout, myname);
}

//**********************************************************************

int ConfigurableAdcChannelDataCopyService::
copy(const AdcChannelData& acdin, AdcChannelData& acdout) const {
  const string myname = "ConfigurableAdcChannelDataCopyService:copy: ";
  if ( m_LogLevel >= 2 ) cout << myname << "Copying channel " << acdin.channel << endl;
  if ( m_CopyChannel )    acdout.channel    = acdin.channel;
  if ( m_CopyPedestal )   acdout.pedestal   = acdin.pedestal;
  if ( m_CopyRaw )        acdout.raw        = acdin.raw;
  if ( m_CopySamples )    acdout.samples    = acdin.samples;
  if ( m_CopyFlags )      acdout.flags      = acdin.flags;
  if ( m_CopySignal )     acdout.signal     = acdin.signal;
  if ( m_CopyRois )       acdout.rois       = acdin.rois;
  if ( m_CopyDigit )      acdout.digit      = acdin.digit;
  if ( m_CopyWire )       acdout.wire       = acdin.wire;
  if ( m_CopyDigitIndex ) acdout.digitIndex = acdin.digitIndex;
  if ( m_CopyWireIndex )  acdout.wireIndex  = acdin.wireIndex;
  return 0;
}

//**********************************************************************

std::ostream& ConfigurableAdcChannelDataCopyService::
print(std::ostream& out, std::string prefix) const {
  out << prefix << "ConfigurableAdcChannelDataCopyService:" << endl;
  out << prefix << "        LogLevel: " << m_LogLevel       << endl;
  out << prefix << "     CopyChannel: " << m_CopyChannel    << endl;
  out << prefix << "    CopyPedestal: " << m_CopyPedestal   << endl;
  out << prefix << "         CopyRaw: " << m_CopyRaw        << endl;
  out << prefix << "     CopySamples: " << m_CopySamples    << endl;
  out << prefix << "       CopyFlags: " << m_CopyFlags      << endl;
  out << prefix << "      CopySignal: " << m_CopySignal     << endl;
  out << prefix << "        CopyRois: " << m_CopyRois       << endl;
  out << prefix << "       CopyDigit: " << m_CopyDigit      << endl;
  out << prefix << "        CopyWire: " << m_CopyWire       << endl;
  out << prefix << "  CopyDigitIndex: " << m_CopyDigitIndex << endl;
  out << prefix << "   CopyWireIndex: " << m_CopyWireIndex  << endl;
  return out;
}

//**********************************************************************

DEFINE_ART_SERVICE_INTERFACE_IMPL(ConfigurableAdcChannelDataCopyService, AdcChannelDataCopyService)

//**********************************************************************
