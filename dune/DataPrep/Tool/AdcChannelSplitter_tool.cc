// AdcChannelSplitter_tool.cc

#include "AdcChannelSplitter.h"
#include <iostream>
#include <vector>

using std::string;
using std::cout;
using std::endl;

//**********************************************************************
// Class methods.
//**********************************************************************

AdcChannelSplitter::AdcChannelSplitter(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")), 
  m_Length(ps.get<Index>("Length")), 
  m_DataPath(ps.get<Name>("DataPath")),
  m_DataView(ps.get<Name>("DataView")) {
  const string myname = "AdcChannelSplitter::ctor: ";
  // Display the configuration.
  if ( m_LogLevel >= 1 ) {
    cout << myname << "Configuration: " << endl;
    cout << myname << "         LogLevel: " << m_LogLevel << endl;
    cout << myname << "           Length: " << m_Length << endl;
    cout << myname << "         DataPath: " << m_DataPath << endl;
    cout << myname << "         DataView: " << m_DataView << endl;
  }
}

//**********************************************************************

DataMap AdcChannelSplitter::update(AdcChannelData& acd) const {
  const string myname = "AdcChannelSplitter::update: ";
  DataMap ret;
  Index dtck = m_Length;
  if ( dtck == 0 ) {
    cout << "ERROR: Length is zero." << endl;
    return ret.setStatus(1);
  }
  Index nvie = acd.viewSize(m_DataPath);
  Index nobj = 0;
  if ( m_LogLevel >= 2 ) cout << myname << "Channel " << acd.channel()
                              << " input object count: " << nvie << endl;
  Index nrawCopied = 0;
  Index nsamCopied = 0;
  // Loop over input objects.
  for ( Index ivie=0; ivie<nvie; ++ivie ) {
    AdcChannelData* pacd = acd.mutableViewEntry(m_DataPath, ivie);
    Index nraw = pacd->raw.size();
    Index nsam = pacd->samples.size();
    // Create the new view. It might be left empty.
    AdcChannelData::View& acds = pacd->updateView(m_DataView);
    Index itck = 0;                   // Start copy at this tick
    while ( true ) {
      Index jtck = itck + dtck;         // End copy at this tick
      bool copyRaw = nraw >= jtck;
      bool copySam = nsam >= jtck;
      if ( !copyRaw && !copySam ) break;
      acds.push_back(*pacd);
      AdcChannelData& acdNew = acds.back();
      acdNew.viewParent = pacd;
      acdNew.tick0 = pacd->tick0 + itck;
      if ( copyRaw ) {
        for ( Index iraw=itck; iraw<jtck; ++iraw ) {
          acdNew.raw.push_back(pacd->raw[iraw]);
          ++nrawCopied;
        }
      }
      if ( copySam ) {
        for ( Index isam=itck; isam<jtck; ++isam ) {
          acdNew.samples.push_back(pacd->samples[isam]);
          ++nsamCopied;
        }
      }
      if ( acdNew.samples.size() ) acdNew.sampleUnit = pacd->sampleUnit;
      ++nobj;
      itck = jtck;
    }
  }
  if ( m_LogLevel >= 2 ) cout << myname << "Channel " << acd.channel()
                              << " output object count: " << nobj << endl;
  ret.setInt("splitInputCount", nvie);
  ret.setInt("splitOutputCount", nobj);
  ret.setInt("splitRawCopyCount", nrawCopied);
  ret.setInt("splitSampleCopyCount", nsamCopied);
  return ret;
}

//**********************************************************************

DEFINE_ART_CLASS_TOOL(AdcChannelSplitter)
