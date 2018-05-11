// AdcChannelDumper_tool.cc

#include "AdcChannelDumper.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

using std::string;
using std::cout;
using std::endl;
using std::ofstream;
using std::ostream;
using std::ostringstream;
using std::setw;
using std::fixed;
using std::setprecision;

using Index = unsigned int;

//**********************************************************************
// Class methods.
//**********************************************************************

AdcChannelDumper::AdcChannelDumper(fhicl::ParameterSet const& ps)
: m_FileName(ps.get<string>("FileName")),
  m_Prefix(ps.get<string>("Prefix")),
  m_NewFile(ps.get<bool>("NewFile")),
  m_MaxSample(ps.get<int>("MaxSample")),
  m_pout(nullptr) {
  if ( m_FileName.size() == 0 ) m_pout = &cout;
  else if ( ! m_NewFile ) m_pout = new ofstream(m_FileName.c_str());
}

//**********************************************************************

AdcChannelDumper::~AdcChannelDumper() {
  if ( m_pout != nullptr && m_pout != &cout ) {
    delete m_pout;
    m_pout = nullptr;
  }
}

//**********************************************************************

DataMap AdcChannelDumper::view(const AdcChannelData& acd) const {
  DataMap res;
  ostream* pout = m_pout;
  bool newfile = pout == nullptr;
  // If file is not already set, build the file name and open it.
  if ( newfile ) {
    if ( ! m_NewFile ) return res.setStatus(1);
    string fname = m_FileName;
    pout = new ofstream(fname.c_str());
  }
  if ( pout == nullptr ) return res.setStatus(2);
  ostream& out = *pout;
  string pre = m_Prefix;
  string sbad = "<Unknown>";
  out << pre << "     Run: ";
  if ( acd.run == AdcChannelData::badIndex ) out << sbad;
  else out << acd.run;
  out << "-";
  if ( acd.subRun == AdcChannelData::badIndex ) out << sbad;
  else out << acd.subRun;
  out << endl;
  out << pre << "   Event: ";
  if ( acd.event == AdcChannelData::badIndex ) out << sbad;
  else out << acd.event;
  out << endl;
  out << pre << " Channel: ";
  if ( acd.channel == AdcChannelData::badChannel ) out << sbad;
  else out << acd.channel;
  out << endl;
  out << pre << " Pedestal: ";
  if ( acd.pedestal == AdcChannelData::badSignal ) out << sbad;
  else out << acd.pedestal;
  out << endl;
  Index nraw = acd.raw.size();
  Index nflg = acd.flags.size();
  Index nprp = acd.samples.size();
  Index nsig = acd.signal.size();
  Index nkeep = 0;
  for ( bool keep : acd.signal ) if ( keep ) ++nkeep;
  Index nsam = nraw;
  if ( nprp > nsam ) nsam = nprp;
  if ( m_MaxSample >= 0 ) {
    Index maxsam = m_MaxSample;
    if ( maxsam < nsam ) nsam = maxsam;
  }
  Index widx = 2 + log10(nsam);
  if ( widx < 6 ) widx = 6;
  Index wraw = 6;
  Index wflg = 4;
  Index wprp = 9;
  Index wsig = 2;
  Index wcnt = wprp;
  out << pre << "     Nraw: " << setw(wcnt) << nraw << endl;
  out << pre << "    Nprep: " << setw(wcnt) << nprp << endl;
  out << pre << "     Nsig: " << setw(wcnt) << acd.signal.size() << endl;
  out << pre << "    Nkeep: " << setw(wcnt) << nkeep << endl;
  out << pre << "     Nroi: " << setw(wcnt) << acd.rois.size() << endl;
  out << pre << "    Ndftm: " << setw(wcnt) << acd.dftmags.size() << endl;
  out << pre << "    Ndftp: " << setw(wcnt) << acd.dftphases.size() << endl;
  if ( nsam > 0 ) {
    out << pre
        << setw(widx) << "Data:"
        << setw(wraw) << "Raw"
        << setw(wflg) << "Flg"
        << setw(wprp) << "Prepared"
        << setw(wsig) << "S" << endl;
  }
  for ( Index isam=0; isam<nsam; ++isam ) {
    ostringstream ssline;
    ssline << setw(widx) << isam;
    if ( isam < nraw ) ssline << setw(wraw) << acd.raw[isam];
    else ssline << setw(wraw) << "";
    if ( isam < nflg ) ssline << setw(wflg) << acd.flags[isam];
    else ssline << setw(wflg) << "";
    if ( isam < nprp ) ssline << fixed << setprecision(2) << setw(wprp) << acd.samples[isam];
    else ssline << setw(wprp) << "";
    if ( isam < nsig ) ssline << setw(wsig) << acd.signal[isam];
    else ssline << setw(wsig) << "";
    out << pre << ssline.str() << endl;
  }
  if ( newfile ) delete pout;
  return res;
}

//**********************************************************************
