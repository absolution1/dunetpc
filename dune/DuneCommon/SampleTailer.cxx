// SampleTailer.cxx

#include "SampleTailer.h"
#include <cmath>
#include "TH1F.h"

//**********************************************************************

SampleTailer::SampleTailer(float a_decayTime, float a_alpha)
: m_decayTime(a_decayTime),
  m_beta(exp(-1.0/m_decayTime)),
  m_alpha(a_alpha),
  m_pedestal(0.0),
  m_tail0(0.0) { }

//**********************************************************************

SampleTailer::SampleTailer(float a_decayTime)
: SampleTailer(a_decayTime, 0) {
  m_alpha = m_beta > 0.0 ? 1.0 - 1.0/m_beta : 0.0;
}

//**********************************************************************

SampleTailer::~SampleTailer() {
  clear();
}

//**********************************************************************

void SampleTailer::clear() {
  m_d.clear();
  m_s.clear();
  m_t.clear();
  delete m_phd;
  m_phd = nullptr;
  delete m_phs;
  m_phs = nullptr;
  delete m_pht;
  m_pht = nullptr;
}

//**********************************************************************

int SampleTailer::setDecayTime(float val, bool cancelSignal) {
  clear();
  if ( val > 0.0 ) {
    m_decayTime = val;
    m_beta = exp(-m_decayTime);
    if ( cancelSignal ) m_alpha = 1.0/m_beta - 1.0;
  } else {
    m_decayTime = 0.0;
    m_beta = 0.0;
    if ( cancelSignal ) m_alpha = 0.0;
    return 1;
  }
  return 0;
}

//**********************************************************************

int SampleTailer::setBeta(float val, bool cancelSignal) {
  clear();
  if ( val > 0.0 ) {
    m_beta = val;
    m_decayTime = -log(m_beta);
    if ( cancelSignal ) m_alpha = 1.0/m_beta - 1.0;
  } else {
    m_decayTime = 0.0;
    m_beta = 0.0;
    if ( cancelSignal ) m_alpha = 0.0;
    return 1;
  }
  return 0;
}

//**********************************************************************

int SampleTailer::setAlpha(float val) {
  clear();
  m_alpha = val;
  return 0;
}

//**********************************************************************

int SampleTailer::setPedestal(float val) {
  clear();
  m_pedestal = val;
  return 0;
}

//**********************************************************************

int SampleTailer::setTail0(float val) {
  clear();
  m_tail0 = val;
  return 0;
}

//**********************************************************************

int SampleTailer::setData(const FloatVector& inData) {
  clear();
  if ( ! isValid() ) return 1;
  m_d = inData;
  Index nsam = m_d.size();
  m_s.resize(nsam, 0.0);
  m_t.resize(nsam, 0.0);
  double t = m_tail0;
  double s = 0.0;
  for ( Index isam=0; isam<nsam; ++isam ) {
    if ( isam > 0 ) t = m_beta*t + m_alpha*s;
    s = m_d[isam] - t - m_pedestal;
    m_t[isam] = t;
    m_s[isam] = s;
  }
  return 0;
}

//**********************************************************************

int SampleTailer::setDataZero(Index nsam) {
  FloatVector dats(nsam, 0.0);
  return setData(dats);
}

//**********************************************************************

int SampleTailer::setSignal(const FloatVector& inSignal) {
  clear();
  if ( ! isValid() ) return 1;
  m_s = inSignal;
  Index nsam = m_s.size();
  m_d.resize(nsam, 0.0);
  m_t.resize(nsam, 0.0);
  double t = m_tail0;
  double s = 0.0;
  for ( Index isam=0; isam<nsam; ++isam ) {
    if ( isam > 0 ) t = m_beta*t + m_alpha*s;
    s = m_s[isam];
    m_t[isam] = t;
    m_d[isam] = s + t + m_pedestal;
  }
  return 0;
}

//**********************************************************************

TH1* SampleTailer::dataHist() {
  if ( m_phd != nullptr ) return m_phd;
  Index nsam = m_d.size();
  if ( nsam == 0 ) return nullptr;
  Name hnam = "hst_data";
  Name httl = "Sample data; Tick;" + unit();
  TH1* ph = new TH1F(hnam.c_str(), httl.c_str(), nsam, 0, nsam);
  ph->SetDirectory(nullptr);
  ph->SetLineWidth(2);
  ph->SetStats(0);
  for ( Index isam=0; isam<nsam; ++isam ) {
    ph->SetBinContent(isam+1, m_d[isam]);
  }
  m_phd = ph;
  return ph;
}

//**********************************************************************

TH1* SampleTailer::signalHist() {
  if ( m_phs != nullptr ) return m_phs;
  Index nsam = m_s.size();
  if ( nsam == 0 ) return nullptr;
  Name hnam = "hst_signal";
  Name httl = "Sample signal; Tick;" + unit();
  TH1* ph = new TH1F(hnam.c_str(), httl.c_str(), nsam, 0, nsam);
  ph->SetDirectory(nullptr);
  ph->SetLineWidth(2);
  ph->SetStats(0);
  for ( Index isam=0; isam<nsam; ++isam ) {
    ph->SetBinContent(isam+1, m_s[isam]);
  }
  m_phs = ph;
  return ph;
}

//**********************************************************************

TH1* SampleTailer::tailHist() {
  if ( m_pht != nullptr ) return m_pht;
  Index nsam = m_t.size();
  if ( nsam == 0 ) return nullptr;
  Name hnam = "hst_tail";
  Name httl = "Sample tail; Tick;" + unit();
  TH1* ph = new TH1F(hnam.c_str(), httl.c_str(), nsam, 0, nsam);
  ph->SetDirectory(nullptr);
  ph->SetLineWidth(2);
  ph->SetStats(0);
  for ( Index isam=0; isam<nsam; ++isam ) {
    ph->SetBinContent(isam+1, m_t[isam]);
  }
  m_pht = ph;
  return ph;
}

//**********************************************************************

void SampleTailer::makeHists() {
  dataHist();
  tailHist();
  signalHist();
}

//**********************************************************************

