// RootParFormula.h

#include "dune/DuneCommon/Utility/RootParFormula.h"

using Names = ParFormula::Names;

//**********************************************************************

RootParFormula::RootParFormula(Name snam, Name sform) :
m_ptf(new TFormula(snam.c_str(), sform.c_str())),
m_sform(sform),
m_defval(0.0) {
  for ( int ipar=0; ipar<m_ptf->GetNpar(); ++ipar ) {
    m_parNames.push_back(m_ptf->GetParName(ipar));
  }
  m_setCounts.resize(npar(), 0);
  m_parValues.resize(npar(), 0.0);
}

//**********************************************************************

bool RootParFormula::isPar(Name parnam) const {
  return std::find(m_parNames.begin(), m_parNames.end(), parnam) != m_parNames.end();
}

//**********************************************************************

Names RootParFormula::setPars() const {
  Names nams;
  for ( Index ipar=0; ipar<npar(); ++ipar ) {
    if ( m_setCounts[ipar] ) nams.push_back(m_parNames[ipar]);
  }
  return nams;
}
  
//**********************************************************************

Names RootParFormula::unsetPars() const {
  Names nams;
  for ( Index ipar=0; ipar<npar(); ++ipar ) {
    if ( m_setCounts[ipar] == 0 ) nams.push_back(m_parNames[ipar]);
  }
  return nams;
}
  
//**********************************************************************

Names RootParFormula::resetPars() const {
  Names nams;
  for ( Index ipar=0; ipar<npar(); ++ipar ) {
    if ( m_setCounts[ipar] > 1 ) nams.push_back(m_parNames[ipar]);
  }
  return nams;
}
  
//**********************************************************************

int RootParFormula::setParValue(Name parnam, Value parval) {
  Index ipar = 0;
  for ( ; ipar<npar(); ++ipar ) {
    if ( m_parNames[ipar] == parnam ) break;
  }
  if ( ipar >= npar() ) return 1;
  m_parValues[ipar] = parval;
  ++m_setCounts[ipar];
  m_ptf->SetParameter(parnam.c_str(), parval);
  return 0;
}

//**********************************************************************

int RootParFormula::unsetParValues() {
  for ( Index ipar=0; ipar<npar(); ++ipar ) {
    m_setCounts[ipar] = 0;
  }
  return 0;
}

//**********************************************************************
