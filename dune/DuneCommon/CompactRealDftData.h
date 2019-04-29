// CompactRealDftData.h
//
// David Adams
// April 2019
//
// Concrete class that holds a compact representation of DFT data.

#ifndef CompactRealDftData_H
#define CompactRealDftData_H

#include "dune/DuneCommon/RealDftData.h"

//**********************************************************************

template<typename F>
class CompactRealDftData : public RealDftData<F> {

public:

  using Index = RealDftNormalization::Index;
  using FloatVector = std::vector<F>;
  using GlobalNormalization = RealDftNormalization::GlobalNormalization;
  using TermNormalization = RealDftNormalization::TermNormalization;

  // Default ctor.
  CompactRealDftData() =default;

  // Ctor from normalization.
  CompactRealDftData(GlobalNormalization gnorm, TermNormalization tnorm)
  : m_gnorm(gnorm), m_tnorm(tnorm) { }

  // Ctor from normalization and data.
  // Must have amps.size() - phas.size() = 0 or 1 for valid DFT.
  // If not, the object is left empty.
  CompactRealDftData(GlobalNormalization gnorm, TermNormalization tnorm,
                     const FloatVector& amps, const FloatVector& phas)
  : m_gnorm(gnorm), m_tnorm(tnorm), m_amps(amps), m_phas(phas) {
    Index namp = amps.size();
    Index npha = phas.size();
    if ( namp < npha || namp > npha + 1 ) clear();
  }
  
  // Ctor from data with default normalization.
  // Must have amps.size() - phas.size() = 0 or 1 for valid DFT.
  // If not, the object is left empty.
  CompactRealDftData(const FloatVector& amps, const FloatVector& phas)
  : m_amps(amps), m_phas(phas) {
    Index namp = amps.size();
    Index npha = phas.size();
    if ( namp < npha || namp > npha + 1 ) clear();
  }
  
  // Copy ctor. Changes representation and normalization if needed.
  CompactRealDftData(const RealDftData<F>& rhs) {
    if ( &rhs == this ) return;
    if ( ! rhs.isValid() ) return;
    Index nsam = rhs.nSample();
    if ( nsam == 0 ) return;
    Index namp = rhs.nCompact();
    Index npha = nsam - namp + 1;
    if ( namp < npha || namp > npha + 1 ) { clear(); return; }
    m_amps.resize(namp);
    m_phas.resize(npha);
    if ( nsam == 0 ) return;
    F gfac = 1.0;
    if ( rhs.isPower() ) gfac *= sqrt(nsam);
    if ( rhs.isBin() ) gfac *= nsam;
    if ( this->isPower() ) gfac /= sqrt(nsam);
    if ( this->isBin() ) gfac /= nsam;
    F tfac = 1.0;
    if ( rhs.isUnit() && ! this->isUnit() ) tfac = 1/sqrt(2.0);
    if ( ! rhs.isUnit() && this->isUnit() ) tfac = sqrt(2.0);
    for ( Index ifrq=0; ifrq<namp; ++ifrq ) {
      F amp = gfac*rhs.amplitude(ifrq);
      if ( this->isAliased(ifrq) ) amp *= tfac;
      m_amps[ifrq] = amp;
      if ( ifrq < npha ) m_phas[ifrq] = rhs.phase(ifrq);
    }
  }

  // Clear data.
  void clear() {
    m_amps.clear();
    m_phas.clear();
  }

  // Move data in.
  int moveIn(FloatVector& amps, FloatVector& phas) {
    Index namp = amps.size();
    Index npha = phas.size();
    if ( namp < npha || namp > npha + 1 ) {
      clear();
      return 1;
    }
    m_amps = std::move(amps);
    m_phas = std::move(phas);
    return 0;
  }

  // Copy data in.
  int copyIn(const FloatVector& amps, const FloatVector& phas) {
    Index namp = amps.size();
    Index npha = phas.size();
    if ( namp < npha || namp > npha + 1 ) {
      clear();
      return 1;
    }
    m_amps = amps;
    m_phas = phas;
    return 0;
  }

  // Move data out.
  int moveOut(FloatVector& amps, FloatVector& phas) {
    amps = std::move(m_amps);
    phas = std::move(m_phas);
    return 0;
  }

  // Copy the data out.
  int copyOut(FloatVector& amps, FloatVector& phas) const {
    amps = m_amps;
    phas = m_phas;
    return 0;
  }

  // Global normalization.
  Index globalNormalization() const override { return m_gnorm; }

  // Term normalization.
  Index termNormalization() const override { return m_tnorm; }

  // Check and return dimension information.
  Index nSample() const override { return nCompact() ? nCompact() + nPhase() - 1 : 0; }
  Index nCompact() const override { return nAmplitude(); }

  // Overide methods that return DFT terms for any representation.
  F amplitude(Index ifrq) const {
    return ifrq < this->nAmplitude() ? this->compactAmplitude(ifrq) :
           ifrq < this->size() ? this->compactAmplitude(this->size()-ifrq) :
           this->badValue();
  }
  F phase(Index ifrq) const {
    return ifrq < this->nPhase() ? this->compactPhase(ifrq) :
           2*ifrq == this->size() ? 0.0 :
           ifrq < this->size() ? this->compactPhase(this->size()-ifrq) :
           this->badValue();
  }
  F real(Index ifrq) const {
    return ifrq < this->size() ? amplitude(ifrq)*cos(phase(ifrq)) :
    this->badValue();
  }
  F imag(Index ifrq) const {
    return ifrq < this->nPhase() ? amplitude(ifrq)*sin(phase(ifrq)) :
           2*ifrq == this->size() ? 0.0 :
           ifrq < this->size() ? -amplitude(ifrq)*sin(phase(ifrq)) :
           this->badValue();
  }

  // Return the compact data.
  Index nAmplitude() const { return m_amps.size(); }
  Index nPhase() const { return m_phas.size(); }
  F compactAmplitude(Index ifrq) const {
    return this->inRange(ifrq) ? m_amps[ifrq] : this->badValue();
  }
  F compactPhase(Index ifrq) const {
    return this->inRange(ifrq) ? m_phas[ifrq] : this->badValue();
  }

private:

  // Data.
  GlobalNormalization m_gnorm = RealDftNormalization::defaultGlobalNormalization();
  TermNormalization m_tnorm = RealDftNormalization::defaultTermNormalization();
  FloatVector m_amps;
  FloatVector m_phas;

};

//**********************************************************************

#endif
