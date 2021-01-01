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
  using Norm = RealDftNormalization;

  // Default ctor.
  CompactRealDftData() =default;

  // Ctor from normalization. Leaves object empty and so invalid.
  explicit CompactRealDftData(Norm norm)
  : m_norm(norm) { }

  // Ctor from normalization.
  CompactRealDftData(Norm norm, Index nsam)
  : m_norm(norm) {
    reset(nsam);
  }

  // Ctor from normalization and data.
  // Must have amps.size() - phas.size() = 0 or 1 for valid DFT.
  // If not, the object is left empty.
  CompactRealDftData(Norm norm, const FloatVector& amps, const FloatVector& phas)
  : m_norm(norm), m_amps(amps), m_phas(phas) {
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

  // Normalization.
  const Norm& normalization() const override { return m_norm; }

  // Clear data.
  void clear() override {
    m_amps.clear();
    m_phas.clear();
  }

  // Reset the DFT data.
  void reset(Index nsam) override {
    Index namp = nsam > 0 ? nsam/2 + 1 : 0;
    Index npha = (nsam + 1)/2;
    m_amps.resize(namp);
    m_phas.resize(npha);
    for ( F& amp : m_amps ) amp = 0.0;
    for ( F& pha : m_phas ) pha = 0.0;
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

  // Set the data.
  int setAmplitude(Index ifrq, F val) {
    if ( ifrq >= nAmplitude() ) return 1;
    m_amps[ifrq] = val;
    return 0;
  }
  int setPhase(Index ifrq, F val) {
    if ( ifrq >= nPhase() ) return 1;
    m_phas[ifrq] = val;
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

  // Check and return dimension information.
  Index nSample() const override { return nCompact() ? nCompact() + nPhase() - 1 : 0; }
  Index nCompact() const override { return nAmplitude(); }

  // Overide methods that return DFT terms for any representation.
  F amplitude(Index ifrq) const override {
    return ifrq < this->nAmplitude() ? this->compactAmplitude(ifrq) :
           ifrq < this->size() ? this->compactAmplitude(this->size()-ifrq) :
           this->badValue();
  }
  F phase(Index ifrq) const override {
    return ifrq < this->nPhase() ? this->compactPhase(ifrq) :
           2*ifrq == this->size() ? 0.0 :
           ifrq < this->size() ? -this->compactPhase(this->size()-ifrq) :
           this->badValue();
  }
  F real(Index ifrq) const override {
    return ifrq < this->size() ? amplitude(ifrq)*cos(phase(ifrq)) :
    this->badValue();
  }
  F imag(Index ifrq) const override {
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
  Norm m_norm;
  FloatVector m_amps;
  FloatVector m_phas;

};

//**********************************************************************

#endif
