// RealDftDataCompact.h
//
// David Adams
// April 2019
//
// Concrete class that holds a compact representation of DFT data.

#ifndef RealDftDataCompact_H
#define RealDftDataCompact_H

#include "dune/DuneCommon/RealDftData.h"

//**********************************************************************

// Interface for a compact representation of DFT data.

template<typename F>
class IRealDftDataCompact : public RealDftData<F> {

public:

  using Index = RealDftNormalization::Index;

  // Add interface that provides access to compact terms.
  virtual Index nAmplitude() const =0;
  virtual Index nPhase() const =0;
  virtual F compactAmplitude(Index ifrq) const =0;
  virtual F compactPhase(Index ifrq) const =0;

  // Overide methods that return DFT terms for either representation.
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

};

//**********************************************************************

// Concrete compact representation of DFT data.

template<typename F,
         RealDftNormalization::GlobalNormalization GNORM,
         RealDftNormalization::TermNormalization TNORM>
class RealDftDataCompact : public IRealDftDataCompact<F> {

public:

  using Index = typename RealDftData<F>::Index;
  using FloatVector = std::vector<F>;

  // Ctor for no data.
  RealDftDataCompact() =default;

  // Ctor from data.
  // Must have amps.size() - phas.size() = 0 or 1 for valid DFT.
  // If not, the object is left empty.
  RealDftDataCompact(const FloatVector& amps, FloatVector& phas)
  : m_amps(amps), m_phas(phas) {
    Index namp = amps.size();
    Index npha = phas.size();
    if ( namp < npha || namp > npha + 1 ) clear();
  }
  
  // Copy ctor. Changes representation and normalization if needed.
  RealDftDataCompact(const RealDftData<F>& rhs) {
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
      if ( isAliased(ifrq) ) amp *= tfac;
      m_amps[ifrq] = amp;
      if ( ifrq < npha ) m_phas[ifrq] = rhs.phase(ifrq);
    }
  }

  // Clear data.
  void clear() {
    m_amps.clear();
    m_phas.clear();
  }

  // Global normalization.
  Index globalNormalization() const override { return GNORM; }

  // Term normalization.
  Index termNormalization() const override { return TNORM; }

  // Check and return dimension information.
  Index nAmplitude() const override { return m_amps.size(); }
  Index nPhase() const override { return m_phas.size(); }
  Index nSample() const override { return nCompact() ? nCompact() + nPhase() - 1 : 0; }
  Index nCompact() const override { return nAmplitude(); }

  // Return the data.
  F compactAmplitude(Index ifrq) const override {
    return this->inRange(ifrq) ? m_amps[ifrq] : this->badValue();
  }
  F compactPhase(Index ifrq) const override {
    return this->inRange(ifrq) ? m_phas[ifrq] : this->badValue();
  }

private:

  // Data.
  FloatVector m_amps;
  FloatVector m_phas;

};

//**********************************************************************

//using ReadDftDataFloat = RealDftData<float>;
//using ReadDftDataDouble = RealDftData<double>;

#endif
