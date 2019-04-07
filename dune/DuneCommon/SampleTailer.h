// SampleTailer.h
//
// David Adams
// April 2019
//
// Class that adds or remove exponential tails from a sampled
// distribution, i.e. a vector of floats.
//
// The data is a sum of signal, tail and pedestal:
//   d[i] = s[i] + t[i] + p
// Or solving for the signal:
//   s[i] = d[i] - t[i] - p
// The tail decays and has contribution from signal as follows:
//   t[i] = beta t[i-1] - alpha s[i-1]
// Tails are added or removed from a distibution by combining this
// in iteration with the repective preceding expression.
//
// The tails are characterized by four parameters:
//   tdec (or beta) - Decay time in ticks. By definition beta = exp(-1/tdec).
//   alpha - Starting tail fraction of the signal (see above).
//   pedestal - Data pedestal (p in the expression above).
//   tail0 - Tail level in the first sample, i.e. t[0] above.
// Note that if
//   alpha = 1/beta - 1
// then the tail precisely cancels the signal:
//   SUM_i (s[i] + t[i]) = 0

#ifndef SampleTailer_H
#define SampleTailer_H

#include <vector>
#include <string>

class TH1;

class SampleTailer {

public:

  using Index = unsigned int;
  using FloatVector = std::vector<float>;
  using Name = std::string;

  // Ctor from the decay time.
  // The tail fraction is chosen so signal cancels tail.
  // The pedestal and t0 are set to zero.
  SampleTailer(float a_decayTime);

  // Ctor from decay time and alpha.
  // The pedestal and t0 are set to zero.
  SampleTailer(float a_decayTime, float a_alpha);

  // Dtor.
  ~SampleTailer();

  // Clear the data, signal and tail vectors.
  void clear();

  // Parameter setters.
  // Data, signal and tail are cleared.
  int setDecayTime(float val, bool cancelSignal);
  int setBeta(float val, bool cancelSignal);
  int setAlpha(float val);
  int setPedestal(float val);
  int setTail0(float val);
  int setUnit(Name val) { m_unit = val; return 0; }

  // Evaluate the signal and tail for provided data.
  int setData(const FloatVector& inData);

  // Evaluate the signal and tail for data value zero for nsam smples..
  int setDataZero(Index nsam);

  // Evaluate the data and tail for provided signal.
  int setSignal(const FloatVector& inSignal);

  // Return if the parameters are valid, e.g. beta > 0.
  bool isValid() const { return m_beta > 0.0; }

  // Getters.
  float decayTime() const { return m_decayTime; }
  float beta() const { return m_beta; }
  float alpha() const { return m_alpha; }
  float pedestal() const { return m_pedestal; }
  float tail0() const { return m_tail0; }
  const FloatVector& data() const { return m_d; }
  float data(Index isam) const { return isam < m_d.size() ? m_d[isam] : 0.0; }
  const FloatVector& signal() const { return m_s; }
  float signal(Index isam) const { return isam < m_s.size() ? m_s[isam] : 0.0; }
  const FloatVector& tail() const { return m_t; }
  float tail(Index isam) const { return isam < m_t.size() ? m_t[isam] : 0.0; }
  Index size() const { return m_d.size(); }
  Name unit() const { return m_unit; }
  TH1* dataHist() const { return m_phd; }
  TH1* signalHist() const { return m_phs; }
  TH1* tailHist() const { return m_pht; }

  // Return histograms.
  TH1* dataHist();
  TH1* signalHist();
  TH1* tailHist();
  void makeHists();

private:

  float m_decayTime;
  float m_beta;
  float m_alpha;
  float m_pedestal;
  float m_tail0;
  Name m_unit;
  FloatVector m_d;
  FloatVector m_s;
  FloatVector m_t;
  TH1* m_phd =nullptr;
  TH1* m_phs =nullptr;
  TH1* m_pht =nullptr;

};

#endif
