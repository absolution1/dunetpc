// RunData.h
//
// David Adams
// May 2018
//
// Detector conditions that are typically constant over the
// course of a run.
//
// Schema are documented at
//  https://wiki.dunescience.org/wiki/ProtoDUNE_run_configuration

#ifndef RunData_H
#define RunData_H

#include <string>
#include <vector>
#include <iostream>

class RunData {

public:

  using Name = std::string;
  using Index = unsigned int;
  using IndexVector = std::vector<Index>;

  // Ctor from a run.
  explicit RunData(Index a_run =0) : m_run(a_run) { };

  // Return a value.
  Index run() const { return m_run; }
  Name cryostat() const { return m_cryostat; }
  const IndexVector& apas() const { return m_apas; }
  float gain() const { return m_gain; }
  float shaping() const { return m_shaping; }
  float leakage() const { return m_leakage; }
  float hvfrac() const { return m_hvfrac; }
  Index pulserAmplitude() const { return m_pulserAmplitude; }
  Index pulserSource() const { return m_pulserSource; }
  Index pulserPeriod() const { return m_pulserPeriod; }
  Name phaseGroup() const { return m_phaseGroup; }
  const IndexVector& phases() const { return m_phases; }

  // Return if value is defined.
  bool isValid() const { return run(); }
  bool haveRun() const { return run(); }
  bool haveCryostat() const { return cryostat().size(); }
  bool haveApas() const { return apas().size(); }
  bool haveGain() const { return m_gain; }
  bool haveShaping() const { return m_shaping; }
  bool haveLeakage() const { return m_leakage; }
  bool haveHvfrac() const { return m_hvfrac; }
  bool havePulserAmplitude() const { return m_pulserAmplitude != 999; }
  bool havePulserSource() const { return m_pulserSource != 0; }
  bool havePulserPeriod() const { return m_pulserPeriod != 0; }
  bool havePhaseGroup() const { return m_phaseGroup.size(); }
  bool havePhases() const { return m_phaseGroup.size(); }

  // Setters.
  void setRun(Index val) { m_run = val; }
  void setCryostat(Name val) { m_cryostat = val; }
  void setApas(const IndexVector& val) { m_apas = val; }
  void setGain(float val) { m_gain = val; }
  void setShaping(float val) { m_shaping = val; }
  void setLeakage(float val) { m_leakage = val; }
  void setHvfrac(float val) { m_hvfrac = val; }
  void setPulserAmplitude(Index val) { m_pulserAmplitude = val; }
  void setPulserSource(Index val) { m_pulserSource = val; }
  void setPulserPeriod(Index val) { m_pulserPeriod = val; }
  void setPhaseGroup(Name val) { m_phaseGroup = val; }
  void setPhases(const IndexVector& val) { m_phases = val; }

  // Accessors.
  Index&       accessRun()             { return m_run; }
  Name&        accessCryostat()        { return m_cryostat; }
  IndexVector& accessApas()            { return m_apas; }
  float&       accessGain()            { return m_gain; }
  float&       accessShaping()         { return m_shaping; }
  float&       accessLeakage()         { return m_leakage; }
  float&       accessHvfrac()          { return m_hvfrac; }
  Index&       accessPulserAmplitude() { return m_pulserAmplitude; }
  Index&       accessPulserSource()    { return m_pulserSource; }
  Index&       accessPulserPeriod()    { return m_pulserPeriod; }
  Name&        accessPhaseGroup()      { return m_phaseGroup; }
  IndexVector& accessPhases()          { return m_phases; }

  std::ostream& print(std::ostream& lhs =std::cout) const {
    const std::string sep = "\n";
    if ( ! isValid() ) {
      lhs << "Invalid run data.";
      return lhs;
    }
    lhs << "Data for run " << run();
    if ( haveCryostat() )        lhs << sep << "      Cryostat: " << cryostat();
    if ( haveApas() )            lhs << sep << "          APAs: ";
    int doComma = 0;
    for ( Index iapa : apas() ) lhs << (doComma++ ? ", " : "") << iapa;
    if ( haveGain() )            lhs << sep << "          Gain: " << gain() << " mV/fC";
    if ( haveShaping() )         lhs << sep << "  Shaping time: " << shaping() << " us";
    if ( haveLeakage() )         lhs << sep << "  Leakage cur.: " << leakage() << " pA";
    if ( haveHvfrac() )          lhs << sep << "      HV frac.: " << hvfrac();
    if ( havePulserAmplitude() ) lhs << sep << "  Pulser ampl.: " << pulserAmplitude();
    if ( havePulserSource() )    lhs << sep << " Pulser source: " << pulserSource();
    if ( havePulserPeriod() )    lhs << sep << " Pulser period: " << pulserPeriod();
    if ( havePhaseGroup()   ) {  lhs << sep << " Phases grouped by " << phaseGroup() << ": [";
      bool first = true;
      for ( Index ipha : phases() ) {
        if ( first ) first = false;
        else lhs << ", ";
        lhs << ipha;
      }
      lhs << "]";
    }
    return lhs;
  }

private:

  // Data.
  Index m_run = 0;
  Name m_cryostat;
  IndexVector m_apas;
  float m_gain = 0.0;
  float m_shaping = 0.0;
  float m_leakage = 0.0;
  float m_hvfrac = 0.0;
  Index m_pulserAmplitude = 999;
  Index m_pulserSource = 0;  // 1=preamp, 2=FEMB
  Index m_pulserPeriod = 0;
  Name  m_phaseGroup;
  IndexVector m_phases;

};

std::ostream& operator<<(std::ostream& lhs, const RunData& rhs) {
  return rhs.print(lhs);
}

#endif
