#ifndef RUNNINGSUMTPFINDERPASS4_H
#define RUNNINGSUMTPFINDERPASS4_H

#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"

#include "dune/DAQSimAna/RunningSumHitFinder/RunningSumTPFinderTool.h"

class RunningSumTPFinderPass4 : public RunningSumTPFinderTool {
 public:
  explicit RunningSumTPFinderPass4(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  virtual std::vector<RunningSumTPFinderTool::Hit>
    findHits(const std::vector<unsigned int>& channel_numbers, 
             const std::vector<std::vector<short>>& collection_samples);
    

 protected:
  std::vector<short> downSample  (const std::vector<short>& orig);
  std::vector<short> findPedestal(const std::vector<short>& orig);
  std::vector<short> filter      (const std::vector<short>& orig);
  std::vector<short> runSum      (const std::vector<short>& orig);
  void hitFinding(const std::vector<short>& waveform, 
		  const std::vector<short>& iqr,
		  std::vector<RunningSumTPFinderTool::Hit>& hits, 
		  int channel);

  bool               m_useSignalKill;
  short              m_signalKillLookahead;
  short              m_signalKillThreshold;
  short              m_signalKillNContig;
  short              m_frugalNContig;
  bool               m_doFiltering;
  unsigned int       m_downsampleFactor;
  std::vector<short> m_filterTaps;
  int                m_multiplier;
  float              m_sigmaThreshold;
};

#endif
