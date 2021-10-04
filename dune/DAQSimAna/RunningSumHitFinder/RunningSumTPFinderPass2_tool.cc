////////////////////////////////////////////////////////////////////////
// Class:       RunningSumTPFinderPass2
// Plugin Type: service (art v2_10_03)
// File:        RunningSumTPFinderPass2_service.cc
//
// Generated at Tue Jun  5 07:51:38 2018 by Philip Rodrigues using cetskelgen
// from cetlib version v3_02_00.
////////////////////////////////////////////////////////////////////////

#include "dune/DAQSimAna/RunningSumHitFinder/RunningSumTPFinderPass2.h"

#include "dune/DAQSimAna/AlgParts.h"

#include <algorithm> // for std::transform
#include <numeric> // for std::accumulate


RunningSumTPFinderPass2::RunningSumTPFinderPass2(fhicl::ParameterSet const & p)
  : m_threshold          (p.get<unsigned int>      ("Threshold"            ,                          10)),
    m_useSignalKill      (p.get<bool>              ("UseSignalKill"        ,                        true)),
    m_signalKillLookahead(p.get<short>             ("SignalKillLookahead"  ,                           5)),
    m_signalKillThreshold(p.get<short>             ("SignalKillThreshold"  ,                          15)),
    m_signalKillNContig  (p.get<short>             ("SignalKillNContig"    ,                           1)),
    m_frugalNContig      (p.get<short>             ("FrugalPedestalNContig",                          10)),
    m_doFiltering        (p.get<bool>              ("DoFiltering"          ,                        true)),
    m_downsampleFactor   (p.get<unsigned int>      ("DownsampleFactor"     ,                           1)),
    m_filterTaps         (p.get<std::vector<short>>("FilterCoeffs"         , {2,  9, 23, 31, 23,  9,  2})),
    m_multiplier         (std::accumulate(m_filterTaps.begin(), m_filterTaps.end(), 0))
    // Default filter taps calculated by:
    // np.round(scipy.signal.firwin(7, 0.1)*100)
    // Initialize member data here.
{
  std::cout << "Threshold is " << m_threshold << std::endl;
}

std::vector<short> RunningSumTPFinderPass2::downSample(const std::vector<short>& orig) {

  //---------------------------------------------
  // Do the downsampling
  //---------------------------------------------
  if (m_downsampleFactor==1) {
    return orig;
  }
  else {
    std::vector<short> waveform;
    for(size_t i=0; i<orig.size(); i+=m_downsampleFactor) {
      waveform.push_back(orig[i]);
    }
    return waveform;
  }
}

std::vector<short> RunningSumTPFinderPass2::findPedestal(const std::vector<short>& waveform)
{
  //---------------------------------------------
  // Pedestal subtraction
  //---------------------------------------------
  const std::vector<short>& pedestal=m_useSignalKill ?
    frugal_pedestal_sigkill(waveform,
			    m_signalKillLookahead,
			    m_signalKillThreshold,
			    m_signalKillNContig) :
    frugal_pedestal(waveform, m_frugalNContig);
  return pedestal;
}

std::vector<short> RunningSumTPFinderPass2::filter(const std::vector<short>& pedsub) {

  //---------------------------------------------
  // Filtering
  //---------------------------------------------
  const size_t ntaps = m_filterTaps.size();
  const short*  taps = m_filterTaps.data();

  std::vector<short> filtered(m_doFiltering ? 
			      apply_fir_filter(pedsub, ntaps, taps) :
			      pedsub);
  if (!m_doFiltering) {
    std::transform(filtered.begin(), filtered.end(),
		   filtered.begin(), 
		   [=](short a) { return a*m_multiplier; });
  }
  return filtered;
}

void
RunningSumTPFinderPass2::hitFinding(const std::vector<short>& waveform,
				    std::vector<RunningSumTPFinderTool::Hit>& hits,
				    int channel) {

  //---------------------------------------------
  // Hit finding
  //---------------------------------------------
  bool is_hit  = false;
  bool was_hit = false;
  RunningSumTPFinderTool::Hit hit(channel, 0, 0, 0);
  for(size_t isample=0; isample<waveform.size()-1; ++isample){
    int   sample_time = isample * m_downsampleFactor;
    short adc         = waveform[isample];
    is_hit = adc >  (short)m_threshold*10;
    if(is_hit && !was_hit) {
      hit.startTime         = sample_time;
      hit.charge            = adc;
      hit.timeOverThreshold = m_downsampleFactor;
    }
    if(is_hit && was_hit) {
      hit.charge            += adc*m_downsampleFactor;
      hit.timeOverThreshold += m_downsampleFactor;
    }
    if(!is_hit && was_hit) {
      hit.charge /= m_multiplier;
      hits.push_back(hit);
    }
    was_hit = is_hit;
  }
}

std::vector<RunningSumTPFinderTool::Hit>
RunningSumTPFinderPass2::findHits(const std::vector<unsigned int>& channel_numbers, 
				  const std::vector<std::vector<short>>& collection_samples) {

  auto hits = std::vector<RunningSumTPFinderTool::Hit>();
  std::cout << "findHits called with "      << collection_samples.size()
	    << " channels. First chan has " << collection_samples[0].size() << " samples" << std::endl;

  for(size_t ich=0; ich<collection_samples.size(); ++ich){
    const std::vector<short>& waveformOrig = collection_samples[ich];

    std::vector<short> waveform  = downSample  (waveformOrig);
    std::vector<short> pedestal  = findPedestal(waveform    );
    std::vector<short> pedsub(waveform.size(), 0);
    for(size_t i=0; i<pedsub.size(); ++i)
      pedsub[i]=waveform[i]-pedestal[i];
    std::vector<short> filtered  = filter(pedsub);
    hitFinding(filtered, hits, channel_numbers[ich]);
  }
  std::cout << "Returning " << hits.size() << " hits" << std::endl;
  std::cout << "hits/channel=" << float(hits.size())/collection_samples   .size() << std::endl;
  std::cout << "hits/tick="    << float(hits.size())/collection_samples[0].size() << std::endl;
  return hits;
}

DEFINE_ART_CLASS_TOOL(RunningSumTPFinderPass2)
