////////////////////////////////////////////////////////////////////////
// Class:       TriggerPrimitiveFinderPass1
// Plugin Type: service (art v2_10_03)
// File:        TriggerPrimitiveFinderPass1_service.cc
//
// Generated at Tue Jun  5 07:51:38 2018 by Philip Rodrigues using cetskelgen
// from cetlib version v3_02_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"

#include "dune/DAQSimAna/TriggerPrimitiveFinderTool.h"

#include "dune/DAQSimAna/AlgParts.h"

#include <algorithm> // for std::transform

class TriggerPrimitiveFinderPass1 : public TriggerPrimitiveFinderTool {
public:
  explicit TriggerPrimitiveFinderPass1(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

    virtual std::vector<TriggerPrimitiveFinderTool::Hit>
    findHits(const std::vector<unsigned int>& channel_numbers, 
             const std::vector<std::vector<short>>& collection_samples);
    

protected:
    std::vector<short> downSample(const std::vector<short>& orig);
    std::vector<short> pedestalSubtract(const std::vector<short>& orig);
    std::vector<short> filter(const std::vector<short>& orig);
    void hitFinding(const std::vector<short>& waveform, std::vector<TriggerPrimitiveFinderTool::Hit>& hits, int channel);

private:
    unsigned int m_threshold;

    bool m_useSignalKill;
    short m_signalKillLookahead;
    short m_signalKillThreshold;
    short m_signalKillNContig;

    short m_frugalNContig;

    bool m_doFiltering;

    unsigned int m_downsampleFactor;
    std::vector<short> m_filterTaps;
};


TriggerPrimitiveFinderPass1::TriggerPrimitiveFinderPass1(fhicl::ParameterSet const & p)
    : m_threshold(p.get<unsigned int>("Threshold", 10)),
      m_useSignalKill(p.get<bool>("UseSignalKill", true)),
      m_signalKillLookahead(p.get<short>("SignalKillLookahead", 5)),
      m_signalKillThreshold(p.get<short>("SignalKillThreshold", 15)),
      m_signalKillNContig(p.get<short>("SignalKillNContig", 1)),
      m_frugalNContig(p.get<short>("FrugalPedestalNContig", 10)),
      m_doFiltering(p.get<bool>("DoFiltering", true)),
      m_downsampleFactor(p.get<unsigned int>("DownsampleFactor", 1)),
      // Default filter taps calculated by:
      //  np.round(scipy.signal.firwin(7, 0.1)*100)
      m_filterTaps(p.get<std::vector<short>>("FilterCoeffs", {2,  9, 23, 31, 23,  9,  2}))

// Initialize member data here.
{
    std::cout << "Threshold is " << m_threshold << std::endl;
}

std::vector<short> TriggerPrimitiveFinderPass1::downSample(const std::vector<short>& orig)
{
    //---------------------------------------------
    // Do the downsampling
    //---------------------------------------------
    if(m_downsampleFactor==1){
        return orig;
    }
    else{
        std::vector<short> waveform;
        for(size_t i=0; i<waveformOrig.size(); i+=m_downsampleFactor){
            waveform.push_back(waveformOrig[i]);
        }
        return waveform;
    }
}

std::vector<short> TriggerPrimitiveFinderPass1::pedestalSubtract(const std::vector<short>& waveform)
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

    std::vector<short> pedsub(waveform.size(), 0);
    for(size_t i=0; i<pedsub.size(); ++i){
        pedsub[i]=waveform[i]-pedestal[i];
    }
    return pedsub;
}

std::vector<short> TriggerPrimitiveFinderPass1::filter(const std::vector<short>& pedsub)
{
    //---------------------------------------------
    // Filtering
    //---------------------------------------------
    std::vector<short> filtered(m_doFiltering ? 
                                apply_fir_filter(pedsub, ntaps, taps) :
                                pedsub);
    if(!m_doFiltering){
        std::transform(filtered.begin(), filtered.end(),
                       filtered.begin(), 
                       [=](short a) { return a*multiplier; });
    }
    return filtered;
}

void
TriggerPrimitiveFinderPass1::hitFinding(const std::vector<short>& waveform,
                                        std::vector<TriggerPrimitiveFinderTool::Hit>& hits,
                                        int channel)
{
    //---------------------------------------------
    // Hit finding
    //---------------------------------------------
    bool is_hit=false;
    bool was_hit=false;
    TriggerPrimitiveFinderTool::Hit hit(channel, 0, 0, 0);
    for(size_t isample=0; isample<filtered.size()-1; ++isample){
        // if(ich>11510) std::cout << isample << " " << std::flush;
        int sample_time=isample*m_downsampleFactor;
        short adc=filtered[isample];
        is_hit=adc>(short)m_threshold;
        if(is_hit && !was_hit){
            // We just started a hit. Set the start time
            hit.startTime=sample_time;
            hit.charge=adc;
            hit.timeOverThreshold=m_downsampleFactor;
        }
        if(is_hit && was_hit){
            hit.charge+=adc*m_downsampleFactor;
            hit.timeOverThreshold+=m_downsampleFactor;
        }
        if(!is_hit && was_hit){
            // The hit is over. Push it to the output vector
            hit.charge/=multiplier;
            hits.push_back(hit);
        }
        was_hit=is_hit;
    }
}


std::vector<TriggerPrimitiveFinderTool::Hit>
TriggerPrimitiveFinderPass1::findHits(const std::vector<unsigned int>& channel_numbers, 
                                      const std::vector<std::vector<short>>& collection_samples)
{
    auto hits=std::vector<TriggerPrimitiveFinderTool::Hit>();
    std::cout << "findHits called with " << collection_samples.size()
              << " channels. First chan has " << collection_samples[0].size() << " samples" << std::endl;
    // std::cout << "First few samples: ";
    // for(int i=0; i<10; ++i) std::cout << collection_samples[0][i] << " ";
    // std::cout << std::endl;

    const size_t ntaps=m_filterTaps.size();
    const short* taps=m_filterTaps.data();
    const int multiplier=std::accumulate(m_filterTaps.begin(), m_filterTaps.end(), 0);

    for(size_t ich=0; ich<collection_samples.size(); ++ich){
        const std::vector<short>& waveformOrig=collection_samples[ich];

        std::vector<short> waveform=downSample(waveformOrig);
        std::vector<short> pedsub=pedestalSubtract(waveform);
        std::vector<short> filtered=filter(pedsub);
        hitFinding(filtered, hits, channel_numbers[ich]);
    }
    std::cout << "Returning " << hits.size() << " hits" << std::endl;
    std::cout << "hits/channel=" << float(hits.size())/collection_samples.size() << std::endl;
    std::cout << "hits/tick=" << float(hits.size())/collection_samples[0].size() << std::endl;
    return hits;
}

DEFINE_ART_CLASS_TOOL(TriggerPrimitiveFinderPass1)
