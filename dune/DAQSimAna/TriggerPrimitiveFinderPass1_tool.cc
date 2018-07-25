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
    

private:
    unsigned int m_threshold;

    bool m_useSignalKill;
    short m_signalKillLookahead;
    short m_signalKillThreshold;
    short m_signalKillNContig;

    short m_frugalNContig;

    bool m_doFiltering;
};


TriggerPrimitiveFinderPass1::TriggerPrimitiveFinderPass1(fhicl::ParameterSet const & p)
    : m_threshold(p.get<unsigned int>("Threshold", 10)),
      m_useSignalKill(p.get<bool>("UseSignalKill", true)),
      m_signalKillLookahead(p.get<short>("SignalKillLookahead", 5)),
      m_signalKillThreshold(p.get<short>("SignalKillThreshold", 15)),
      m_signalKillNContig(p.get<short>("SignalKillNContig", 1)),
      m_frugalNContig(p.get<short>("FrugalPedestalNContig", 10)),
      m_doFiltering(p.get<bool>("DoFiltering", true))

// Initialize member data here.
{
    std::cout << "Threshold is " << m_threshold << std::endl;
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

    // Taps calculated by:
    //  np.round(scipy.signal.firwin(7, 0.1)*100)
    const size_t ntaps=7;
    const short taps[ntaps]={2,  9, 23, 31, 23,  9,  2};
    const int multiplier=100;

    for(size_t ich=0; ich<collection_samples.size(); ++ich){
        const std::vector<short>& waveform=collection_samples[ich];

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

        //---------------------------------------------
        // Filtering
        //---------------------------------------------
        std::vector<short> filtered(m_doFiltering ? 
                                    apply_fir_filter(pedsub, ntaps, taps) :
                                    pedsub);
        if(!m_doFiltering){
            std::transform(filtered.begin(), filtered.end(),
                           filtered.begin(), 
                           [](short a) { return a*multiplier; });
        }
        // Print out the waveforms on one channel for debugging

        // if(channel_numbers[ich]==1600){
        //     for(auto s: waveform){ std::cout << s << " ";}
        //     std::cout << std::endl;
        //     for(auto s: pedsub){ std::cout << s << " ";}
        //     std::cout << std::endl;
        //     for(auto s: filtered){ std::cout << s << " ";}
        //     std::cout << std::endl;
        // }

        // if(ich>10) exit(0);
        //---------------------------------------------
        // Hit finding
        //---------------------------------------------
        bool is_hit=false;
        bool was_hit=false;
        TriggerPrimitiveFinderTool::Hit hit(channel_numbers[ich], 0, 0, 0);
        for(size_t isample=0; isample<filtered.size()-1; ++isample){
            // if(ich>11510) std::cout << isample << " " << std::flush;
            short adc=filtered[isample];
            is_hit=adc>(short)m_threshold;
            if(is_hit && !was_hit){
                // We just started a hit. Set the start time
                hit.startTime=isample;
                hit.charge=adc;
                hit.timeOverThreshold=1;
            }
            if(is_hit && was_hit){
                hit.charge+=adc;
                hit.timeOverThreshold++;
            }
            if(!is_hit && was_hit){
                // The hit is over. Push it to the output vector
                hit.charge/=multiplier;
                hits.push_back(hit);
            }
            was_hit=is_hit;
        }
        // std::cout << std::endl;
    }
    std::cout << "Returning " << hits.size() << " hits" << std::endl;
    std::cout << "hits/channel=" << float(hits.size())/collection_samples.size() << std::endl;
    std::cout << "hits/tick=" << float(hits.size())/collection_samples[0].size() << std::endl;
    return hits;
}

DEFINE_ART_CLASS_TOOL(TriggerPrimitiveFinderPass1)
