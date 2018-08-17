////////////////////////////////////////////////////////////////////////
//
// An offline-friendly packaging of the Central Trigger Board (CTB)'s data for ProtoDUNE-SP
// Tom Junk, August 10, 2018
//
////////////////////////////////////////////////////////////////////////

#ifndef  pdspctb_H
#define  pdspctb_H

#include "RtypesCore.h"
#include <stdint.h>

namespace raw {

  namespace ctb {

    struct Trigger {
      uint32_t word_type;
      ULong64_t trigger_word;
      ULong64_t timestamp;
    };

    struct ChStatus {
      uint32_t word_type;
      uint32_t pds;
      uint32_t crt;
      uint32_t beam_hi;
      uint32_t beam_lo;
      ULong64_t timestamp;
    };

    struct Feedback {
      uint32_t word_type;
      uint32_t padding;
      uint32_t source;
      uint32_t code;
      ULong64_t timestamp;
    };

    struct Misc {
      uint32_t word_type;
      ULong64_t payload;
      ULong64_t timestamp;
    };

    struct WordIndex {
      uint32_t word_type;
      uint32_t index;
    };

    class pdspctb
    {

    public:

      pdspctb() {}; // Constructor of an emtpy data product

      // constructor with all of the vectors set
    pdspctb(std::vector<raw::ctb::Trigger> &trigs,
	    std::vector<raw::ctb::ChStatus> &chstats,
	    std::vector<raw::ctb::Feedback> &fbs,
	    std::vector<raw::ctb::Misc> &m,
	    std::vector<raw::ctb::WordIndex> &wordindexes) : 
      fTriggers(trigs), fChStatuses(chstats), fFeedbacks(fbs), fMiscs(m), fIndexes(wordindexes) {};


      const std::vector<raw::ctb::Trigger>&     GetTriggers() const;   
      const std::vector<raw::ctb::ChStatus>&    GetChStatuses() const; 
      const std::vector<raw::ctb::Feedback>&    GetFeedbacks() const;  
      const std::vector<raw::ctb::Misc>&        GetMiscs() const;
      const std::vector<raw::ctb::WordIndex>&   GetIndexes() const;

      size_t  GetNTriggers() const;   
      size_t  GetNChStatuses() const; 
      size_t  GetNFeedbacks() const;  
      size_t  GetNMiscs() const;      
      size_t  GetNIndexes() const;      

      const raw::ctb::Trigger&    GetTrigger(size_t i) const;   
      const raw::ctb::ChStatus&   GetChStatuse(size_t i) const; 
      const raw::ctb::Feedback&   GetFeedback(size_t i) const;  
      const raw::ctb::Misc&       GetMisc(size_t i) const;      
      const raw::ctb::WordIndex&  GetIndex(size_t i) const;      

    private:

      std::vector<raw::ctb::Trigger> fTriggers;
      std::vector<raw::ctb::ChStatus> fChStatuses;
      std::vector<raw::ctb::Feedback> fFeedbacks;
      std::vector<raw::ctb::Misc> fMiscs;
      std::vector<raw::ctb::WordIndex> fIndexes;

    };


  } // namespace ctb
} // namespace raw

// accessors

const std::vector<raw::ctb::Trigger>&       raw::ctb::pdspctb::GetTriggers()   const { return fTriggers; }
const std::vector<raw::ctb::ChStatus>&      raw::ctb::pdspctb::GetChStatuses() const { return fChStatuses; }
const std::vector<raw::ctb::Feedback>&      raw::ctb::pdspctb::GetFeedbacks()  const { return fFeedbacks; }
const std::vector<raw::ctb::Misc>&          raw::ctb::pdspctb::GetMiscs()      const { return fMiscs; }
const std::vector<raw::ctb::WordIndex>&     raw::ctb::pdspctb::GetIndexes()    const { return fIndexes; }

size_t  raw::ctb::pdspctb::GetNTriggers()   const { return fTriggers.size(); }
size_t  raw::ctb::pdspctb::GetNChStatuses() const { return fChStatuses.size(); }
size_t  raw::ctb::pdspctb::GetNFeedbacks()  const { return fFeedbacks.size(); }
size_t  raw::ctb::pdspctb::GetNMiscs()      const { return fMiscs.size(); }
size_t  raw::ctb::pdspctb::GetNIndexes()    const { return fIndexes.size(); }

const raw::ctb::Trigger&     raw::ctb::pdspctb::GetTrigger(size_t i)   const { return fTriggers.at(i); }
const raw::ctb::ChStatus&    raw::ctb::pdspctb::GetChStatuse(size_t i) const { return fChStatuses.at(i); }
const raw::ctb::Feedback&    raw::ctb::pdspctb::GetFeedback(size_t i)  const { return fFeedbacks.at(i); }
const raw::ctb::Misc&        raw::ctb::pdspctb::GetMisc(size_t i)      const { return fMiscs.at(i); }
const raw::ctb::WordIndex&   raw::ctb::pdspctb::GetIndex(size_t i)     const { return fIndexes.at(i); }


#endif // pdspctb_H
