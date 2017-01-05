/** ****************************************************************************
 * @file EventMixingSummary.h
 * @brief Definition of basic mixing information
 * @author wketchum@fnal.gov
 * 
 * ****************************************************************************/

#ifndef DATAOVERLAYPRODUCTS_EVENTMIXINGSUMMARY_DUNE_H
#define DATAOVERLAYPRODUCTS_EVENTMIXINGSUMMARY_DUNE_H

#include <stdint.h>

namespace dunemix {

  //This is just gonna be a stupid copy of the basic
  //info contained in the art event. Nothing special,
  //though more could be added as necessary
  class EventMixingSummary {

  public:
  EventMixingSummary():
    fEvent(0),fSubrun(0),fRun(0) {}
    
#ifndef __GCCXML__
  public:
      
    EventMixingSummary(uint32_t e,uint32_t s,uint32_t r)
      {
	fEvent  = e;
	fSubrun = s;
	fRun    = r;
      }    
    
    uint32_t Event()  const ;
    uint32_t SubRun() const ;
    uint32_t Run()    const ;

#endif // !__GCCXML__
  private:
    uint32_t fEvent;
    uint32_t fSubrun;
    uint32_t fRun;
  }; // class EventMixingSummary()
  
#ifndef __GCCXML__
  inline uint32_t dunemix::EventMixingSummary::Event()  const { return fEvent;  }
  inline uint32_t dunemix::EventMixingSummary::SubRun() const { return fSubrun; }
  inline uint32_t dunemix::EventMixingSummary::Run()    const { return fRun;    }
#endif // !__GCCXML__
  
} // namespace dunemix


#endif // DATAOVERLAYPRODUCTS_EVENTMIXINGSUMMARY_DUNE_H

////////////////////////////////////////////////////////////////////////
