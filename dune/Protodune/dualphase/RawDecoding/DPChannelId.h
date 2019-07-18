////////////////////////////////////////////////////////////////////////
// Class:   DPChannelId
//          Encapsulates various indecies for channel mapping for DP
//          detector
//
// Created: vgalymov, Thu Jul 18 10:36:39 CEST 2019
//
////////////////////////////////////////////////////////////////////////

#ifndef __DP_CHANNEL_ID__
#define __DP_CHANNEL_ID__

//
namespace dune
{
  class DPChannelId
  {
  private:
    unsigned seqn_;
    unsigned short crate_, card_, cardch_;
    unsigned short crp_, view_, viewch_;
    //bool exists_;
    unsigned short state_;
    
  public:
    DPChannelId() {;}
  DPChannelId( unsigned seqn, 
		 unsigned short crate, unsigned short card, unsigned short cardch, 
		 unsigned short crp, unsigned short view, unsigned short viewch, 
	       unsigned short state = 0 ) : 
    seqn_(seqn), crate_(crate), card_(card), cardch_(cardch), crp_(crp), view_(view), viewch_(viewch), state_(state) {}
    
    bool operator<(const DPChannelId &rhs) const { return seqn_ < rhs.seqn_; }
    
    const unsigned seqn() const { return seqn_; }
    const unsigned short crate() const { return crate_; }
    const unsigned short card() const { return card_; }
    const unsigned short cardch() const { return cardch_; }
    const unsigned short crp() const { return crp_; }
    const unsigned short view() const { return view_; }
    const unsigned short viewch() const { return viewch_; }
    const unsigned short state() const { return state_; }
    const bool exists() const { return (state_ == 0); }
  };
}


#endif
