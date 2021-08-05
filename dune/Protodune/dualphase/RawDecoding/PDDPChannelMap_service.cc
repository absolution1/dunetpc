////////////////////////////////////////////////////////////////////////
// Class:       PDDPChannelMap
// Plugin Type: service (art v3_02_06)
// File:        PDDPChannelMap_service.cc
//
// Generated at Thu Jul 18 10:29:03 2019 by Vyacheslav Galymov using cetskelgen
// from cetlib version v3_07_02.
//
// 
// Classes to facility channel order translation between different 
// representations
// 
// seqn   - unique counter to give channel data position in the record
// crate  - utca crate number starting
// card   - AMC card number
// cardch - AMC card channel number 
// crp    - CRP number
// view   - view number (0/1)
// viewch - view channel number 
// 
// Boost multi_index_container provides interface to search and order various 
// indicies. The container structure is defined in DPChannelTable, which has
// the following interfaces:
//  - raw sequence index, tag IndexRawSeqn
//  - crate number, tag IndexCrate, to get all channels to a given crate
//  - crate number and card number, tag IndexCrateCard, 
//    to get all channels for a given card in a crate
//  - crate, card, and channel number, tag IndexCrateCardChan, to access
//    a given channel of a given card in a given crate
//  - CRP index, tag IndexCrp, to access channels assigned to a given CRP
//  - CRP index and view index, tag IndexCrpView, to access channels assigned to 
//    a given view in a specified CRP
//  - CRP index, view index, and channel number, tag IndexCrpViewChan, to access 
//    a given view channel in a given CRP
// 
////////////////////////////////////////////////////////////////////////

#include <boost/range.hpp>
#include <boost/range/adaptors.hpp>

#include <iostream>
#include <iomanip>
#include <utility>

#include "PDDPChannelMap.h"

#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"

using namespace dune;

//class PDDPChannelMap;



//
// ctor
PDDPChannelMap::PDDPChannelMap(fhicl::ParameterSet const& p, art::ActivityRegistry& areg)
{
  mapname_  = "";
  ntot_     = 0;
  ncrates_  = 0;
  ncrps_    = 0;
  nch_      = 64;

  std::string MapName  = p.get<std::string>("MapName", "pddp2crp");
  // for test maps one can also specify Nb of crates, cards per crate, and views
  unsigned ncrateInMap = p.get<unsigned>("MapCrateNb", 1);
  unsigned ncardsInMap = p.get<unsigned>("MapCardNb", 10);
  unsigned nviewsInMap = p.get<unsigned>("MapViewNb",  1);
  
  //initialize channel map
  initMap( MapName, ncrateInMap, ncardsInMap, nviewsInMap );
}


// 
// initMap
void PDDPChannelMap::initMap(  std::string mapname, unsigned ncrates,
		       unsigned ncards, unsigned nviews )
{
  // already defined?
  if( mapname_.compare( mapname ) == 0 )
    {
      //cout<<"nothing to do"<<endl;
      return;
    }
  
  clearMap();
  mapname_ = mapname;
  if( mapname.compare("pddp2crp") == 0 ) 
    pddp2crpMap();
  else if( mapname.compare("pddp4crp") == 0 ) 
    pddp4crpMap();
  else 
    simpleMap( ncrates, ncards, nviews );
}


//
// clearMap 
void PDDPChannelMap::clearMap()
{
  DPChannelTable().swap( chanTable );
  ncrates_ = 0;
  ncrps_   = 0;
  ntot_    = 0;
  mapname_ = "";
  crateidx_.clear();
  crpidx_.clear();
}


//
// a simple channel map for testing purposes
void PDDPChannelMap::simpleMap( unsigned ncrates, unsigned ncards, unsigned nviews )
{
  unsigned nctot  = ncards * ncrates; // total number of cards
  unsigned ncview = nctot / nviews;   // allocate the same for each view
  unsigned nch    = nch_;             // number of ch per card (fixed)
  
  unsigned seqn  = 0;
  unsigned crate = 0;
  unsigned crp   = 0;
  unsigned view  = 0;
  unsigned vch   = 0;
  for( unsigned card = 0; card < nctot; card++ )
    {
      if( card > 0 ) 
	{
	  if( card % ncards == 0 ) crate++;
	  if( card % ncview == 0 ) {view++; vch=0;}
	}
      for( unsigned ch = 0; ch < nch; ch++ )
	{
	  add( seqn++, crate, card % ncards, ch, crp, view, vch++);
	  //cout<<seqn<<" "<<crate<<" "<<card%ncards
	}
    }
  //
}


//
// The channel map for 2 CRP configuration
void PDDPChannelMap::pddp2crpMap()
{
  // all idices run from 0
  std::vector<unsigned> cards_per_crate_real{5, 5, 5, 10, 5, 5, 10, 10, 10, 5};
  unsigned ncrates = 12;  // only 10 are active
  std::vector<unsigned> cards_per_crate(ncrates, 10);
  unsigned nch  = 64;
  unsigned nchc = nch/2;  // logical data channels per KEL/VHDCI connector 
  unsigned ncrp = 4;      // only 2 are active
  std::vector<unsigned> crpv(2*ncrp, 0);
  
  // map to TPC # in larsoft 
  std::vector<unsigned> dunetpcno = {3, 1, 0, 2};
  
  // all connector mappings should include ADC channel inversion on AMC 
  // this inversion is in groups of 8ch, i.e., AMC ch 0 -> 7 should be remapped to 7 -> 0
  // the connector mappings below should be (re)generated with the python script card2crp.py
  
  // kel connector orientation in a given view, chans 0 -> 31
  std::vector<unsigned> kel_nor = { 7,  6,  5,  4,  3,  2,  1,  0, 15, 14, 13, 12, 11, 10,  9,  8, 23, 22, 21, 20, 19, 18, 17, 16, 31, 30, 29, 28, 27, 26, 25, 24 };
  // kel connector orientation in a given view, chans 31 -> 0
  std::vector<unsigned> kel_inv = {24, 25, 26, 27, 28, 29, 30, 31, 16, 17, 18, 19, 20, 21, 22, 23,  8,  9, 10, 11, 12, 13, 14, 15,  0,  1,  2,  3,  4,  5,  6,  7 };
  
  //
  unsigned seqn = 0;
  unsigned crp, view;
  for(unsigned crate = 0;crate<ncrates;crate++)
    {
      // number of cards in this crate
      unsigned ncards = cards_per_crate[ crate ];
      for( unsigned card = 0;card<ncards;card++ )
	{
	  // which view ?
	  if( crate < 6 ) view = 0;
	  else view = 1;
	  
	  // which CRP ?
	  if( crate < 3 ) 
	    {
	      crp = 1;
	      if( card >= 5 ) crp = 2;
	    }
	  else if( crate >= 3 && crate < 6 )
	    {
	      crp = 0;
	      if( card >= 5 ) crp = 3;
	    }
	  else if( crate >=6 && crate < 9 )
	    {
	      crp  = 0;
	      if( card >= 5 ) crp = 1;
	    }
	  else
	    {
	      crp = 3;
	      if( card >= 5 ) crp = 2;
	    }
	  
	  
	  // get iterator to view channel numbering
	  auto kel = kel_inv; // inverted order
	  bool topAmcFirst = true;
	  if( view == 0 )
	    {
	      kel = kel_nor; //normal order
	      topAmcFirst = false;
	    }

	  auto vchIt      = crpv.begin() + view * ncrp + crp;

	  if( topAmcFirst ) // bot connector goes second
	    *vchIt = *vchIt + nchc;
	  
	  unsigned cstart = 0; // bottom AMC connector
	  for( unsigned ch = cstart; ch<cstart+nchc; ch++ )
	    {
	      int idx = ch - cstart;
	      if( idx < 0 || idx >= (int)kel.size() )
		{
		  // should not happen
		  std::cerr<<"I screwed this up\n";
		  continue;
		}
	      // view channel
	      unsigned vch = *vchIt + kel[ idx ];
	      if( view == 1 )
		{
		  int ival = -vch + 959;
		  if( ival < 0 )
		    {
		      std::cerr<<"Bad view channel number\n";
		      continue;
		    }
		  vch = (unsigned)ival;
		}
	      
	      // check if the channel actually exits
	      // since l1evb writes in the same format
	      unsigned state = 0;
	      if( crate >= cards_per_crate_real.size() )
		{
		  state = 1;
		}
	      else
		{
		  if( card >= cards_per_crate_real[ crate ] )
		    state = 1;
		}

	      //add( seqn++, crate, card, ch, crp, view, vch, state );
	      unsigned tpc = dunetpcno[ crp ];
	      unsigned kview = (view == 0 ) ? 1:0;
	      add( seqn++, crate, card, ch, tpc, kview, vch, state );
	    }

	  if( topAmcFirst )  // first half: top  connector is first
	    *vchIt = *vchIt - nchc;
	  else
	    *vchIt = *vchIt + nchc;
	  
	  // move to the second connector
	  cstart = nchc; // top amc connector
	  for( unsigned ch = cstart; ch<cstart+nchc; ch++ )
	    {
	      int idx = ch - cstart;
	      if( idx < 0 || idx >= (int)kel.size() )
		{
		  // should not happen
		  std::cerr<<"I screwed this up\n";
		  continue;
		}
	      // view channel
	      unsigned vch = *vchIt + kel[ idx ];
	      if( view == 1 )
		{
		  int ival = -vch + 959;
		  if( ival < 0 )
		    {
		      std::cerr<<"Bad view channel number\n";
		      continue;
		    }
		  vch = (unsigned)ival;
		}
	      
	      // check if the channel actually exits
	      // since l1evb writes in the same format
	      unsigned state = 0;
	      if( crate >= cards_per_crate_real.size() )
		{
		  state = 1;
		}
	      else
		{
		  if( card >= cards_per_crate_real[ crate ] )
		    state = 1;
		}
	      
	      //add( seqn++, crate, card, ch, crp, view, vch, state );
	      unsigned tpc = dunetpcno[ crp ];
	      unsigned kview = (view == 0 ) ? 1:0;
	      add( seqn++, crate, card, ch, tpc, kview, vch, state );
	      
	    }
	  
	  // all done
	  if( topAmcFirst )
	    *vchIt = *vchIt + nch;
	  else
	    *vchIt = *vchIt + nchc;
	  
	  //
	}
    }
}

//
// The channel map for 4 CRP configuration 
// DO NOT USE until this configuration actually exists
void PDDPChannelMap::pddp4crpMap()
{
  // NEEDS TO BE FIXED

  // all idices run from 0
  unsigned ncrates = 12;
  std::vector<unsigned> cards_per_crate(ncrates, 10);
  unsigned nch  = 64;
  unsigned ncrp = 4;
  std::vector<unsigned> crpv(2*ncrp, 0);

  //
  unsigned seqn = 0;
  unsigned crp, view;
  for(unsigned crate = 0;crate<ncrates;crate++)
    {
      // number of cards in this crate
      unsigned ncards = cards_per_crate[ crate ];
      for( unsigned card = 0;card<ncards;card++ )
	{
	  // which view ?
	  if( crate < 6 ) view = 0;
	  else view = 1;
	  
	  // which CRP ?
	  if( crate < 3 ) 
	    {
	      crp = 1;
	      if( card >= 5 ) crp = 2;
	    }
	  else if( crate >= 3 && crate < 6 )
	    {
	      crp = 0;
	      if( card >= 5 ) crp = 3;
	    }
	  else if( crate >=6 && crate < 9 )
	    {
	      crp  = 0;
	      if( card >= 5 ) crp = 1;
	    }
	  else
	    {
	      crp = 3;
	      if( card >= 5 ) crp = 2;
	    }
	  
	  
	  // get iterator to view channel numbering
	  auto vchIt = crpv.begin() + view * ncrp + crp;
	  
	  int apara = 1;
	  int bpara = 0;
	  if( view == 0 )
	    {
	      // order 31 -> 0 for these KEL connectors
	      apara = -1;
	      bpara = 31;
	    }
	  // card channels
	  for( unsigned ch = 0; ch<nch; ch++ )
	    {
	      if( ch == nch/2 ) *vchIt = *vchIt + nch/2;
	      int tmp = apara * (ch % 32) + bpara;
	      if( tmp < 0 ) 
		{
		  std::cerr<<"oh oh I screwed up\n";
		  continue;
		}
	      
	      unsigned vch = *vchIt + tmp;
	      add( seqn++, crate, card, ch, crp, view, vch );
	    }
	  // add the second connector
	  *vchIt = *vchIt + nch/2;

	  //
	}
    }
}

//
// add channl ID to map
void PDDPChannelMap::add( unsigned seq, unsigned crate, unsigned card, unsigned cch,
			  unsigned crp, unsigned view, unsigned vch, unsigned short state )
{
  chanTable.insert( DPChannelId(seq, crate, card, cch, crp, view, vch, state) );
  //
  ntot_    = chanTable.size();

  //if( state == 0 ) //existing channels only
      //{
  crateidx_.insert( crate );
  ncrates_ = crateidx_.size();
      
  crpidx_.insert( crp );
  ncrps_ = crpidx_.size();
  //}
  
  //
  //ncrates_ = cdistinct( chanTable.get<IndexCrateCardChan>(), 
  //[](const DPChannelId& c){ return c.crate();} );
  
  //
  //ncrps_   = cdistinct( chanTable.get<IndexCrpViewChan>(), 
  //[](const DPChannelId& c){ return c.crp();} );
}

//
// 
//
boost::optional<DPChannelId> PDDPChannelMap::find_by_seqn( unsigned seqn ) const
{
  auto it = chanTable.get<IndexRawSeqnHash>().find( seqn );
  if( it != chanTable.get<IndexRawSeqnHash>().end() )
    return *it;
  
  return boost::optional<DPChannelId>();
}

//
// the most low level info
std::vector<DPChannelId> PDDPChannelMap::find_by_seqn( unsigned from, unsigned to ) const
{
  if( to < from ) std::swap( from, to );
  
  std::vector<DPChannelId> res;
  
  if( from == to )
    {
      if( boost::optional<DPChannelId> id = find_by_seqn( from ) ) 
	res.push_back( *id );
      //auto it = chanTable.find(from);
      //if( it != chanTable.end() )
      //res.push_back( *it );
    }
  else
    {
      auto first = chanTable.get<IndexRawSeqn>().lower_bound( from );
      auto last  = chanTable.get<IndexRawSeqn>().upper_bound( to );
      res.insert( res.begin(), first, last );
    }
  
  return res;
}

//
//
std::vector<DPChannelId> PDDPChannelMap::find_by_crate( unsigned crate, bool ordered ) const
{
  if( not ordered ) // get from hashed index
    {
      const auto r = chanTable.get<IndexCrate>().equal_range( crate );
      std::vector<DPChannelId> res(r.first, r.second);
      return res;
    }

  const auto r = chanTable.get<IndexCrateCardChan>().equal_range( boost::make_tuple(crate) );
  std::vector<DPChannelId> res(r.first, r.second);

  return res;
}

//
//
std::vector<DPChannelId> PDDPChannelMap::find_by_crate_card( unsigned crate, unsigned card, bool ordered ) const
{
  if( not ordered ) // get from hashed index
    {
      const auto r = chanTable.get<IndexCrateCard>().equal_range( boost::make_tuple(crate, card) );
      std::vector<DPChannelId> res(r.first, r.second);
      return res;
    }
  
  // ordered accodring to channel number
  const auto r = chanTable.get<IndexCrateCardChan>().equal_range( boost::make_tuple( crate, card) );
  std::vector<DPChannelId> res(r.first, r.second);
  
  return res;
}

//
//
boost::optional<DPChannelId> PDDPChannelMap::find_by_crate_card_chan( unsigned crate,
								unsigned card, unsigned chan ) const
{
  auto it = chanTable.get<IndexCrateCardChanHash>().find( boost::make_tuple(crate, card, chan) );
  if( it != chanTable.get<IndexCrateCardChanHash>().end() )
    return *it;
  
  return boost::optional<DPChannelId>();
}

//
//
std::vector<DPChannelId> PDDPChannelMap::find_by_crp( unsigned crp, bool ordered ) const
{
  if( not ordered ) // get from hashed index
    {
      const auto r = chanTable.get<IndexCrp>().equal_range( crp );
      std::vector<DPChannelId> res(r.first, r.second);
      return res;
    }

  const auto r = chanTable.get<IndexCrpViewChan>().equal_range( crp );
  std::vector<DPChannelId> res(r.first, r.second);
  //return res;
  return res;
}

//
//
std::vector<DPChannelId> PDDPChannelMap::find_by_crp_view( unsigned crp, unsigned view, bool ordered ) const
{
  if( not ordered ) // get from hashed index
    {
      const auto r = chanTable.get<IndexCrpView>().equal_range( boost::make_tuple(crp, view) );
      std::vector<DPChannelId> res(r.first, r.second);
      return res;
    }
  
  // ordered accodring to channel number
  const auto r = chanTable.get<IndexCrpViewChan>().equal_range( boost::make_tuple(crp, view) );
  std::vector<DPChannelId> res(r.first, r.second);
  return res;
}

//
//
boost::optional<DPChannelId> PDDPChannelMap::find_by_crp_view_chan( unsigned crp,
							       unsigned view, unsigned chan ) const
{
  auto it = chanTable.get<IndexCrpViewChanHash>().find( boost::make_tuple(crp, view, chan) );
  if( it != chanTable.get<IndexCrpViewChanHash>().end() )
    return *it;
  
  return boost::optional<DPChannelId>();
}
  
//
// Map to CRP channels
int PDDPChannelMap::MapToCRP(int seqch, int &crp, int &view, int &chv) const
{
  crp = view = chv = -1;
  if( boost::optional<DPChannelId> id = find_by_seqn( (unsigned)seqch ) ) 
    {
      if( !id->exists() ) return -1;
      crp  = id->crp();
      view = id->view();
      chv  = id->viewch();
      return 1;
    }
  
  return -1;
}

//
// Map to DAQ channel sequence
int PDDPChannelMap::MapToDAQ(int crp, int view, int chv, int &seqch) const
{
  seqch = -1;
  if( boost::optional<DPChannelId> id = find_by_crp_view_chan( (unsigned)crp, (unsigned)view, (unsigned)chv ) )
    {
      if( !id->exists() ) return -1;
      seqch = id->seqn();
      return 1;
    }
  return -1;
}

//
// number of cards assigned to a given crate
unsigned PDDPChannelMap::ncards( unsigned crate ) const
{
  // assumes it is sorted according to card number
  unsigned count  = 0;
  auto r = chanTable.get<IndexCrateCardChan>().equal_range( crate );
  ssize_t last    = -1;
  for( DPChannelId const &ch : boost::make_iterator_range( r ) )
    {
      if( !ch.exists() ) continue;
      unsigned val = ch.card();
      if( val != last ) 
	{
	  count++;
	  last = val;
	}
    }
  
  return count;
}

//
// number of views assigned to a given CRP
unsigned PDDPChannelMap::nviews( unsigned crp ) const 
{
  // assumes it is sorted according to crp number
  unsigned count  = 0;
  auto r = chanTable.get<IndexCrpViewChan>().equal_range( crp );
  ssize_t last    = -1;
  for( DPChannelId const &ch : boost::make_iterator_range( r ) )
    {
      if( !ch.exists() ) continue;
      unsigned val = ch.view();
      if( val != last ) 
	{
	  count++;
	  last = val;
	}
    }
  
  return count;
}

//
//
void PDDPChannelMap::print( std::vector<DPChannelId> &vec )
{
  for( auto it = vec.begin();it!=vec.end();++it )
    {
      unsigned seqn   = it->seqn();
      unsigned crate  = it->crate();
      unsigned card   = it->card();
      unsigned cch    = it->cardch();
      unsigned crp    = it->crp();
      unsigned view   = it->view();
      unsigned vch    = it->viewch();
      unsigned state  = it->state();
      bool     exists = it->exists();
      std::cout<<std::setw(7)<<seqn
	       <<std::setw(4)<<crate
	       <<std::setw(3)<<card
	       <<std::setw(3)<<cch
	       <<std::setw(3)<<crp
	       <<std::setw(2)<<view
	       <<std::setw(4)<<vch
	       <<std::setw(2)<<state
	       <<std::setw(2)<<exists<<std::endl;
    }
}


 
 DEFINE_ART_SERVICE(dune::PDDPChannelMap)
