//////////////////////////////////////////////////////////////////////////////////////
//
//  Decoder for the DAQ events for charg readout
//  The general idea is to try to make this object threadsave when possible
//
//  TODO: test compression
// 
//////////////////////////////////////////////////////////////////////////////////////

//#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstdlib>

#include "EventDecoder.h"
#include "LogMsg.h"
#include "Timer.h"

using namespace std;
using namespace dlardaq;

//
// ctor
//
EventDecoder::EventDecoder(size_t nch, size_t nsample)
{
  // this is needed to read compressed data
  m_nadc    = dlardaq::BitsADC;  // bits for ADC resolution
  m_nch     = nch;               // total number of channels
  m_nsample = nsample;           // number of samples in each channel

  m_RunHeadBuf.resize( dlardaq::RunHeadSz );
  m_EveHeadBuf.resize( dlardaq::EveHeadSz  );
  m_EndFootBuf.resize( dlardaq::FileFootSz );
  
  //
  m_bytes_left = 0;

  // data mutex initialization
  pthread_mutex_init(&m_data_mutex, NULL);
}


//
// dtor
//
EventDecoder::~EventDecoder()
{
  Close();
  pthread_mutex_destroy(&m_data_mutex);
}


//
// close input file
//
void EventDecoder::Close()
{
  // lock mutex
  lock( m_data_mutex );

  if(m_file.is_open())
    m_file.close();  
  
  m_totev = 0;
  
  // unlock mutex
  unlock( m_data_mutex );
}



//
//
//
ssize_t EventDecoder::Open(std::string finname)
{
  // attempt to close any previously opened files
  Close();
  
  lock( m_data_mutex );

  m_EveData.clear(); // not used when reading from file

  // clear event bookmarks
  m_events.clear();

  //
  m_file.open(finname.c_str(), ios::in | ios::binary);
  
  if( !m_file.is_open() )
    {
      msg_err<<"Could not open "<<finname<<" for reading"<<endl;
      unlock( m_data_mutex );
      return -1;
    }
  
  // read run header bytes
  ReadBytes( m_RunHeadBuf );
  dlardaq::decode_runhead(&m_RunHeadBuf[0], m_rnh);
  
  
  m_pstart = m_file.tellg();
  
  // read footer
  // fast forward to the end
  m_file.seekg(-dlardaq::FileFootSz, m_file.end );
  m_pend = m_file.tellg();
  
  ReadBytes( m_EndFootBuf );
  dlardaq::decode_filefoot(&m_EndFootBuf[0], m_flf);

  // rewind back to the begining
  m_file.seekg( m_pstart );

  
  // 

  m_totev = m_flf.num_events;
  
  if( m_totev == 0 )
    {
      msg_err<<"File "<<finname<<" is empty"<<endl;
      unlock(m_data_mutex);
      Close();
      return 0;
    }

  unlock( m_data_mutex );
  
  return m_totev;
}

//
// read a given number of bytes from file
//
void EventDecoder::ReadBytes( std::vector<BYTE> &bytes )
{
  m_file.read( &bytes[0], bytes.size() );
}


//
//
//
void EventDecoder::ReadEvent( std::vector<adc16_t> &adc, bool headonly )
{
  // first read the header
  ReadBytes( m_EveHeadBuf );
  decode_evehead(&m_EveHeadBuf[0], m_evh);  
  //dlardaq::msg_info<<">>> Event number "<<m_evh.ev_num<<endl;
  //dlardaq::msg_info<<">>> Event size   "<<m_evh.ev_size<<endl;
  //dlardaq::msg_info<<">>> Trig number  "<<m_evh.trig_info.num<<endl;
  //dlardaq::msg_info<<">>> Time stamp   "<<m_evh.trig_info.ts.tv_sec<<" s "
  //<<m_evh.trig_info.ts.tv_nsec<<" ns"<<endl;
  //dlardaq::msg_info<<">>> Flags "<<bitset<8>(m_evh.dq_flag)<<endl;

  
  if(headonly) // read only header and skip to next event
    {
      // get event size
      size_t evsz = m_evh.ev_size;

      // move to next event
      m_file.seekg(evsz, std::ios::cur);
      return;
    }

  //
  // decode event data
  if(!m_EveDataBuf.empty()) m_EveDataBuf.clear();
  
  m_EveDataBuf.resize(m_evh.ev_size, '\0');
  
  // read content of the file
  ReadBytes(m_EveDataBuf);
  Decode( &m_EveDataBuf[0], m_EveDataBuf.size(), GETDCFLAG(m_evh.dq_flag), adc );
}

//
// get a specific event from file
// 
ssize_t EventDecoder::GetEvent(size_t evnum, dlardaq::evheader_t &eh, 
			       std::vector<adc16_t> &adc) 
{
  adc.clear();
  if(!m_file.is_open()) return -1;
  //if(evnum >= m_totev)  return -1; //no-can-do
  lock(m_data_mutex);

  if(evnum < m_events.size())  // already know where it is
    {
      m_file.seekg(m_events[evnum]);
      ReadEvent( adc );
      eh = m_evh;
    }
  else
    {
      size_t curev = m_events.size();
      while( curev <= evnum )
	{
	  // our event begins at this position
	  streampos pos = m_file.tellg();
	  m_events.push_back( pos );
	  
	  // readonly header to get event data size
	  if( curev < evnum ) 
	    ReadEvent(adc, true);
	  else // event we want
	    ReadEvent(adc, false);
	  
	  curev++;
	}
      
      // last event is the one with want ...
      eh = m_evh;
    }

  unlock(m_data_mutex);

  return evnum;
}

//
//
//
ssize_t EventDecoder::GetEvent( dlardaq::evheader_t &eh,
				std::vector<adc16_t> &adc )
{
  if(m_file.is_open())
    {
      // return first event
      return GetEvent( 0, eh, adc);
    }

  ssize_t rval = -1;
  //else // we have an event buffer
  if( m_EveData.empty() ) return -1;

  // lock mutex
  lock(m_data_mutex);

  eh  = m_evh;  
  adc = m_EveData;

  // clear the internal data buffer
  m_EveData.clear();
  
  rval = eh.ev_num;
  
  // unlock mutex
  unlock(m_data_mutex);
  return rval;
}

//
//
//
ssize_t EventDecoder::Decode( const char *buf, size_t nb, bool cflag, 
			      std::vector<adc16_t> &adc )
{
  adc.clear();

  if( !cflag ) //not compressed
    {
      adc.resize( (nb*8)/m_nadc, 0 );
      unpack12into16( buf, &adc[0], nb );
    }
  else
    {
      size_t byteidx = 0;
      HuffDataCompressor::Instance().DecompressEventData(m_nadc, m_nch, m_nsample, buf, nb, byteidx, adc);
      
    }  

  return adc.size();
}


//
// we assume that the buffer is 
// run header + event header + event data but it can come in different packets
//
void EventDecoder::ReadBuffer(const char *buf, size_t nb)
{
  if(m_file.is_open())
    {
      msg_err<<"Cannot use this function while reading data from a file"<<endl;
      return;
    }
  
  // lock mutex
  lock(m_data_mutex);
  m_totev = 0;  
  if(!m_EveData.empty()) m_EveData.clear();
  
  // we have finished reading previous event
  if( m_bytes_left <= 0 || 
      Timer::GetTimer().splittime(false, false) > TMAXWAIT)
    {
      m_EveDataBuf.clear();
      if(IsFirstPacket(buf, nb)) //check we have a correct packet
	{
	  // decode header
	  size_t  nb_read = 0;
	  const char *pdata = buf;
	  
	  // decode run header
	  ssize_t ret;
	  ret = dlardaq::decode_runhead(pdata, m_rnh);
	  if( ret <= 0) 
	    {
	      // unlock mutex
	      unlock(m_data_mutex);
	      return;
	    }
	  pdata   += ret;
	  nb_read += ret;

	  // decode event header
	  ret = dlardaq::decode_evehead(pdata, m_evh);
	  if( ret <= 0) 
	    {
	      // unlock mutex
	      unlock(m_data_mutex);
	      return;
	    }
	  pdata   += ret;
	  nb_read += ret;
	  
	  m_bytes_left = m_evh.ev_size;
	  if(m_bytes_left <= 0) 
	    {
	      // unlock mutex
	      unlock(m_data_mutex);
	      return;
	    }
	  
	  //msg_info<<"Start recieving event "<<m_evh.ev_num<<endl;

	  // copy the data
	  m_EveDataBuf.insert(m_EveDataBuf.end(), pdata, pdata + nb-nb_read);
	  m_bytes_left -= (nb - nb_read);
	  	  
	  // start timer 
	  Timer::GetTimer().start();

	  //
	  unlock(m_data_mutex);
	  return;
	}
      else
	{
	  unlock(m_data_mutex);
	  return;
	}
    }
  
  // append to buffer
  m_EveDataBuf.insert(m_EveDataBuf.end(), buf, buf+nb);
  m_bytes_left -= nb;
  if( m_bytes_left < 0 )
    {
      msg_err<<"Byte packet mismatch"<<endl;
      
      // unlock mutex
      unlock(m_data_mutex);
      return;
    }

  // we have full event
  if( m_bytes_left == 0 )
    {
      // decode data
      bool cflag =  GETDCFLAG(m_evh.dq_flag);
      Decode( &m_EveDataBuf[0], m_EveDataBuf.size(), cflag, m_EveData );
      m_totev = 1;
      
      msg_info<<"Recieved: run "<<m_rnh.run_num<<" event "<<m_evh.ev_num<<endl;
    }
    
  // unlock mutex
  unlock(m_data_mutex);
}


//
// detect new event key after run header
//
bool EventDecoder::IsFirstPacket(const char *buf, size_t nb)
{
  BYTE k0 = buf[dlardaq::RunHeadSz];
  BYTE k1 = buf[dlardaq::RunHeadSz+1];
  bool rval = ( ((k0 & 0xFF) == EVSKEY) && ((k1 & 0xFF) == EVSKEY) );
  
  return rval;
}
