/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
//
//  Class for compressing ADC data internally represented
//  as 16 bit unsigned int 
//  Max ADC resolution, however, cannot be greater than 15bits as
//  1bit is reserved to signal whether the packet is compressed or not 
//
//  The Huffman encoding scheme used here is the one developed by uBooNE
//
//  The basic idea is to handle binnary codes as strings --> maybe not 
//  the most performant, but maybe easier to understand
//
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstdlib>

#include "LogMsg.h"
#include "HuffDataCompressor.h"

using namespace std;
using namespace dlardaq;

//
//
//
//
// ctor (private)
HuffDataCompressor::HuffDataCompressor()
{
  // max ADC supported

  // note that we will use "short" to store either ADC or diff
  // for 16 bit short goes between -/+ 32767 or 15 bits
  m_MaxAdcBits = 15; //+1 bit is reserved for compression flag
  
  // number of reserved bits
  m_NbitsHead  = 1; 
  
  // we use 8 bit words on the output
  m_NbitsByte  = 8;  //sizeof(BYTE)*8;
  
  m_Verbosity  = 0;

  // define encoding map
  SetEncoding();
}


//
// The Huffman encoding scheme 
// Needless to say it should not change 
// between comression/decompression steps
//
// set encoding system
void HuffDataCompressor::SetEncoding()
{

  // basic huffman encoding scheme
  // same as used by uBooNE
  std::map<short, std::string> codeMap;
  codeMap[0]  = "1";
  codeMap[-1] = "01";
  codeMap[1]  = "001";
  codeMap[-2] = "0001";
  codeMap[2]  = "00001";
  codeMap[-3] = "000001";
  codeMap[3]  = "0000001";

  // max difference 
  m_MaxDiff     = 3;

  // 0 no repetion encoding >MaxDiff repetion encoding
  m_NSeqRep    = 4; 
  m_NSeqRepVal = 0;
  m_SeqEnable  = (m_NSeqRep > 0); // && (m_NSeqRep > m_MaxDiff);

  m_CmMap.clear();
  m_UCmMap.clear();
  
  //
  m_MinCodeSize  = 999;
  m_MaxCodeSize  = 0;
  
  std::map<short, std::string>::iterator it;
  for(it = codeMap.begin();it!=codeMap.end();it++)
    {
      
      short deltaval = it->first;
      string bincode = it->second;
      
      if(m_SeqEnable) 
	{
	  // add zeros to nominal huffman codes
	  bincode = "0" + bincode;
	}
      
      size_t codesize = (bincode).size();
      if(codesize > m_MaxCodeSize)
	m_MaxCodeSize = codesize;
      if(codesize < m_MinCodeSize)
	m_MinCodeSize = codesize;
      
      // set code map
      m_CmMap[deltaval] = bincode;
    }

  // run-length encoding enabled
  if(m_SeqEnable)
    {
      m_NSeqRepVal          = m_MaxDiff + m_NSeqRep;
      m_CmMap[m_NSeqRepVal] = "1";
      //m_UCmMap["1"]         = m_NSeqRepVal;
      m_MinCodeSize         = 1;
    }

  // !!! m_MinCodeSize should be 1 !!!
  m_UCmMap.resize( m_MaxCodeSize );
  for(it = m_CmMap.begin();it!=m_CmMap.end();it++)
    {
      short deltaval = it->first;
      string bincode = it->second;
      m_UCmMap[bincode.size()-1] = std::make_pair(bincode, deltaval);
    }
}

//
//
//
// Print encoding
void HuffDataCompressor::PrintEncoding()
{
  cout<<endl<<"Huffman coding scheme: "<<endl;

  
  for(size_t i=0;i<m_UCmMap.size();i++)
    {
      if(i == 0 && m_SeqEnable)
	{
	  cout<<setw(10)<<m_UCmMap[i].first
	      <<setw(5)<<m_NSeqRep<<" x no change"
	      <<setw(4)<<i+1<<" bit code length"<<endl;
	  
	}
      else
	{
	  cout<<setw(10)<<m_UCmMap[i].first<<setw(5)<<m_UCmMap[i].second
	      <<setw(4)<<i+1<<" bit code length"<<endl;
	}
      
    }

  cout<<endl;
}




//
//
//
// access map to compress
std::string HuffDataCompressor::GetCodeFromValue(short val)
{
  string emptystr = "";
  
  if(std::abs(val) > m_MaxDiff && val != m_NSeqRepVal)
    {
      msg_err<<__FILE__<<", "<<__LINE__<<" the value '"<<val<<"' exceeds range"<<endl;
      return emptystr;
    }

  std::map<short, std::string>::iterator it;
  it = m_CmMap.find( val );
  if(it != m_CmMap.end())
    return it->second;
  
  //
  msg_err<<__FILE__<<", "<<__LINE__<<" the value '"<<val<<"' could not be found"<<endl;
  return emptystr;
}


//
// Internal buffer to store raw adc / adc differences
//
// Add to huff buffer
void HuffDataCompressor::AddToHuffBuffer(short val, bool rawadc)
{
  HuffAccum_t tmp;
  tmp.value   = val;
  tmp.nrep    = 0;
  tmp.codestr = "";

  if(!rawadc) 
    tmp.codestr = GetCodeFromValue(val);
  else
    {
      HuffBuffer.push_back( tmp );
      return;
    }
  
  // check if the last value was the same
  if( HuffBuffer.back().value == tmp.value && 
      !HuffBuffer.back().codestr.empty() )
    {
      // just increment repetition counter
      HuffBuffer.back().nrep += 1;
    }
  else
    {
      // add new entry
      HuffBuffer.push_back( tmp );
    }
}

//
//
//
// get binary string
std::string HuffDataCompressor::GetBinaryString(short val, size_t strsize)
{
  std::string s = bitset<16> ( val ).to_string();

  if(strsize < 16)
    return s.substr(s.size() - strsize); //crop to required size
    
  // ok we have a long string. Why?
  msg_warn<<__FILE__<<", "<<__LINE__<<" the string looks too long"<<endl;
  while(s.size() < strsize) s = "0" + s; //pad 0 to front
  return s;
}

//
//
//
// set number of ADC bits
bool HuffDataCompressor::SetNbitsAdc( short nbadc )
{
  if( nbadc > m_MaxAdcBits ) return false;
  
  // this is max bits that can be occupied by data
  m_NbitsHC    = nbadc;
  
  // our basic packet size (should not exceed m_MaxAdcBits)
  m_PacketSize = m_NbitsHC + m_NbitsHead;

  return true;
}

//
//
//
// add words to byte buffer keeping track of the remainder
// partbyte ("partial byte") contains the remainder string
void HuffDataCompressor::AddWordsToByteBuffer(std::string words, 
					      std::vector<BYTE> &buf, 
					      std::string &partbyte)
{
  if( words.empty() && partbyte.empty() ) return;

  // nothing to add yet
  if( words.size() < m_NbitsByte && partbyte.empty() )
    {
      partbyte = words;
      return;
    }
  
  // just pad with 0 and finish
  if(words.empty() && !partbyte.empty())
    {
      // pad partial byte until we get to a full byte
      while(partbyte.size() < m_NbitsByte) partbyte += "0";
      buf.push_back( (strtoul(partbyte.c_str(), 0, 2) & 0xff) );
      partbyte.clear();
      return;
    }
  
  for(size_t i=0;i<words.size();i++)
    {
      if(partbyte.size() == m_NbitsByte)
	{
	  buf.push_back( (strtoul(partbyte.c_str(), 0, 2) & 0xff) );
	  partbyte.clear();
	}
      
      partbyte += words[i];
    }
}



//
//
//
// compress ch data using huffman codes
void HuffDataCompressor::CompressChData( short nbadc, std::vector<adc16_t> &raw_in,
					 std::vector<BYTE> &bin_out )
{
  if(!SetNbitsAdc( nbadc ))
    {
      bin_out.clear();
      msg_err<<"ADC "<<nbadc<<" bits exceeds max supported "<<m_MaxAdcBits<<endl;
      return;
    }
  
  
  HuffBuffer.clear();
  
  // calculate differences and fill our HC buffer
  short delta;
  for(size_t i=0;i<raw_in.size();i++)
    {
      if(i==0) delta = m_MaxDiff + 1;
      else delta = raw_in[i] - raw_in[i-1];
      
      if(std::abs(delta) > m_MaxDiff)  // store uncompressed raw data
	AddToHuffBuffer( (short)raw_in[i], true );
      else
	AddToHuffBuffer( delta, false );
    }
  
  
  // not defined unless SeqEnable is set
  string bitrepcode = "";
  if(m_SeqEnable) 
    {
      bitrepcode = GetCodeFromValue( m_NSeqRepVal ); //NSeqRep );
      if(bitrepcode.empty())
	{
	  msg_err<<__FILE__<<", "<<__LINE__
		 <<" encoding repetion failed. Aborting ..."<<endl;
	  return;
	}
    }  
  
  // now compress stream
  short nhcrep;  //repetition 
  string bitword, bitcode, partbyte;
  bitword  = "";
  partbyte = ""; // string to hold ramainder of byte

  for(size_t i=0;i<HuffBuffer.size();i++)
    {
      
      //cout<<i
      //<<" "<<HuffBuffer[i].value
      //<<" "<<HuffBuffer[i].nrep;
      //if( !HuffBuffer[i].codestr.empty() ) cout<<" "<<HuffBuffer[i].codestr<<endl;
      //else cout<<" rawadc "<<endl;
	
      if(HuffBuffer[i].codestr.empty()) //raw adc
	{
	  // write out whatever we have in the buffer
	  if(bitword.size() > m_NbitsHead) // more than just a header
	    {
	      // pad with zero's to buffer length
	      while(bitword.size() < m_PacketSize) bitword += "0"; 
	      AddWordsToByteBuffer( bitword, bin_out, partbyte );
	    }
	  
	  // uncompressed data
	  bitword = "0" + GetBinaryString( HuffBuffer[i].value, m_NbitsHC );
	  
	  //
	  AddWordsToByteBuffer( bitword, bin_out, partbyte );
	  bitword = "1";
	  continue;
	}

      // compress data
      bitcode  = HuffBuffer[i].codestr;
      nhcrep   = HuffBuffer[i].nrep;

      for(short jj=0;jj<=nhcrep;jj++)
	{
	  //
	  if(jj==0)
	    bitword += bitcode;
	  else if(jj+m_NSeqRep <= nhcrep && m_SeqEnable)
	    {
	      bitword += bitrepcode;
	      jj += (m_NSeqRep-1);
	    }
	  else
	    bitword += bitcode;
	  
	  
	  // attempt to write to compresssed stream
	  if(bitword.size() >= m_PacketSize) // check the size
	    {
	      if(bitword.size() == m_PacketSize) //write it out
		{
		  AddWordsToByteBuffer( bitword, bin_out, partbyte );
		  bitword = "1"; //next iteration of compression
		}
	      else
		{
		  // too large for one packet 
		  string towrite = bitword.substr(0, m_PacketSize);
		  string tosave  = bitword.substr(m_PacketSize);
		  
		  // write out
		  AddWordsToByteBuffer( towrite, bin_out, partbyte );
		  
		  // carry over for next iteration
		  bitword = "1" + tosave; 
		}
	    }
	} // end for

    } // end i-loop
  
  //cout<<"Last bitword and partial byte"<<bitword<<" "<<partbyte<<endl;
  
  // write whatever is left
  if(bitword.size() > m_NbitsHead) // more than just a header
    AddWordsToByteBuffer( bitword, bin_out, partbyte );
  
  //cout<<"Last partial byte word "<<partbyte<<endl;
  
  // write out remaning byte padded with 0s to next byte boundary
  AddWordsToByteBuffer( "", bin_out, partbyte );
}


//
//
//
// compress event data using huffman codes
void HuffDataCompressor::CompressEventData( short nbadc, size_t nch, size_t seqlen,
					    std::vector<adc16_t> &raw_in,
					    std::vector<BYTE> &bin_out )
{
  bin_out.clear();
  if(!SetNbitsAdc( nbadc ))
    {
      msg_err<<"ADC "<<nbadc<<" bits exceeds max supported "<<m_MaxAdcBits<<endl;
      return;
    }

  if( seqlen * nch != raw_in.size() )
    {
      msg_err<<"Length of raw data vector does not match with expected "<<endl;
      return;
    }
  
  for(size_t i=0;i<nch;i++)
    {
      size_t istart = i*seqlen;
      // get ch data
      std::vector<adc16_t> chdata(&raw_in[istart], &raw_in[istart+seqlen]);
      CompressChData( nbadc, chdata, bin_out );
    }
}


//
//
// compress event data using huffman codes
void HuffDataCompressor::CompressEventData( short nbadc, size_t nch, size_t seqlen,
					    std::vector< std::vector<adc16_t> > &raw_in,
					    std::vector<BYTE> &bin_out )
{
  bin_out.clear();
  if(!SetNbitsAdc( nbadc ))
    {
      msg_err<<"ADC "<<nbadc<<" bits exceeds max supported "<<m_MaxAdcBits<<endl;
      return;
    }

  if( nch != raw_in.size() )
    {
      msg_err<<"Number of ch in raw data vector does not match with expected "<<endl;
      return;
    }
  
  for(size_t i=0;i<nch;i++)
    {
      if( raw_in[i].size() != seqlen ) 
	{
	  msg_err<<"No support for compression of unequal sequence lenghts at the moment"<<endl;
	  bin_out.clear();
	  return;
	}

      CompressChData( nbadc, raw_in[i], bin_out );
    }
}


//
//
//
// decode individual bits from a byte
void HuffDataCompressor::ReadNextByte(size_t &byteidx, const char *buf,
				      std::deque<bitset<1> > &bits)
{
  unsigned short bitmask = 1 << (m_NbitsByte - 1);
  
  // decode individual bits
  for(size_t i=0;i<m_NbitsByte;i++)
    bits.push_back( bitset<1>( (buf[byteidx] & (bitmask >> i)) != 0) );
  
  byteidx++;
}


//
//
//
// decode individual bits from a byte
void HuffDataCompressor::ReadNextByte(std::ifstream &fin, std::deque< std::bitset<1> > &bits, 
				      bool &status)
{
  status = true;
  char byteword;
  
  if(!fin.get(byteword))
    {
      status = false;
      return;
    }

  unsigned short bitmask = 1 << (m_NbitsByte - 1);
  
  // decode individual bits
  for(size_t i=0;i<m_NbitsByte;i++)
    bits.push_back( bitset<1>( (byteword & (bitmask >> i)) != 0) );
}


//
//
// 
// Decompress event from binary sequence (only single event)
//
void HuffDataCompressor::DecompressEventData( short nbadc,
					      size_t nch,
					      size_t seqlen, 
					      const char *buf, size_t bufsize, size_t &byteidx,
					      std::vector<adc16_t> &adc )
{
  if(!SetNbitsAdc( nbadc ))
    {
      msg_err<<"ADC "<<nbadc<<" bits exceeds max supported "<<m_MaxAdcBits<<endl;
      return;
    }
  adc.clear();

  vector< adc16_t > chdata;
  deque< bitset<1> > bitqueue;

  stringstream ss;
  short lastdelta  = -999;
  size_t bitsread  = 0;
  int noread       = 0;

  size_t chread  = 0;

  byteidx   = 0;
  // basic idea is to read bits one byte at a time from the input buffer
  // the first entry in the bitqueue deck should always be aligned on the 
  // compressed / raw flag and the deck should contain at least m_PacketSize bits
  while( adc.size() != nch*seqlen )
    {
      if(byteidx<bufsize) ReadNextByte(byteidx, buf, bitqueue);
      else noread++;
      
      if(noread > 1) //loopcounter == 1, we process last bits but not read new byte
	{
	  msg_err<<"There seems to be a problem with decoding"<<endl
		 <<"Bytes read "<<byteidx<<" out of "<<bufsize<<endl
		 <<"Samples accumulated in this channel "<<chdata.size()<<endl
		 <<"Remaining bits are "<<bitqueue.size()<<" "<<bitqueue.empty()<<endl;
	  break;
	}
      
      //read more bytes if possible
      if( bitqueue.size() < m_PacketSize && byteidx != bufsize) 
	continue;

      bool iscomp = bitqueue.front().test(0);
      bitsread++;
      bitqueue.pop_front(); // compressed / raw flag      

      if(!iscomp) // uncompressed
	{
	  ss.str(""); //clear our code string (could carry over from compressed)
	  
	  if( bitqueue.size() < m_NbitsHC ) 
	    {
	      msg_err<<"Fatal decoding error has been encountered : "<<endl
		     <<" Number of bits in the uncompressed stream should be at least "
		     <<m_NbitsHC<<" the current value is "<<bitqueue.size()<<endl;
	      abort();
	    }
	  
	  size_t bitcounter = 0;
	  while(bitcounter < m_NbitsHC) 
	    {
	      ss << bitqueue.front();
	      bitqueue.pop_front();
	      bitcounter++;
	    }
	  bitsread += bitcounter;

	  adc16_t val = ( strtoul(ss.str().c_str(), 0, 2) & 0x7FFF );
	  chdata.push_back( val );
	  ss.str("");
	}
      else // handle compressed bits
	{
	  size_t bitcounter = 0;
	  size_t bitstoread = m_NbitsHC;
	  if( bitqueue.size() < m_NbitsHC ) 
	    bitstoread = bitqueue.size();
	  
	  while( bitcounter < bitstoread )
	    {	      
	      ss << bitqueue.front();
	      bitqueue.pop_front();
	      bitcounter++;
	      bitsread++;
	      // check if we have a good code
	      size_t strsize = ss.str().size();
	      if( strsize >= m_MinCodeSize && strsize <= m_MaxCodeSize ) 
		{
		  bool ok = (m_UCmMap[strsize-1].first == ss.str());
		  if(!ok) continue; //not a good code
		  short val = m_UCmMap[strsize-1].second;
		  
		  ss.str(""); //reset string code

		  
		  // add to our adc vector
		  if( std::abs(val) <= m_MaxDiff)
		    {
		      adc16_t newval = chdata.back() + val;
		      chdata.push_back( newval );
		      lastdelta = val;
		    }
		  else if(val == m_NSeqRepVal)
		    {
		      for(short j=0;j<m_NSeqRep;j++)
			{
			  adc16_t newval = chdata.back() + lastdelta; //add same value
			  chdata.push_back( newval );
			}
		    }		
		}// end of if
	      
	      // we got the whole sequence now
	      if(chdata.size() == seqlen) break;

	    } // end while bitcounter
	} // end else 
      
      // check if we finished decoding data for a given channel
      if(chdata.size() == seqlen)
	{
	  
	  // calculate number of padded bits to the byte boundary
	  size_t padbits = bitsread % m_NbitsByte;
	  if(padbits > 0) padbits = m_NbitsByte - padbits;
	  
	  if(m_Verbosity > 1)
	    {
	      msg_info<<"Bytes read : "<<byteidx<<endl
		      <<"Bits read  : "<<bitsread<<endl
		      <<"Bits to boundary : "<<padbits<<endl
		      <<"Last value : "<<chdata.back()<<" ADC "<<endl;
	    }
	  
	  // error checking
	  if( bitqueue.size() < padbits )
	    {
	      msg_err<<"Fatal decoding error has been encountered : "<<endl
		     <<"Byte boundary does not appear to be valid"<<endl
		     <<"Check that the codes are padded with 0 to the next byte boundary"<<endl;
	      abort();
	    }
	  
	  // remove padded bits
	  while(padbits > 0)
	    {
	      bitqueue.pop_front();
	      padbits--;
	    }
	  
	  // some basic check
	  if(!bitqueue.empty())
	    {
	      if(bitqueue.front().test(0) && chread <= (nch-1) )
		{
		  msg_err<<"Fatal decoding error had been encounter : "<<endl
			 <<"The first bit of the next ch sequence should always be 0 and not "
			 <<bitqueue.front()<<endl;
		  abort();
		}
	    }

	  // reset the bitsread counter
	  bitsread = 0;
	  
	  // save channel data
	  adc.insert(adc.end(), chdata.begin(), chdata.end() );
	  chread++;
	  chdata.clear();
	  
	  if(m_Verbosity > 1)
	    cout<<"Decoded "<<adc.size()<<" samples"<<endl<<endl;
	}
      
      // basic consistency check
      if(chread == nch)
	{
	  // let's see if we grabbed more bytes than needed
	  if(bitqueue.size() > 0) // should be no more than 8
	    {
	      msg_info<<"Bits remaning in the queue "<<bitqueue.size()<<endl;
	      byteidx--;
	    }
	  
	  
	  if(m_Verbosity > 0)
	    {
	      if(m_Verbosity > 1) 	      
		msg_info<<"Bits remaning in the queue "<<bitqueue.size()<<endl;
	      msg_info<<"DECOMPRESSED EVENT STATUS: OK "<<endl;
	    }

	  // ok we are done
	  break;
	}
    } // end while loop
}

//
//
// 
// Decompress event data from binary file
// It IS MANDATORY that the position in the file is 
// set to the beginning of the event data before calling this function
//
void HuffDataCompressor::DecompressEventData( std::ifstream &fin, 
					      short nbadc, 
					      size_t nch,
					      size_t seqlen,
					      std::vector< adc16_t > &adc)
{
  if(!SetNbitsAdc( nbadc ))
    {
      msg_err<<"ADC "<<nbadc<<" bits exceeds max supported "<<m_MaxAdcBits<<endl;
      return;
    }
  adc.clear();
  
  vector< adc16_t > chdata;
  deque< bitset<1> > bitqueue;

  stringstream ss;
  short lastdelta  = -999;
  size_t byteidx   = 0;
  size_t bitsread  = 0;

  int noread = 0;  
  bool flag;
  
  size_t chread  = 0;
  
  // basic idea is to read bits one byte at a time from the input buffer
  // the first entry in the bitqueue deck should always be aligned on the 
  // compressed / raw flag and the deck should contain at least m_PacketSize bits
  while( adc.size() != nch*seqlen )
    {
      ReadNextByte(fin, bitqueue, flag);
      if(!flag) noread++;
      
      if(noread > 1) //noread == 1, we process last bits but not read new byte
	{
	  msg_err<<"There seems to be a problem with decoding"<<endl
		 <<"Samples accumulated in this channel "<<chdata.size()<<endl
		 <<"Remaining bits are "<<bitqueue.size()<<" "<<bitqueue.empty()<<endl;
	  break;
	}
      
      //read more bytes if possible
      if( bitqueue.size() < m_PacketSize && noread < 1) 
	continue;

      bool iscomp = bitqueue.front().test(0);
      bitsread++;
      bitqueue.pop_front(); // compressed / raw flag      

      if(!iscomp) // uncompressed
	{
	  ss.str(""); //clear our code string (could carry over from compressed)
	  
	  if( bitqueue.size() < m_NbitsHC ) 
	    {
	      msg_err<<"Fatal decoding error has been encountered : "<<endl
		     <<"Number of bits in the uncompressed stream should be at least "
		     <<m_NbitsHC<<" the current value is "<<bitqueue.size()<<endl;
	      abort();
	    }
	  
	  size_t bitcounter = 0;
	  while(bitcounter < m_NbitsHC) 
	    {
	      ss << bitqueue.front();
	      bitqueue.pop_front();
	      bitcounter++;
	    }
	  bitsread += bitcounter;

	  adc16_t val = ( strtoul(ss.str().c_str(), 0, 2) & 0x7FFF );
	  chdata.push_back( val );
	  ss.str("");
	}
      else // handle compressed bits
	{
	  size_t bitcounter = 0;
	  size_t bitstoread = m_NbitsHC;
	  if( bitqueue.size() < m_NbitsHC ) 
	    bitstoread = bitqueue.size();
	  
	  while( bitcounter < bitstoread )
	    {	      
	      ss << bitqueue.front();
	      bitqueue.pop_front();
	      bitcounter++;
	      bitsread++;
	      // check if we have a good code
	      size_t strsize = ss.str().size();
	      if( strsize >= m_MinCodeSize && strsize <= m_MaxCodeSize ) 
		{
		  bool ok = (m_UCmMap[strsize-1].first == ss.str());
		  if(!ok) continue; //not a good code
		  short val = m_UCmMap[strsize-1].second;
		  
		  // otherwise reset the code
		  ss.str(""); //reset string code
		  
		  // add to our adc vector
		  if( std::abs(val) <= m_MaxDiff)
		    {
		      adc16_t newval = chdata.back() + val;
		      chdata.push_back( newval );
		      lastdelta = val;
		    }
		  else if(val == m_NSeqRepVal)
		    {
		      for(short j=0;j<m_NSeqRep;j++)
			{
			  adc16_t newval = chdata.back() + lastdelta; //add same value
			  chdata.push_back( newval );
			}
		    }		
		}// end of if
	      
	      // we got the whole sequence now
	      if(chdata.size() == seqlen) break;

	    } // end while bitcounter
	} // end else 
      
      // check if we finished decoding data for a given channel
      if(chdata.size() == seqlen)
	{
	  
	  // calculate number of padded bits to the byte boundary
	  size_t padbits = bitsread % m_NbitsByte;
	  if(padbits > 0) padbits = m_NbitsByte - padbits;
	  
	  if(m_Verbosity > 1)
	    {
	      msg_info<<"Bytes read : "<<byteidx<<endl
		      <<"Bits read  : "<<bitsread<<endl
		      <<"Bits to boundary : "<<padbits<<endl
		      <<"Last value : "<<chdata.back()<<" ADC "<<endl;
	    }
	  
	  // error checking
	  if( bitqueue.size() < padbits )
	    {
	      msg_err<<"Fatal decoding error has been encountered : "<<endl
		     <<" Byte boundary does not appear to be valid"<<endl
		     <<" Check that the codes are padded with 0 to the next byte boundary"<<endl;
	      abort();
	    }
	  
	  // remove padded bits
	  while(padbits > 0)
	    {
	      bitqueue.pop_front();
	      padbits--;
	    }
	  
	  // some basic check
	  if(!bitqueue.empty())
	    {
	      if(bitqueue.front().test(0) && chread <= (nch-1) )
		{
		  msg_err<<"Fatal decoding error had been encounter : "<<endl
			 <<" The first bit of the next ch sequence should always be 0 and not "
			 <<bitqueue.front()<<endl;
		  abort();
		}
	    }

	  // reset the bitsread counter
	  bitsread = 0;
	  
	  // save channel data
	  adc.insert(adc.end(), chdata.begin(), chdata.end() );
	  chread++;
	  chdata.clear();
	  
	  if(m_Verbosity > 1)
	    msg_info<<"Decoded "<<adc.size()<<" samples"<<endl<<endl;
	}
      
      // basic consistency check
      if(chread == nch)
	{
	  // let's see if we grabbed more bytes than needed
	  if(bitqueue.size() > 0) // should be no more than 8
	    {
	      msg_info<<"Bits remaning in the queue "<<bitqueue.size()<<endl;
	      msg_info<<"Current position in file "<<fin.tellg();
	      fin.seekg( -1, std::ios::cur );
	      msg_info<<" Prv position in the file "<<fin.tellg()<<endl;
	    }
	  
	  
	  if(m_Verbosity > 0)
	    {
	      if(m_Verbosity > 1) 	      
		msg_info<<"Bits remaning in the queue "<<bitqueue.size()<<endl;
	      msg_info<<"DECOMPRESSED EVENT STATUS: OK "<<endl;
	    }

	  // ok we are done
	  break;
	}
      
    } // end while 
}
