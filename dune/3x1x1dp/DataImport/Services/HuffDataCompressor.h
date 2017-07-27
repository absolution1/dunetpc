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
//  the most performant, but hopefully easier to understand ;-)
//
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

#ifndef __HUFFDATACOMPRESSOR_H__
#define __HUFFDATACOMPRESSOR_H__

#include <map>
#include <vector>
#include <deque>
#include <bitset>
#include <fstream>
#include <utility>

#include "dlardaq.h"

namespace dlardaq
{
  //
  // compress/decompress raw ADC data using Huffman encoding
  //
  class HuffDataCompressor
  {
  public:
    static HuffDataCompressor& Instance()
    {
      static HuffDataCompressor inst;
      return inst;
    }
    
    // compress single sequence
    // need to provide number of bits of the ADC
    // and the raw data vector
    // outputs compressed stream packed into 8 bit vector
    // 
    void CompressChData( short nbadc, std::vector<adc16_t> &raw_in,
			 std::vector<BYTE> &bin_out );
    
    // the raw adc should just be adc values without event header
    // these data are in 1D array: [tdc + ch x ntdc]
    void CompressEventData( short nbadc, size_t nch, size_t seqlen,
			    std::vector<adc16_t> &raw_in,
			    std::vector<BYTE> &bin_out );

    // the raw adc should just be adc values without event header
    // these data are in 2D array: [ch, tdc]
    void CompressEventData( short nbadc, size_t nch, size_t seqlen,
			    std::vector< std::vector<adc16_t> > &raw_in,
			    std::vector<BYTE> &bin_out );

    
    // decompress event from binary sequence in memory
    // 
    void DecompressEventData( short nbadc,       
			      size_t nch,
			      size_t seqlen, 
			      const char *buf, size_t bufsize, size_t &byteidx,
			      std::vector<adc16_t>  &adc );

    // NOTE: Before calling this function 
    //       one must ensure that the position of the input file 
    //       is always set at the begining of the event data
    //       This function reads the file byte by byte so it could 
    //       be rather slow
    void DecompressEventData( std::ifstream &fin, 
			      short nbadc, 
			      size_t nch,
			      size_t seqlen,
			      std::vector< adc16_t > &adc);
			      

    //
    void PrintEncoding();
    
    void SetVerbosity(int val){ m_Verbosity = val; }

  private:

    // ctor
    HuffDataCompressor();
    HuffDataCompressor(const HuffDataCompressor&);
    HuffDataCompressor& operator=(const HuffDataCompressor&);
    // Prevent unwanted destruction
    ~HuffDataCompressor(){;}
    
    // define encoding scheme
    void SetEncoding();

    // structure to collect Huffman codes
    struct HuffAccum_t
    {
      std::string codestr;  // binary code; leave empty for raw data
      short value;            // value: diff / or raw adc
      short nrep;           // number of repetitions of this code
    };  

    //
    std::vector<HuffAccum_t> HuffBuffer;  // buffer for Huffman codes

    // add to buffer
    void AddToHuffBuffer(short val, bool rawadc);
    
    
    //
    bool SetNbitsAdc( short nbadc );

    //
    void AddWordsToByteBuffer(std::string words, 
			      std::vector<BYTE> &buf, 
			      std::string &partbyte);
    
    //
    void ReadNextByte( size_t &byteidx, const char *buf,
		       std::deque< std::bitset<1> > &bits );

    // read next byte from file stream
    void ReadNextByte( std::ifstream &fin, std::deque< std::bitset<1> > &bits, bool &status );
    
    
    // get code from value
    std::string GetCodeFromValue(short val); 
    
    // get value from code
    //short GetValueFromCode(std::string code, bool &ok);
    
    std::string GetBinaryString(short val, size_t strsize);
        
    //
    //
    std::map<short, std::string> m_CmMap;     // map to compresss
    //std::map<std::string, short> m_UCmMap;    // map to decompress

    // map to decompress is in the sorted vector of codes
    std::vector< std::pair<std::string, short> >  m_UCmMap;    

    int    m_Verbosity;                       //

    short  m_MaxAdcBits;                      // max ADC bits
    short  m_MaxDiff;                         // max difference
    size_t m_PacketSize;                      // size of packet in bits
    size_t m_MinCodeSize;                     // shortest code length
    size_t m_MaxCodeSize;                     // longest code length

    size_t m_NbitsByte;    // size in bits of buffer
    size_t m_NbitsHead;   // number of header bits
    size_t m_NbitsHC;     // number of bits for huffman codes

    short  m_NSeqRep;     // encode sequence repetition
    short  m_NSeqRepVal;  // dummy value
    bool   m_SeqEnable;   // flag to encode repetitions 
  };
}


#endif
