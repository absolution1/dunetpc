/*
    Input source driver for ProtoDUNE-DP raw data
    

 */
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/ProductRegistryHelper.h"
#include "art/Framework/IO/Sources/SourceHelper.h"
#include "art/Framework/Core/FileBlock.h"
#include "art/Framework/Principal/RunPrincipal.h"
#include "art/Framework/Principal/SubRunPrincipal.h"
#include "art/Framework/Principal/EventPrincipal.h"

#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <vector>
#include <stdint.h>

// this namespace lris seems to be used for other raw data converters
namespace lris  
{
  //typedef struct timespec trigstamp_t; // time_t tv_sec, long tv_nsec
  typedef std::vector<uint16_t> adcdata_t;
  typedef std::vector<uint8_t>  eflags_t;
  typedef char BYTE;

  class PDDPRawInputDriver 
  {
  public: 
    // Required constructor
    PDDPRawInputDriver( fhicl::ParameterSet const &pset,
			art::ProductRegistryHelper &helper,
			art::SourceHelper const &pm );

    // Required by FileReaderSource
    void closeCurrentFile();
    void readFile( std::string const &name,
		   art::FileBlock* &fb);
    bool readNext( art::RunPrincipal* const &inR,
		   art::SubRunPrincipal* const &inSR,
		   art::RunPrincipal* &outR,
		   art::SubRunPrincipal* &outSR,
		   art::EventPrincipal* &outE );

  private: 
    art::SourceHelper const&	__sourceHelper;
    art::SubRunID 		__currentSubRunID;
    uint16_t 			__eventCnt; 
    uint16_t                    __eventNum;


    void __open(std::string finname);
    void __close();

    // read a chunk of binary data
    void __readChunk( std::vector<BYTE> &bytes, size_t sz );
    bool __unpackEvent( std::vector<BYTE> &buf );

    // file locations
    std::vector<std::streampos> __events;
    std::vector<uint32_t> __evsz;
    
    // input file
    size_t __filesz;
    std::ifstream __file;

    // header sizes
    enum { evinfoSz = 44 };
    
    // trigger structure from WR trigserver
    typedef struct triginfo_t
    {
      uint8_t type;
      uint32_t num;
      struct timespec ts; //{ time_t ts.tv_sec, long tv_nsec }
    } triginfo_t;
    
    // strucutre to hold decoded event header
    typedef struct eveinfo_t
    {
      uint32_t runnum;
      uint8_t  runflags; 
      triginfo_t ti;       // trigger info
      uint8_t  evflag;     // data quality flag
      uint32_t evnum;      // event number
      uint32_t evszlro;    // size of event in bytes
      uint32_t evszcro;    // size of event in bytes
    } eveinfo_t;
    
    // fragment from each L1 event builder
    typedef struct fragment_t
    {
      eveinfo_t ei;
      const BYTE* bytes;
      adcdata_t crodata;
      adcdata_t lrodata;
    } fragment_t;

    //
    unsigned __unpack_evtable();
    unsigned __unpack_eve_info( const char *buf, eveinfo_t &ei );
 
  };
}
