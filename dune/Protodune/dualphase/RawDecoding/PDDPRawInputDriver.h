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

#include "lardataobj/RawData/RawDigit.h"

#include <fstream>
#include <string>


// anonymous namespace 
namespace 
{
  typedef struct timespec trigstamp_t; // time_t tv_sec, long tv_nsec
  //typedef std::vector<uint16_t> adcdata_t;
  typedef std::vector< raw::RawDigit::ADCvector_t > adcbuf_t;
  typedef std::vector<uint8_t>  eflags_t;
  typedef char BYTE;

  //
  struct DaqEvent
  {
    // event has been correctly loaded
    //bool valid;

    // global event quality flag
    bool good;
    
    // run info
    uint32_t runnum;       // run number
    uint8_t  runflags;
    
    // event info
    uint32_t evnum;
    eflags_t evflags;
    
    // trigger info
    uint8_t     trigtype;
    uint32_t    trignum;
    trigstamp_t trigstamp;
    
    // unpacked CRO ADC buffer
    adcbuf_t crodata;

    // unpacked CRO ADC data
    //adcdata_t lrodata;
      
    // number of decoded channels
    unsigned chcro() const { return crodata.size(); }

    // CRO data compression
    raw::Compress_t compression;
  };
}

//
// this namespace lris seems to be used for other raw data converters
namespace lris  
{

  ///
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
    //std::string                 __output_label;
    std::string                 __outlbl_digits;
    std::string                 __outlbl_status;
    std::string                 __outlbl_rdtime;
    std::string                 __output_inst;

    uint32_t 			__eventCtr; 
    uint32_t                    __eventNum;

    // number of uncompressed samples per channel
    size_t __nsacro;
    
    // ped inversion to deal with the inverted signal polarity
    unsigned __invped; 

    // close binary file
    void __close();

    // read a chunk of binary data
    void __readChunk( std::vector<BYTE> &bytes, size_t sz );

    // unpack binary data written by each L1 evb builder
    bool __unpackEvent( std::vector<BYTE> &buf, DaqEvent &event );

    // crp to daq mapping
    std::vector<unsigned> __daqch;

    // file locations
    std::vector<std::streampos> __events;
    std::vector<uint32_t> __evsz;
    
    // input file
    size_t __filesz;
    std::ifstream __file;
    unsigned __file_seqno;

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
      adcbuf_t crodata;
      //adcbuf_t lrodata;
    } fragment_t;

    //
    unsigned __unpack_evtable();
    unsigned __unpack_eve_info( const char *buf, eveinfo_t &ei );
    unsigned __get_file_seqno( std::string s);
  };
}

