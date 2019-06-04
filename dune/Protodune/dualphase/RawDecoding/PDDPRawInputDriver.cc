/*
    Input source driver for ProtoDUNE-DP raw data
    

 */

#include "canvas/Persistency/Provenance/SubRunID.h"	
#include "lardataobj/RawData/RawDigit.h"


#include "PDDPRawInputDriver.h"


//
//
//
namespace 
{
  void unpackCroData( const char *buf, size_t nb, bool cflag, unsigned nsa, adcdata_t &data )
  {
    data.clear();
    if( !cflag ) 
      {
	const BYTE* start = buf;
	const BYTE* stop  = start + nb;
	while(start!=stop)
	  {
	    BYTE v1 = *start++;
	    BYTE v2 = *start++;
	    BYTE v3 = *start++;
	    
	    uint16_t tmp1 = ((v1 << 4) + ((v2 >> 4) & 0xf)) & 0xfff;
	    uint16_t tmp2 = (((v2 & 0xf) << 8 ) + (v3 & 0xff)) & 0xfff;
	    data.push_back( tmp1 );
	    data.push_back( tmp2 );
	  }
      }
    else//TODO finalize the format of compressed data
      {
	msg_err<<"The compression should not be set yet"<<endl;
	//DecompressEventDataFast( buf, nb, nsa, data );
      }  
    
    //return data.size();
  }

  // get byte content for a given data type
  // NOTE: cast assumes host byte order
  template<typename T> T ConvertToValue(const void *in)
    {
      const T *ptr = static_cast<const T*>(in);
      return *ptr;
    }
}



//
namespace lris
{
  
  //
  PDDPRawInputDriver::PDDPRawInputDriver( fhicl::ParameterSet const &pset,
					  art::ProductRegistryHelper &helper,
					  art::SourceHelper const &pm ) :
    __sourceHelper( pm ),
    __currentSubRunID(),
    __eventCnt( 0 ),
    __eventNum( 0 )
  {
    helper.reconstitutes<std::vector<raw::RawDigit>, art::InEvent>("daq");
  }

  //
  void PDDPRawInputDriver::closeCurrentFile()
  {
    //mf::LogInfo(__FUNCTION__)<<"File boundary: processed " <<fEventCounter <<" events out of " <<fNEvents <<"\n";
    __close();
  }
   
  //
  void PDDPRawInputDriver::readFile( std::string const &name,
				     art::FileBlock* &fb )
  {
    
  }


  //
  bool PDDPRawInputDriver::readNext(  art::RunPrincipal* const &inR,
				      art::SubRunPrincipal* const &inSR,
				      art::RunPrincipal* &outR,
				      art::SubRunPrincipal* &outSR,
				      art::EventPrincipal* &outE )
  {
    
    return true;
  }
  
}
