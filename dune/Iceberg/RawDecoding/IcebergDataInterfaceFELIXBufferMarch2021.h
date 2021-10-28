////////////////////////////////////////////////////////////////////////
// IcebergDataInterfaceFELIXBufferMarch2021.h
//
// Tool to unpack RCE and FELIX fragments.  A restructuring from the IcebergRawDecoder module
//
// These methods take references to vectors of raw::RawDigit, raw::RDTimeStamp, and raw::RDStatus data products as arguments.
// These vectors are not cleared on input, and so when data are retrieved, they are appended to any existing data already
// in the vectors.  The RDStatus vector is an exception, where just one RDStatus instance will be in the vector.  Previously
// accumulated RDStatus values from previous calls will be logically ORed into the single RDStatus instances contents.
//
//  Methods are provided to retrieve all data from fragments on an input label, or by specified APA list.  In cases where
//  data from a specified APA are requested but no labels are provided by the caller, labels are input via FCL parameters.
//  This is true because data from an APA may appear with different labels during the course of the ProtoDUNE-SP run.
//
/////////////////////////////////////////////////////////////////////////
#ifndef IcebergDataInterfaceFELIXBufferMarch2021_H
#define IcebergDataInterfaceFELIXBufferMarch2021_H

#include <vector>

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
//#include "art/Framework/Principal/Run.h"
//#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/RDTimeStamp.h"
#include "dune/DuneObj/PDSPTPCDataInterfaceParent.h"

class IcebergDataInterfaceFELIXBufferMarch2021 : public PDSPTPCDataInterfaceParent {

 public:

  IcebergDataInterfaceFELIXBufferMarch2021(fhicl::ParameterSet const& ps);

  int retrieveData(art::Event &evt, std::string inputlabel, std::vector<raw::RawDigit> &raw_digits, std::vector<raw::RDTimeStamp> &rd_timestamps,
                   std::vector<raw::RDStatus> &rdstatuses );

  // method to get raw digits, RDTimeStamps, RDStatuses from input files
  // specified in the fcl configuration

  int retrieveDataAPAListWithLabels(art::Event &evt, std::string inputlabel, std::vector<raw::RawDigit> &raw_digits, std::vector<raw::RDTimeStamp> &rd_timestamps,
                                    std::vector<raw::RDStatus> &rdstatuses, 
                                    std::vector<int> &apalist);

  // method to get raw digits, RDTimeStamps, RDStatuses for a specified list of APA's.  The list of possible labels on which to find
  // APA data is proved by fcl configuration.

  int retrieveDataForSpecifiedAPAs(art::Event &evt, std::vector<raw::RawDigit> &raw_digits, std::vector<raw::RDTimeStamp> &rd_timestamps,
                                   std::vector<raw::RDStatus> &rdstatuses,  
                                   std::vector<int> &apalist);

  // inputLabel examples:  "daq:TPC" or "daq::ContainerTPC" for RCE, "daq:FELIX" or "daq::ContainerFELIX" for FELIX
  // returns:  0:  success, or   1: discarded corrupted data, or 2: kept some corrupted data

 private:

  // open files

  std::vector<FILE*> fInputFilePointers;

  // configuration parameters

  std::vector<std::string>   fInputFiles; 
  size_t                     fNSamples;
  bool                       fCompressHuffman;
  ULong64_t                  fDesiredStartTimestamp;
  bool                       fFirstRead;

  // some convenience typedefs for porting old code

  typedef std::vector<raw::RawDigit> RawDigits;
  typedef std::vector<raw::RDTimeStamp> RDTimeStamps;
  typedef std::vector<raw::RDStatus> RDStatuses;

  // private methods

  void unpack14(const uint32_t *packed, uint16_t *unpacked);

  void computeMedianSigma(raw::RawDigit::ADCvector_t &v_adc, 
                          float &median, 
                          float &sigma);

};


#endif
