////////////////////////////////////////////////////////////////////////////
// \file RawData311InputDriver.h
// brief Source to convert raw binary file from 311 to root file useful to LArSoft
// Adapted from LarRawInputDriverShortBo.h
//
// Created: April 20th 2017, Last Modified: 
// Author: Kevin Fusshoeller, kevin.fusshoeller@cern.ch
////////////////////////////////////////////////////////////////////////////


#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/ProductRegistryHelper.h"
#include "art/Framework/IO/Sources/SourceHelper.h"
#include "art/Framework/Core/FileBlock.h"
#include "art/Framework/Principal/RunPrincipal.h"
#include "art/Framework/Principal/SubRunPrincipal.h"
#include "art/Framework/Principal/EventPrincipal.h"
#include "canvas/Persistency/Provenance/SubRunID.h"	
#include "lardataobj/RawData/RawDigit.h"

#include "dlardaq.h"
#include "EventDecoder.h"

#include <fstream>
#include <string>
#include <vector>

// Conversion of binary data to root files
namespace lris 
{

class RawData311InputDriver 
{
  // Class to fill constraints on a template argument to the class FileReaderSource
  public: 
    // Required constructor
    RawData311InputDriver(fhicl::ParameterSet const &pset,
			art::ProductRegistryHelper &helper,
			art::SourceHelper const &pm);

    // Functions required by FileReaderSource
    void closeCurrentFile();
    void readFile(std::string const &name,
		  art::FileBlock* &fb);
    bool readNext(art::RunPrincipal* const &inR,
		  art::SubRunPrincipal* const &inSR,
		  art::RunPrincipal* &outR,
		  art::SubRunPrincipal* &outSR,
		  art::EventPrincipal* &outE);

  private: 
    art::SourceHelper const&	fSourceHelper;
    art::SubRunID 		fCurrentSubRunID;
    //std::ifstream		infile;
    uint16_t 			fEventCounter;
    uint16_t 			fNEvents;
    std::string 		filename;
   
    short 			nchannels = 1280;
    short 			nsamples = 1667;
    dlardaq::EventDecoder	DataDecode;  

    dlardaq::runheader_t 	file_head;
    dlardaq::footer_t 		file_foot;

    std::vector< std::pair<double, double> > fPedMap;
    std::string 		fPedestalFile;

    void process_Event311(std::vector<raw::RawDigit>& digitList,
			     dlardaq::evheader_t &event_head,
			     uint16_t evt_num);

    double GetPedMean(size_t LAr_chan, std::vector< std::pair<double, double> > *fPedMap){ return fPedMap->at(LAr_chan).first; }

    double GetPedRMS(size_t LAr_chan, std::vector< std::pair<double, double> > *fPedMap){ return fPedMap->at(LAr_chan).second; }

     
}; //RawData311InputDriver

} //namespace lris
