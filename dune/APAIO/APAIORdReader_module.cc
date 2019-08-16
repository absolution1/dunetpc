////////////////////////////////////////////////////////////////////////
// Class:       APAIORdReader
// Plugin Type: producer (art v3_02_06)
// File:        APAIORdReader_module.cc
//
// Generated at Fri Aug  9 10:57:30 2019 by Thomas Junk using cetskelgen
// from cetlib version v3_07_02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
class APAIORdReader;
#include "TSystem.h"
#include "TROOT.h"
#include "TString.h"
#include "lardataobj/RawData/RawDigit.h"
#include <memory>


class APAIORdReader : public art::EDProducer {
public:
  explicit APAIORdReader(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  APAIORdReader(APAIORdReader const&) = delete;
  APAIORdReader(APAIORdReader&&) = delete;
  APAIORdReader& operator=(APAIORdReader const&) = delete;
  APAIORdReader& operator=(APAIORdReader&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  // Declare member data here.
  std::string fInputDirBaseName;
  uint32_t fAPAmin;
  uint32_t fAPAmax;

};


APAIORdReader::APAIORdReader(fhicl::ParameterSet const& p)
  : EDProducer{p}  // ,
// More initializers here.
{
  fInputDirBaseName = p.get<std::string>("InputDirBaseName","APAIO");
  fAPAmin = p.get<uint32_t>("APAMin",0);
  fAPAmax = p.get<uint32_t>("APAMax",150);
  produces<std::vector<raw::RawDigit>>();

}

void APAIORdReader::produce(art::Event& e)
{
  // Implementation of required member function here.

  std::unique_ptr<std::vector<raw::RawDigit>>  digcol(new std::vector<raw::RawDigit>);

  auto runno = e.run();
  auto subrunno = e.subRun();
  auto eventno = e.event();

  TString inputdir=fInputDirBaseName;
  inputdir += "r";
  TString rnst="";
  rnst.Form("%08d",runno);
  inputdir += rnst;
  inputdir += "s";
  TString srnst="";
  srnst.Form("%05d",subrunno);
  inputdir += srnst;
  inputdir += "e";
  TString est="";
  est.Form("%08d",eventno);
  inputdir += est;
  gSystem->ExpandPathName(inputdir);

  void* dirp = gSystem->OpenDirectory(inputdir);
  const char* entry;

  size_t version=0;

  TString fullfilename = "";
  while((entry = (char*)gSystem->GetDirEntry(dirp))) 
    {
      TString str=entry;
      if (str.Contains("apa"))
	{
	  fullfilename = inputdir;
	  fullfilename += "/";
	  fullfilename += entry;
	  //std::cout << " APAIORdReader found file: " << fullfilename << std::endl;
	  FILE *infile = fopen(fullfilename,"r");
	  while (infile != NULL)
	    {
	      uint32_t chan=0;
	      uint32_t samples=0;
	      raw::Compress_t compression = raw::kNone;
	      float pedestal=0;
	      float sigma=0;
	      uint32_t nadc=0;

	      fread(&version,sizeof(version),1,infile);
	      if (feof(infile)) break;

	      fread(&chan,sizeof(chan),1,infile);
	      //std::cout << "chan: " << chan << std::endl;
	      uint32_t apanum = chan/2560;
	      if (apanum > fAPAmax || apanum < fAPAmin) break;

	      fread(&samples,sizeof(samples),1,infile);
	      fread(&compression,sizeof(compression),1,infile);
	      fread (&pedestal,sizeof(pedestal),1,infile);
	      fread (&sigma,sizeof(sigma),1,infile);
	      fread(&nadc,sizeof(nadc),1,infile);
	      raw::RawDigit::ADCvector_t adcs(nadc);
	      fread(adcs.data(),sizeof(adcs.front()),adcs.size(),infile);
	      digcol->emplace_back(chan,samples,adcs,compression);
	      digcol->back().SetPedestal(pedestal,sigma);
	    }
	  if (infile != NULL) fclose(infile);
	}
    }
  gSystem->FreeDirectory(dirp);
  e.put(std::move(digcol));
}

DEFINE_ART_MODULE(APAIORdReader)
