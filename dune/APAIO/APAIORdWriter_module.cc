////////////////////////////////////////////////////////////////////////
// Class:       APAIORdWriter
// Plugin Type: analyzer (art v3_02_06)
// File:        APAIORdWriter_module.cc
//
// Generated at Thu Aug  8 16:29:17 2019 by Thomas Junk using cetskelgen
// from cetlib version v3_07_02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TString.h"
#include "lardataobj/RawData/RawDigit.h"

class APAIORdWriter;


class APAIORdWriter : public art::EDAnalyzer {
public:
  explicit APAIORdWriter(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  APAIORdWriter(APAIORdWriter const&) = delete;
  APAIORdWriter(APAIORdWriter&&) = delete;
  APAIORdWriter& operator=(APAIORdWriter const&) = delete;
  APAIORdWriter& operator=(APAIORdWriter&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:

  // Declare member data here.

  std::string fRawDigitLabel;
  std::string fOutputDirBaseName;

};


APAIORdWriter::APAIORdWriter(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.

  fRawDigitLabel = p.get<std::string>("RawDigitLabel","daq");
  fOutputDirBaseName = p.get<std::string>("OutputDirBaseName","APAIO");
  consumes<std::vector<raw::RawDigit>>(fRawDigitLabel);
}

void APAIORdWriter::analyze(art::Event const& e)
{
  auto runno = e.run();
  auto subrunno = e.subRun();
  auto eventno = e.event();

  auto rdighandle = e.getValidHandle<std::vector<raw::RawDigit>>(fRawDigitLabel);
  if (rdighandle->size()==0) return;

  TString outdir=fOutputDirBaseName;
  outdir += "r";
  TString rnst="";
  rnst.Form("%08d",runno);
  outdir += rnst;
  outdir += "s";
  TString srnst="";
  srnst.Form("%05d",subrunno);
  outdir += srnst;
  outdir += "e";
  TString est="";
  est.Form("%08d",eventno);
  outdir += est;

  TString cmd;
  cmd = "touch " + outdir;
  gSystem->Exec(cmd);  // todo -- test these system commands for failure
  cmd = "rm -rf " + outdir;
  gSystem->Exec(cmd);
  cmd = "mkdir " + outdir;
  gSystem->Exec(cmd);  

  // todo -- get these from the geometry
  size_t napas = 150;
  size_t nchans_per_apa = 2560;

  size_t version = 1;  // data format version

  // todo -- make an index of channel ID's found in the raw digits and sort them.  Then look in the sorted vector using lower_bound
  // so we don't have to do a double loop of APAs*channels.  This may be fast enough as is.

  for (uint32_t iapaindex = 0; iapaindex < napas; ++iapaindex)
    {
      //std::cout << "Writing output for apa: " << iapaindex << std::endl;
      FILE *ofile=NULL;
      TString ofilename=outdir;
      ofilename += "/apa";
      TString apast="";
      apast.Form("%03d",iapaindex);
      ofilename += apast;  
      for (auto const& rawdigit : *rdighandle)
        {
          uint32_t chan = rawdigit.Channel();
          uint32_t apa = chan / nchans_per_apa;
	  if (apa == iapaindex)
	    {
	      if (ofile == NULL)
		{
		  ofile = fopen(ofilename,"w");
		}
	      fwrite(&version,sizeof(version),1,ofile);
	      fwrite(&chan,sizeof(chan),1,ofile);
	      uint32_t samples = rawdigit.Samples();
	      fwrite(&samples,sizeof(samples),1,ofile);
	      raw::Compress_t compression = rawdigit.Compression();	      
	      fwrite(&compression,sizeof(compression),1,ofile);
	      float pedestal = rawdigit.GetPedestal();
	      float sigma = rawdigit.GetSigma();
	      fwrite (&pedestal,sizeof(pedestal),1,ofile);
	      fwrite (&sigma,sizeof(sigma),1,ofile);

	      // todo: optionally compress data further
	      auto adcs = rawdigit.ADCs();
	      uint32_t nadc = adcs.size();
	      fwrite(&nadc,sizeof(nadc),1,ofile);
	      //std::cout << "writing " << adcs.size() << " adcs with size: " << sizeof(adcs.front()) << " for channel: " << chan << std::endl;
	      fwrite(adcs.data(),sizeof(adcs.front()),adcs.size(),ofile);
	    }
        }
      if (ofile != NULL) fclose(ofile);
    }
}

DEFINE_ART_MODULE(APAIORdWriter)
