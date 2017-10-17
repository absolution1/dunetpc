#include <functional>
#include <iostream>
#include <string>
#include <vector>

#include "canvas/Utilities/InputTag.h"
#include "gallery/Event.h"

#include "TFile.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TGraph.h"
#include "TH2I.h"
#include "TStyle.h"
#include "TColor.h"
#include "lardataobj/RawData/RawDigit.h"

using namespace art;
using namespace std;

// make a poor-man's event display of raw::RawDigits for the ievcount'th event in the file using gallery
// limitation:  does not call the raw::Uncompress method, and so zero-suppressed rawdigits cannot be used.

// Tom Junk, September 2017
// Oct 2017 -- use channel indexes from raw::RawDigits instead of which rawdigit it is
// arguments:  filename -- input file, larsoft formatted
// ievcount:  which event to display.  This is the tree index in the file and not the event number
// autoped:  true if you want to subtract the average of the adc values of a channel before displaying
// minval, maxval -- in order to color the plot well.  Minval: white; Maxval: black
// inputtag:  use "daq" for 3x1x1 data or MC, and SplitterInput:TPC for split 35t data

// Example invocation for a 3x1x1 imported rootfile, make an event display for the first event in the file.
// root [0] .L evd_rawdigits.C++
// root [1] evd_rawdigits("/pnfs/dune/tape_backed/dunepro/test-data/dune/raw/01/85/12/09/wa105_r842_s32_1501156823.root",0);


void
evd_rawdigits(std::string const& filename, size_t ievcount, bool autoped=true, int minval=0, int maxval=30, std::string const& inputtag="daq")
{
  //gStyle->SetPalette(kGreyScale);
  //gStyle->SetPalette(kInvertedDarkBodyRadiator);
  //TColor::InvertPalette();


  gStyle->SetOptStat(0);
  Int_t MyPalette[100];
  Double_t Red[]    = {1., 0.0};
  Double_t Green[]  = {1., 0.0};
  Double_t Blue[]   = {1., 0.0};
  Double_t Length[] = {0., 1.0};
  Int_t FI = TColor::CreateGradientColorTable(2, Length, Red, Green, Blue, 100);
  for (int i=0;i<100;i++) MyPalette[i] = FI+i;
  gStyle->SetPalette(100, MyPalette);

  size_t evcounter=0;

  InputTag rawdigit_tag(inputtag);

  // InputTag rawdigit_tag{ "daq" }; for 3x1x1 data and MC
  //InputTag rawdigit_tag{ "SplitterInput:TPC" };   // for split 35t data
  // Create a vector of length 1, containing the given filename.
  vector<string> filenames(1, filename);

  for (gallery::Event ev(filenames); !ev.atEnd(); ev.next()) {
    if (evcounter == ievcount)
      {
	auto const& rawdigits = *ev.getValidHandle<vector<raw::RawDigit>>(rawdigit_tag);
	if (!rawdigits.empty())
	  {
	    size_t nrawdigits = rawdigits.size();
	    size_t nchans=0;
	    for (size_t i=0; i<nrawdigits; ++i)
	      {
		size_t ic = rawdigits[i].Channel();
		if (nchans<ic) nchans=ic;
	      }
	    nchans++;  // set plots to go from channel 0 to the maximum channel we find.

	    size_t nticks = rawdigits[0].Samples();  // assume uncompressed, all channels have the same
                                                    // number of samples.
	    TH2I *rdh = (TH2I*) new TH2I("rdh","Raw Digits",nchans,-0.5,nchans-0.5,nticks,-0.5,nticks-0.5);
	    rdh->SetDirectory(0);
	    for (size_t ichan=0;ichan<nrawdigits; ++ichan)
	      {
		size_t ic = rawdigits[ichan].Channel();
		// cout << "filling for channel: " << ic << " " << nticks << endl;
		int iped = 0;  // an integer for raw digits
		if (autoped)
		  {
		    double csum=0;
		    for (size_t itick=0;itick<nticks;++itick)
		      {
			csum += rawdigits[ichan].ADC(itick);
		      }
		    csum /= ( (double) nticks );
		    iped = (int) csum;
		  }
		for (size_t itick=0;itick<nticks;++itick)
		  {
		    //cout << "channel: " << ichan << " tick: " << itick << endl;
		    rdh->SetBinContent(ic+1,itick+1,rawdigits[ichan].ADC(itick)-iped);
		  }
	      }
            rdh->SetMinimum(minval);
	    rdh->SetMaximum(maxval);
	    rdh->GetXaxis()->SetTitle("Channel");
	    rdh->GetYaxis()->SetTitle("Tick");
	    rdh->Draw("colz");
	    //gDirectory->ls();
	  }
      }
    ++evcounter;
  }
}
