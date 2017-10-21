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
#include "TH2D.h"
#include "TStyle.h"
#include "TColor.h"
#include "TMath.h"
#include "TVirtualFFT.h"
#include "TPad.h"
#include "lardataobj/RawData/RawDigit.h"

using namespace art;
using namespace std;

// ROOT script using gallery to make a FFT plot by channel the ievtcount'th event in the input file using raw::RawDigits.  Produces the plot interactively on a TCanvas
// Tom Junk, Sep. 2017
// Oct 2017 -- use the channel index from raw::RawDigits instead of the index of which rawdigit it is.
// Example invocation for a 3x1x1 imported rootfile, make an FFT plot for the first event in the file.
// root [0] .L chanfft.C++
// root [1] chanfft("/pnfs/dune/tape_backed/dunepro/test-data/dune/raw/01/85/12/09/wa105_r842_s32_1501156823.root",0,0,10000000,"daq",2.5);

// limitation:  does not call the raw::Uncompress method, and so zero-suppressed rawdigits cannot be used.
// Needs gallery to be setup.

// arguments:  filename -- input file, larsoft formatted
// ievcount:  which event to display the FFT for.  This is the tree index in the file and not the event number
// tickmin, tickmax -- to truncate the part of the event to run the FFT on.  Set to big and small numbers for no truncation.
// digifreq:  ADC Sampling frequency in MHz.  Used to label plots.  It's 2 MHz for 35t, DUNE FD MC, and 2.5 MHz for 3x1x1
// plotmin, plotmax:  used to aid in coloring the plot.
// inputtag: use "daq" for MC and 3x1x1 imported data.  It's SplitterInput:TPC for split 35t data

void chanfft(std::string const& filename, 
             size_t ievcount=0, 
             size_t tickmin=0,
             size_t tickmax=100000, 
	     std::string const& inputtag="daq", 
             double digifreq=2.0,
             double plotmin=0, 
             double plotmax=1E21)
{
  //gStyle->SetPalette(kInvertedDarkBodyRadiator);
  //gStyle->SetPalette(kVisibleSpectrum);
  gStyle->SetPalette(kRainBow);
  //gStyle->SetPalette(kColorPrintableOnGrey);

  gStyle->SetOptStat(0);

  size_t evcounter=0;

  InputTag rawdigit_tag{ inputtag };
  //InputTag rawdigit_tag{ "daq" };
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
	    nchans++;

	    size_t tlow = TMath::Max(tickmin, (size_t) 0);
	    size_t thigh = TMath::Min(tickmax, (size_t) rawdigits[0].Samples()-1); // assume uncompressed; all channels have the same number of samples
	    size_t nticks = thigh - tlow + 1;

	    double x[nticks];
	    double re=0;
	    double im=0;
	    double mag=0;

	    TH2D *cf = (TH2D*) new TH2D("cf","FFT by channel",nchans,-0.5,nchans-0.5,nticks/2,0,digifreq/2.0);

	    Int_t nti = nticks;
	    TVirtualFFT *fftr2c = TVirtualFFT::FFT(1,&nti,"R2C ES K");
            for (size_t ichan=0;ichan<nrawdigits;++ichan)
	      {
		for (size_t itick=tlow; itick <= thigh; ++itick) x[itick-tlow] = rawdigits[ichan].ADC(itick); 
		//cout << x[0] << " " << x[1] << " " << x[2] << endl;
		fftr2c->SetPoints(x);
		fftr2c->Transform();
		size_t ic = rawdigits[ichan].Channel(); 
		for (size_t i=0;i<nticks/2;++i)
		  {
		    fftr2c->GetPointComplex(i, re, im);
		    mag = TMath::Sqrt(re*re + im*im);
		    // cout << mag << endl;
		    cf->SetBinContent(ic+1,i,mag);
		  }
	      }
	    
	    cout << "Finished fft'ing channels" << endl;

	    cf->SetDirectory(0);
	    cf->SetMaximum(cf->GetMaximum(plotmax));
	    cf->SetMinimum(cf->GetMinimum(plotmin));
	    cf->GetXaxis()->SetTitle("Channel");
	    cf->GetYaxis()->SetTitle("Frequency (MHz)");
	    cf->Draw("colz");
	    gPad->SetLogz();
	  }
      }
    ++evcounter;
  }
}
