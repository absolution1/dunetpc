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
#include "TH1D.h"
#include "TStyle.h"
#include "TColor.h"
#include "TMath.h"
#include "TVectorT.h"
#include "TCanvas.h"
#include "lardataobj/RawData/RawDigit.h"
#include "TLegend.h"

using namespace art;
using namespace std;

//  ROOT script using gallery to make mean, rms, and drms plots using raw::RawDigits for the ievcount'th event in the file
// Tom Junk, Fermilab, Sep. 2017
// Oct 2017 -- update to use the channel number from the rawdigit and not just the index of the rawdigit
// Example invocation for a 3x1x1 imported rootfile, make a correlation plot for the first event in the file.
// root [0] .L meanrms.C++
// root [1] meanrms("/pnfs/dune/tape_backed/dunepro/test-data/dune/raw/01/85/12/09/wa105_r842_s32_1501156823.root");

// arguments:  filename -- input file, larsoft formatted
// ievcount:  which event to display the mean and RMS for.  This is the tree index in the file and not the event number
// tickmin, tickmax -- to truncate the part of the event run on.  Set to big and small numbers for no truncation.
// inputtag: use "daq" for MC and 3x1x1 imported data.  It's SplitterInput:TPC for split 35t data


void meanrms(std::string const& filename, 
             size_t ievcount=0, 
             size_t tickmin=0, 
             size_t tickmax=100000,  
             std::string const& inputtag="daq")
{

  gStyle->SetOptStat(0);

  size_t evcounter=0;

  double s2i = 1.0/TMath::Sqrt(2.0);

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
	    const size_t nrawdigits = rawdigits.size();
	    size_t nchans=0;
	    for (size_t i=0; i<nrawdigits; ++i)
	      {
		size_t ic = rawdigits[i].Channel();
		if (nchans<ic) nchans=ic;
	      }
	    nchans++;  // set plots to go from channel 0 to the maximum channel we find.

	    size_t tlow = TMath::Max(tickmin, (size_t) 0);
	    size_t thigh = TMath::Min(tickmax, (size_t) rawdigits[0].Samples()-1); // assume uncompressed; all channels have the same number of samples
	    size_t nticks = thigh - tlow + 1;

	    TVectorT<double> x(nchans);
	    TVectorT<double> dx(nchans);
	    TVectorT<double> sumx(nchans);
	    TVectorT<double> sumxx(nchans);
	    TVectorT<double> sumdx(nchans);
	    TVectorT<double> sumdxx(nchans);

	    TH1D *mean = (TH1D*) new TH1D("mean","Mean",nchans,-0.5,nchans-0.5);
	    TH1D *mean2 = (TH1D*) new TH1D("mean2","mean2",nchans,-0.5,nchans-0.5); // declare with nchans so it can be overplotted
                                                                                 // even though we have one fewer channel
	    TH1D *rmshist = (TH1D*) new TH1D("rmshist","RMS",nchans,-0.5,nchans-0.5);
	    TH1D *drmshist = (TH1D*) new TH1D("drmshist","dRMS",nchans,-0.5,nchans-0.5);

            for (size_t itick=tlow;itick<=thigh;++itick)
	      {
		for (size_t ichan=0;ichan<nrawdigits; ++ichan) 
		  { 
		    size_t ic = rawdigits[ichan].Channel();
		    x[ic] = rawdigits[ichan].ADC(itick);
		  } 
		 for (size_t ichan=1;ichan<nchans; ++ichan) 
		   {
		     dx[ichan] = (x[ichan]-x[ichan-1])*s2i;
		   }
		 dx[0] = 0;
		 sumx += x;
		 sumdx += dx;
		 sumxx += x.Sqr();  // nb -- this overwrites x with its square so do it last
		 sumdxx += dx.Sqr(); // nb -- this overwrites dx with its square so do it last
	      }
	    
	    cout << "Finished filling vectors" << endl;

	    for (size_t ichan=0;ichan<nchans;++ichan)
	      {
		  {
		    double avg = sumx[ichan]/nticks;
		    double avgsquare = sumxx[ichan]/nticks;
		    mean->SetBinContent(ichan+1,avg);
		    mean2->SetBinContent(ichan+1,avg);
		    double rms = 0;
		    double tval = avgsquare - avg*avg;
		    if (tval>0) rms = TMath::Sqrt(tval);
		    mean->SetBinError(ichan+1,rms);
		    rmshist->SetBinContent(ichan+1,rms);

		    avg = sumdx[ichan]/nticks;
		    avgsquare = sumdxx[ichan]/nticks;
		    tval = avgsquare - avg*avg;
		    rms = 0;
		    if (tval>0) rms = TMath::Sqrt(tval);
		    mean2->SetBinError(ichan+1,avg);
		    drmshist->SetBinContent(ichan+1,rms);

		  }
	      }

	    TCanvas *mycanvas = new TCanvas("c","c",700,600);
	    mycanvas->Divide(1,2);
            mycanvas->cd(1);
	    mean->SetDirectory(0);
	    mean2->SetDirectory(0);
	    mean->GetXaxis()->SetTitle("Channel");
	    mean->GetYaxis()->SetTitle("Mean ADC");
	    mean->Draw("HIST");

	    mycanvas->cd(2);
	    rmshist->SetDirectory(0);
	    rmshist->GetXaxis()->SetTitle("Channel");
	    rmshist->GetYaxis()->SetTitle("RMS and dRMS");
	    drmshist->SetDirectory(0);
	    rmshist->SetLineColor(kBlack);
	    drmshist->SetLineColor(kGreen);
	    rmshist->Draw("HIST");
	    drmshist->Draw("SAME,HIST");

	  }
      }
    ++evcounter;
  }
}
