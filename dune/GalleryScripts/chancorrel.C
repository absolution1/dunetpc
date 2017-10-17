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
#include "TVectorT.h"
#include "lardataobj/RawData/RawDigit.h"

using namespace art;
using namespace std;

// ROOT script using gallery to make a correlation plot of raw::RawDigits for the ievcount'th event in the file, using ticks tickmin to tickmax
// Tom Junk, Sep. 2017.  
// Oct 2017, change to plot the channel ID's from raw::RawDigit intead of just the index of which one it is.
// Example invocation for a 3x1x1 imported rootfile, make a correlation plot for the first event in the file, for ticks 800 to 1000
// root [0] .L chancorrel.C++
// root [1] chancorrel("/pnfs/dune/tape_backed/dunepro/test-data/dune/raw/01/85/12/09/wa105_r842_s32_1501156823.root",0,800,1000,"daq");

// also computes the slope of a least-squares straight-line fit to the ADC's each channel vs. every other channel.  This is computed using
// the same intermediate data and so incurs almost no extra time.  Results not plotted by default -- draw the "slope" histo if you want.

// arguments:  filename -- input file, larsoft formatted
// ievcount:  which to use to display the correlation plot.  This is the tree index in the file and not the event number
// tickmin, tickmax -- to truncate the part of the event to run the correlation calc on.  Set to big and small numbers for no truncation.
// inputtag: use "daq" for MC and 3x1x1 imported data.  It's SplitterInput:TPC for split 35t data

void chancorrel(std::string const& filename,
             size_t ievcount=0, 
             size_t tickmin=0,
             size_t tickmax=100000, 
	     std::string const& inputtag="daq")

{
  // black = +1 correlation, red=-1 correlation, white=no correlation

  gStyle->SetOptStat(0);
  Int_t MyPalette[100];
  Double_t Red[]    = {1.0, 1., 0.0};
  Double_t Green[]  = {0.0, 1., 0.0};
  Double_t Blue[]   = {0.0, 1., 0.0};
  Double_t Length[] = {0.0, 0.5, 1.0};
  Int_t FI = TColor::CreateGradientColorTable(3, Length, Red, Green, Blue, 100);
  for (int i=0;i<100;i++) MyPalette[i] = FI+i;
  gStyle->SetPalette(100, MyPalette);

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
	    const size_t nrawdigits = rawdigits.size();
	    size_t tlow = TMath::Max(tickmin, (size_t) 0);
	    size_t thigh = TMath::Min(tickmax, (size_t) rawdigits[0].Samples()-1); // assume uncompressed; all channels have the same number of samples
	    size_t nticks = thigh - tlow + 1;

	    // find the maximum channel; make the plot from channel zero to the maximum channel

	    size_t nchans=0;
	    for (size_t i=0; i<nrawdigits; ++i)
	      {
		size_t ic = rawdigits[i].Channel();
		if (nchans<ic) nchans=ic;
	      }
	    nchans++;

	    TVectorT<double> x(nchans);
	    TVectorT<double> sumx(nchans);
	    TVectorT<double> sumxx(nchans);
	    TVectorT<double> sumxy[nchans];
	    for (size_t ichan=0; ichan < nchans; ++ichan)
	      {
		sumxy[ichan].ResizeTo(nchans);
	      }

	    TH2D *correl = (TH2D*) new TH2D("correl","Correl",nchans,-0.5,nchans-0.5,nchans,-0.5,nchans-0.5);
	    TH2D *slope = (TH2D*) new TH2D("slope","Slope",nchans,-0.5,nchans-0.5,nchans,-0.5,nchans-0.5);

            for (size_t itick=tlow;itick<=thigh;++itick)
	      {
	         for (size_t ichan=0;ichan<nrawdigits; ++ichan) 
		   { 
		     size_t ic = rawdigits[ichan].Channel(); 
		     //cout << ic << endl;
		     if ( ic > 2047) cout << ichan << " " << ic << " " << nrawdigits << endl;
		     x[ic] = rawdigits[ichan].ADC(itick); 
		   }
		 sumx += x;
		 for (size_t ic=0;ic<nchans;++ic) 
		   { 
		     sumxy[ic] += x[ic]*x;
		   }
		 sumxx += x.Sqr();  // nb -- this overwrites x with its square so do it last
	      }
	    
	    cout << "Finished filling vectors" << endl;

	    for (size_t ichan=0;ichan<nchans;++ichan)
	      {
		for (size_t ichan2=0;ichan2<nchans;++ichan2)
		  {
		    double cnum = ( nticks*sumxy[ichan][ichan2] -
				    sumx[ichan]*sumx[ichan2] );
		    double cden = 
		      ( TMath::Sqrt( nticks*sumxx[ichan] - 
				     sumx[ichan]*sumx[ichan] ) *
		        TMath::Sqrt( nticks*sumxx[ichan2] - 
				     sumx[ichan2]*sumx[ichan2] ) );
		    double c = 0;
		    if (cden>0) c=cnum/cden;

		    correl->SetBinContent(ichan+1,ichan2+1,c);

		    double s = ( nticks*sumxy[ichan][ichan2] -
				 sumx[ichan]*sumx[ichan2] ) /
		      ( nticks*sumxx[ichan] - sumx[ichan]*sumx[ichan] );
		    if (c>0.8) slope->SetBinContent(ichan+1,ichan2+1,s);
		    else slope->SetBinContent(ichan+1,ichan2+1,0);
		  }
	      }

	    correl->SetDirectory(0);
            correl->SetMinimum(-1);
	    correl->SetMaximum(1);
	    correl->GetXaxis()->SetTitle("Channel");
	    correl->GetYaxis()->SetTitle("Channel");
	    correl->Draw("colz");
	    slope->SetDirectory(0);
	    slope->GetXaxis()->SetTitle("Channel");
	    slope->GetYaxis()->SetTitle("Channel");
	    //slope->Draw("colz");
	  }
      }
    ++evcounter;
  }
}
