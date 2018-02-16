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
#include "TCanvas.h"
#include "TPad.h"
#include "TGaxis.h"

using namespace art;
using namespace std;

// This is a slighlty modified version of Tom Junk's event display "evd_rawdigits.C" tuned for the 3x1x1.
// Many things are hard-coded, to be changed in the future.
// Chirstoph Alt, November 2017
//
// arguments:  filename -- input file, larsoft formatted
// ievcount:  which event to display.  This is the tree index in the file and not the event number
// inputtag:  use "daq" for 3x1x1 data or MC
// autoped:  true if you want to subtract the average of the adc values of a channel before displaying
// minval, maxval -- in order to color the plot well.  Minval: white; Maxval: black

// Example invocation for a 3x1x1 imported rootfile, make an event display for the first event in the file.
// root [0] .L evd_311.C++
// root [1] evd_311("/pnfs/dune/tape_backed/dunepro/test-data/dune/raw/01/85/12/09/wa105_r842_s32_1501156823.root",0);


void ReverseYAxis(TH2I *h);

void
evd_311(std::string const& filename, size_t ievcount, std::string const& inputtag="daq", bool autoped=true, int minval=-10, int maxval=30)
{

  gStyle->SetOptStat(0);

  size_t evcounter=0;

  InputTag rawdigit_tag(inputtag);

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
	    TCanvas *c = new TCanvas("c","c",1800,500);
	    TPad *pview0 = new TPad("pview0","",0.00,0,0.265,1);
	    TPad *pview1 = new TPad("pview1","",0.265,0,1,1);
	    pview0->Draw();
	    pview0->SetLeftMargin(0.25);
	    pview0->SetRightMargin(0.02);
	    pview0->SetTopMargin(0.07);
	    pview0->SetBottomMargin(0.17);
	    pview1->Draw();
	    pview1->SetLeftMargin(0.09);
	    pview1->SetRightMargin(0.12);
	    pview1->SetTopMargin(0.07);
	    pview1->SetBottomMargin(0.17);
	    pview0->cd();

	    TH2I *hview0 = (TH2I*) new TH2I("hview0","View 0",320,-0.5,320-0.5,nticks,-1667+0.5,0.5);
	    hview0->SetDirectory(0);
	    for (size_t ichan=0;ichan<320; ++ichan)
	      {
		size_t ic = rawdigits[ichan].Channel();
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
		    hview0->SetBinContent(ic+1,itick+1,rawdigits[ichan].ADC(nticks-itick-1)-iped);
		  }
	      }

	    TH2I *hview1 = (TH2I*) new TH2I("hview1","View 1",960,321-0.5,1280-0.5,nticks,-1667+0.5,0.5);
	    hview1->SetDirectory(0);
	    for (size_t ichan=320;ichan<nrawdigits; ++ichan)
	      {
		size_t ic = rawdigits[ichan].Channel();
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

		    hview1->SetBinContent(ic+1-320,itick+1,rawdigits[ichan].ADC(nticks-itick-1)-iped);
		  }
	      }

	    hview0->GetXaxis()->SetLabelSize(0.075);
	    hview0->GetXaxis()->SetLabelOffset(0.001);
	    hview0->GetXaxis()->SetTitleOffset(999);
   	    hview0->GetXaxis()->SetTitleSize(0);
	    hview0->GetXaxis()->SetNdivisions(4);

	    hview0->GetYaxis()->SetLabelSize(0.08);
	    hview0->GetYaxis()->SetLabelOffset(0.01);
	    hview0->GetYaxis()->SetTitleSize(0.075);
	    hview0->GetYaxis()->SetTitleOffset(1.75);

            hview0->SetMinimum(minval);
	    hview0->SetMaximum(maxval);
	    hview0->GetXaxis()->SetTitle("Channel");
	    hview0->GetYaxis()->SetTitle("Tick");
	    hview0->Draw("col");

//	    ReverseYAxis(hview0);

	    pview1->cd();

	    hview1->GetXaxis()->SetLabelSize(0.075);
	    hview1->GetXaxis()->SetTitleSize(0.075);

//	    hview1->GetYaxis()->SetLabelColor(kWhite);
//   	    hview1->GetYaxis()->SetLabelSize(0);
	    hview1->GetYaxis()->SetLabelSize(0.08);
	    hview1->GetYaxis()->SetLabelOffset(0.005);
	    hview1->GetYaxis()->SetTitleOffset(999);
   	    hview1->GetYaxis()->SetTitleSize(0);

	    hview1->GetZaxis()->SetLabelSize(0.075);
	    hview1->GetZaxis()->SetTitle("ADC");
//	    hview1->GetZaxis()->SetTitleOffset(0.5);
	    hview1->GetZaxis()->SetTitleSize(1);

            hview1->SetMinimum(minval);
	    hview1->SetMaximum(maxval);
	    hview1->GetXaxis()->SetTitle("Channel");
	    hview1->GetYaxis()->SetTitle("Tick");
	    hview1->Draw("colz");

	  }
      }
    ++evcounter;
  }
}

void ReverseYAxis(TH2I *h)
{
       // Remove the current axis
       h->GetYaxis()->SetLabelOffset(999);
       h->GetYaxis()->SetTickLength(0);
//       h->GetYaxis()->SetTitleOffset(999);
//       h->GetYaxis()->SetTitleSize(0);

       // Redraw the new axis
       gPad->Update();
       TGaxis *newaxis = new TGaxis(gPad->GetUxmin(),
                                    gPad->GetUymax(),
                                    gPad->GetUxmin()-0.001,
                                    gPad->GetUymin(),
                                    h->GetYaxis()->GetXmin(),
                                    h->GetYaxis()->GetXmax(),
                                    510,"R");
       newaxis->SetLabelFont(42);
       newaxis->SetLabelSize(0.075);
       newaxis->SetLabelOffset(-0.01);
//       newaxis->SetTitle("Ticks");
//       newaxis->SetTitleFont(42);
//       newaxis->SetTitleSize(0.075);
//       newaxis->SetTitleOffset(-0.03);
       newaxis->Draw();
}
