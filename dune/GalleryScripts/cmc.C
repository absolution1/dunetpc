#include <functional>
#include <iostream>
#include <string>
#include <vector>

#include "canvas/Utilities/InputTag.h"
#include "gallery/Event.h"

#include "TFile.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TGraph2D.h"
#include "TMath.h"
//#include "lardataobj/RecoBase/SpacePoint.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "nusimdata/SimulationBase/MCParticle.h"

using namespace art;
using namespace std;

// make a poor-man's event display of mc -- for the ievcount'th event in the file
// Tom Junk, Fermilab, Oct 5, 2017.  Based on gallery demos by Marc Paterno.
// Usage:  setup dunetpc and gallery, run root.
//  .L cmc.C++
// cmc(0,"inputfile.root");
// the clip region is defined with mins and maxes.  A margin can be added around the outside.

void
cmc(size_t ievcount=0, std::string const& filename="prodgenie_nu_dune10kt_1x2x6_gen_g4_lar.root",  std::string tagstring="largeant",
    double xclipmin=-99999, double xclipmax=99999, double yclipmin=-99999, double yclipmax=99999, double zclipmin=-99999, double zclipmax=99999,
    double margin=200, bool showneutrons=false)
{

  size_t evcounter=0;

  InputTag mcptag{ tagstring };
  // Create a vector of length 1, containing the given filename.
  vector<string> filenames(1, filename);

  for (gallery::Event ev(filenames); !ev.atEnd(); ev.next()) {
    if (evcounter == ievcount)
      {
	auto const& mcparticles = *ev.getValidHandle<vector<simb::MCParticle>>(mcptag);
	if (!mcparticles.empty())
	  {

	    bool first=true;
	    size_t icolor=1;

            TCanvas *c = new TCanvas("c","TGraph2D Event Display",0,0,800,800);

	    // find the right size for the graph and make a graph with points in the corners

	    double minx=9999;
	    double maxx=-9999;
	    double miny=9999;
	    double maxy=-9999;
	    double minz=9999;
	    double maxz=-9999;
      	    TGraph2D *grt = new TGraph2D();
	    for (size_t imcp=0;imcp<mcparticles.size(); ++imcp)
	      {
		int ipdg = std::abs(mcparticles[imcp].PdgCode());
		// plot only leptons, p, k, pi
		if (ipdg == 11 || ipdg == 13 || ipdg == 15 || ipdg == 211 || ipdg == 321 || ipdg == 2212)
		  {
		    size_t numtpts = mcparticles[imcp].NumberTrajectoryPoints();
		    for (size_t ipt=0; ipt<numtpts; ++ipt)
		      {
			TLorentzVector pos = mcparticles[imcp].Position(ipt);
			double x = pos.X();
			double y = pos.Y();
			double z = pos.Z();
			if (x<xclipmin || x>xclipmax) continue;
			if (y<yclipmin || y>yclipmax) continue;
			if (z<zclipmin || z>zclipmax) continue;
			maxx = max(x,maxx);
			maxy = max(y,maxy);
			maxz = max(z,maxz);
			minx = min(x,minx);
			miny = min(y,miny);
			minz = min(z,minz);
		      }
		  }
	      }
	    grt->SetPoint(0,minz-margin,minx-margin,miny-margin);
	    grt->SetPoint(1,maxz+margin,maxx+margin,maxy+margin);
	    grt->SetMarkerColor(0);
	    grt->SetLineColor(0);
	    grt->SetMarkerStyle(1);
	    std::string titlestring=tagstring;
	    titlestring += " MCParticle Display";
	    grt->SetTitle(titlestring.c_str());
	    grt->Draw("P");
	    grt->GetXaxis()->SetTitle("Z");
	    grt->GetXaxis()->SetTitleColor(4);
	    grt->GetYaxis()->SetTitle("X");
	    grt->GetYaxis()->SetTitleColor(4);
	    grt->GetZaxis()->SetTitle("Y");
	    grt->GetZaxis()->SetTitleColor(4);
	    //cout << "boundaries: " << minx << " " << miny << " " << minz << " " << maxx << " " << maxy << " " << maxz << endl;

	    for (size_t imcp=0;imcp<mcparticles.size(); ++imcp)
	      {
		int ipdg = std::abs(mcparticles[imcp].PdgCode());
		// plot only leptons, p, k, pi
		int icolor=0;
		if (ipdg == 11) icolor = 6;  // magenta for electrons
		if (ipdg == 13) icolor = 1;  // black for muons
		if (ipdg == 15) icolor = 5;  // taus
		if (ipdg == 211) icolor = 2;  // red for pions
		if (ipdg == 321) icolor = 7; // cyan for kaons
		if (ipdg == 2212) icolor = 4; // protons are blue
		if (ipdg == 2112 && showneutrons) icolor = 3; // neutrons are green

		if (icolor>0)
		  {
		    size_t numtpts = mcparticles[imcp].NumberTrajectoryPoints();
		    TGraph2D *gri = new TGraph2D();
	            int ipc=0;
		    for (size_t ipt=0; ipt<numtpts; ++ipt)
		      {
			TLorentzVector pos = mcparticles[imcp].Position(ipt);
			double x = pos.X();
			double y = pos.Y();
			double z = pos.Z();
			if (x<xclipmin || x>xclipmax) continue;
			if (y<yclipmin || y>yclipmax) continue;
			if (z<zclipmin || z>zclipmax) continue;
			gri->SetPoint(ipc,z,x,y);
			++ipc;
			if (ipc>49)
			  {
			    gri->SetMarkerColor(icolor);
			    gri->SetLineColor(icolor);
			    gri->SetMarkerStyle(1);
			    gri->Draw("LINE,SAME");
			    gri = new TGraph2D();
			    ipc = 0;
			    gri->SetPoint(ipc,z,x,y);
			    ++ipc;
			  }
			//cout << ipt << " " << x << " " << y << " " << z << endl;
		      }
		    gri->SetMarkerColor(icolor);
		    gri->SetLineColor(icolor);
		    gri->SetMarkerStyle(1);
		    if (gri->GetN()) gri->Draw("LINE,SAME");
		  }
	      }
	    //c->Update();
	  }
      }
    ++evcounter;
  }
}
