#define MakeHistograms_cxx
#include "MakeHistograms.h"
#include <TStyle.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TLatex.h>

#include <iostream>
#include <sstream>
#include <vector>

void MakeHistograms::PrintPlots(UInt_t eventnb)
{
//   In a ROOT session, you can do:
//      Root > .L MakeHistograms.C
//      Root > MakeHistograms t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.PrintPlots(eventnb);       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

	unsigned int nch;
	unsigned int nADC;
	unsigned int pn;
	std::stringstream  rootfn;
	fEvent=eventnb;

	for(int l=0;l<4;l++) {
   	fTimeChanU[l]->Reset();
   	fTimeChanV[l]->Reset();
   	fTimeChanZ0[l]->Reset();
   	fTimeChanZ1[l]->Reset();
	}

   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
	
	if(Event!=eventnb) continue;

	nch=Chan->size();
	for(unsigned int j=0;j<nch;j++) {
		nADC=(*ADC)[j].size();
		pn=Plane->at(j);
		if(nADC<1) continue;
		for(unsigned int k=0;k<nADC;k++) {
				switch(pn){
					case 0: fTimeChanU[APA->at(j)]->Fill(Chan->at(j),(*TDC)[j].at(k),(*ADC)[j].at(k)); break;
					case 1: fTimeChanV[APA->at(j)]->Fill(Chan->at(j),(*TDC)[j].at(k),(*ADC)[j].at(k)); break;
					case 2: fTimeChanZ0[APA->at(j)]->Fill(Chan->at(j),(*TDC)[j].at(k),(*ADC)[j].at(k)); break;
					case 3: fTimeChanZ1[APA->at(j)]->Fill(Chan->at(j),(*TDC)[j].at(k),(*ADC)[j].at(k)); break;
					default: break;
				}
			}
		}
   }

	// Save histograms to a root file
	rootfn.str("");
	rootfn<<"TimevsChanHistos_E"<<fEvent<<".root";
   TFile f(rootfn.str().c_str(),"recreate");
	for(int l=0;l<4;l++) {
   	fTimeChanU[l]->Write();
   	fTimeChanV[l]->Write();
   	fTimeChanZ0[l]->Write();
   	fTimeChanZ1[l]->Write();
	}
	DrawUplane();
	DrawVplane();
	DrawZplane();
	PrintHistos();
	PrintGraphs(fTimeChanU);
	PrintGraphs(fTimeChanV);
	PrintGraphs(fTimeChanZ0);
	PrintGraphs(fTimeChanZ1);
	f.Close();
}

void MakeHistograms::HistDef()
{
	//Histogram names and titles
	std::stringstream  name, title;
   TH2I* TempHisto;

	UInt_t NChPerAPA;
   unsigned int UChMin;
   unsigned int UChMax;
   unsigned int VChMin;
   unsigned int VChMax;
   unsigned int Z0ChMin;
   unsigned int Z0ChMax;
   unsigned int Z1ChMin;
   unsigned int Z1ChMax;
   unsigned int minT=0, maxT=3200;
   unsigned int binT = (maxT-minT);

	LoadTree(0); 	GetEntry(0);
	maxT=Nticks;
	NChPerAPA=NofUChan+NofVChan+NofZ0Chan+NofZ1Chan;

 for(unsigned int i=0;i<4;i++){
    UChMin=i*NChPerAPA;
    UChMax=UChMin + NofUChan-1;
    VChMin=UChMax+1;
    VChMax=VChMin + NofVChan-1;
    Z0ChMin=VChMax+1;
    Z0ChMax=Z0ChMin+NofZ0Chan-1;
    Z1ChMin=Z0ChMax+1;
    Z1ChMax=Z1ChMin+NofZ1Chan-1;

    // construct the histograms; TH2 constructors: ("Name", "Title", NxBin, xMin, xMax, NyBin, yMax, yMin)
    name.str("");
    name << "TimeChanU";
    name <<  i;
    title.str("");
    title << "Time vs Channel(Plane U, APA";
    title << i<<")";
	 TempHisto = new TH2I(name.str().c_str(),title.str().c_str(), NofUChan, UChMin, UChMax, binT, minT, maxT);
	 fTimeChanU.push_back(TempHisto);

    name.str("");
    name << "TimeChanV";
    name << i;
    title.str("");
    title << "Time vs Channel(Plane V, APA";
    title << i<<")";
    TempHisto = new TH2I(name.str().c_str(),title.str().c_str(), NofVChan, VChMin, VChMax, binT, minT, maxT);
	 fTimeChanV.push_back(TempHisto);

    name.str("");
    name << "TimeChanZ0";
    name << i;
    title.str("");
    title << "Time vs Channel(Plane Z0, APA";
    title <<i<<")";
    TempHisto = new TH2I(name.str().c_str(),title.str().c_str(), NofZ0Chan, Z0ChMin, Z0ChMax, binT, minT, maxT);
	 fTimeChanZ0.push_back(TempHisto);

    name.str("");
    name << "TimeChanZ1";
    name << i;
    title.str("");
    title << "Time vs Channel(Plane Z1, APA";
    title << i<<")";
    TempHisto = new TH2I(name.str().c_str(),title.str().c_str(), NofZ1Chan, Z1ChMin, Z1ChMax, binT, minT, maxT);
	 fTimeChanZ1.push_back(TempHisto);

    fTimeChanU[i]->SetStats(0);    fTimeChanV[i]->SetStats(0);    
	 fTimeChanZ0[i]->SetStats(0);    fTimeChanZ1[i]->SetStats(0);

    fTimeChanU[i]->GetXaxis()->SetTitle("Channel"); fTimeChanU[i]->GetYaxis()->SetTitle("TDC");
    fTimeChanV[i]->GetXaxis()->SetTitle("Channel"); fTimeChanV[i]->GetYaxis()->SetTitle("TDC");
    fTimeChanZ0[i]->GetXaxis()->SetTitle("Channel"); fTimeChanZ0[i]->GetYaxis()->SetTitle("TDC");
    fTimeChanZ1[i]->GetXaxis()->SetTitle("Channel"); fTimeChanZ1[i]->GetYaxis()->SetTitle("TDC");
	}

}

void MakeHistograms::Pal1() {
   const Int_t NRGBs = 5;
   const Int_t NCont = 255;

   Double_t Red[5]    = {0.2, 0.0, 1.0, 0.0, 1.0};
   Double_t Green[5]  = {0.2, 0.0, 1.0, 1.0, 0.0};
   Double_t Blue[5]   = {0.2, 1.0, 1.0, 0.0, 0.0};
   Double_t Length[5] = {0.0, .47, .50, .53, 1.0};

   TColor::CreateGradientColorTable(NRGBs, Length, Red, Green, Blue, NCont);
   gStyle->SetNumberContours(NCont);
}

void MakeHistograms::Pal2() {
   const Int_t NRGBs = 3;
   const Int_t NCont = 255;
   Double_t Red[3]    = { 0.00, 0.00, 1.00};
   Double_t Green[3]  = { 0.00, 1.00, 0.00};
   Double_t Blue[3]   = { 1.00, 0.00, 0.00};
   Double_t Length[3] = { 0.00, 0.50, 1.00};
   TColor::CreateGradientColorTable(NRGBs, Length, Red, Green, Blue, NCont);
   gStyle->SetNumberContours(NCont);
}

void MakeHistograms::PrintHistos() {
	TCanvas *c1 = new TCanvas("c1","DUNE 35t detector event display",800,600);
	std::stringstream  filename;
	c1->SetTickx(1);c1->SetTicky(1);

	for(int l=0;l<4;l++) {
		Pal1();
   	fTimeChanU[l]->SetMinimum(-500);
   	fTimeChanU[l]->SetMaximum(500);
		fTimeChanU[l]->GetXaxis()->CenterTitle(kTRUE);
		fTimeChanU[l]->GetYaxis()->CenterTitle(kTRUE);
		fTimeChanU[l]->GetXaxis()->SetTitleOffset(1.0);
		fTimeChanU[l]->GetYaxis()->SetTitleOffset(1.3);
		fTimeChanU[l]->Draw("colz");
		filename.str("");
		filename<<"TimevsChan_E"<<fEvent<<"_A"<<l<<"_U.png";
		gPad->RedrawAxis();
		c1->Print(filename.str().c_str());

//		Pal1();
   	fTimeChanV[l]->SetMinimum(-500);
   	fTimeChanV[l]->SetMaximum(500);
		fTimeChanV[l]->GetXaxis()->CenterTitle(kTRUE);
		fTimeChanV[l]->GetYaxis()->CenterTitle(kTRUE);
		fTimeChanV[l]->GetXaxis()->SetTitleOffset(1.0);
		fTimeChanV[l]->GetYaxis()->SetTitleOffset(1.3);
		fTimeChanV[l]->Draw("colz");
		filename.str("");
		filename<<"TimevsChan_E"<<fEvent<<"_A"<<l<<"_V.png";
		gPad->RedrawAxis();
		c1->Print(filename.str().c_str());

		Pal2();
   	fTimeChanZ0[l]->SetMinimum(0);
   	fTimeChanZ0[l]->SetMaximum(500);
		fTimeChanZ0[l]->GetXaxis()->CenterTitle(kTRUE);
		fTimeChanZ0[l]->GetYaxis()->CenterTitle(kTRUE);
		fTimeChanZ0[l]->GetXaxis()->SetTitleOffset(1.0);
		fTimeChanZ0[l]->GetYaxis()->SetTitleOffset(1.3);
		fTimeChanZ0[l]->Draw("colz");
		filename.str("");
		filename<<"TimevsChan_E"<<fEvent<<"_A"<<l<<"_Z0.png";
		gPad->RedrawAxis();
		c1->Print(filename.str().c_str());

//		Pal2();
   	fTimeChanZ1[l]->SetMinimum(0);
   	fTimeChanZ1[l]->SetMaximum(500);
		fTimeChanZ1[l]->GetXaxis()->CenterTitle(kTRUE);
		fTimeChanZ1[l]->GetYaxis()->CenterTitle(kTRUE);
		fTimeChanZ1[l]->GetXaxis()->SetTitleOffset(1.0);
		fTimeChanZ1[l]->GetYaxis()->SetTitleOffset(1.3);
		fTimeChanZ1[l]->Draw("colz");
		filename.str("");
		filename<<"TimevsChan_E"<<fEvent<<"_A"<<l<<"_Z1.png";
		gPad->RedrawAxis();
		c1->Print(filename.str().c_str());
	}
		c1->Close();
}

void MakeHistograms::PrintGraphs(std::vector<TH2I*> h) {
	TCanvas *c2 = new TCanvas("c2","DUNE 35t detector event display",800,600);
	std::stringstream  filename,ytitle,title;
	Int_t i, j;
	std::string temp, pn;

	pn=h[0]->GetName(); 
	if(pn.find("U")!=std::string::npos) pn="U";
	if(pn.find("V")!=std::string::npos) pn="V";
	if(pn.find("Z0")!=std::string::npos) pn="Z0";
	if(pn.find("Z1")!=std::string::npos) pn="Z1";

   c2->SetLeftMargin(0.14);
   c2->SetRightMargin(0.06);
	c2->SetTickx(1);c2->SetTicky(1);

	for(int l=0;l<4;l++) {
		filename.str("");
		filename<<"ADCvsChan_E"<<fEvent<<"_A"<<l<<"_"<<pn<<".png";
		TMultiGraph *mg = new TMultiGraph();
		const Int_t nchannels=h[l]->GetNbinsX();
		Int_t nticks=h[l]->GetNbinsY(); 
		std::vector<Double_t> tick;
		std::vector<Double_t> chadc;
		std::vector< std::vector<Double_t> > adc;
		for(i=0;i<nticks;i++) {tick.push_back(i);} 

		for(i=0;i<nchannels;i++){
			for(j=0;j<nticks;j++) chadc.push_back(50*i+h[l]->GetBinContent(i+1,j+1));
			adc.push_back(chadc);
			chadc.clear();
		}
		TGraph *gr;
		 for(i=0;i<nchannels;i++){
			gr= new TGraph(nticks,&tick[0], &adc[i][0]); 
			mg->Add(gr);
			}
		title.str("");
		title<<"ADC vs Time for all channels, E"<<fEvent<<"_A"<<l<<"_"<<pn;
		mg->Draw("al");
		ytitle.str("");
		ytitle<<"Channels ("<<h[l]->GetBinLowEdge(1)<<"~"<<(h[l]->GetBinLowEdge(nchannels)+h[l]->GetBinWidth(nchannels))<<", "<<nchannels<<" channels)";
		mg->SetTitle(title.str().c_str());
	//   mg->GetXaxis()->SetNdivisions(32, kTRUE);
		mg->GetXaxis()->SetTitle("Time (ticks)");
		mg->GetXaxis()->SetTitleSize(0.03);
		mg->GetXaxis()->SetTitleOffset(1.0);
		mg->GetXaxis()->CenterTitle(kTRUE);
	//   mg->GetXaxis()->SetRangeUser(1,64);
		mg->GetXaxis()->SetLabelSize(0.03);
	//   mg->GetXaxis()->SetLabelOffset(0.05);

		mg->GetYaxis()->SetTitle(ytitle.str().c_str());
		mg->GetYaxis()->SetTitleSize(0.03);
		mg->GetYaxis()->SetTitleOffset(1.5);
		mg->GetYaxis()->CenterTitle(kTRUE);
	//   mg->GetYaxis()->SetRangeUser(0,500); 
		mg->GetYaxis()->SetLabelSize(0.03);
	//   mg->GetYaxis()->SetLabelOffset(0.01);
	 	c2->Update();
		c2->Print(filename.str().c_str());
		mg->Delete();
	}
		c2->Close();
}

void MakeHistograms::DrawUplane() {
	std::stringstream  filename;
	int i;

	fCanvas = new TCanvas("Ecanvas","DUNE 35t detector event display",700,700);
   fCanvas->SetFillColor(42);
//   fCanvas->SetFrameFillColor(33);

	fCanvas->Divide(3,1,0.01,0.05); 
	TPad* Ecanvas_1 = (TPad*)(fCanvas->GetPrimitive("Ecanvas_1"));
	TPad* Ecanvas_2 = (TPad*)(fCanvas->GetPrimitive("Ecanvas_2"));
	TPad* Ecanvas_3 = (TPad*)(fCanvas->GetPrimitive("Ecanvas_3"));
	Ecanvas_2->Divide(1,2,0.001,0.001); 

	Pal1();

	for(i=0;i<4;i++){
		if(i==0) {
			Ecanvas_1->cd(); 
		}
		if(i==1||i==2) {
			Ecanvas_2->cd(i); 
		}
		if(i==3) {
			Ecanvas_3->cd(); 
		}
   	fTimeChanU[i]->SetMinimum(0);
   	fTimeChanU[i]->SetMaximum(500);
		fTimeChanU[i]->SetTitle("");
		fTimeChanU[i]->SetTitle("");				
		fTimeChanU[i]->Draw("ahcol");
		gPad->SetEditable(false);
		gPad->RedrawAxis();
	}
   fCanvas->cd();
   TLatex *tex = new TLatex(0.45,0.03,"APA2/1");
   tex->SetTextSize(0.03);
   tex->SetTextAngle(0);
   tex->SetLineWidth(2);
   tex->Draw();
   tex = new TLatex(0.13,0.03,"APA0");
   tex->SetTextSize(0.03);
   tex->SetLineWidth(2);
   tex->Draw();
   tex = new TLatex(0.8,0.03,"APA3");
   tex->SetTextSize(0.03);
   tex->SetLineWidth(2);
   tex->Draw();

   fCanvas->Modified();
	fCanvas->Update();

	filename.str("");
	filename<<"UPlaneView_E"<<fEvent<<".png";
	fCanvas->Print(filename.str().c_str());
}

void MakeHistograms::DrawVplane() {
	std::stringstream  filename;
	int i;

	fCanvas = new TCanvas("Ecanvas","DUNE 35t detector event display",700,700);
   fCanvas->SetFillColor(42);
//   fCanvas->SetFrameFillColor(33);

	fCanvas->Divide(3,1,0.01,0.05); 
	TPad* Ecanvas_1 = (TPad*)(fCanvas->GetPrimitive("Ecanvas_1"));
	TPad* Ecanvas_2 = (TPad*)(fCanvas->GetPrimitive("Ecanvas_2"));
	TPad* Ecanvas_3 = (TPad*)(fCanvas->GetPrimitive("Ecanvas_3"));
	Ecanvas_2->Divide(1,2,0.001,0.001); 

	Pal1();

	for(i=0;i<4;i++){
		if(i==0) {
			Ecanvas_1->cd(); 
		}
		if(i==1||i==2) {
			Ecanvas_2->cd(i); 
		}
		if(i==3) {
			Ecanvas_3->cd(); 
		}
   	fTimeChanV[i]->SetMinimum(0);
   	fTimeChanV[i]->SetMaximum(500);
		fTimeChanV[i]->SetTitle("");
		fTimeChanV[i]->SetTitle("");				
		fTimeChanV[i]->Draw("ahcol");
		gPad->SetEditable(false);
		gPad->RedrawAxis();
	}
   fCanvas->cd();
   TLatex *tex = new TLatex(0.45,0.03,"APA2/1");
   tex->SetTextSize(0.03);
   tex->SetTextAngle(0);
   tex->SetLineWidth(2);
   tex->Draw();
   tex = new TLatex(0.13,0.03,"APA0");
   tex->SetTextSize(0.03);
   tex->SetLineWidth(2);
   tex->Draw();
   tex = new TLatex(0.8,0.03,"APA3");
   tex->SetTextSize(0.03);
   tex->SetLineWidth(2);
   tex->Draw();

   fCanvas->Modified();
	fCanvas->Update();

	filename.str("");
	filename<<"VPlaneView_E"<<fEvent<<".png";
	fCanvas->Print(filename.str().c_str());
}


void MakeHistograms::DrawZplane() {
	std::stringstream  filename;
	int i,j;

	fCanvas = new TCanvas("Ecanvas","DUNE 35t detector event display",700,700);
   fCanvas->SetFillColor(42);
//   fCanvas->SetFrameFillColor(33);

	fCanvas->Divide(3,1,0.001,0.05); 
	TPad* Ecanvas_1 = (TPad*)(fCanvas->GetPrimitive("Ecanvas_1"));
	TPad* Ecanvas_2 = (TPad*)(fCanvas->GetPrimitive("Ecanvas_2"));
	TPad* Ecanvas_3 = (TPad*)(fCanvas->GetPrimitive("Ecanvas_3"));

	Ecanvas_1->Divide(2,1,0.001,0.001); 
	Ecanvas_2->Divide(2,2,0.001,0.001); 
	Ecanvas_3->Divide(2,1,0.001,0.001); 

	Pal2();
	filename.str("");

	for(i=0;i<4;i++){
		if(i==0) {
			for(j=1;j<=2;j++){
				Ecanvas_1->cd(j); 
				if(j==1) fTimeChanZ1[i]->Draw("ahcol");
					else fTimeChanZ0[i]->Draw("ahcol");
				fTimeChanZ0[i]->SetTitle("");
				fTimeChanZ1[i]->SetTitle("");				
				gPad->SetEditable(false);
			}
		}
		if(i==1) {
			for(j=1;j<=2;j++){
				Ecanvas_2->cd(j); 
				if(j==1) fTimeChanZ1[i]->Draw("ahcol");
					else fTimeChanZ0[i]->Draw("ahcol");
				fTimeChanZ0[i]->SetTitle("");
				fTimeChanZ1[i]->SetTitle("");				
				gPad->SetEditable(false);
			}
		}
		if(i==2) {
			for(j=1;j<=2;j++){
				Ecanvas_2->cd(i+j); 
				if(j==1) fTimeChanZ1[i]->Draw("ahcol");
					else fTimeChanZ0[i]->Draw("ahcol");
				fTimeChanZ0[i]->SetTitle("");
				fTimeChanZ1[i]->SetTitle("");				
				gPad->SetEditable(false);
			}
		}
		if(i==3) {
			for(j=1;j<=2;j++){
				Ecanvas_3->cd(j); 
				if(j==1) fTimeChanZ1[i]->Draw("ahcol");
					else fTimeChanZ0[i]->Draw("ahcol");
				fTimeChanZ0[i]->SetTitle("");
				fTimeChanZ1[i]->SetTitle("");				
				gPad->SetEditable(false);
			}
		}
   	fTimeChanZ0[i]->SetMinimum(0);
   	fTimeChanZ0[i]->SetMaximum(500);
   	fTimeChanZ1[i]->SetMinimum(0);
   	fTimeChanZ1[i]->SetMaximum(500);
		gPad->RedrawAxis();
	}

   fCanvas->cd();
   TLatex *tex = new TLatex(0.45,0.03,"APA2/1");
   tex->SetTextSize(0.03);
   tex->SetTextAngle(0);
   tex->SetLineWidth(2);
   tex->Draw();
   tex = new TLatex(0.13,0.03,"APA0");
   tex->SetTextSize(0.03);
   tex->SetLineWidth(2);
   tex->Draw();
   tex = new TLatex(0.8,0.03,"APA3");
   tex->SetTextSize(0.03);
   tex->SetLineWidth(2);
   tex->Draw();

   fCanvas->Modified();
	fCanvas->Update();
	filename.str("");
	filename<<"ZPlaneView_E"<<fEvent<<".png";
	fCanvas->Print(filename.str().c_str());
}



