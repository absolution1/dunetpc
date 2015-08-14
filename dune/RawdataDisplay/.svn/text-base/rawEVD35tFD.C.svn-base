//Event display for DUNE 10kt FD
//A histogram data file 'duneEVDraw.root' is required
//to run this script, do this in the root session
//>.x rawEVD35tFD.C++ 
//Click one of the buttons on the bottom to display desired anode plane histograms,
//and click one of the cells on the main dislpay to see the corresponding time vs channel plot
//Nov. 11, 2013, Seongtae Park(UTA)

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <TROOT.h>
#include <TRint.h>
#include <TStyle.h>
#include <TFile.h>
#include <TColor.h>
#include "TExec.h"
#include <TH1.h>
#include <TH2.h>
#include <TGClient.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TSystem.h>
#include <TRandom.h>
#include <TGButton.h>
#include <TGLabel.h>
#include <TLine.h>
#include <TArrow.h>
#include <TLatex.h>
#include <TGFrame.h>
#include <TRootEmbeddedCanvas.h>
#include <RQ_OBJECT.h>
#include "TGFileDialog.h"

unsigned int APA[4]={3,2,1,0};

class MyMainFrame {
	RQ_OBJECT("MyMainFrame")
private:
	TGMainFrame *fMain;
	TRootEmbeddedCanvas *fEcanvas;
	Int_t fpAPA; //remember the previously displayed APA(histograms)
	Int_t fcAPA; //current display APA(histograms)
	std::string fname;
	TFile *fFile;
public:
	MyMainFrame(const TGWindow *p,UInt_t w,UInt_t h);
	virtual ~MyMainFrame();
	string file_open();
	int DrawUplane();
	int DrawVplane();
	int DrawZplane();
	int StackedADCvsTime(TVirtualPad*, TObject*, Int_t);
	void showHits(TVirtualPad* , TObject*, Int_t);
};

MyMainFrame::MyMainFrame(const TGWindow *p,UInt_t w,UInt_t h) {
	fpAPA=-1;
	// Create a main frame
	fMain = new TGMainFrame(p,w,h);
	// Create canvas widget
	fEcanvas = new TRootEmbeddedCanvas("Ecanvas",fMain,700,700);
	fMain->AddFrame(fEcanvas, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 5,5,5,5));
	// Create a horizontal frame widget with buttons
	TGHorizontalFrame *hframe = new TGHorizontalFrame(fMain,200,40);

	TGTextButton *fileOpen = new TGTextButton(hframe,"  &Open  ");
	fileOpen->Connect("Clicked()","MyMainFrame",this,"file_open()");
	hframe->AddFrame(fileOpen, new TGLayoutHints(kLHintsCenterX,	2,15,3,3));

	TGTextButton *drawUplane = new TGTextButton(hframe,"  plot&Uplane  ");
	drawUplane->Connect("Clicked()","MyMainFrame",this,"DrawUplane()");
	hframe->AddFrame(drawUplane, new TGLayoutHints(kLHintsCenterX,	2,2,3,3));

	TGTextButton *drawVplane = new TGTextButton(hframe,"  plot&Vplane  ");
	drawVplane->Connect("Clicked()","MyMainFrame",this,"DrawVplane()");
	hframe->AddFrame(drawVplane, new TGLayoutHints(kLHintsCenterX,	2,2,3,3));

	TGTextButton *drawZplane = new TGTextButton(hframe,"  plot&Zplane  ");
	drawZplane->Connect("Clicked()","MyMainFrame",this,"DrawZplane()");
	hframe->AddFrame(drawZplane, new TGLayoutHints(kLHintsCenterX,	2,2,3,3));

	TGTextButton *exit = new TGTextButton(hframe,"   &Exit   ", "gApplication->Terminate(0)");
	hframe->AddFrame(exit, new TGLayoutHints(kLHintsCenterX,	15,5,3,4));
	fMain->AddFrame(hframe, new TGLayoutHints(kLHintsCenterX, 2,2,2,2));
	// Set a name to the main frame
	fMain->SetWindowName("DUNE Far Detector Event Display");
	// Map all subwindows of main frame
	fMain->MapSubwindows();
	// Initialize the layout algorithm
	fMain->Resize(fMain->GetDefaultSize());
	// Map main frame
	fMain->MapWindow();
}

string MyMainFrame::file_open() {
   TGFileInfo file_info_;
   const char *filetypes[] = {"Root Files", "*.root", 0, 0};
   file_info_.fFileTypes = filetypes;
   file_info_.fIniDir = StrDup(".");
   new TGFileDialog(gClient->GetDefaultRoot(), gClient->GetDefaultRoot(),
                    kFDOpen, &file_info_);
   if( file_info_.fFilename ){
      cout << "'" << file_info_.fFilename << "' selected." << endl;
		fname=file_info_.fFilename;
		fFile=new TFile(fname.c_str());
		DrawUplane();
		return file_info_.fFilename; 
	}
return "";
}

int MyMainFrame::StackedADCvsTime(TVirtualPad*  selpad, TObject* selected, Int_t event) {
	Int_t i, j;
	std::stringstream ytitle, hname;
	std::string temp(""), title;

   if (event != 2||string(selected->IsA()->GetName())!="TH2I") return -1;
   TCanvas *c1 = new TCanvas("c1","ADC vs Time",750, 200,800,600);
   c1->SetLeftMargin(0.14);
   c1->SetRightMargin(0.06);
	temp=selected->GetName(); 
	cout<<"Selected histogram= "<<temp<<endl;

	if ( !fFile->IsOpen() ) return -1; 
	fFile->cd("rawdraw");
	hname.str(""); 
	hname << "rawdraw/"<<temp;
	TH2I* h=(TH2I*)fFile->Get(hname.str().c_str()); 

   TMultiGraph *mg = new TMultiGraph();
	const Int_t nchannels=h->GetNbinsX();
	Int_t nticks=h->GetNbinsY(); 
	std::vector<Double_t> tick;
	std::vector<Double_t> chadc;
	std::vector< std::vector<Double_t> > adc;
	for(i=0;i<nticks;i++) {tick.push_back(i);} 

	for(i=0;i<nchannels;i++){
		for(j=0;j<nticks;j++) chadc.push_back(50*i+h->GetBinContent(i+1,j+1));
		adc.push_back(chadc);
		chadc.clear();
	}
   TGraph *gr;
	 for(i=0;i<nchannels;i++){
		gr= new TGraph(nticks,&tick[0], &adc[i][0]); 
		mg->Add(gr);
		}
 
	title="ADC vs Time for all channels, "+temp;
   mg->Draw("al");
	ytitle<<"Channels ("<<h->GetBinLowEdge(1)<<"~"<<(h->GetBinLowEdge(nchannels)+h->GetBinWidth(nchannels))<<", "<<nchannels<<" channels)";
   mg->SetTitle(title.c_str());
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

 	c1->Update();
return 0;
}

void Pal1() {
   static Int_t  colors[50];
   static Bool_t initialized = kFALSE;

   Double_t Red[5]    = {0.2, 0.0, 1.0, 0.0, 1.0};
   Double_t Green[5]  = {0.2, 0.0, 1.0, 1.0, 0.0};
   Double_t Blue[5]   = {0.2, 1.0, 1.0, 0.0, 0.0};
   Double_t Length[5] = {0.0, .47, .50, .53, 1.0};

   if(!initialized){
      Int_t FI = TColor::CreateGradientColorTable(5,Length,Red,Green,Blue,50);
      for (int i=0; i<50; i++) colors[i] = FI+i;
      initialized = kTRUE;
   }
   gStyle->SetPalette(50,colors);
}

void Pal2() {
   static Int_t  colors[50];
   static Bool_t initialized = kFALSE;

   Double_t Red[3]    = { 1.00, 0.00, 1.00};
   Double_t Green[3]  = { 1.00, 0.00, 0.00};
   Double_t Blue[3]   = { 1.00, 1.00, 0.00};
   Double_t Length[3] = { 0.00, 0.07, 1.00};

   if(!initialized){
      Int_t FI = TColor::CreateGradientColorTable(3,Length,Red,Green,Blue,50);
      for (int i=0; i<50; i++) colors[i] = FI+i;
      initialized = kTRUE;
   }
   gStyle->SetPalette(50,colors);
}

int MyMainFrame::DrawUplane() {
	std::stringstream  hname;
	int i;
	TH2I *htemp;

   TCanvas *fCanvas = (TCanvas*)gROOT->GetListOfCanvases()->FindObject("Ecanvas");
	if(fCanvas) fCanvas->Clear();
	fCanvas = fEcanvas->GetCanvas();
   fCanvas->SetFillColor(42);
//   fCanvas->SetFrameFillColor(33);

	fCanvas->cd();
	fCanvas->Divide(3,1,0.01,0.05); 
	TPad* Ecanvas_1 = (TPad*)(fCanvas->GetPrimitive("Ecanvas_1"));
	TPad* Ecanvas_2 = (TPad*)(fCanvas->GetPrimitive("Ecanvas_2"));
	TPad* Ecanvas_3 = (TPad*)(fCanvas->GetPrimitive("Ecanvas_3"));

	Ecanvas_2->Divide(1,2,0.001,0.001); 

	gStyle->SetNumberContours(64);
//   Int_t MyPalette[20];
   Double_t r[]    = {0.2, 0.0, 1.0, 0.0, 1.0};
   Double_t g[]    = {0.2, 0.0, 1.0, 1.0, 0.0};
   Double_t b[]    = {0.2, 1.0, 1.0, 0.0, 0.0};
   Double_t stop[] = {0.0, .25, .50, .75, 1.0};
//   Int_t FI = TColor::CreateGradientColorTable(5, stop, r, g, b, 20);
   TColor::CreateGradientColorTable(5, stop, r, g, b, 64);

	if ( !fFile->IsOpen() ) return -1;
	fFile->cd("rawdraw");
	hname.str(""); 
	fCanvas->Connect("Selected(TVirtualPad*, TObject*, Int_t)", "MyMainFrame", this,"showHits(TVirtualPad*, TObject*, Int_t)");
	for(i=0;i<4;i++){
		if(i==0) {
			Ecanvas_3->cd(); 
		}
		if(i==1||i==2) {
			Ecanvas_2->cd(i); 
		}
		if(i==3) {
			Ecanvas_1->cd(); 
		}
		hname << "rawdraw/fTimeChanThumbU"<<APA[i];
		htemp=(TH2I*)fFile->Get(hname.str().c_str()); 
		htemp->SetLabelOffset();
		htemp->SetTitle("");
		htemp->Draw("ahcol");
		hname.str("");

		gPad->SetEditable(false);
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
	fMain->SetWindowName("U-Plane View"); 

	return 0;
}

int MyMainFrame::DrawVplane() {
	std::stringstream  hname;
	int i;
	TH2I *htemp;

   TCanvas *fCanvas = (TCanvas*)gROOT->GetListOfCanvases()->FindObject("Ecanvas");
	if(fCanvas) fCanvas->Clear();
	fCanvas = fEcanvas->GetCanvas();
   fCanvas->SetFillColor(42);
//   fCanvas->SetFrameFillColor(33);

	fCanvas->cd();
	fCanvas->Divide(3,1,0.01,0.05); 
	TPad* Ecanvas_1 = (TPad*)(fCanvas->GetPrimitive("Ecanvas_1"));
	TPad* Ecanvas_2 = (TPad*)(fCanvas->GetPrimitive("Ecanvas_2"));
	TPad* Ecanvas_3 = (TPad*)(fCanvas->GetPrimitive("Ecanvas_3"));

	Ecanvas_2->Divide(1,2,0.001,0.001); 

//   Int_t MyPalette[20];
	gStyle->SetNumberContours(64);
   Double_t r[]    = {0.2, 0.0, 1.0, 0.0, 1.0};
   Double_t g[]    = {0.2, 0.0, 1.0, 1.0, 0.0};
   Double_t b[]    = {0.2, 1.0, 1.0, 0.0, 0.0};
   Double_t stop[] = {0.0, .25, .50, .75, 1.0};
//   Int_t FI = TColor::CreateGradientColorTable(5, stop, r, g, b, 20);
   TColor::CreateGradientColorTable(5, stop, r, g, b, 64);

	if ( !fFile->IsOpen() ) return -1;
	fFile->cd("rawdraw");

	hname.str(""); 
	fCanvas->Connect("Selected(TVirtualPad*, TObject*, Int_t)", "MyMainFrame", this,"showHits(TVirtualPad*, TObject*, Int_t)");
	for(i=0;i<4;i++){
		if(i==0) {
			Ecanvas_3->cd(); 
		}
		if(i==1||i==2) {
			Ecanvas_2->cd(i); 
		}
		if(i==3) {
			Ecanvas_1->cd(); 
		}
		hname << "rawdraw/fTimeChanThumbV"<<APA[i];
		htemp=(TH2I*)fFile->Get(hname.str().c_str()); 
		htemp->SetLabelOffset();
		htemp->SetTitle("");
		htemp->Draw("ahcol");
		hname.str("");

		gPad->SetEditable(false);
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
	fMain->SetWindowName("V-Plane View");
	return 0;
}

int MyMainFrame::DrawZplane() {
	std::stringstream  hname;
	int i,j;
	TH2I *htemp;

   TCanvas *fCanvas = (TCanvas*)gROOT->GetListOfCanvases()->FindObject("Ecanvas");
	if(fCanvas) fCanvas->Clear();
	fCanvas = fEcanvas->GetCanvas();
   fCanvas->SetFillColor(42);
//   fCanvas->SetFrameFillColor(33);

	fCanvas->cd();
	fCanvas->Divide(3,1,0.001,0.05); 
	TPad* Ecanvas_1 = (TPad*)(fCanvas->GetPrimitive("Ecanvas_1"));
	TPad* Ecanvas_2 = (TPad*)(fCanvas->GetPrimitive("Ecanvas_2"));
	TPad* Ecanvas_3 = (TPad*)(fCanvas->GetPrimitive("Ecanvas_3"));
/*	Ecanvas_1->SetLeftMargin(0.01);
	Ecanvas_1->SetRightMargin(0.01);
	Ecanvas_1->SetTopMargin(0.05);
	Ecanvas_1->SetBottomMargin(0.05);
	Ecanvas_2->SetLeftMargin(0.01);
	Ecanvas_2->SetRightMargin(0.01);
	Ecanvas_2->SetTopMargin(0.05);
	Ecanvas_2->SetBottomMargin(0.05);
	Ecanvas_3->SetLeftMargin(0.01);
	Ecanvas_3->SetRightMargin(0.01);
	Ecanvas_3->SetTopMargin(0.05);
	Ecanvas_3->SetBottomMargin(0.05);
*/
	Ecanvas_1->Divide(2,1,0.001,0.001); 
	Ecanvas_2->Divide(2,2,0.001,0.001); 
	Ecanvas_3->Divide(2,1,0.001,0.001); 

	gStyle->SetNumberContours(64);
   Double_t r[]    = {0.2, 0.0, 1.0, 0.0, 1.0};
   Double_t g[]    = {0.2, 0.0, 1.0, 1.0, 0.0};
   Double_t b[]    = {0.2, 1.0, 1.0, 0.0, 0.0};
   Double_t stop[] = {0.0, .25, .50, .75, 1.0};
   TColor::CreateGradientColorTable(5, stop, r, g, b, 64);

	if ( !fFile->IsOpen() ) return -1;
	fFile->cd("rawdraw");

	hname.str(""); 
	fCanvas->Connect("Selected(TVirtualPad*, TObject*, Int_t)", "MyMainFrame", this,"showHits(TVirtualPad*, TObject*, Int_t)");
	for(i=0;i<4;i++){
		if(i==0) {
			for(j=1;j<=2;j++){
				Ecanvas_3->cd(j); 
				if(j==1) hname << "rawdraw/fTimeChanThumbZ1"<<APA[i];
					else hname << "rawdraw/fTimeChanThumbZ0"<<APA[i];
				htemp=(TH2I*)fFile->Get(hname.str().c_str()); 
				htemp->SetLabelOffset();
				htemp->SetTitle("");
				htemp->Draw("ahcol");
				hname.str("");
				gPad->SetEditable(false);
			}
		}
		if(i==1) {
			for(j=1;j<=2;j++){
				Ecanvas_2->cd(j); 
				if(j==1) hname << "rawdraw/fTimeChanThumbZ1"<<APA[i];
					else hname << "rawdraw/fTimeChanThumbZ0"<<APA[i];
				htemp=(TH2I*)fFile->Get(hname.str().c_str()); 
				htemp->SetLabelOffset();
				htemp->SetTitle("");
				htemp->Draw("ahcol");
				hname.str("");
				gPad->SetEditable(false);
			}
		}
		if(i==2) {
			for(j=1;j<=2;j++){
				Ecanvas_2->cd(i+j); 
				if(j==1) hname << "rawdraw/fTimeChanThumbZ1"<<APA[i];
					else hname << "rawdraw/fTimeChanThumbZ0"<<APA[i];
				htemp=(TH2I*)fFile->Get(hname.str().c_str()); 
				htemp->SetLabelOffset();
				htemp->SetTitle("");
				htemp->Draw("ahcol");
				hname.str("");
				gPad->SetEditable(false);
			}
		}
		if(i==3) {
			for(j=1;j<=2;j++){
				Ecanvas_1->cd(j); 
				if(j==1) hname << "rawdraw/fTimeChanThumbZ1"<<APA[i];
					else hname << "rawdraw/fTimeChanThumbZ0"<<APA[i];
				htemp=(TH2I*)fFile->Get(hname.str().c_str()); 
				htemp->SetLabelOffset();
				htemp->SetTitle("");
				htemp->Draw("ahcol");
				hname.str("");
				gPad->SetEditable(false);
			}
		}

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
	fMain->SetWindowName("Z-Plane View");
	return 0;
}

MyMainFrame::~MyMainFrame() {
// Clean up used widgets: frames, buttons, layout hints
fMain->Cleanup();
delete fMain;
}

void MyMainFrame::showHits(TVirtualPad*  selpad, TObject* selected, Int_t event)
{
	std::stringstream  cmd;
	std::string apa("");
	Int_t found;

//   int event = gPad->GetEvent();
   if (event != 1||string(selected->IsA()->GetName())!="TH2I") return;
	apa=selected->GetName(); 
	found=apa.find("U",0); if(found!=-1) {apa.erase(0,found+1); fcAPA=atoi(apa.c_str()); cout<<endl<<"Selected APA= "<<apa<<endl;}
	found=apa.find("V",0); if(found!=-1) {apa.erase(0,found+1); fcAPA=atoi(apa.c_str()); cout<<endl<<"Selected APA= "<<apa<<endl;}
	found=apa.find("Z0",0); if(found!=-1) {apa.erase(0,found+2); fcAPA=atoi(apa.c_str()); cout<<endl<<"Selected APA= "<<apa<<endl;}
	found=apa.find("Z1",0); if(found!=-1) {apa.erase(0,found+2); fcAPA=atoi(apa.c_str()); cout<<endl<<"Selected APA= "<<apa<<endl;}
	
	if(fpAPA!=-1){
		cmd << "gDirectory"<<"->Delete(\"fTimeChanU"<<fpAPA<<"\");";
		gROOT->ProcessLine(cmd.str().c_str()); cmd.str("");
		cmd << "gDirectory"<<"->Delete(\"fTimeChanV"<<fpAPA<<"\");";
		gROOT->ProcessLine(cmd.str().c_str()); cmd.str("");
		cmd << "gDirectory"<<"->Delete(\"fTimeChanZ0"<<fpAPA<<"\");";
		gROOT->ProcessLine(cmd.str().c_str()); cmd.str("");
		cmd << "gDirectory"<<"->Delete(\"fTimeChanZ1"<<fpAPA<<"\");";
		gROOT->ProcessLine(cmd.str().c_str()); cmd.str("");
	}
   TCanvas *c2 = (TCanvas*)gROOT->GetListOfCanvases()->FindObject("c2");
	if(c2) c2->Close();
   c2 = new TCanvas("c2","Time vs Channel",715,50,730,740);
	c2->Divide(1,4,0.015,0.015);

	cout<<"Ploted histograms: "<<endl;
   c2->cd(1);
	c2->cd(1)->SetTickx(1);c2->cd(1)->SetTicky(1);
	cmd <<"fTimeChanU"<<fcAPA<<"->SetMinimum(-500);";
	gROOT->ProcessLine(cmd.str().c_str()); cmd.str("");
	cmd <<"fTimeChanU"<<fcAPA<<"->SetMaximum(500);";
	gROOT->ProcessLine(cmd.str().c_str()); cmd.str("");
	cmd <<"fTimeChanU"<<fcAPA<<"->Draw(\"colz\")";
	cout<<"fTimeChanU"<<fcAPA<<endl;
	gROOT->ProcessLine(cmd.str().c_str()); cmd.str("");
   TExec *ex1 = new TExec("ex1","Pal1();");
   ex1->Draw();
	cmd <<"fTimeChanU"<<fcAPA<<"->Draw(\"colz same\")";
	gROOT->ProcessLine(cmd.str().c_str()); cmd.str("");
	gPad->RedrawAxis();

   c2->cd(2);
	c2->cd(2)->SetTickx(1);c2->cd(2)->SetTicky(1);
	cmd <<"fTimeChanV"<<fcAPA<<"->SetMinimum(-500);";
	gROOT->ProcessLine(cmd.str().c_str()); cmd.str("");
	cmd <<"fTimeChanV"<<fcAPA<<"->SetMaximum(500);";
	gROOT->ProcessLine(cmd.str().c_str()); cmd.str("");
	cmd <<"fTimeChanV"<<fcAPA<<"->Draw(\"colz\")";
	cout<<"fTimeChanV"<<fcAPA<<endl;
	gROOT->ProcessLine(cmd.str().c_str()); cmd.str("");
   TExec *ex2 = new TExec("ex2","Pal1();");
   ex2->Draw();
	cmd <<"fTimeChanV"<<fcAPA<<"->Draw(\"colz same\")";
	gROOT->ProcessLine(cmd.str().c_str()); cmd.str("");
	gPad->RedrawAxis();

   c2->cd(3);
	c2->cd(3)->SetTickx(1);c2->cd(3)->SetTicky(1);
	cmd <<"fTimeChanZ0"<<fcAPA<<"->SetMinimum(0);";
	gROOT->ProcessLine(cmd.str().c_str()); cmd.str("");
	cmd <<"fTimeChanZ0"<<fcAPA<<"->SetMaximum(500);";
	gROOT->ProcessLine(cmd.str().c_str()); cmd.str("");
	cmd <<"fTimeChanZ0"<<fcAPA<<"->Draw(\"colz\")";
	cout<<"fTimeChanZ0"<<fcAPA<<endl;
	gROOT->ProcessLine(cmd.str().c_str()); cmd.str("");
   TExec *ex3 = new TExec("ex3","Pal2();");
   ex3->Draw();
	cmd <<"fTimeChanZ0"<<fcAPA<<"->Draw(\"colz same\")";
	gROOT->ProcessLine(cmd.str().c_str()); cmd.str("");
	gPad->RedrawAxis();

   c2->cd(4);
	c2->cd(4)->SetTickx(1);c2->cd(4)->SetTicky(1);
	cmd <<"fTimeChanZ1"<<fcAPA<<"->SetMinimum(0);";
	gROOT->ProcessLine(cmd.str().c_str()); cmd.str("");
	cmd <<"fTimeChanZ1"<<fcAPA<<"->SetMaximum(500);";
	gROOT->ProcessLine(cmd.str().c_str()); cmd.str("");
	cmd <<"fTimeChanZ1"<<fcAPA<<"->Draw(\"colz\")";
	cout<<"fTimeChanZ1"<<fcAPA<<endl;
	gROOT->ProcessLine(cmd.str().c_str()); cmd.str("");
   TExec *ex4 = new TExec("ex4","Pal2();");
   ex4->Draw();
	cmd <<"fTimeChanZ1"<<fcAPA<<"->Draw(\"colz same\")";
	gROOT->ProcessLine(cmd.str().c_str()); cmd.str("");
	gPad->RedrawAxis();

	fpAPA=fcAPA;
	c2->Connect("Selected(TVirtualPad*, TObject*, Int_t)", "MyMainFrame", this,"StackedADCvsTime(TVirtualPad*, TObject*, Int_t)");
}

void rawEVD35tFD() {
// Popup the GUI...
new MyMainFrame(gClient->GetRoot(),200,200);


}


