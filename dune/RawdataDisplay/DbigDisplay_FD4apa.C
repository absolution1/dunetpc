//Event display for DUNE FD Workspace disambiguation
//Dec. 23, 2013, Tyler Alion
// tylerdalion@gmail.com


#include <iostream>
#include <sstream>
#include <string>
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
#include "TGStatusBar.h"
#include "TGNumberEntry.h"
#include "TGButtonGroup.h"


class MyMainFrame {
  RQ_OBJECT("MyMainFrame")
  
private:

  TFile               *fFile;
  TGMainFrame         *fMain;
  TRootEmbeddedCanvas *fEcanvas;
  TGHorizontalFrame   *hframe;
  TGStatusBar         *fStatusBar, *fPlotLabels;
  TGNumberEntry       *fEVT;
  TGButtonGroup       *APAbgr;

  TCanvas *fCanvas;
  TPad* fEvanvas_1;
  TPad* fEvanvas_2;
  TPad* fEvanvas_3;
  TPad* fEvanvas_4;

  unsigned int         CanvWidth = 1400;
  int                  fEvent = 1;
  int                  fAPA   = 0;
  std::string          fname;

public:

  MyMainFrame(const TGWindow *p,UInt_t w,UInt_t h);
  virtual ~MyMainFrame();
  void pickAPA0();
  void pickAPA1();
  void pickAPA2();
  void pickAPA3();
  int DrawDisambig();
  int DrawCheated();
  int DrawIncorrect();
  int DrawZBot();
  int DrawZTop();
  int NextEvt();
  int ChooseEvt();
  int PrevEvt();
  void ZoomDraw(TH2I& h);
  std::string file_open();
  
};




//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
MyMainFrame::MyMainFrame(const TGWindow *p,UInt_t w,UInt_t h) {
  
  fMain    = new TGMainFrame(p,w,h);
  fEcanvas = new TRootEmbeddedCanvas("Ecanvas",fMain,CanvWidth,700);
  fMain->AddFrame(fEcanvas, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 5,5,5,5));
  hframe = new TGHorizontalFrame(fMain,CanvWidth*0.8,0);
  
  TGRadioButton* APA[4];
  APAbgr = new TGButtonGroup( hframe,"APA",kHorizontalFrame); 
  APA[0] = new TGRadioButton( APAbgr, new TGHotString("&0") ); 
  APA[0]->Connect("Clicked()","MyMainFrame",this,"pickAPA0()");
  APA[1] = new TGRadioButton( APAbgr, new TGHotString("&1") ); 
  APA[1]->Connect("Clicked()","MyMainFrame",this,"pickAPA1()");
  APA[2] = new TGRadioButton( APAbgr, new TGHotString("&2") ); 
  APA[2]->Connect("Clicked()","MyMainFrame",this,"pickAPA2()");
  APA[3] = new TGRadioButton( APAbgr, new TGHotString("&3") ); 
  APA[3]->Connect("Clicked()","MyMainFrame",this,"pickAPA3()");
  APA[0]->SetState(kButtonDown);

  TGTextButton *fileOpen = new TGTextButton(hframe," &Open ");
  fileOpen->Connect("Clicked()","MyMainFrame",this,"file_open()");
  hframe->AddFrame(fileOpen, new TGLayoutHints(kLHintsCenterX,10,10,3,3));
  
  TGTextButton *drawD = new TGTextButton(hframe," Show &Disambiguation ");
  drawD->Connect("Clicked()","MyMainFrame",this,"DrawDisambig()");
  hframe->AddFrame(drawD, new TGLayoutHints(kLHintsCenterX, 5,2,3,3));

  TGTextButton *drawI = new TGTextButton(hframe," Show &Incorrect ");
  drawI->Connect("Clicked()","MyMainFrame",this,"DrawIncorrect()");
  hframe->AddFrame(drawI, new TGLayoutHints(kLHintsCenterX, 5,2,3,3));

  /*
  TGTextButton *drawAPA0 = new TGTextButton(hframe," APA 0 ");
  drawAPA0->Connect("Clicked()","MyMainFrame",this,"pickAPA0()");
  hframe->AddFrame(drawAPA0, new TGLayoutHints(kLHintsCenterX, 5,2,3,3));

  TGTextButton *drawAPA1 = new TGTextButton(hframe," APA 1 ");
  drawAPA1->Connect("Clicked()","MyMainFrame",this,"pickAPA1()");
  hframe->AddFrame(drawAPA1, new TGLayoutHints(kLHintsCenterX, 5,2,3,3));

  TGTextButton *drawAPA2 = new TGTextButton(hframe," APA 2 ");
  drawAPA2->Connect("Clicked()","MyMainFrame",this,"pickAPA2()");
  hframe->AddFrame(drawAPA2, new TGLayoutHints(kLHintsCenterX, 5,2,3,3));

  TGTextButton *drawAPA3 = new TGTextButton(hframe," APA 3 ");
  drawAPA3->Connect("Clicked()","MyMainFrame",this,"pickAPA3()");
  hframe->AddFrame(drawAPA3, new TGLayoutHints(kLHintsCenterX, 5,2,3,3));
  */

  hframe->AddFrame(APAbgr, new TGLayoutHints(kLHintsCenterX,10,10,0,7));

  TGTextButton *drawZBot = new TGTextButton(hframe,"  &Bot Z  ");
  drawZBot->Connect("Clicked()","MyMainFrame",this,"DrawZBot()");
  hframe->AddFrame(drawZBot, new TGLayoutHints(kLHintsCenterX,	2,2,3,3));

  TGTextButton *drawZTop = new TGTextButton(hframe,"  &Top Z  ");
  drawZTop->Connect("Clicked()","MyMainFrame",this,"DrawZTop()");
  hframe->AddFrame(drawZTop, new TGLayoutHints(kLHintsCenterX,	2,2,3,3));


  TGTextButton *prevEvt = new TGTextButton(hframe," &Prev Evt ");
  prevEvt->Connect("Clicked()","MyMainFrame",this,"PrevEvt()");
  hframe->AddFrame(prevEvt, new TGLayoutHints(kLHintsCenterX, 5,2,3,3));

  fEVT = new TGNumberEntry(hframe, 0, 5, -1, 
			   TGNumberFormat::kNESInteger, 
			   TGNumberFormat::kNEAPositive);
  fEVT->Connect("ValueSet(Int_t)","MyMainFrame",this,"ChooseEvt()");
  (fEVT->GetNumberEntry())->Connect("ReturnPressed()", "MyMainFrame", this,"ChooseEvt()");
  hframe->AddFrame(fEVT, new TGLayoutHints(kLHintsCenterX, 5,2,3,3));
  fEVT->SetIntNumber(1);  

  TGTextButton *nextEvt = new TGTextButton(hframe," &Next Evt ");
  nextEvt->Connect("Clicked()","MyMainFrame",this,"NextEvt()");
  hframe->AddFrame(nextEvt, new TGLayoutHints(kLHintsCenterX, 2,5,3,3));
  
  TGTextButton *exit = new TGTextButton(hframe,"   &Exit   ", "gApplication->Terminate(0)");
  hframe->AddFrame(exit, new TGLayoutHints(kLHintsCenterX,	10,10,3,4));
  
  
  Int_t parts[] = {80, 10, 10}; 
  fStatusBar = new TGStatusBar(fMain,600,40,kHorizontalFrame); 
  fStatusBar->SetParts(parts,3); 
  fMain->AddFrame( fStatusBar,new TGLayoutHints(kLHintsCenterX, 2,2,0,0) );
  fStatusBar->SetText("OPEN A FILE",0); 
  fStatusBar->SetText("",1); 
  fStatusBar->SetText("",2); 
 
  Int_t Lparts[] = {25, 25, 25, 25}; 
  fPlotLabels = new TGStatusBar( fMain, 1100, 40, kHorizontalFrame); 
  fPlotLabels->SetParts(Lparts,4); 
  fMain->AddFrame( fPlotLabels,new TGLayoutHints(kLHintsCenterX, 2,2,0,0) );
  fPlotLabels->SetText("",0);  // label top left plot
  fPlotLabels->SetText("",1);  // label bottom left plot 
  fPlotLabels->SetText("",2);  // label top right plot 
  fPlotLabels->SetText("",3);  // label bottom right plot 
  
  
  fMain->AddFrame( fPlotLabels, new TGLayoutHints(kLHintsCenterX, 2,2,2,2) );
  fMain->AddFrame( fStatusBar,  new TGLayoutHints(kLHintsCenterX, 2,2,2,2) );
  fMain->AddFrame( hframe,      new TGLayoutHints(kLHintsCenterX, 2,2,2,2) );
  APAbgr->Show();

  fMain->SetWindowName("FD Workspace Disambiguation Display");
  fMain->MapSubwindows();
  fMain->Resize(fMain->GetDefaultSize());
  fMain->MapWindow();

  fCanvas = (TCanvas*)gROOT->GetListOfCanvases()->FindObject("Ecanvas");
  if(fCanvas) fCanvas->Clear();
  fCanvas = fEcanvas->GetCanvas();  
  fCanvas->cd();
  fCanvas->Divide(2,2,0,0);
  fEvanvas_1 = (TPad*)(fCanvas->GetPrimitive("Ecanvas_1"));
  fEvanvas_2 = (TPad*)(fCanvas->GetPrimitive("Ecanvas_2"));
  fEvanvas_3 = (TPad*)(fCanvas->GetPrimitive("Ecanvas_3"));
  fEvanvas_4 = (TPad*)(fCanvas->GetPrimitive("Ecanvas_4"));

  this->file_open();

}


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

int MyMainFrame::NextEvt(){
  fEvent++;
  fEVT->SetIntNumber(fEvent);
  std::stringstream newLabel;
  newLabel.str(""); newLabel << "Event "<<fEvent;
  fStatusBar->SetText(newLabel.str().c_str(),1);
  fAPA = 0;
  this->DrawZBot();
  //this->DrawCheated();
  return fEvent;
}
int MyMainFrame::PrevEvt(){
  if(fEvent==1) return 1;
  fEvent--;
  fEVT->SetIntNumber(fEvent);
  std::stringstream newLabel;
  newLabel.str(""); newLabel << "Event "<<fEvent;
  fStatusBar->SetText(newLabel.str().c_str(),1);
  fAPA = 0;
  this->DrawZBot();
  //this->DrawCheated();
  return fEvent;
}
int MyMainFrame::ChooseEvt(){
  fEvent = fEVT->GetNumberEntry()->GetIntNumber();
  std::stringstream newLabel;
  newLabel.str(""); newLabel << "Event "<<fEvent;
  fStatusBar->SetText(newLabel.str().c_str(),1);
  fAPA = 0;
  this->DrawZBot();
  //this->DrawCheated();
  return fEvent;
}


void MyMainFrame::pickAPA0(){
  fAPA = 0;
  //this->DrawCheated();
}
void MyMainFrame::pickAPA1(){
  fAPA = 1;
  //this->DrawCheated();
}
void MyMainFrame::pickAPA2(){
  fAPA = 2;
  //this->DrawCheated();
}
void MyMainFrame::pickAPA3(){
  fAPA = 3;
  //this->DrawCheated();
}



//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
string MyMainFrame::file_open() {
  fEvent = 1;
  fAPA   = 0;
  TGFileInfo file_info_;
  const char *filetypes[] = {"Root Files", "*.root", 0, 0};
  file_info_.fFileTypes = filetypes;
  file_info_.fIniDir = StrDup(".");
  new TGFileDialog(gClient->GetDefaultRoot(), gClient->GetDefaultRoot(),
		   kFDOpen, &file_info_);
  if( file_info_.fFilename ){
    cout << "'" << file_info_.fFilename << "' selected." << endl;

    fname=file_info_.fFilename;
    fStatusBar->SetText(fname.c_str(),0);
    fStatusBar->SetText("Event 1",1);

    fFile = new TFile(fname.c_str());

    //this->DrawCheated();
    this->DrawZBot();
    return file_info_.fFilename; 

  }
  return "";
}

/*
//////////////////////////////////////////////////////////
int MyMainFrame::DrawCheated() {

  this->DrawDisambig();

  return fEvent;

}
*/

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
int MyMainFrame::DrawCheated() {

  if ( !fFile->IsOpen() ) return -1;
  fFile->cd("dbigevd");
   
  //fPlotLabels->SetText("(top left): ....WORKING....",0); 
  //fPlotLabels->SetText("",1); 
  //fPlotLabels->SetText("",2); 
  //fPlotLabels->SetText("",3); 
  //fPlotLabels->DoRedraw();
  fPlotLabels->SetText("Top Left: U1",0); 
  fPlotLabels->SetText("Bottom Left: V1",1); 
  fPlotLabels->SetText("Top Right: U0",2); 
  fPlotLabels->SetText("Bottom Right: V0",3); 

  std::stringstream  hname;
  TH2I *h1;
  hname.str(""); 
  float margin = 0.05;


  fEvanvas_1->cd();
  fEvanvas_1->SetLeftMargin(margin);      fEvanvas_1->SetRightMargin(margin);
  fEvanvas_1->SetTopMargin(margin);       fEvanvas_1->SetBottomMargin(margin);
  hname << "dbigevd/fCheatWireU1_" << fAPA << "_" << fEvent;
  h1=(TH2I*)fFile->Get(hname.str().c_str());
  hname.str("");
  h1->SetTitle("U Side 1");
  this->ZoomDraw(*h1);

  fCanvas->cd();  
  fCanvas->Modified();
  fCanvas->Update();

  fEvanvas_2->cd();
  fEvanvas_2->SetLeftMargin(margin);      fEvanvas_2->SetRightMargin(margin);
  fEvanvas_2->SetTopMargin(margin);       fEvanvas_2->SetBottomMargin(margin);
  hname << "dbigevd/fCheatWireV1_" << fAPA << "_" << fEvent;
  h1=(TH2I*)fFile->Get(hname.str().c_str());
  hname.str("");
  h1->SetTitle("V Side 1");
  this->ZoomDraw(*h1);
  
  fCanvas->cd();  
  fCanvas->Modified();
  fCanvas->Update();

  fEvanvas_3->cd();
  fEvanvas_3->SetLeftMargin(margin);      fEvanvas_3->SetRightMargin(margin);
  fEvanvas_3->SetTopMargin(margin);       fEvanvas_3->SetBottomMargin(margin);
  hname << "dbigevd/fCheatWireU0_" << fAPA << "_" << fEvent;
  h1=(TH2I*)fFile->Get(hname.str().c_str());
  hname.str("");
  h1->SetTitle("U Side 0");
  this->ZoomDraw(*h1);

  fCanvas->cd();  
  fCanvas->Modified();
  fCanvas->Update();

  fEvanvas_4->cd();
  fEvanvas_4->SetLeftMargin(margin);      fEvanvas_4->SetRightMargin(margin);
  fEvanvas_4->SetTopMargin(margin);       fEvanvas_4->SetBottomMargin(margin);
  hname << "dbigevd/fCheatWireV0_" << fAPA << "_" << fEvent;
  h1=(TH2I*)fFile->Get(hname.str().c_str());
  hname.str("");
  h1->SetTitle("V Side 0");
  this->ZoomDraw(*h1);

  fCanvas->cd();  
  fCanvas->Modified();
  fCanvas->Update();

  hname.str(""); hname << "APA "<< fAPA << " Cheated disambiguation; evt " << fEvent;
  fMain->SetWindowName(hname.str().c_str());
  hname.str("");

  delete h1;

  return 0;
}



//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
int MyMainFrame::DrawDisambig() {

  if ( !fFile->IsOpen() ) return -1;
  fFile->cd("dbigevd");
  
  fPlotLabels->SetText("Top Left: U1",0); 
  fPlotLabels->SetText("Bottom Left: V1",1); 
  fPlotLabels->SetText("Top Right: U0",2); 
  fPlotLabels->SetText("Bottom Right: V0",3); 

  std::stringstream  hname;
  TH2I *h1, *h2;
  hname.str(""); 
  float margin = 0.05;

  fEvanvas_1->cd();
  fEvanvas_1->SetLeftMargin(margin);      fEvanvas_1->SetRightMargin(margin);
  fEvanvas_1->SetTopMargin(margin);       fEvanvas_1->SetBottomMargin(margin);
  hname << "dbigevd/fCheatWireU1_" << fAPA << "_" << fEvent;
  h1=(TH2I*)fFile->Get(hname.str().c_str());
  hname.str("");
  hname << "dbigevd/fDisambigWireU1_" << fAPA << "_" << fEvent;
  h2=(TH2I*)fFile->Get(hname.str().c_str());
  hname.str("");
  h1->SetTitle("U Side 1");
  this->ZoomDraw(*h1);
  h2->Draw("same");

  fCanvas->cd();  
  fCanvas->Modified();
  fCanvas->Update();

  fEvanvas_2->cd();
  fEvanvas_2->SetLeftMargin(margin);      fEvanvas_2->SetRightMargin(margin);
  fEvanvas_2->SetTopMargin(margin);       fEvanvas_2->SetBottomMargin(margin);
  hname << "dbigevd/fCheatWireV1_" << fAPA << "_" << fEvent;
  h1=(TH2I*)fFile->Get(hname.str().c_str());
  hname.str("");
  hname << "dbigevd/fDisambigWireV1_" << fAPA << "_" << fEvent;
  h2=(TH2I*)fFile->Get(hname.str().c_str());
  hname.str("");
  h1->SetTitle("V Side 1");
  this->ZoomDraw(*h1);
  h2->Draw("same");
  
  fCanvas->cd();  
  fCanvas->Modified();
  fCanvas->Update();

  fEvanvas_3->cd();
  fEvanvas_3->SetLeftMargin(margin);      fEvanvas_3->SetRightMargin(margin);
  fEvanvas_3->SetTopMargin(margin);       fEvanvas_3->SetBottomMargin(margin);
  hname << "dbigevd/fCheatWireU0_" << fAPA << "_" << fEvent;
  h1=(TH2I*)fFile->Get(hname.str().c_str());
  hname.str("");
  hname << "dbigevd/fDisambigWireU0_" << fAPA << "_" << fEvent;
  h2=(TH2I*)fFile->Get(hname.str().c_str());
  hname.str("");
  h1->SetTitle("U Side 0");
  this->ZoomDraw(*h1);
  h2->Draw("same");

  fCanvas->cd();  
  fCanvas->Modified();
  fCanvas->Update();

  fEvanvas_4->cd();
  fEvanvas_4->SetLeftMargin(margin);      fEvanvas_4->SetRightMargin(margin);
  fEvanvas_4->SetTopMargin(margin);       fEvanvas_4->SetBottomMargin(margin);
  hname << "dbigevd/fCheatWireV0_" << fAPA << "_" << fEvent;
  h1=(TH2I*)fFile->Get(hname.str().c_str());
  hname.str("");
  hname << "dbigevd/fDisambigWireV0_" << fAPA << "_" << fEvent;
  h2=(TH2I*)fFile->Get(hname.str().c_str());
  hname.str("");
  h1->SetTitle("V Side 0");
  this->ZoomDraw(*h1);
  h2->Draw("same");
  
  
  //gPad->SetEditable(false);
  
  fCanvas->cd();  
  fCanvas->Modified();
  fCanvas->Update();
  hname.str(""); hname << "APA "<< fAPA << " Attempted disambiguations; evt " << fEvent;
  fMain->SetWindowName(hname.str().c_str());
  hname.str("");

  delete h1; delete h2;

  return 0;
}


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
int MyMainFrame::DrawIncorrect() {

  if ( !fFile->IsOpen() ) return -1;
  fFile->cd("dbigevd");

  fPlotLabels->SetText("Top Left: U1",0); 
  fPlotLabels->SetText("Bottom Left: V1",1); 
  fPlotLabels->SetText("Top Right: U0",2); 
  fPlotLabels->SetText("Bottom Right: V0",3);

  std::stringstream  hname;
  TH2I *h1, *h2;
  hname.str(""); 
  float margin = 0.05;

  fEvanvas_1->cd();
  fEvanvas_1->SetLeftMargin(margin);      fEvanvas_1->SetRightMargin(margin);
  fEvanvas_1->SetTopMargin(margin);       fEvanvas_1->SetBottomMargin(margin);
  hname << "dbigevd/fCheatWireU1_" << fAPA << "_" << fEvent;
  h1=(TH2I*)fFile->Get(hname.str().c_str());
  hname.str("");
  hname << "dbigevd/fDisambigWireU1_" << fAPA << "_" << fEvent;
  h2=(TH2I*)fFile->Get(hname.str().c_str());
  hname.str("");
  h1->SetTitle("U Side 1");
  this->ZoomDraw(*h1); // needs to be the cheated for correct zooming?
  h2->Draw("same");
  h1->Draw("same");
  // ^^^ order of h1, h2 draws are the only thing different from DrawDisambig()

  
  fEvanvas_2->cd();
  fEvanvas_2->SetLeftMargin(margin);      fEvanvas_2->SetRightMargin(margin);
  fEvanvas_2->SetTopMargin(margin);       fEvanvas_2->SetBottomMargin(margin);
  hname << "dbigevd/fCheatWireV1_" << fAPA << "_" << fEvent;
  h1=(TH2I*)fFile->Get(hname.str().c_str());
  hname.str("");
  hname << "dbigevd/fDisambigWireV1_" << fAPA << "_" << fEvent;
  h2=(TH2I*)fFile->Get(hname.str().c_str());
  hname.str("");
  h1->SetTitle("V Side 1");
  this->ZoomDraw(*h1);
  h2->Draw("same");
  h1->Draw("same");
  

  fEvanvas_3->cd();
  fEvanvas_3->SetLeftMargin(margin);      fEvanvas_3->SetRightMargin(margin);
  fEvanvas_3->SetTopMargin(margin);       fEvanvas_3->SetBottomMargin(margin);
  hname << "dbigevd/fCheatWireU0_" << fAPA << "_" << fEvent;
  h1=(TH2I*)fFile->Get(hname.str().c_str());
  hname.str("");
  hname << "dbigevd/fDisambigWireU0_" << fAPA << "_" << fEvent;
  h2=(TH2I*)fFile->Get(hname.str().c_str());
  hname.str("");
  h1->SetTitle("U Side 0");
  this->ZoomDraw(*h1);
  h2->Draw("same");
  h1->Draw("same");

  
  fEvanvas_4->cd();
  fEvanvas_4->SetLeftMargin(margin);      fEvanvas_4->SetRightMargin(margin);
  fEvanvas_4->SetTopMargin(margin);       fEvanvas_4->SetBottomMargin(margin);
  hname << "dbigevd/fCheatWireV0_" << fAPA << "_" << fEvent;
  h1=(TH2I*)fFile->Get(hname.str().c_str());
  hname.str("");
  hname << "dbigevd/fDisambigWireV0_" << fAPA << "_" << fEvent;
  h2=(TH2I*)fFile->Get(hname.str().c_str());
  hname.str("");
  h1->SetTitle("V Side 0");
  this->ZoomDraw(*h1);
  h2->Draw("same");
  h1->Draw("same");
  
  
  //gPad->SetEditable(false);
  
  
  fCanvas->cd();  
  fCanvas->Modified();
  fCanvas->Update();
  hname.str(""); hname << "APA "<< fAPA << " Incorrect disambiguations; evt " << fEvent;
  fMain->SetWindowName(hname.str().c_str());
  hname.str("");

  delete h1; delete h2;

  return 0;
}


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
int MyMainFrame::DrawZBot() {

  if ( !fFile->IsOpen() ) return -1;
  fFile->cd("dbigevd");

  std::stringstream  hname;
  TH2I *htemp;
  fStatusBar->SetText("View Z",2);
  fPlotLabels->SetText("Top Left: APA 0, Z1 side",0); 
  fPlotLabels->SetText("Bot Left: APA 0, Z0 side",1); 
  fPlotLabels->SetText("Top Right: APA 2, Z1 side",2); 
  fPlotLabels->SetText("Bot Right: APA 2, Z1 side",3); 
  
  
  //APA 0, Z1
  fEvanvas_1->cd();
  fEvanvas_1->SetLeftMargin(0.01);  fEvanvas_1->SetRightMargin(0.01);
  fEvanvas_1->SetTopMargin(0.01);   fEvanvas_1->SetBottomMargin(0.01);
  hname.str(""); hname << "dbigevd/fTimeChanZ1_0_" << fEvent;
  htemp=(TH2I*)fFile->Get(hname.str().c_str()); 
  htemp->SetLabelOffset();
  htemp->SetTitle("");
  htemp->Draw("ahcol");
  
  //APA 2, Z1
  fEvanvas_2->cd();
  fEvanvas_2->SetLeftMargin(0.01);  fEvanvas_2->SetRightMargin(0.01);
  fEvanvas_2->SetTopMargin(0.01);   fEvanvas_2->SetBottomMargin(0.01);
  hname.str(""); hname << "dbigevd/fTimeChanZ1_2_" << fEvent;
  htemp=(TH2I*)fFile->Get(hname.str().c_str()); 
  htemp->SetLabelOffset();
  htemp->SetTitle("");
  htemp->Draw("ahcol");

  //APA 0, Z0
  fEvanvas_3->cd();
  fEvanvas_3->SetLeftMargin(0.01);  fEvanvas_3->SetRightMargin(0.01);
  fEvanvas_3->SetTopMargin(0.01);   fEvanvas_3->SetBottomMargin(0.01);
  hname.str(""); hname << "dbigevd/fTimeChanZ0_0_" << fEvent;
  htemp=(TH2I*)fFile->Get(hname.str().c_str()); 
  htemp->SetLabelOffset();
  htemp->SetTitle("");
  htemp->Draw("ahcol");

  //APA 2, Z0
  fEvanvas_4->cd();
  fEvanvas_4->SetLeftMargin(0.01);  fEvanvas_4->SetRightMargin(0.01);
  fEvanvas_4->SetTopMargin(0.01);   fEvanvas_4->SetBottomMargin(0.01);  
  hname.str(""); hname << "dbigevd/fTimeChanZ0_2_" << fEvent;
  htemp=(TH2I*)fFile->Get(hname.str().c_str()); 
  htemp->SetLabelOffset();
  htemp->SetTitle("");
  htemp->Draw("ahcol");
 

  //gPad->SetEditable(false);
  
  fCanvas->cd();  
  fCanvas->Modified();
  fCanvas->Update();
  hname.str(""); hname << "Z-Plane; Bottom APAs (0,2); evt " << fEvent;
  fMain->SetWindowName(hname.str().c_str());
  hname.str("");

  delete htemp;

  return 0;
}



//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
int MyMainFrame::DrawZTop() {

  if ( !fFile->IsOpen() ) return -1;
  fFile->cd("dbigevd");

  std::stringstream  hname;
  TH2I *htemp;
  fStatusBar->SetText("View Z",2);
  fPlotLabels->SetText("Top Left: APA 1, Z1 side",0); 
  fPlotLabels->SetText("Bot Left: APA 1, Z0 side",1); 
  fPlotLabels->SetText("Top Right: APA 3, Z1 side",2); 
  fPlotLabels->SetText("Bot Right: APA 3, Z1 side",3);   
  
  //APA 1, Z1
  fEvanvas_1->cd();
  fEvanvas_1->SetLeftMargin(0.01);  fEvanvas_1->SetRightMargin(0.01);
  fEvanvas_1->SetTopMargin(0.01);   fEvanvas_1->SetBottomMargin(0.01);
  hname.str(""); hname << "dbigevd/fTimeChanZ1_1_" << fEvent;
  htemp=(TH2I*)fFile->Get(hname.str().c_str()); 
  htemp->SetLabelOffset();
  htemp->SetTitle("");
  htemp->Draw("ahcol");
  
  //APA 3, Z1
  fEvanvas_2->cd();
  fEvanvas_2->SetLeftMargin(0.01);  fEvanvas_2->SetRightMargin(0.01);
  fEvanvas_2->SetTopMargin(0.01);   fEvanvas_2->SetBottomMargin(0.01);
  hname.str(""); hname << "dbigevd/fTimeChanZ1_3_" << fEvent;
  htemp=(TH2I*)fFile->Get(hname.str().c_str()); 
  htemp->SetLabelOffset();
  htemp->SetTitle("");
  htemp->Draw("ahcol");

  //APA 1, Z0
  fEvanvas_3->cd();
  fEvanvas_3->SetLeftMargin(0.01);  fEvanvas_3->SetRightMargin(0.01);
  fEvanvas_3->SetTopMargin(0.01);   fEvanvas_3->SetBottomMargin(0.01);
  hname.str(""); hname << "dbigevd/fTimeChanZ0_1_" << fEvent;
  htemp=(TH2I*)fFile->Get(hname.str().c_str()); 
  htemp->SetLabelOffset();
  htemp->SetTitle("");
  htemp->Draw("ahcol");

  //APA 3, Z0
  fEvanvas_4->cd();
  fEvanvas_4->SetLeftMargin(0.01);  fEvanvas_4->SetRightMargin(0.01);
  fEvanvas_4->SetTopMargin(0.01);   fEvanvas_4->SetBottomMargin(0.01);  
  hname.str(""); hname << "dbigevd/fTimeChanZ0_3_" << fEvent;
  htemp=(TH2I*)fFile->Get(hname.str().c_str()); 
  htemp->SetLabelOffset();
  htemp->SetTitle("");
  htemp->Draw("ahcol");


  //gPad->SetEditable(false);
  
  fCanvas->cd();  
  fCanvas->Modified();
  fCanvas->Update(); 
  hname.str(""); hname << "Z-Plane; Top APAs (1,3); evt " << fEvent;
  fMain->SetWindowName(hname.str().c_str());
  hname.str("");

  delete htemp;

  return 0;
}



//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
MyMainFrame::~MyMainFrame() {
// Clean up used widgets: frames, buttons, layout hints
fMain->Cleanup();
delete fMain;
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
void DbigDisplay_FD4apa() {
// Popup the GUI...
new MyMainFrame(gClient->GetRoot(),200,200);
}



//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
void MyMainFrame::ZoomDraw(TH2I& h){

  int firstXbin =h.GetXaxis()->GetFirst();
  int lastXbin  =h.GetXaxis()->GetLast();
  int firstYbin =h.GetYaxis()->GetFirst();
  int lastYbin  =h.GetYaxis()->GetLast();
  for(int Xbin = firstXbin; Xbin<=lastXbin; Xbin++){
    bool found=false;
    for(int Ybin = firstYbin; Ybin<=lastYbin; Ybin++){
      if(h.GetBinContent(Xbin, Ybin)<=0) continue;
      firstXbin = Xbin;
      found=true;      break;
    }
    if(found) break;
  }
  for(int Ybin = firstYbin; Ybin<=lastYbin; Ybin++){
    bool found=false;
    for(int Xbin = firstXbin; Xbin<=lastXbin; Xbin++){
      if(h.GetBinContent(Xbin, Ybin)<=0) continue;
      firstYbin = Ybin;
      found=true;      break;
    }
    if(found) break;
  }
  for(int Xbin = lastXbin; Xbin>=firstXbin; Xbin--){
    bool found=false;
    for(int Ybin = lastYbin; Ybin>=firstYbin; Ybin--){
      if(h.GetBinContent(Xbin, Ybin)<=0) continue;
      lastXbin = Xbin;
      found=true;      break;
    }
    if(found) break;
  }
  for(int Ybin = lastYbin; Ybin>=firstYbin; Ybin--){
    bool found=false;
    for(int Xbin = lastXbin; Xbin>=firstXbin; Xbin--){
      if(h.GetBinContent(Xbin, Ybin)<=0) continue;
      lastYbin = Ybin;
      found=true;      break;
    }
    if(found) break;
  }

  //double xmin = h.GetXaxis()->GetBinLowEdge(firstXbin);
  //double xmax = h.GetXaxis()->GetBinUpEdge(lastXbin);
  //double ymin = h.GetYaxis()->GetBinLowEdge(firstYbin);
  //double ymax = h.GetYaxis()->GetBinUpEdge(lastYbin);
  //h.GetXaxis()->SetLimits(xmin,xmax);
  //h.GetYaxis()->SetLimits(ymin,ymax);
  h.GetXaxis()->SetRange(firstXbin,lastXbin);
  h.GetYaxis()->SetRange(firstYbin-10,lastYbin+10);
  h.Draw();

}






/*

void MyMainFrame::MakeLogRainbowPalette()
{

  this->MakeRainbow();

  for(unsigned int x=1; x<=256; x++){
    float f = (log(x)-log(1))/(log(256)-log(1));
    int indx = (int)(f*(float)fNcolors);
    if (indx<0)        indx = 0;
    if (indx>=fNcolors) indx = fNcolors-1;
    fLogColors[x-1] = fColors[indx];
  }

  gStyle->SetPalette(fNcolors,fLogColors);

  return;
}




void MyMainFrame::MakeRainbow(){

    fColors[  0] = TColor::GetColor(  45,   0,  36 );
    fColors[  1] = TColor::GetColor(  45,   0,  36 );
    fColors[  2] = TColor::GetColor(  56,   0,  46 );
    fColors[  3] = TColor::GetColor(  60,   0,  49 );
    fColors[  4] = TColor::GetColor(  67,   0,  54 );
    fColors[  5] = TColor::GetColor(  70,   0,  59 );
    fColors[  6] = TColor::GetColor(  71,   0,  61 );
    fColors[  7] = TColor::GetColor(  75,   0,  68 );
    fColors[  8] = TColor::GetColor(  74,   0,  73 );
    fColors[  9] = TColor::GetColor(  74,   0,  77 );
    fColors[ 10] = TColor::GetColor(  73,   0,  81 );
    fColors[ 11] = TColor::GetColor(  71,   0,  87 );
    fColors[ 12] = TColor::GetColor(  69,   1,  90 );
    fColors[ 13] = TColor::GetColor(  68,   2,  94 );
    fColors[ 14] = TColor::GetColor(  66,   3,  97 );
    fColors[ 15] = TColor::GetColor(  63,   6, 102 );
    fColors[ 16] = TColor::GetColor(  61,   7, 106 );
    fColors[ 17] = TColor::GetColor(  58,  10, 109 );
    fColors[ 18] = TColor::GetColor(  56,  12, 113 );
    fColors[ 19] = TColor::GetColor(  53,  15, 116 );
    fColors[ 20] = TColor::GetColor(  48,  18, 119 );
    fColors[ 21] = TColor::GetColor(  47,  20, 121 );
    fColors[ 22] = TColor::GetColor(  44,  23, 124 );
    fColors[ 23] = TColor::GetColor(  41,  27, 128 );
    fColors[ 24] = TColor::GetColor(  40,  28, 129 );
    fColors[ 25] = TColor::GetColor(  37,  32, 132 );
    fColors[ 26] = TColor::GetColor(  34,  36, 134 );
    fColors[ 27] = TColor::GetColor(  29,  43, 137 );
    fColors[ 28] = TColor::GetColor(  25,  52, 138 );
    fColors[ 29] = TColor::GetColor(  24,  57, 139 );
    fColors[ 30] = TColor::GetColor(  24,  62, 141 );
    fColors[ 31] = TColor::GetColor(  24,  64, 142 );
    fColors[ 32] = TColor::GetColor(  23,  65, 142 );
    fColors[ 33] = TColor::GetColor(  23,  69, 143 );
    fColors[ 34] = TColor::GetColor(  23,  71, 142 );
    fColors[ 35] = TColor::GetColor(  23,  71, 142 );
    fColors[ 36] = TColor::GetColor(  23,  73, 142 );
    fColors[ 37] = TColor::GetColor(  23,  75, 142 );
    fColors[ 38] = TColor::GetColor(  23,  75, 142 );
    fColors[ 39] = TColor::GetColor(  23,  78, 142 );
    fColors[ 40] = TColor::GetColor(  23,  80, 142 );
    fColors[ 41] = TColor::GetColor(  23,  80, 142 );
    fColors[ 42] = TColor::GetColor(  23,  82, 141 );
    fColors[ 43] = TColor::GetColor(  23,  85, 141 );
    fColors[ 44] = TColor::GetColor(  23,  85, 141 );
    fColors[ 45] = TColor::GetColor(  23,  87, 140 );
    fColors[ 46] = TColor::GetColor(  23,  87, 140 );
    fColors[ 47] = TColor::GetColor(  24,  90, 140 );
    fColors[ 48] = TColor::GetColor(  24,  90, 140 );
    fColors[ 49] = TColor::GetColor(  24,  93, 139 );
    fColors[ 50] = TColor::GetColor(  24,  93, 139 );
    fColors[ 51] = TColor::GetColor(  24,  93, 139 );
    fColors[ 52] = TColor::GetColor(  24,  93, 139 );
    fColors[ 53] = TColor::GetColor(  24,  97, 139 );
    fColors[ 54] = TColor::GetColor(  24,  97, 139 );
    fColors[ 55] = TColor::GetColor(  25, 101, 138 );
    fColors[ 56] = TColor::GetColor(  25, 101, 138 );
    fColors[ 57] = TColor::GetColor(  25, 104, 137 );
    fColors[ 58] = TColor::GetColor(  25, 104, 137 );
    fColors[ 59] = TColor::GetColor(  25, 104, 137 );
    fColors[ 60] = TColor::GetColor(  26, 108, 137 );
    fColors[ 61] = TColor::GetColor(  26, 108, 137 );
    fColors[ 62] = TColor::GetColor(  27, 111, 136 );
    fColors[ 63] = TColor::GetColor(  27, 111, 136 );
    fColors[ 64] = TColor::GetColor(  27, 111, 136 );
    fColors[ 65] = TColor::GetColor(  27, 115, 135 );
    fColors[ 66] = TColor::GetColor(  27, 115, 135 );
    fColors[ 67] = TColor::GetColor(  28, 118, 134 );
    fColors[ 68] = TColor::GetColor(  28, 118, 134 );
    fColors[ 69] = TColor::GetColor(  29, 122, 133 );
    fColors[ 70] = TColor::GetColor(  29, 122, 133 );
    fColors[ 71] = TColor::GetColor(  29, 122, 133 );
    fColors[ 72] = TColor::GetColor(  29, 122, 133 );
    fColors[ 73] = TColor::GetColor(  29, 125, 132 );
    fColors[ 74] = TColor::GetColor(  29, 125, 132 );
    fColors[ 75] = TColor::GetColor(  30, 128, 131 );
    fColors[ 76] = TColor::GetColor(  30, 128, 131 );
    fColors[ 77] = TColor::GetColor(  31, 131, 130 );
    fColors[ 78] = TColor::GetColor(  31, 131, 130 );
    fColors[ 79] = TColor::GetColor(  31, 131, 130 );
    fColors[ 80] = TColor::GetColor(  32, 134, 128 );
    fColors[ 81] = TColor::GetColor(  32, 134, 128 );
    fColors[ 82] = TColor::GetColor(  33, 137, 127 );
    fColors[ 83] = TColor::GetColor(  33, 137, 127 );
    fColors[ 84] = TColor::GetColor(  33, 137, 127 );
    fColors[ 85] = TColor::GetColor(  34, 140, 125 );
    fColors[ 86] = TColor::GetColor(  34, 140, 125 );
    fColors[ 87] = TColor::GetColor(  35, 142, 123 );
    fColors[ 88] = TColor::GetColor(  35, 142, 123 );
    fColors[ 89] = TColor::GetColor(  36, 145, 121 );
    fColors[ 90] = TColor::GetColor(  36, 145, 121 );
    fColors[ 91] = TColor::GetColor(  36, 145, 121 );
    fColors[ 92] = TColor::GetColor(  37, 147, 118 );
    fColors[ 93] = TColor::GetColor(  37, 147, 118 );
    fColors[ 94] = TColor::GetColor(  38, 150, 116 );
    fColors[ 95] = TColor::GetColor(  38, 150, 116 );
    fColors[ 96] = TColor::GetColor(  40, 152, 113 );
    fColors[ 97] = TColor::GetColor(  40, 152, 113 );
    fColors[ 98] = TColor::GetColor(  41, 154, 111 );
    fColors[ 99] = TColor::GetColor(  41, 154, 111 );
    fColors[100] = TColor::GetColor(  42, 156, 108 );
    fColors[101] = TColor::GetColor(  42, 156, 108 );
    fColors[102] = TColor::GetColor(  43, 158, 106 );
    fColors[103] = TColor::GetColor(  43, 158, 106 );
    fColors[104] = TColor::GetColor(  43, 158, 106 );
    fColors[105] = TColor::GetColor(  45, 160, 104 );
    fColors[106] = TColor::GetColor(  45, 160, 104 );
    fColors[107] = TColor::GetColor(  46, 162, 101 );
    fColors[108] = TColor::GetColor(  46, 162, 101 );
    fColors[109] = TColor::GetColor(  48, 164,  99 );
    fColors[110] = TColor::GetColor(  48, 164,  99 );
    fColors[111] = TColor::GetColor(  50, 166,  97 );
    fColors[112] = TColor::GetColor(  50, 166,  97 );
    fColors[113] = TColor::GetColor(  51, 168,  95 );
    fColors[114] = TColor::GetColor(  53, 170,  93 );
    fColors[115] = TColor::GetColor(  53, 170,  93 );
    fColors[116] = TColor::GetColor(  53, 170,  93 );
    fColors[117] = TColor::GetColor(  55, 172,  91 );
    fColors[118] = TColor::GetColor(  55, 172,  91 );
    fColors[119] = TColor::GetColor(  57, 174,  88 );
    fColors[120] = TColor::GetColor(  57, 174,  88 );
    fColors[121] = TColor::GetColor(  59, 175,  86 );
    fColors[122] = TColor::GetColor(  62, 177,  84 );
    fColors[123] = TColor::GetColor(  64, 178,  82 );
    fColors[124] = TColor::GetColor(  64, 178,  82 );
    fColors[125] = TColor::GetColor(  67, 180,  80 );
    fColors[126] = TColor::GetColor(  67, 180,  80 );
    fColors[127] = TColor::GetColor(  69, 181,  79 );
    fColors[128] = TColor::GetColor(  72, 183,  77 );
    fColors[129] = TColor::GetColor(  72, 183,  77 );
    fColors[130] = TColor::GetColor(  72, 183,  77 );
    fColors[131] = TColor::GetColor(  75, 184,  76 );
    fColors[132] = TColor::GetColor(  77, 186,  74 );
    fColors[133] = TColor::GetColor(  80, 187,  73 );
    fColors[134] = TColor::GetColor(  83, 189,  72 );
    fColors[135] = TColor::GetColor(  87, 190,  72 );
    fColors[136] = TColor::GetColor(  91, 191,  71 );
    fColors[137] = TColor::GetColor(  95, 192,  70 );
    fColors[138] = TColor::GetColor(  99, 193,  70 );
    fColors[139] = TColor::GetColor( 103, 194,  70 );
    fColors[140] = TColor::GetColor( 107, 195,  70 );
    fColors[141] = TColor::GetColor( 111, 196,  70 );
    fColors[142] = TColor::GetColor( 111, 196,  70 );
    fColors[143] = TColor::GetColor( 115, 196,  70 );
    fColors[144] = TColor::GetColor( 119, 197,  70 );
    fColors[145] = TColor::GetColor( 123, 197,  70 );
    fColors[146] = TColor::GetColor( 130, 198,  71 );
    fColors[147] = TColor::GetColor( 133, 199,  71 );
    fColors[148] = TColor::GetColor( 137, 199,  72 );
    fColors[149] = TColor::GetColor( 140, 199,  72 );
    fColors[150] = TColor::GetColor( 143, 199,  73 );
    fColors[151] = TColor::GetColor( 143, 199,  73 );
    fColors[152] = TColor::GetColor( 147, 199,  73 );
    fColors[153] = TColor::GetColor( 150, 199,  74 );
    fColors[154] = TColor::GetColor( 153, 199,  74 );
    fColors[155] = TColor::GetColor( 156, 199,  75 );
    fColors[156] = TColor::GetColor( 160, 200,  76 );
    fColors[157] = TColor::GetColor( 167, 200,  78 );
    fColors[158] = TColor::GetColor( 170, 200,  79 );
    fColors[159] = TColor::GetColor( 173, 200,  79 );
    fColors[160] = TColor::GetColor( 173, 200,  79 );
    fColors[161] = TColor::GetColor( 177, 200,  80 );
    fColors[162] = TColor::GetColor( 180, 200,  81 );
    fColors[163] = TColor::GetColor( 183, 199,  82 );
    fColors[164] = TColor::GetColor( 186, 199,  82 );
    fColors[165] = TColor::GetColor( 190, 199,  83 );
    fColors[166] = TColor::GetColor( 196, 199,  85 );
    fColors[167] = TColor::GetColor( 199, 198,  85 );
    fColors[168] = TColor::GetColor( 199, 198,  85 );
    fColors[169] = TColor::GetColor( 203, 198,  86 );
    fColors[170] = TColor::GetColor( 206, 197,  87 );
    fColors[171] = TColor::GetColor( 212, 197,  89 );
    fColors[172] = TColor::GetColor( 215, 196,  90 );
    fColors[173] = TColor::GetColor( 218, 195,  91 );
    fColors[174] = TColor::GetColor( 224, 194,  94 );
    fColors[175] = TColor::GetColor( 224, 194,  94 );
    fColors[176] = TColor::GetColor( 230, 193,  96 );
    fColors[177] = TColor::GetColor( 233, 192,  98 );
    fColors[178] = TColor::GetColor( 236, 190, 100 );
    fColors[179] = TColor::GetColor( 238, 189, 104 );
    fColors[180] = TColor::GetColor( 240, 188, 106 );
    fColors[181] = TColor::GetColor( 240, 188, 106 );
    fColors[182] = TColor::GetColor( 242, 187, 110 );
    fColors[183] = TColor::GetColor( 244, 185, 114 );
    fColors[184] = TColor::GetColor( 245, 184, 116 );
    fColors[185] = TColor::GetColor( 247, 183, 120 );
    fColors[186] = TColor::GetColor( 248, 182, 123 );
    fColors[187] = TColor::GetColor( 248, 182, 123 );
    fColors[188] = TColor::GetColor( 250, 181, 125 );
    fColors[189] = TColor::GetColor( 251, 180, 128 );
    fColors[190] = TColor::GetColor( 252, 180, 130 );
    fColors[191] = TColor::GetColor( 253, 180, 133 );
    fColors[192] = TColor::GetColor( 253, 180, 133 );
    fColors[193] = TColor::GetColor( 254, 180, 134 );
    fColors[194] = TColor::GetColor( 254, 179, 138 );
    fColors[195] = TColor::GetColor( 255, 179, 142 );
    fColors[196] = TColor::GetColor( 255, 179, 145 );
    fColors[197] = TColor::GetColor( 255, 179, 145 );
    fColors[198] = TColor::GetColor( 255, 179, 152 );
    fColors[199] = TColor::GetColor( 255, 180, 161 );
    fColors[200] = TColor::GetColor( 255, 180, 164 );
    fColors[201] = TColor::GetColor( 255, 180, 167 );
    fColors[202] = TColor::GetColor( 255, 180, 167 );
    fColors[203] = TColor::GetColor( 255, 181, 169 );
    fColors[204] = TColor::GetColor( 255, 181, 170 );
    fColors[205] = TColor::GetColor( 255, 182, 173 );
    fColors[206] = TColor::GetColor( 255, 183, 176 );
    fColors[207] = TColor::GetColor( 255, 183, 176 );
    fColors[208] = TColor::GetColor( 255, 184, 179 );
    fColors[209] = TColor::GetColor( 255, 185, 179 );
    fColors[210] = TColor::GetColor( 255, 185, 182 );
    fColors[211] = TColor::GetColor( 255, 186, 182 );
    fColors[212] = TColor::GetColor( 255, 186, 182 );
    fColors[213] = TColor::GetColor( 255, 187, 185 );
    fColors[214] = TColor::GetColor( 255, 188, 185 );
    fColors[215] = TColor::GetColor( 255, 189, 188 );
    fColors[216] = TColor::GetColor( 255, 189, 188 );
    fColors[217] = TColor::GetColor( 255, 190, 188 );
    fColors[218] = TColor::GetColor( 255, 191, 191 );
    fColors[219] = TColor::GetColor( 255, 192, 191 );
    fColors[220] = TColor::GetColor( 255, 194, 194 );
    fColors[221] = TColor::GetColor( 255, 194, 194 );
    fColors[222] = TColor::GetColor( 255, 197, 197 );
    fColors[223] = TColor::GetColor( 255, 198, 198 );
    fColors[224] = TColor::GetColor( 255, 200, 200 );
    fColors[225] = TColor::GetColor( 255, 201, 201 );
    fColors[226] = TColor::GetColor( 255, 201, 201 );
    fColors[227] = TColor::GetColor( 255, 202, 202 );
    fColors[228] = TColor::GetColor( 255, 203, 203 );
    fColors[229] = TColor::GetColor( 255, 205, 205 );
    fColors[230] = TColor::GetColor( 255, 206, 206 );
    fColors[231] = TColor::GetColor( 255, 206, 206 );
    fColors[232] = TColor::GetColor( 255, 208, 208 );
    fColors[233] = TColor::GetColor( 255, 209, 209 );
    fColors[234] = TColor::GetColor( 255, 211, 211 );
    fColors[235] = TColor::GetColor( 255, 215, 215 );
    fColors[236] = TColor::GetColor( 255, 216, 216 );
    fColors[237] = TColor::GetColor( 255, 216, 216 );
    fColors[238] = TColor::GetColor( 255, 218, 218 );
    fColors[239] = TColor::GetColor( 255, 219, 219 );
    fColors[240] = TColor::GetColor( 255, 221, 221 );
    fColors[241] = TColor::GetColor( 255, 223, 223 );
    fColors[242] = TColor::GetColor( 255, 226, 226 );
    fColors[243] = TColor::GetColor( 255, 228, 228 );
    fColors[244] = TColor::GetColor( 255, 230, 230 );
    fColors[245] = TColor::GetColor( 255, 230, 230 );
    fColors[246] = TColor::GetColor( 255, 232, 232 );
    fColors[247] = TColor::GetColor( 255, 235, 235 );
    fColors[248] = TColor::GetColor( 255, 237, 237 );
    fColors[249] = TColor::GetColor( 255, 240, 240 );
    fColors[250] = TColor::GetColor( 255, 243, 243 );
    fColors[251] = TColor::GetColor( 255, 246, 246 );
    fColors[252] = TColor::GetColor( 255, 249, 249 );
    fColors[253] = TColor::GetColor( 255, 251, 251 );
    fColors[254] = TColor::GetColor( 255, 253, 253 );
    fColors[255] = TColor::GetColor( 255, 255, 255 );

}




*/

