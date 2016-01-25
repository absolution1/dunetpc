#ifndef _NLPlotMkr
#define _NLPlotMkr

#include <fstream>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TPaveText.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TDatime.h>
#include <TMultiGraph.h>
#include <TGaxis.h>
#include <TPaveStats.h>

// REMOVE includes below to be removed when I add partition to the OnMon header
#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TSQLStatement.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cstring>

// NOTE:  Go through ALL of the above and remove what isn't necessary!!!


using namespace std;

//
// Define global variables to be used in NearlinePlotMaker.C
//

ifstream inFile;

TDatime *TNowGMT; // current GMT time
TDatime *TNow;    // current local time
TDatime *SRtime;
Int_t GMToffset;

unsigned int run;
unsigned int subrun;

UInt_t LastRun;
UInt_t LastSR;



struct NearlinePlot{
  
  TH1F* fHistogram;
  TGraph* fGraphTime;
  std::string fHistName;
  int fNumPoints;
  std::vector<float> fMetricVec;
  std::vector<float> fTimeVec;
  int fPlotCount;
  bool fRMSPlot;

  NearlinePlot(int Npoint, std::string this_hist_name, std::string hist_title, int num_bins, int min_x, int max_x, bool RMS=false);
  bool AddHistogram(TFile const & file, TTree* header, int Xstrtime, int Xsrtime, int XNow, int GMToffset);
  TCanvas* makeHistoCanvas(std::string can_name, std::string can_title, int width, int height, TPaveText* updateText);
  TCanvas* makeGraphTimeCanvas(std::string can_name, std::string can_title, int width, int height);


};


NearlinePlot::NearlinePlot(int Npoint, std::string this_hist_name, std::string hist_title, int num_bins, int min_x, int max_x, bool RMS)
{
  fNumPoints = Npoint;
  fHistName = this_hist_name;
  fMetricVec = std::vector<float>(Npoint, 0);
  fTimeVec = std::vector<float>(Npoint, 0);
  fPlotCount = 0;
  fRMSPlot=RMS;
  fHistogram = new TH1F(fHistName.c_str(), hist_title.c_str(), num_bins, min_x, max_x);
  fGraphTime = 0;
}

bool NearlinePlot::AddHistogram(TFile const & file, TTree* header, int Xstrtime, int Xsrtime, int XNow, int GMToffset){
  TH1F *hist_temp = (TH1F*)file.FindObjectAny(fHistName.c_str());    
  if(hist_temp != 0){
    fHistogram->Add(hist_temp,1.0);
  }    
  else return false;

  if(hist_temp != 0 && header != 0 && Xsrtime != XNow - GMToffset) {
    fTimeVec.at(fPlotCount) = Xsrtime;
    if(fRMSPlot) fMetricVec.at(fPlotCount) = hist_temp->GetRMS(1);
    else fMetricVec.at(fPlotCount) = hist_temp->GetMean(1);
    fPlotCount++;
    return true;
  }

  return false;

}

TCanvas* NearlinePlot::makeHistoCanvas(std::string can_name, std::string can_title, int width, int height, TPaveText* updateText){
  TCanvas *can = new TCanvas(can_name.c_str(), can_title.c_str(), width, height);
  can->cd();
  can->SetLogy();
  gStyle->SetOptStat(111111);
  fHistogram->SetLineWidth(2);
  fHistogram->SetLineColor(kRed);
  updateText->Draw();

  return can;
}

TCanvas* NearlinePlot::makeGraphTimeCanvas(std::string can_name, std::string can_title, int width, int height){
  TCanvas* can = new TCanvas(can_name.c_str(), can_title.c_str(), width, height);
  can->cd();
  gPad->SetGridx();
  fGraphTime = new TGraph(fPlotCount);
  for(int i=0;i<fPlotCount;i++) fGraphTime->SetPoint(i, fTimeVec.at(i), fMetricVec.at(i));
  return can;
  
}


#endif
