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

void histogramZoom(TH1* hist, double n_sigma);
void graphZoom(TGraph* gr, double n_sigma);




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
  TGraph* fGraphMetricTime;
  TGraph* fGraphMetricRmsTime;
  std::string fHistName;
  std::string fHistOutputName;
  std::string fGraphOutputName;
  std::string fGraphRmsOutputName;
  std::string fGraphOutputCanName;
  std::string fGraphRmsOutputCanName;
  int fNumPoints;
  std::vector<float> fMetricVec;
  std::vector<float> fMetricRmsVec;
  std::vector<float> fTimeVec;
  int fPlotCount;

  bool fMakeHistogram;

  NearlinePlot(int Npoint, std::string this_hist_name, std::string this_hist_output_name, std::string this_graph_output_name, std::string this_graph_rms_output_name, std::string this_can_metric_name, std::string this_can_metric_rms_name, std::string hist_title, int num_bins, int min_x, int max_x);
  bool AddHistogram(TFile const & file, TTree* header, int Xsrtime, int XNow, int GMToffset);
  TCanvas* makeHistoCanvas(std::string can_name, std::string can_title, int width, int height, TPaveText* updateText);
  TCanvas* makeGraphMetricTimeCanvas(int width, int height, TPaveText* updateText, std::string taxis_labels, int time_ago, int XNow, bool zoom=false);

  TCanvas* makeGraphMetricRmsTimeCanvas(int width, int height, TPaveText* updateText, std::string taxis_labels, int time_ago, int XNow, bool zoom=false);


};


NearlinePlot::NearlinePlot(int Npoint, std::string this_hist_name, std::string this_hist_output_name, std::string this_graph_output_name, std::string this_graph_rms_output_name, std::string this_can_metric_name, std::string this_can_metric_rms_name, std::string hist_title, int num_bins, int min_x, int max_x)
{
  fNumPoints = Npoint;
  fHistName = this_hist_name;
  fHistOutputName = this_hist_output_name;
  fGraphOutputName = this_graph_output_name;
  fGraphRmsOutputName = this_graph_rms_output_name;
  fGraphOutputCanName = this_can_metric_name;
  fGraphRmsOutputCanName = this_can_metric_rms_name;
  fMetricVec = std::vector<float>(Npoint, 0);
  fMetricRmsVec = std::vector<float>(Npoint, 0);
  fTimeVec = std::vector<float>(Npoint, 0);
  fPlotCount = 0;
   fHistogram = new TH1F(fHistName.c_str(), hist_title.c_str(), num_bins, min_x, max_x);
  fGraphMetricTime = 0;
}

bool NearlinePlot::AddHistogram(TFile const & file, TTree* header, int Xsrtime, int XNow, int GMToffset){
  TH1F *hist_temp = (TH1F*)file.FindObjectAny(fHistName.c_str());    
  if(hist_temp != 0){
    fHistogram->Add(hist_temp,1.0);
  }    
  else return false;

  if(hist_temp != 0 && header != 0 && Xsrtime != XNow - GMToffset) {
    fTimeVec.at(fPlotCount) = Xsrtime;
    fMetricVec.at(fPlotCount) = hist_temp->GetMean(1);
    fMetricRmsVec.at(fPlotCount) = hist_temp->GetRMS(1);
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
  fHistogram->Draw();
  updateText->Draw();

  return can;
}

TCanvas* NearlinePlot::makeGraphMetricTimeCanvas(int width, int height, TPaveText* updateText, std::string taxis_labels, int time_ago, int XNow, bool zoom){
  TCanvas* can = new TCanvas("can_metric_time", fGraphOutputCanName.c_str(), width, height);
  can->cd();
  gPad->SetGridx();
  fGraphMetricTime = new TGraph(fPlotCount);
  for(int i=0;i<fPlotCount;i++) fGraphMetricTime->SetPoint(i, fTimeVec.at(i), fMetricVec.at(i));
  
  fGraphMetricTime->SetTitle(fGraphOutputCanName.c_str());
  fGraphMetricTime->SetMarkerColor(kBlue);
  fGraphMetricTime->GetXaxis()->SetTimeDisplay(1);
  fGraphMetricTime->GetXaxis()->SetLabelSize(0.03);
  fGraphMetricTime->GetXaxis()->SetTimeFormat(taxis_labels.c_str());
  fGraphMetricTime->GetXaxis()->SetLimits(time_ago,XNow);
  fGraphMetricTime->GetXaxis()->SetTitle("(central time)");
  fGraphMetricTime->Draw("A*");
  if(zoom){
    graphZoom(fGraphMetricTime, 2.0);
    fGraphMetricTime->SetTitle(std::string(fGraphOutputCanName+" - Zoom").c_str());
  }
  updateText->Draw();

  int maxtime = 0;
  double max  = 0.0, ave = 0.0;
  TPaveText *LastPoint = new TPaveText(0.3,0.88,0.93,0.93,"NDC");
  LastPoint->SetLineColor(1);
  LastPoint->SetFillColor(0);
  LastPoint->SetBorderSize(1);
  char lptext[128];
  maxtime = 0;
  max = 0.0;
  ave = 0.0;
  for(int i = 0; i < fPlotCount; ++i) {
    ave += (double)fMetricVec.at(i);
    if(fTimeVec.at(i) > maxtime) {
      maxtime = fTimeVec.at(i);
      max     = fMetricVec.at(i);
    }
  }
  if(fPlotCount > 0) ave = ave/(double)fPlotCount;
  sprintf(lptext,"Last Point = %f  /  Average = %f",max,ave);
  LastPoint->Clear();
  LastPoint->AddText(lptext);


  LastPoint->Draw();

  return can;
  
}

TCanvas* NearlinePlot::makeGraphMetricRmsTimeCanvas(int width, int height, TPaveText* updateText, std::string taxis_labels, int time_ago, int XNow, bool zoom){
  TCanvas* can = new TCanvas("can_metric_rms_time", fGraphRmsOutputCanName.c_str(), width, height);
  can->cd();
  gPad->SetGridx();
  fGraphMetricRmsTime = new TGraph(fPlotCount);
  for(int i=0;i<fPlotCount;i++) fGraphMetricRmsTime->SetPoint(i, fTimeVec.at(i), fMetricRmsVec.at(i));
  
  fGraphMetricRmsTime->SetTitle(fGraphRmsOutputCanName.c_str());
  fGraphMetricRmsTime->SetMarkerColor(kBlue);
  fGraphMetricRmsTime->GetXaxis()->SetTimeDisplay(1);
  fGraphMetricRmsTime->GetXaxis()->SetLabelSize(0.03);
  fGraphMetricRmsTime->GetXaxis()->SetTimeFormat(taxis_labels.c_str());
  fGraphMetricRmsTime->GetXaxis()->SetLimits(time_ago,XNow);
  fGraphMetricRmsTime->GetXaxis()->SetTitle("(central time)");
  fGraphMetricRmsTime->Draw("A*");

  if(zoom){
    graphZoom(fGraphMetricRmsTime, 2.0);
    fGraphMetricRmsTime->SetTitle(std::string(fGraphRmsOutputCanName+" - Zoom").c_str());
  }

  updateText->Draw();

  int maxtime = 0;
  double max  = 0.0, ave = 0.0;
  TPaveText *LastPoint = new TPaveText(0.3,0.88,0.93,0.93,"NDC");
  LastPoint->SetLineColor(1);
  LastPoint->SetFillColor(0);
  LastPoint->SetBorderSize(1);
  char lptext[128];
  maxtime = 0;
  max = 0.0;
  ave = 0.0;
  for(int i = 0; i < fPlotCount; ++i) {
    ave += (double)fMetricVec.at(i);
    if(fTimeVec.at(i) > maxtime) {
      maxtime = fTimeVec.at(i);
      max     = fMetricVec.at(i);
    }
  }
  if(fPlotCount > 0) ave = ave/(double)fPlotCount;
  sprintf(lptext,"Last Point = %f  /  Average = %f",max,ave);
  LastPoint->Clear();
  LastPoint->AddText(lptext);


  LastPoint->Draw();

  return can;
  
}


#endif
