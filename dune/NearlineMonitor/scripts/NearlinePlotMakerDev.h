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


struct NearlinePlotInfo{

  std::string fMetricName;
  int fChannel;
  int fNumDays;
  std::string fFileExtension;
  bool fGotData;

  NearlinePlotInfo(){
    fGotData=false;
  }


  NearlinePlotInfo(std::string metric_name, int channel, int num_days, std::string file_extension){
    fMetricName = metric_name;
    fChannel = channel;
    fNumDays = num_days;
    fFileExtension = file_extension;
    fGotData = true;
  }
  
  std::string GetHistOutputName(){
    char name[256];
    sprintf(name, "%sChan%04i_%.3i_days.%s", fMetricName.c_str(), fChannel, fNumDays, fFileExtension.c_str());
    return std::string(name);
  }
  
  std::string GetMetricMeanTimeGraphName(bool zoom=false){
    char name[256];
    if(zoom) sprintf(name, "%sMeanTimeChan%04i_%.3i_days_zoom.%s", fMetricName.c_str(), fChannel, fNumDays, fFileExtension.c_str());
    else sprintf(name, "%sMeanTimeChan%04i_%.3i_days.%s", fMetricName.c_str(), fChannel, fNumDays, fFileExtension.c_str());
    return std::string(name);
  }

  std::string GetMetricRmsTimeGraphName(bool zoom=false){
    char name[256];
    if(zoom) sprintf(name, "%sRmsTimeChan%04i_%.3i_days_zoom.%s", fMetricName.c_str(), fChannel, fNumDays, fFileExtension.c_str());
    else sprintf(name, "%sRmsTimeChan%04i_%.3i_days.%s", fMetricName.c_str(), fChannel, fNumDays, fFileExtension.c_str());
    return std::string(name);
  }
};

struct NearlinePlot{
  
  TH1F* fHistogram;
  TGraph* fGraphMetricTime;
  TGraph* fGraphMetricRmsTime;
  std::string fHistName;
  std::string fHistTitle;
  std::string fHistOutputName;
  std::string fGraphOutputName;
  std::string fGraphRmsOutputName;
  std::string fGraphOutputCanName;
  std::string fGraphRmsOutputCanName;

  NearlinePlotInfo fPlotInfo;


  int fNumPoints;
  std::vector<float> fMetricVec;
  std::vector<float> fMetricRmsVec;
  std::vector<float> fTimeVec;
  int fPlotCount;
  

  bool fMakeHistogram;



  NearlinePlot(std::string this_hist_name, std::string this_hist_output_name, std::string this_graph_output_name, std::string this_graph_rms_output_name, std::string this_can_metric_name, std::string this_can_metric_rms_name);


  bool AddHistogram(TFile const & file, TTree* header, int Xsrtime, int XNow, int GMToffset);
  TCanvas* makeHistoCanvas(std::string can_name, std::string can_title, int width, int height, TPaveText* updateText);
  TCanvas* makeGraphMetricTimeCanvas(int width, int height, TPaveText* updateText, std::string taxis_labels, int time_ago, int XNow, bool zoom=false);
  TCanvas* makeGraphMetricRmsTimeCanvas(int width, int height, TPaveText* updateText, std::string taxis_labels, int time_ago, int XNow, bool zoom=false);
  void setHistTitle(std::string hist_title);
  void setPlotInfo(NearlinePlotInfo this_plot_info);

};


NearlinePlot::NearlinePlot(std::string this_hist_name, std::string this_hist_output_name, std::string this_graph_output_name, std::string this_graph_rms_output_name, std::string this_can_metric_name, std::string this_can_metric_rms_name)
{
  fNumPoints = 0;
  fHistName = this_hist_name;
  fHistOutputName = this_hist_output_name;
  fGraphOutputName = this_graph_output_name;
  fGraphRmsOutputName = this_graph_rms_output_name;
  fGraphOutputCanName = this_can_metric_name;
  fGraphRmsOutputCanName = this_can_metric_rms_name;
  fMetricVec.resize(0);
  fMetricRmsVec.resize(0);
  fTimeVec.resize(0);
  fPlotCount = 0;
  fHistogram=0;// = new TH1F(fHistName.c_str(), hist_title.c_str(), this_binning.num_bins_x, this_binning.x_min, this_binning.x_max);
  fGraphMetricTime = 0;

}

void NearlinePlot::setHistTitle(std::string hist_title){
  fHistTitle = hist_title;
}

void NearlinePlot::setPlotInfo(NearlinePlotInfo this_plot_info){
  fPlotInfo = this_plot_info;
}

bool NearlinePlot::AddHistogram(TFile const & file, TTree* header, int Xsrtime, int XNow, int GMToffset){
  TH1F *hist_temp = (TH1F*)file.FindObjectAny(fHistName.c_str());    
  if(hist_temp != 0){
    if(fHistogram == 0){
      fHistogram = (TH1F*) hist_temp->Clone((fHistName + "_temp").c_str());
      std::cerr << "hist name: " << fHistogram->GetName() << std::endl;
      if(fHistTitle!=0) fHistogram->SetTitle(fHistTitle.c_str());
      fHistogram->SetDirectory(0);
    }
    else fHistogram->Add(hist_temp,1.0);
  }    
  else{
    std::cerr << "INFO: Failed to find Histogram " << fHistName << std::endl;
    return false;
  }
  if(hist_temp != 0 && header != 0 && Xsrtime != XNow - GMToffset) {
    fTimeVec.push_back(Xsrtime);
    fMetricVec.push_back(hist_temp->GetMean(1));
    fMetricRmsVec.push_back(hist_temp->GetRMS(1));
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
