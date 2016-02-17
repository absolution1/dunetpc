#ifndef _NLPlotMkr
#define _NLPlotMkr

#include "TDatime.h"
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

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <fstream>

// NOTE:  Go through ALL of the above and remove what isn't necessary!!!
#include <sstream>


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

struct NearlinePlotLogScale{
  bool fMetricTimeGraphLog;
  bool fHisto1DLog;
  bool fHisto2DLog;
  bool fBinByBinLog;
  NearlinePlotLogScale(bool metric_time_graph_log=true, bool histo_1d_log=true, bool histo_2d_log=true, bool bin_by_bin_log=true){
    fMetricTimeGraphLog=metric_time_graph_log;
    fHisto1DLog=histo_1d_log;
    fHisto2DLog=histo_2d_log;
    fBinByBinLog=bin_by_bin_log;
  }
};


struct NearlinePlotEnables{
  bool fMakeMetricTimeGraph;
  bool fMake2DHisto;
  bool fMakeBinByBinPlots;
  bool fNormaliseHisto1D;
  NearlinePlotEnables(bool normalise_histo_1d=false, bool make_metric_time_graph=true, bool make_2d_histo=true, bool make_bin_by_bin_plots=false){
    fMakeMetricTimeGraph=make_metric_time_graph;
    fMake2DHisto=make_2d_histo;
    fMakeBinByBinPlots=make_bin_by_bin_plots;
    fNormaliseHisto1D=normalise_histo_1d;
  }
};

struct NearlinePlotInfo{

  std::string fMetricName;
  int fChannel;
  int fNumDays;
  std::string fFileExtension;
  std::string fMetricDetails;

  void AddMetricDetails(std::string metric_details);
  std::string GetMetricMeanTimeGraphName(bool zoom=false);
  std::string GetMetricRmsTimeGraphName(bool zoom=false);
  std::string GetBinByBinTimeGraphName(unsigned int bin);
  

  NearlinePlotInfo(){};

  NearlinePlotInfo(std::string metric_name, int channel, int num_days, std::string file_extension, std::string metric_details=""){
    fMetricName = metric_name;
    fChannel = channel;
    fNumDays = num_days;
    fFileExtension = file_extension;
    fMetricDetails = metric_details;
  }
  
  std::string GetHistOutputName(){
    char name[256];
    if(fChannel >= 0) sprintf(name, "%sChan%04i_%.3i_days.%s", fMetricName.c_str(), fChannel, fNumDays, fFileExtension.c_str());
    else sprintf(name, "%s_%.3i_days.%s", fMetricName.c_str(), fNumDays, fFileExtension.c_str());
    return std::string(name);
  }

  std::string GetHist2DOutputName(){
    char name[256];
    if(fChannel >=0) sprintf(name, "%sVsTimeChan%04i_%.3i_days.%s", fMetricName.c_str(), fChannel, fNumDays, fFileExtension.c_str());
    else sprintf(name, "%sVsTime_%.3i_days.%s", fMetricName.c_str(), fNumDays, fFileExtension.c_str());
    return std::string(name);
  }

};

void NearlinePlotInfo::AddMetricDetails(std::string metric_details){
  fMetricDetails = metric_details;
}

std::string NearlinePlotInfo::GetMetricMeanTimeGraphName(bool zoom){
  char name[256];

  if(fChannel >= 0) sprintf(name, "%sMeanTimeChan%04i_%.3i_days", fMetricName.c_str(), fChannel, fNumDays);
  else  sprintf(name, "%sMeanTime_%.3i_days", fMetricName.c_str(), fNumDays);

  if(zoom) sprintf(name, "%s_zoom.%s", name, fFileExtension.c_str());
  else sprintf(name, "%s.%s", name, fFileExtension.c_str());

  return std::string(name);
}

std::string NearlinePlotInfo::GetMetricRmsTimeGraphName(bool zoom){
  char name[256];
  if(fChannel >= 0) sprintf(name, "%sRmsTimeChan%04i_%.3i_days", fMetricName.c_str(), fChannel, fNumDays);
  else  sprintf(name, "%sRmsTime_%.3i_days", fMetricName.c_str(), fNumDays);

  if(zoom) sprintf(name, "%s_zoom.%s", name, fFileExtension.c_str());
  else sprintf(name, "%s.%s", name, fFileExtension.c_str());

  return std::string(name);
}

std::string NearlinePlotInfo::GetBinByBinTimeGraphName(unsigned int bin){
  char name[256];
  if(fChannel>=0) sprintf(name, "%sTimeChan%04iBin%u_%.3i_days.%s", fMetricName.c_str(), fChannel, bin, fNumDays, fFileExtension.c_str());
  else sprintf(name, "%sTimeBin%u_%.3i_days.%s", fMetricName.c_str(), bin, fNumDays, fFileExtension.c_str());
  return std::string(name);
}




struct NearlinePlot{

  NearlinePlotEnables fPlotEnables;
  NearlinePlotLogScale fPlotLogScale;

  TH1F* fHistogram;
  std::string fHistName;
  NearlinePlotInfo fPlotInfo;
  std::string fHistTitle;

  TH2F* fHistogram2D;
  TH1F* fHistogram2DNormalisation;
  bool fNormalised;

  TGraph* fGraphMetricTime;
  TGraph* fGraphMetricRmsTime;
  std::vector<float> fMetricVec;
  std::vector<float> fMetricRmsVec;
  std::vector<float> fTimeVec;
  int fPlotCount;

  std::vector<TGraphErrors*> fBinByBinGraphMetricTime;
  std::vector<std::vector<float>> fBinByBinMetricVec;
  std::vector<std::vector<float>> fBinByBinMetricErrorVec;
  std::vector<std::string> fBinByBinLabels;
  std::string fBinByBinYAxisTitle;

  NearlinePlot(std::string this_hist_name, NearlinePlotInfo this_plot_info, NearlinePlotEnables this_plot_enable=NearlinePlotEnables(), NearlinePlotLogScale this_plot_log_scale=NearlinePlotLogScale());



  bool AddHistogram(TFile const & file, TTree* header, int Xsrtime, int XNow, int GMToffset, int time_ago);
  bool AddHistogram1D(TFile const & file, TTree* header, int Xsrtime, int XNow, int GMToffset, int time_ago);
  bool AddHistogram2D(TFile const & file, TTree* header, int Xsrtime, int XNow, int GMToffset, int time_ago);
  void normaliseHisto1D();
  
  TCanvas* makeHistoCanvas(TPaveText* updateText, int width=1200, int height=800);
  TCanvas* makeHisto2DCanvas(TPaveText* updateText, int time_ago, int XNow, int width=1200, int height=800, std::string taxis_labels="");
  TCanvas* makeGraphMetricTimeCanvas(TPaveText* updateText, int time_ago, int XNow, bool rms=false, bool zoom=false, int width=1200, int height=800, std::string taxis_labels="");
  TCanvas* makeBinByBinGraphTime(unsigned int bin, TPaveText* updateText, int time_ago, int XNow, int width=1200, int height=800, std::string taxis_labels="");
  void printPlots(std::string plot_dir, TPaveText* updateText, int time_ago, int XNow, int width=1200, int height=800, std::string taxis_labels="");
  void printHistogram1D(std::string plot_dir, TPaveText* updateText, int time_ago, int XNow, int width=1200, int height=800, std::string taxis_labels="");
  void printHistogram2D(std::string plot_dir, TPaveText* updateText, int time_ago, int XNow, int width=1200, int height=800, std::string taxis_labels="");
  void printGraphs(std::string plot_dir, TPaveText* updateText, int time_ago, int XNow, int width=1200, int height=800, std::string taxis_labels="");
  void printBinByBinGraphs(std::string plot_dir, TPaveText* updateText, int time_ago, int XNow, int width=1200, int height=800, std::string taxis_labels="");
  void setMetricDetails(std::string metric_details);

  void setHistTitle(std::string hist_title);
  void setPlotInfo(NearlinePlotInfo this_plot_info);

};


NearlinePlot::NearlinePlot(std::string this_hist_name, NearlinePlotInfo this_plot_info, NearlinePlotEnables this_plot_enable, NearlinePlotLogScale this_plot_log_scale)
{

  fPlotEnables=this_plot_enable;
  fPlotLogScale=this_plot_log_scale;

  fHistogram=0;
  fHistName = this_hist_name;  
  fPlotInfo = this_plot_info;
  fHistTitle = "";

  fHistogram2D=0;
  fHistogram2DNormalisation=0;
  fNormalised = false;

  fMetricVec.resize(0);
  fMetricRmsVec.resize(0);
  fTimeVec.resize(0);
  fPlotCount = 0;
  fGraphMetricTime = 0;
  fGraphMetricRmsTime = 0;

  fBinByBinGraphMetricTime.resize(0);
  fBinByBinMetricVec.resize(0);
  fBinByBinMetricErrorVec.resize(0);
  fBinByBinLabels.resize(0);
  fBinByBinYAxisTitle = "";
}

void NearlinePlot::setHistTitle(std::string hist_title){
  fHistTitle = hist_title;
}

void NearlinePlot::setPlotInfo(NearlinePlotInfo this_plot_info){
  fPlotInfo = this_plot_info;
}

void NearlinePlot::setMetricDetails(std::string metric_details){
  fPlotInfo.AddMetricDetails(metric_details);
}


bool NearlinePlot::AddHistogram(TFile const & file, TTree* header, int Xsrtime, int XNow, int GMToffset, int time_ago){

  bool result;
  result = AddHistogram1D(file, header, Xsrtime, XNow, GMToffset, time_ago);
  if(fPlotEnables.fMake2DHisto) result = result & AddHistogram2D(file, header, Xsrtime, XNow, GMToffset, time_ago);
  return result;

}

bool NearlinePlot::AddHistogram1D(TFile const & file, TTree* header, int Xsrtime, int XNow, int GMToffset, int time_ago){

  TH1F *hist_temp = (TH1F*)file.FindObjectAny(fHistName.c_str());    
  if(hist_temp != 0){
    if(fHistogram == 0){
      fHistogram = (TH1F*) hist_temp->Clone((fHistName + "_temp").c_str());
      if(fHistTitle!="") fHistogram->SetTitle(fHistTitle.c_str());
      else fHistTitle = fHistogram->GetTitle();
      fHistogram->SetDirectory(0);
    }
    else{
      if(fHistogram->GetNbinsX() != hist_temp->GetNbinsX()){
	//	std::cerr << "ERROR: " << file.GetName() << " " << fHistogram->GetName() << " " << fHistogram->GetNbinsX() << " bins vs " << hist_temp->GetNbinsX() << std::endl;
	return false;
      }
      if(hist_temp->GetEntries() == 0){
	//	std::cerr << "ERROR: " << file.GetName() << " " << fHistogram->GetName() << " Has no entries" << std::endl;
	return false;
      }
      fHistogram->Add(hist_temp,1.0);
    }
  }    
  else{
    //    std::cerr << "NearlinePlot::AddHistogram - Failed to find histogram - " << fHistName << std::endl;
    return false;
  }
  if(hist_temp != 0 && header != 0 && Xsrtime != XNow - GMToffset) {
    fTimeVec.push_back(Xsrtime);
    fMetricVec.push_back(hist_temp->GetMean(1));
    fMetricRmsVec.push_back(hist_temp->GetRMS(1));

    if(fPlotEnables.fMakeBinByBinPlots){
      //Initiliase the binbybin vector of vectors
      if(fBinByBinMetricVec.size()==0){
        fBinByBinMetricVec.resize(hist_temp->GetNbinsX());
        fBinByBinMetricErrorVec.resize(hist_temp->GetNbinsX());
	fBinByBinLabels.resize(hist_temp->GetNbinsX());
	fBinByBinYAxisTitle = std::string(hist_temp->GetYaxis()->GetTitle());
	for(int bin=1; bin<=hist_temp->GetNbinsX();bin++){
	  std::string label = hist_temp->GetXaxis()->GetBinLabel(bin);
	  fBinByBinLabels.at(bin-1) = (label);
	  //	  std::cerr << "ERROR: " << hist_temp->GetName() << " bin " << bin << " " << label << std::endl;
	}	
      }
      if(fHistogram->GetNbinsX() != hist_temp->GetNbinsX()){
	//	std::cerr << "ERROR: " << file.GetName() << " " << fHistogram->GetName() << " " << fHistogram->GetNbinsX() << " bins vs " << hist_temp->GetNbinsX() << std::endl;
	return false;;
      }
      if(hist_temp->GetEntries() == 0){
	//	std::cerr << "ERROR: " << file.GetName() << " " << fHistogram->GetName() << " Has no entries" << std::endl;
	return false;
      }
      for(int bin=1; bin<=hist_temp->GetNbinsX();bin++){
        fBinByBinMetricVec.at(bin-1).push_back(hist_temp->GetBinContent(bin));
        fBinByBinMetricErrorVec.at(bin-1).push_back(hist_temp->GetBinError(bin));
      }
    }//fMakeBinByBinPlots

    fPlotCount++;
    //    std::cerr << "NearlinePlot::AddHistogram - Success" << std::endl;
    return true;
  }

  std::cerr << "NearlinePlot::AddHistogram - Failed - Shouldn't get here" << std::endl;
  return false;

}

bool NearlinePlot::AddHistogram2D(TFile const & file, TTree* header, int Xsrtime, int XNow, int GMToffset, int time_ago){

  TH1F *hist_temp = (TH1F*)file.FindObjectAny(fHistName.c_str());    
  if(hist_temp != 0){
    if(fHistogram2D == 0){
      int nbinsx = hist_temp->GetNbinsX();
      double xmin = hist_temp->GetBinLowEdge(1);
      double xmax = hist_temp->GetBinLowEdge(1+nbinsx);
      int nbinst = 72;//FIXME
      double tmin = time_ago;
      double tmax = XNow;
      fHistogram2D = new TH2F((fHistName + "_vs_time").c_str(), (std::string(hist_temp->GetTitle()) + " - vs Time").c_str(), nbinst, tmin, tmax, nbinsx, xmin, xmax);
      //      fHistogram2D->GetXaxis()->SetTitle("x axis");
      fHistogram2D->GetYaxis()->SetTitle(hist_temp->GetXaxis()->GetTitle());
      fHistogram2D->GetZaxis()->SetTitle(hist_temp->GetYaxis()->GetTitle());

      //Also make the normalisation histogram
      fHistogram2DNormalisation = new TH1F((fHistName + "_normalisation").c_str(), (std::string(hist_temp->GetName()) + " Normalisation").c_str(), nbinst, tmin, tmax);

      //Make sure we keep the histograms in memory
      fHistogram2D->SetDirectory(0);
      fHistogram2DNormalisation->SetDirectory(0);

      //Get bin labels
      int nbins = hist_temp->GetNbinsX();
      for(int i=0;i<=nbins+1;i++){
	std::string label = hist_temp->GetXaxis()->GetBinLabel(i);
	if(label != "") {
	  fHistogram2D->GetYaxis()->SetBinLabel(i, label.c_str());
	  fHistogram2D->GetYaxis()->SetLabelSize(hist_temp->GetXaxis()->GetLabelSize());
	}
      }//bins
    }//make TH2
  }//Got TH1
  else{
    //    std::cerr << "NearlinePlot::AddHistogram2D - Failed to find histogram - " << fHistName << std::endl;
    return false;
  }
  if(hist_temp != 0 && header != 0 && Xsrtime != XNow - GMToffset) {  
    //Now add this th1 to the th2
    int nbins = hist_temp->GetNbinsX();
    if(fHistogram->GetNbinsX() != hist_temp->GetNbinsX()){
      //      std::cerr << "ERROR: " << file.GetName() << " " << fHistogram->GetName() << " " << fHistogram->GetNbinsX() << " bins vs " << hist_temp->GetNbinsX() << std::endl;
      return false;
    }
    if(hist_temp->GetEntries() == 0){
      //      std::cerr << "ERROR: " << file.GetName() << " " << fHistogram->GetName() << " Has no entries" << std::endl;
      return false;
    }

    for(int i=0;i<=nbins+1;i++){
      double content = hist_temp->GetBinContent(i);
      double center = hist_temp->GetBinCenter(i);
      fHistogram2D->Fill(Xsrtime, center, content);
    }//loops over the histogram bins
    
    //Add Entry to our normalisation
    fHistogram2DNormalisation->Fill(Xsrtime);
    return true;
  }
  //  std::cerr << "NearlinePlot::AddHistogram2D - Failed - Shouldn't get here" << std::endl;  
  return false;
}

void NearlinePlot::normaliseHisto1D(){
  if(fPlotCount>0) fHistogram->Scale(1./fPlotCount);
}

TCanvas* NearlinePlot::makeHistoCanvas(TPaveText* updateText, int width, int height){

  if(fPlotEnables.fNormaliseHisto1D) normaliseHisto1D();

  std::string can_name = fHistName + "_can";
  std::string can_title = fHistName + "_can";

  TCanvas *can = new TCanvas(can_name.c_str(), can_title.c_str(), width, height);
  can->cd();
  gPad->SetGridx();
  gPad->SetGridy();
  if(fPlotLogScale.fHisto1DLog) can->SetLogy();
  if(!(fPlotEnables.fNormaliseHisto1D)) gStyle->SetOptStat(111111);
  fHistogram->SetLineWidth(2);
  fHistogram->SetLineColor(kRed);
  fHistogram->SetName(fHistName.c_str());
  fHistogram->Draw();
  updateText->Draw();


  return can;
}

TCanvas* NearlinePlot::makeHisto2DCanvas(TPaveText* updateText, int time_ago, int XNow, int width, int height, std::string taxis_labels){

  std::string can_name = fHistName + "_vs_time_can";
  std::string can_title = fHistName + "_vs_time_can";

  TCanvas *can = new TCanvas(can_name.c_str(), can_title.c_str(), width, height);
  can->cd();
  if(fPlotLogScale.fHisto2DLog)  can->SetLogz();
  can->SetRightMargin(0.15);
  gStyle->SetOptStat(0);
  gPad->SetGridx();
  gPad->SetGridy();
  if(taxis_labels==""){
    if(fPlotInfo.fNumDays <= 2) taxis_labels = "%H:%M";
    else taxis_labels = "%m/%d";
  }

  if(fNormalised!=true){
    for(int xbin=1;xbin<=fHistogram2D->GetNbinsX();xbin++){
      double norm = fHistogram2DNormalisation->GetBinContent(xbin);
      if(norm <= 0.0) continue;
      for(int ybin=1;ybin<=fHistogram2D->GetNbinsY();ybin++){
        int global_bin = fHistogram2D->GetBin(xbin, ybin);
        double content = fHistogram2D->GetBinContent(global_bin);
	if(content <= 0.0) continue;
        fHistogram2D->SetBinContent(global_bin, content/norm);
      }//ybin
    }//xbin
    fNormalised=true;
  }//Normalise


  fHistogram2D->GetXaxis()->SetTimeFormat(taxis_labels.c_str());
  fHistogram2D->GetXaxis()->SetLimits(time_ago,XNow);
  fHistogram2D->GetXaxis()->SetTitle("(central time)");
  fHistogram2D->GetXaxis()->SetTimeDisplay(1);
  fHistogram2D->GetXaxis()->SetLabelSize(0.03);
  fHistogram2D->Draw("COLZ");
  updateText->Draw();


  return can;
}

TCanvas* NearlinePlot::makeBinByBinGraphTime(unsigned int bin, TPaveText* updateText, int time_ago, int XNow, int width, int height, std::string taxis_labels){


  if(bin >= fBinByBinMetricVec.size()) return NULL;

  std::string can_name = fHistName + "_time_bin_";
  std::string can_title = fHistTitle;
  if(fBinByBinLabels.at(bin) == "" ){
    can_name += std::to_string(bin);
    can_title +=  + " Bin " + std::to_string(bin);
  }
  else{
    can_title += " " + fBinByBinLabels.at(bin);
  }

  TCanvas* can = new TCanvas(can_name.c_str(), can_title.c_str(), width, height);
  can->cd();
  gPad->SetGridx();

  TGraphErrors* gr = new TGraphErrors(fPlotCount);
  for(int i=0;i<fPlotCount;i++) gr->SetPoint(i, fTimeVec.at(i), fBinByBinMetricVec.at(bin).at(i));
  for(int i=0;i<fPlotCount;i++) gr->SetPointError(i, 0, fBinByBinMetricErrorVec.at(bin).at(i));

  gr->SetTitle(can_title.c_str());
  gr->SetMarkerColor(kBlue);
  gr->GetYaxis()->SetTitle(fBinByBinYAxisTitle.c_str());
  gr->GetXaxis()->SetTimeDisplay(1);
  gr->GetXaxis()->SetLabelSize(0.03);

  if(taxis_labels==""){
    if(fPlotInfo.fNumDays <= 2) taxis_labels = "%H:%M";
    else taxis_labels = "%m/%d";
  }

  gr->GetXaxis()->SetTimeFormat(taxis_labels.c_str());
  gr->GetXaxis()->SetLimits(time_ago,XNow);
  gr->GetXaxis()->SetTitle("(central time)");
  gr->Draw("A*");

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

TCanvas* NearlinePlot::makeGraphMetricTimeCanvas(TPaveText* updateText, int time_ago, int XNow, bool rms, bool zoom, int width, int height, std::string taxis_labels){

  std::string can_name;
  if(rms) can_name= fHistName + "_rms_time_can";
  else can_name= fHistName + "_metric_time_can";

  std::string can_title;
  if(rms)  can_title = "RMS of " + fHistTitle;
  else  can_title = "Mean of " + fHistTitle;


  if(zoom) can_name += "_zoom";


  TCanvas* can = new TCanvas(can_name.c_str(), can_title.c_str(), width, height);
  can->cd();
  gPad->SetGridx();

  TGraph* gr;
  std::vector<float> yValuesVec;
  if(rms){
    fGraphMetricRmsTime = new TGraph(fPlotCount);
    gr = fGraphMetricRmsTime;
    yValuesVec = fMetricRmsVec;
  }
  else{
    fGraphMetricTime = new TGraph(fPlotCount);
    gr = fGraphMetricTime;
    yValuesVec = fMetricVec;
  }

  
  for(int i=0;i<fPlotCount;i++) gr->SetPoint(i, fTimeVec.at(i), yValuesVec.at(i));
  
  gr->SetMarkerColor(kBlue);
  gr->GetXaxis()->SetTimeDisplay(1);
  gr->GetXaxis()->SetLabelSize(0.03);

  if(taxis_labels==""){
    if(fPlotInfo.fNumDays <= 2) taxis_labels = "%H:%M";
    else taxis_labels = "%m/%d";
  }

  gr->GetXaxis()->SetTimeFormat(taxis_labels.c_str());
  gr->GetXaxis()->SetLimits(time_ago,XNow);
  gr->GetXaxis()->SetTitle("(central time)");
  gr->Draw("A*");

  if(zoom){
    graphZoom(gr, 2.0);
    gr->SetTitle(std::string(can_title + " - Zoom").c_str());
  }
  else   gr->SetTitle(can_title.c_str());

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

void NearlinePlot::printPlots(std::string plot_dir, TPaveText* updateText, int time_ago, int XNow, int width, int height, std::string taxis_labels){

  printHistogram1D(plot_dir, updateText, time_ago, XNow, width, height, taxis_labels);

  if(fPlotEnables.fMake2DHisto) printHistogram2D(plot_dir, updateText, time_ago, XNow, width, height, taxis_labels);
  
  if(fPlotEnables.fMakeMetricTimeGraph) printGraphs(plot_dir, updateText, time_ago, XNow, width, height, taxis_labels);

  if(fPlotEnables.fMakeBinByBinPlots) printBinByBinGraphs(plot_dir, updateText, time_ago, XNow, width, height, taxis_labels);

}

void NearlinePlot::printHistogram1D(std::string plot_dir, TPaveText* updateText, int time_ago, int XNow, int width, int height, std::string taxis_labels){
  TCanvas* can = makeHistoCanvas(updateText,  width,  height);
  std::string can_name = plot_dir + "/" + fPlotInfo.GetHistOutputName();
  can->Print(can_name.c_str());
  delete can;
}

void NearlinePlot::printHistogram2D(std::string plot_dir, TPaveText* updateText, int time_ago, int XNow, int width, int height, std::string taxis_labels){
  TCanvas* can = makeHisto2DCanvas(updateText, time_ago, XNow,  width,  height, taxis_labels);
  std::string can_name = plot_dir + "/" + fPlotInfo.GetHist2DOutputName();
  can->Print(can_name.c_str());
  delete can;
}

void NearlinePlot::printBinByBinGraphs(std::string plot_dir, TPaveText* updateText, int time_ago, int XNow, int width, int height, std::string taxis_labels){

  TCanvas *can;
  std::string can_name;

  for(unsigned int bin=0;bin<fBinByBinMetricVec.size();bin++){ 
    can = makeBinByBinGraphTime(bin, updateText, time_ago, XNow); 
    can_name = plot_dir + "/" + fPlotInfo.GetBinByBinTimeGraphName(bin); 
    can->Print(can_name.c_str()); 
    delete can; 
  }//bin 

}

void NearlinePlot::printGraphs(std::string plot_dir, TPaveText* updateText, int time_ago, int XNow, int width, int height, std::string taxis_labels){

  bool zoom;
  bool rms;
  TCanvas *can;
  std::string can_name;

  rms=false;  

  zoom=false;
  can = makeGraphMetricTimeCanvas(updateText, time_ago, XNow, rms, zoom);
  can_name = plot_dir + "/" + fPlotInfo.GetMetricMeanTimeGraphName(zoom);
  can->Print(can_name.c_str());

  delete can;

  zoom=true;
  can = makeGraphMetricTimeCanvas(updateText, time_ago, XNow, rms, zoom);
  can_name = plot_dir + "/" + fPlotInfo.GetMetricMeanTimeGraphName(zoom);
  can->Print(can_name.c_str());
  
  delete can;

  rms=true;  

  zoom=false;
  can = makeGraphMetricTimeCanvas(updateText, time_ago, XNow, rms, zoom);
  can_name = plot_dir + "/" + fPlotInfo.GetMetricRmsTimeGraphName(zoom);
  can->Print(can_name.c_str());

  delete can;

  zoom=true;
  can = makeGraphMetricTimeCanvas(updateText, time_ago, XNow, rms, zoom);
  can_name = plot_dir + "/" + fPlotInfo.GetMetricRmsTimeGraphName(zoom);
  can->Print(can_name.c_str());
  
  delete can;

}



struct NearlineHTML{

  static std::string MakeStartYearPlot(std::string relative_plot_path){

    std::string output;
    output+="<h3>Start Year of the Run VS Run Number.</h3>\n";
    output+="<figure>\n";
    output+="<img src=\"" + relative_plot_path + "\" width=\"800\">\n";
    output+="</figure>\n";
    output+="<p>\n";
    output+="<b>This plot requires some explanation:</b>\n";
    output+="I noticed that for a large number of runs, some or all of the data appeared to be missing from the events in the files. This resulted in one of two behaviors:\n";
    output+="<BR>\n";
    output+="<b>1)</b> The run number was correct in the art::events but the year cooresponding to the event time was 0 (which defaults to 1969.)\n";
    output+="<BR>\n";
    output+="<b>2)</b> The file from that run was entirely devoid of usable events making the run number and year default to zero (again, 1969 for the year.)\n";
    output+="<BR>\n";
    output+="Both the run number and the year are obtained from the art::events by the nearline analyzer module. So this plot should be interpreted as a very general DAQ-health indicator. Situation 1 can occur if the art::events are present in the file but there is no RCE information. Situation 2 can occur if there are no good events in the entire file. I would imagine that once the DAQ is stable during normal running, his plot will become obsolete.\n";
    output+="</p>\n";
    output+="<BR><BR><BR>\n";

    return output;
  }

  static std::string MakeHistogram(std::string plot_location, NearlinePlot* plot){ 

    NearlinePlotInfo plot_info = plot->fPlotInfo;
    std::string metric_name = plot_info.fMetricName;
    int channel = plot_info.fChannel;
    std::string metric_details = plot_info.fMetricDetails;
    plot_location = plot_location + "/" + plot_info.GetHistOutputName();

    std::string output;
    output += "<hr>\n";
    if(channel >=0) output += "<h3>" + metric_name + " spectrum for channel # " + std::to_string(channel);
    else output += "<h3>" + metric_name + " spectrum";
    if(metric_details!="") output += " " + metric_details; 
    output += ".</h3>\n";
    output += "<figure>\n";
    output += "<img src=\"" + plot_location + "\" width=\"800\">\n"; 
    output += "</figure>\n"; 
    return output;
  }

  static std::string MakeHistogram2D(std::string plot_location, NearlinePlot* plot){

    NearlinePlotInfo plot_info = plot->fPlotInfo;
    std::string metric_name = plot_info.fMetricName;
    int channel = plot_info.fChannel;
    std::string metric_details = plot_info.fMetricDetails;
    plot_location = plot_location + "/" + plot_info.GetHist2DOutputName();

    std::string output;
    output += "<hr>\n";
    if(channel >= 0)    output += "<h3>" + metric_name + " spectrum vs. time for channel # " + std::to_string(channel);
    else     output += "<h3>" + metric_name + " spectrum vs. time";
    if(metric_details!="") output += " " + metric_details; 
    output += ".</h3>\n";
    output += "<figure>\n";
    output += "<img src=\"" + plot_location + "\" width=\"800\">\n"; 
    output += "</figure>\n"; 
    return output;

  }

  static std::string MakeHistogramPair(std::string plot_location, NearlinePlot* plot){

    NearlinePlotInfo plot_info = plot->fPlotInfo;
    std::string metric_name = plot_info.fMetricName;
    int channel = plot_info.fChannel;
    std::string metric_details = plot_info.fMetricDetails;

    std::string plot_location_histo = plot_location + "/" + plot_info.GetHistOutputName();
    std::string plot_location_histo2d = plot_location + "/" + plot_info.GetHist2DOutputName();

    std::string output;
    output += "<hr>\n";
    if(channel >= 0) output += "<h3>" + metric_name + " spectrum  for channel # " + std::to_string(channel);
    output += "<h3>" + metric_name + " spectrum";
    if(metric_details!="") output += " " + metric_details; 
    output += ".</h3>\n";
    output += "<table>\n";
    output += "<tr>\n";
    output += "<td> <img src=\"" + plot_location_histo + "\" width=\"800\"></td>\n"; 
    output += "<td> <img src=\"" + plot_location_histo2d + "\" width=\"800\"></td>\n"; 
    output += "</tr>\n";
    output += "</table>\n";

    return output;

  }

  
  static std::string MakeGraphPair(std::string plot_location, NearlinePlot* plot, bool rms){

    NearlinePlotInfo plot_info = plot->fPlotInfo;
    std::string metric_name = plot_info.fMetricName;
    int channel = plot_info.fChannel;
    std::string metric_details = plot_info.fMetricDetails;

    std::string plot_location_zoom;
    std::string plot_location_unzoom;
    if(rms){
      plot_location_zoom = plot_location + "/" + plot_info.GetMetricRmsTimeGraphName(true);
      plot_location_unzoom = plot_location + "/" + plot_info.GetMetricRmsTimeGraphName();
    }
    else{
      plot_location_zoom = plot_location + "/" + plot_info.GetMetricMeanTimeGraphName(true);
      plot_location_unzoom = plot_location + "/" + plot_info.GetMetricMeanTimeGraphName();
    }
    
    std::string output;
    output+="<h3>";
    if(rms) output += "RMS of the ";
    else output += "Mean of the ";
    if(channel >= 0) output += metric_name + " spectra for channel # " + std::to_string(channel);
    else output += metric_name + " spectra";
    if(metric_details!="") output += " " + metric_details;
    output +=".</h3>\n";
    output +="<table>\n";
    output +="<tr>\n";
    output +="<td> <img src=\"" + plot_location_unzoom + "\" width=\"800\"></td>\n";
    output +="<td> <img src=\"" + plot_location_zoom + "\" width=\"800\"></td>\n";
    output +="</tr>\n";
    output +="</table>\n";

    return output;
    
  }

  static std::string MakePlotSet(std::string plot_location, NearlinePlot* plot){

    NearlinePlotInfo plot_info=plot->fPlotInfo;
    NearlinePlotEnables plot_enables=plot->fPlotEnables;

    std::string metric_name = plot_info.fMetricName;
    std::string output;
    output += "\n\n\n";

    if(plot_enables.fMake2DHisto) output += MakeHistogramPair(plot_location, plot) + "\n";
    else output += MakeHistogram(plot_location, plot) + "\n";
    //    output += "\n";
    if(plot_enables.fMakeMetricTimeGraph) output += MakeGraphPair(plot_location, plot, false) + "\n";
    //    output += "\n";
    if(plot_enables.fMakeMetricTimeGraph) output += MakeGraphPair(plot_location, plot, true) +"\n";
    output += "\n\n";
    output += "<p>\n";
    output += "First plot: " + metric_name + " spectra on a specified channel in the specified time period.\n";
    output += "Second plot: Mean of the top plot over time (each blue dot corresponds to one subrun).\n";
    output += "Third plot: RMS of the top plot over time (each blue dot corresponds to one subrun).\n";
    output += "</p>\n";
    output += "<BR><BR><BR>\n\n";

    return output;

  }

  static std::string MakeBinByBinGraphs(std::string plot_location, NearlinePlot* plot){

    NearlinePlotEnables plot_enables = plot->fPlotEnables;
    NearlinePlotInfo plot_info = plot->fPlotInfo;
    if(!plot_enables.fMakeBinByBinPlots) return "";

    std::string output;
    
    for(size_t bin=0; bin<plot->fBinByBinMetricVec.size();bin++){
      std::string metric_name = plot_info.fMetricName;
      std::string plot_name = plot_info.GetBinByBinTimeGraphName(bin); 


      if(bin+1 == (plot->fBinByBinMetricVec.size())){
	output += "<h3>";
	output += metric_name + " - ";
	if(plot->fBinByBinLabels.at(bin)!="") output += plot->fBinByBinLabels.at(bin);
	else output += "bin " + std::to_string(bin);
	output += ".</h3>\n";
	output += "<figure>\n";
	output += "<img src=\"" + plot_location + "/" + plot_name + "\" width=\"800\">\n";
	output += "</figure>\n";
	output += "\n\n\n";
      }//one bin
      else{
	std::string plot_name_next = plot_info.GetBinByBinTimeGraphName(bin+1); 
	output += "<h3>" + metric_name + " - ";
	
	if(plot->fBinByBinLabels.at(bin)!="") output += plot->fBinByBinLabels.at(bin);	
	else output += "bin " + std::to_string(bin);
	
	output += " and ";
	if(plot->fBinByBinLabels.at(bin+1)!="") output += plot->fBinByBinLabels.at(bin+1);	
	else output += "bin " + std::to_string(bin+1);
	
	output +="<table>\n";
	output +="<tr>\n";
	output +="<td> <img src=\"" + plot_location + "/" + plot_name + "\" width=\"800\"></td>\n";
	output +="<td> <img src=\"" + plot_location + "/" + plot_name_next + "\" width=\"800\"></td>\n";
	output +="</tr>\n";
	output +="</table>\n";
	output += "\n\n\n";
	bin++;
      }//two bins

    }//loop over bins
    return output;

  }

  static std::string MakePageHeader(int Ndays){

    std::string output;
    output+="<head>\n";
    output+="<meta http-equiv=\"refresh\" content=\"60\">\n";
    output+="</head>\n\n\n";
    if(Ndays<=2) output+="<h1>35T Nearline Monitoring - 24 Hour Plots</h1>\n\n";
    else if(Ndays==7) output+="<h1>35T Nearline Monitoring - 7 Days Plots</h1>\n\n";
    else if(Ndays==31) output+="<h1>35T Nearline Monitoring - 31 Days Plots</h1>\n\n";
    else output+="<h1>35T Nearline Monitoring - " + std::to_string(Ndays) + " Days Plots</h1>\n\n";
    output+="<h2><b>THIS PAGE IS CURRENTLY UNDER CONSTRUCTION.</b></h2>\n";
    output+="<BR><BR><BR>\n";

    return output;
  }

};

struct NearlineProcessingVersion{

  std::vector<std::string> fVersionStrings;
  std::vector<int> fVersionInts;
  std::vector<int> fRuns;
  TGraph* fGr;

  NearlineProcessingVersion(){
    fVersionStrings.resize(0);
    fVersionInts.resize(0);
    fRuns.resize(0);
    fGr=0;
  }

  std::string GetNearlineVersionFromFileName(std::string fileName){

    std::string file_location = "/lbne/data2/users/lbnedaq/nearline/";
    size_t file_location_size = file_location.size();
    
    std::string output = fileName.substr(file_location_size);
    output = output.substr(0, output.find("/"));
    
    if(output.size()==0){
      output = "Unknown";
    }
    else if(output.find("v") != 0){
      output = "Unknown";
    }
    return output;
  }//GetNearlineVersionFromFileName

  void AddFile(std::string fileName, int run){
    if(run<=0) return;
    std::string version = GetNearlineVersionFromFileName(fileName);
    //see if this is a new version
    int version_index=-1;
    for(size_t i=0;i<fVersionStrings.size();i++){
      if(version == fVersionStrings.at(i)) version_index = i;
    }//versions
    if(version_index == -1){
      version_index = fVersionStrings.size();
      fVersionStrings.push_back(version);
    }
    
    fVersionInts.push_back(version_index);
    fRuns.push_back(run);    
    std::cerr << "INFO : " << "run " << run << " version " << version << std::endl;
  }

  TGraph* GetGraph(){

    fGr = new TGraph(fRuns.size());
    for(size_t i=0;i<fRuns.size();i++) fGr->SetPoint(i, fRuns.at(i), fVersionInts.at(i));
    
    return fGr;
  }

  TCanvas* GetVersionCanvas(TPaveText* updateText, int time_ago, int XNow, int width=1200, int height=800){
    TCanvas* can = new TCanvas("can_processing_version", "can_processing_version", width, height);
    can->SetRightMargin(0.20);
    TGraph* gr = GetGraph();
    //    can->cd()->SetLogy();
    gr->Draw("A*");
    gr->SetTitle("Nearline Processing Version");   
    gr->GetXaxis()->SetTitle("Run Number");   

    //Set the bin labels and range correctly
    gr->GetYaxis()->Set(fVersionStrings.size(),-0.5,fVersionStrings.size()-0.5);
    gr->GetYaxis()->SetRangeUser(-0.5, fVersionStrings.size()-0.5);
    for(size_t i=0;i<fVersionStrings.size();i++){
      gr->GetYaxis()->SetBinLabel(i+1, fVersionStrings.at(i).c_str());
    }

    gr->Draw("A*");

    gPad->SetGridx();
    gPad->SetGridy();

    updateText->Draw();
    return can;
  }
  
  static std::string GetPlotName(int Ndays){ 
    char name[256];
    sprintf(name, "ProcessingVersion_%.3i_days.png", Ndays);
    return std::string(name);
  }

  
  void PrintVersionPlots(std::string plot_dir, int Ndays, TPaveText* updateText, int time_ago, int XNow, int width=1200, int height=800){
    TCanvas *can = GetVersionCanvas(updateText, time_ago, XNow, width, height);
    can->Print((plot_dir + "/" + GetPlotName(Ndays)).c_str());
    delete can;
    
    for(size_t i=0;i<fVersionStrings.size();i++){
      std::cerr << "INFO : Versions Used: " << fVersionStrings.at(i) << std::endl;
    }

  }


  static std::string MakeVersionPlotsHTML(std::string relative_plot_dir, int Ndays){

    std::string output;
    output+="<h3>Nearline Processing Version as a function of Run.</h3>\n";
    output+="<figure>\n";
    output+="<img src=\"" + relative_plot_dir + "/" + GetPlotName(Ndays) + "\" width=\"800\">\n";
    output+="</figure>\n";
    output+="<BR><BR><BR>\n";

    return output;
  }
};//NearlineProcessingVersion


struct NearlineProcessingPedestal{

  static std::string GetPedestalFileName(std::string done_file_name){
    std::ifstream in_file(done_file_name.c_str());
    std::string line;
    std::string pedestal_file="Unknown";
    while(std::getline(in_file, line)){   
      //    std::cerr << "INFO : " << line << std::endl;
      if(line.find("NEARLINE_PEDESTAL ") != std::string::npos)  pedestal_file = line.substr(std::string("NEARLINE_PEDESTAL ").size());
      
    }
    return pedestal_file;
  }

  std::vector<std::string> fPedestalStrings;
  std::vector<int> fPedestalInts;
  std::vector<int> fRuns;
  TGraph* fGr;

  NearlineProcessingPedestal(){
    fPedestalStrings.resize(0);
    fPedestalInts.resize(0);
    fRuns.resize(0);
    fGr=0;
  }
  void AddFile(std::string fileName, int run){
    if(run<=0) return;
    std::string version = GetPedestalFileName(fileName);
    //see if this is a new version
    int version_index=-1;
    for(size_t i=0;i<fPedestalStrings.size();i++){
      if(version == fPedestalStrings.at(i)) version_index = i;
    }//versions
    if(version_index == -1){
      version_index = fPedestalStrings.size();
      fPedestalStrings.push_back(version);
    }
    
    fPedestalInts.push_back(version_index);
    fRuns.push_back(run);    
    std::cerr << "INFO : " << "run " << run << " version " << version << std::endl;
  }

  TGraph* GetGraph(){

    fGr = new TGraph(fRuns.size());
    for(size_t i=0;i<fRuns.size();i++) fGr->SetPoint(i, fRuns.at(i), fPedestalInts.at(i));
    
    return fGr;
  }

  TCanvas* GetPedestalCanvas(TPaveText* updateText, int time_ago, int XNow, int width=1200, int height=800){
    TCanvas* can = new TCanvas("can_processing_version", "can_processing_version", width, height);
    can->SetRightMargin(0.20);
    TGraph* gr = GetGraph();
    //    can->cd()->SetLogy();
    gr->Draw("A*");
    gr->SetTitle("Nearline Processing Pedestal");   
    gr->GetXaxis()->SetTitle("Run Number");   

    //Set the bin labels and range correctly
    gr->GetYaxis()->Set(fPedestalStrings.size(),-0.5,fPedestalStrings.size()-0.5);
    gr->GetYaxis()->SetRangeUser(-0.5, fPedestalStrings.size()-0.5);
    for(size_t i=0;i<fPedestalStrings.size();i++){
      gr->GetYaxis()->SetBinLabel(i+1, fPedestalStrings.at(i).c_str());
    }

    gr->Draw("A*");

    gPad->SetGridx();
    gPad->SetGridy();

    updateText->Draw();
    return can;
  }
  
  static std::string GetPlotName(int Ndays){ 
    char name[256];
    sprintf(name, "ProcessingPedestal_%.3i_days.png", Ndays);
    return std::string(name);
  }

  
  void PrintPedestalPlots(std::string plot_dir, int Ndays, TPaveText* updateText, int time_ago, int XNow, int width=1200, int height=800){
    TCanvas *can = GetPedestalCanvas(updateText, time_ago, XNow, width, height);
    can->Print((plot_dir + "/" + GetPlotName(Ndays)).c_str());
    delete can;
    
    for(size_t i=0;i<fPedestalStrings.size();i++){
      std::cerr << "INFO : Pedestals Used: " << fPedestalStrings.at(i) << std::endl;
    }

  }

  static std::string MakePedestalPlotsHTML(std::string relative_plot_dir, int Ndays){

    std::string output;
    output+="<h3>Nearline Processing Pedestal as a function of Run.</h3>\n";
    output+="<figure>\n";
    output+="<img src=\"" + relative_plot_dir + "/" + GetPlotName(Ndays) + "\" width=\"800\">\n";
    output+="</figure>\n";
    output+="<BR><BR><BR>\n";

    return output;
  }
};

struct NearlineProcessingTime{
  
  TDatime fStartDate;
  TDatime fEndNearlineAnaDate;
  TDatime fEndNearlineMuonDate;
  TDatime fEndDate;

  NearlineProcessingTime(std::string file_name, int Ndays);
  static TDatime GetDateTime(std::string date, int Ndays);
  static std::string GetDoneFileName(std::string file_name);
  static std::string GetEVDDoneFileName(std::string file_name);
  static TDatime const InvalidDateTime;
};

const TDatime NearlineProcessingTime::InvalidDateTime = TDatime(1995,00,00,00,00,00);

NearlineProcessingTime::NearlineProcessingTime(std::string file_name, int Ndays){

  fStartDate=InvalidDateTime;
  fEndNearlineAnaDate=InvalidDateTime;
  fEndNearlineMuonDate=InvalidDateTime;
  fEndDate=InvalidDateTime;
  
  std::ifstream in_file(file_name.c_str());
  std::string line;

  std::string start = "";
  std::string end_nearline_ana = "";
  std::string end_nearline_muon = "";
  std::string end = "";

  while(std::getline(in_file, line)){   
    if(line.find("START_DATE ") != std::string::npos)  start = line.substr(std::string("START_DATE ").size());
    if(line.find("END_NEARLINE_ANA ") != std::string::npos)  end_nearline_ana = line.substr(std::string("END_NEARLINE_ANA ").size());
    if(line.find("END_NEARLINE_MUON ") != std::string::npos)  end_nearline_muon = line.substr(std::string("END_NEARLINE_MUON ").size());
    if(line.find("END_DATE ") != std::string::npos)  end = line.substr(std::string("END_DATE ").size());
  }

  in_file.close();

  if(start != "") fStartDate = GetDateTime(start, Ndays);
  if(end_nearline_ana != "") fEndNearlineAnaDate = GetDateTime(end_nearline_ana, Ndays);
  if(end_nearline_muon != "") fEndNearlineMuonDate = GetDateTime(end_nearline_muon, Ndays);
  if(end != "") fEndDate = GetDateTime(end, Ndays);
  
}


TDatime NearlineProcessingTime::GetDateTime(std::string date, int Ndays){
  std::string temp_file_name = "/tmp/NearlineDate_" + std::to_string(Ndays) + ".txt";
  std::string command = "date -d \"" + date + "\" +\"%Y-%m-%d %T\" > " + temp_file_name;
  int retval = system(command.c_str());
  std::ifstream temp_file(temp_file_name.c_str());
  std::string new_date;
  std::getline(temp_file, new_date);
  temp_file.close();
  TDatime date_time(new_date.c_str());
  return date_time;
}

std::string NearlineProcessingTime::GetDoneFileName(std::string file_name){

  ///lbne/data2/users/lbnedaq/nearline/v04_36_01/011/011425/lbne_r011425_sr01_20160217T083237_nearline_hists.root
  ///lbne/data2/users/lbnedaq/nearline/v04_36_01/011/011425/lbne_r011425_sr01_20160217T083237.root.DONE

  std::string file_name_done = file_name.substr(0, file_name.find("_nearline")) + ".root.DONE";
  return file_name_done;

}

std::string NearlineProcessingTime::GetEVDDoneFileName(std::string file_name){


  ///lbne/data2/users/lbnedaq/nearline/v04_36_01/011/011425/lbne_r011425_sr01_20160217T083237_nearline_hists.root
  ///lbne/data2/users/lbnedaq/nearline_evd/v04_36_01/011/011425_01/lbne_r011425_sr01_20160217T083237.rootEVD.DONE

  std::string prefix_dir = file_name.substr(0, file_name.find("/nearline/")); // = "/lbne/data2/users/lbnedaq"
  std::string suffix_dir = file_name.substr(file_name.find("/v"), std::string("/v04_36_01/011/011425").size());
  std::string short_file_name = file_name.substr(file_name.find(suffix_dir)+suffix_dir.size()+1);
  std::string file_name_done = prefix_dir + "/nearline_evd" + suffix_dir + "_01/" + short_file_name;
  file_name_done = file_name_done.substr(0, file_name_done.find("_nearline")) + ".rootEVD.DONE";

  return file_name_done;
}

struct NearlineProcessingTimePlot{
  
  std::vector<int> fVecRunTotal;
  std::vector<int> fVecRunNearlineAna;
  std::vector<int> fVecRunNearlineMuon;
  std::vector<int> fVecTotalTime;	
  std::vector<int> fVecNearlineAnaTime;	
  std::vector<int> fVecNearlineMuonTime;
  bool fIsEVD;
  //  TMultiGraph* fMultiGraph;
  TLegend fLegend;

  NearlineProcessingTimePlot(){
    fVecRunTotal.resize(0);		
    fVecRunNearlineAna.resize(0);		
    fVecRunNearlineMuon.resize(0);		
    fVecTotalTime.resize(0);	  
    fVecNearlineAnaTime.resize(0);	  
    fVecNearlineMuonTime.resize(0);
    fIsEVD=false;
  }

  void AddFile(std::string filename, int run, int Ndays){
    NearlineProcessingTime this_processing_time(filename, Ndays);
    std::cerr << "INFO : " << "run " << run << " StartDate: " << this_processing_time.fStartDate.AsString() << std::endl;
    std::cerr << "INFO : " << "run " << run << " EndDate: " << this_processing_time.fEndDate.AsString() << std::endl;
    if(run<=0) return;
    if(this_processing_time.fStartDate == NearlineProcessingTime::InvalidDateTime) return;
    if(this_processing_time.fEndDate == NearlineProcessingTime::InvalidDateTime) return;

    int total_time = this_processing_time.fEndDate.Get() - this_processing_time.fStartDate.Get();
    fVecRunTotal.push_back(run);
    fVecTotalTime.push_back(total_time);

    //    std::cerr << "ERROR: run " << run << " time " << total_time << std::endl;
    
    int nearline_ana_time = this_processing_time.fEndNearlineAnaDate.Get() - this_processing_time.fStartDate.Get();
    int nearline_muon_time = this_processing_time.fEndNearlineMuonDate.Get() - this_processing_time.fEndNearlineAnaDate.Get();

    if(this_processing_time.fEndNearlineAnaDate == NearlineProcessingTime::InvalidDateTime){ 
      nearline_ana_time = 0;
      nearline_muon_time = 0;
    }
    if(this_processing_time.fEndNearlineMuonDate == NearlineProcessingTime::InvalidDateTime){ 
      nearline_muon_time = 0;
    }
    if(nearline_ana_time!=0){
      fVecNearlineAnaTime.push_back(nearline_ana_time);
      fVecRunNearlineAna.push_back(run);
    }
    if(nearline_muon_time!=0){
      fVecNearlineMuonTime.push_back(nearline_ana_time);
      fVecRunNearlineMuon.push_back(run);
    }
  }
  
  TGraph* GetTotalTimeGraph(){
    TGraph* gr = new TGraph(fVecRunTotal.size());
    for(unsigned int i=0;i<fVecRunTotal.size();i++) gr->SetPoint(i, fVecRunTotal.at(i), fVecTotalTime.at(i));
    if(fIsEVD) gr->SetTitle("Nearline Event Display Processing Time");
    else gr->SetTitle("Nearline Processing Time");
    gr->SetMarkerColor(kBlue);
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(1.5);
    gr->GetXaxis()->SetTitle("Run Number");
    gr->GetYaxis()->SetTitle("Processing Time in Seconds");
    return gr;
  }
  TGraph* GetNearlineAnaTimeGraph(){
    TGraph* gr = new TGraph(fVecRunNearlineAna.size());
    for(unsigned int i=0;i<fVecRunNearlineAna.size();i++) gr->SetPoint(i, fVecRunNearlineAna.at(i), fVecNearlineAnaTime.at(i));
    gr->SetTitle("Nearline Ana Processing Time");
    gr->SetMarkerColor(kRed);
    gr->SetMarkerStyle(21);
    gr->SetMarkerSize(1.5);
    gr->GetXaxis()->SetTitle("Run Number");
    gr->GetYaxis()->SetTitle("Processing Time in Seconds");
    return gr;
  }
  TGraph* GetNearlineMuonTimeGraph(){
    TGraph* gr = new TGraph(fVecRunNearlineMuon.size());
    for(unsigned int i=0;i<fVecRunNearlineMuon.size();i++) gr->SetPoint(i, fVecRunNearlineMuon.at(i), fVecNearlineMuonTime.at(i));
    gr->SetTitle("Nearline Muon Processing Time");
    gr->SetMarkerColor(kBlack);
    gr->SetMarkerStyle(22);
    gr->SetMarkerSize(1.5);
    gr->GetXaxis()->SetTitle("Run Number");
    gr->GetYaxis()->SetTitle("Processing Time in Seconds");
    return gr;
  }
  TMultiGraph* GetTimeMultiGraph(){

    TMultiGraph* mgr = new TMultiGraph();
    TGraph* gr;
    if(fIsEVD) fLegend = TLegend(0.75,0.85,0.95,0.95);
    else fLegend = TLegend(0.75,0.75,0.95,0.95);

    gr = GetTotalTimeGraph();
    mgr->Add(gr, "P");
    fLegend.AddEntry(gr, gr->GetTitle(), "p");

    gr = GetNearlineAnaTimeGraph();
    if(gr->GetN()!=0){
      mgr->Add(gr, "P");
      fLegend.AddEntry(gr, gr->GetTitle(), "p");
    }
    gr = GetNearlineMuonTimeGraph();
    if(gr->GetN()!=0){
      mgr->Add(gr, "P");
      fLegend.AddEntry(gr, gr->GetTitle(), "p");
    }
    return mgr;
  }

  TCanvas* GetTimeCanvas(TPaveText* updateText, int time_ago, int XNow, int width=1200, int height=800){
    TCanvas* can = new TCanvas("can_processing_time", "can_processing_time", width, height);
    can->SetRightMargin(0.20);
    TMultiGraph* gr = GetTimeMultiGraph();
    can->cd()->SetLogy();
    gr->Draw("A*");
    if(fIsEVD) gr->SetTitle("Nearline Event Display Processing Time");   
    else gr->SetTitle("Nearline Processing Time");   
    gr->GetXaxis()->SetTitle("Run Number");   
    gr->GetYaxis()->SetTitle("Processing Time in Seconds");   

    gPad->SetGridx();
    gPad->SetGridy();

    updateText->Draw();
    fLegend.Draw();
    //    fMultiGraph = gr;
    return can;
  }
  
  std::string GetPlotName(int Ndays){ 
    char name[256];
    if(fIsEVD) sprintf(name, "ProcessingTimeEvd_%.3i_days.png", Ndays);
    else sprintf(name, "ProcessingTime_%.3i_days.png", Ndays);
    return std::string(name);
  }

  void PrintTimePlots(std::string plot_dir, int Ndays, TPaveText* updateText, int time_ago, int XNow, int width=1200, int height=800){
    TCanvas *can = GetTimeCanvas(updateText, time_ago, XNow, width, height);
    can->Print((plot_dir + "/" + GetPlotName(Ndays)).c_str());
    delete can;

    //    delete fMultiGraph;
  }

  std::string MakeTimePlotsHTML(std::string relative_plot_dir, int Ndays){

    std::string output;
    if(fIsEVD) output+="<h3>Nearline Event Display Processing Time as a function of Run.</h3>\n";
    else output+="<h3>Nearline Processing Time as a function of Run.</h3>\n";
    output+="<figure>\n";
    output+="<img src=\"" + relative_plot_dir + "/" + GetPlotName(Ndays) + "\" width=\"800\">\n";
    output+="</figure>\n";
    output+="<p>\n";
    if(fIsEVD)       output+="<b>Nearline Event Display Processing Time</b> - Total time taken to run this Nearline Processing Job<BR>\n";
    else{
      output+="<b>Nearline Processing Time</b> - Total time taken to run this Nearline Processing Job<BR>\n";
      output+="<b>Nearline Ana Processing Time</b> - Time taken to run the NearlineAna part of the Job<BR>\n";
      output+="<b>Nearline Muon Processing Time</b> - Time taken to run the Nearline Muon Counter part of the Job<BR>\n";
    }
    output+="</p>\n";
    output+="<BR><BR><BR>\n";

    return output;
  }


};


#endif
