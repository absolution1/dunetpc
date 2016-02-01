#include "NearlinePlotMaker.h"
#include "TMath.h"

const std::string PLOT_DIR = "/web/sites/lbne-dqm.fnal.gov/htdocs/NearlineMonitoring/plots";
const std::string PLOT_DIR_DEBUG = "/web/sites/lbne-dqm.fnal.gov/htdocs/NearlineMonitoring/plots_testing";

void histogramZoom(TH1* hist, double n_sigma){
  
  //Get the mean and RMS of the y-axis
  double mean = hist->GetMean(2);
  double rms = hist->GetRMS(2);

  double y_low = mean - rms*n_sigma;
  double y_high = mean + rms*n_sigma;
  
  if(y_low<0) y_low=0;

  TAxis* y_axis = hist->GetYaxis();
  y_axis->SetRangeUser(y_low, y_high);
      
}

void graphZoom(TGraph* gr, double n_sigma){
  
  //We need to make sure that the TGraph has already had draw called
  //Remember where we were cd'd into and go back there at the end
  TVirtualPad* old_pad = gPad;
  
  TCanvas* can_temp = new TCanvas("can_temp", "can_temp");
  gr->Draw("AL");
  
  //Get the mean and RMS excluding points that are less than 1 on the y-axis
  double mean = 0;
  double rms = 0;
  int num_entries = 0;
  double* y_values = gr->GetY();
  for(int i=0;i<gr->GetN();i++){
    if(y_values[i] < 1) continue;
    rms += TMath::Power(y_values[i],2);
    mean += y_values[i];
    num_entries++;
  }
  mean /= num_entries;
  rms = TMath::Abs(rms / num_entries - mean*mean);
  rms = TMath::Sqrt(rms);

  double y_low = mean - rms*n_sigma;
  double y_high = mean + rms*n_sigma;
  
  if(y_low<0) y_low=0;

  TAxis* y_axis = gr->GetYaxis();
  y_axis->SetRangeUser(y_low, y_high);
  
  //Return to whence we came
  old_pad->cd();
  delete can_temp;

}


//
// Make the plots for the nearline webpage from the nearline output files.
//

Long64_t NearlinePlotMaker(int Ndays, bool debug=false){

  std::cout << "\n\nMaking 35t Nearline plots for " << Ndays << " days...\n\n";
  


  // Define varibales about time...
  TNowGMT = new TDatime;  // current GMT time
  TNow    = new TDatime;  // current local time
  SRtime  = new TDatime;
  int XNowGMT  = TNowGMT->Convert(kTRUE);
  int XNow     = TNow   ->Convert(kFALSE);
  GMToffset    = XNowGMT - XNow;
  XNow         = XNow - GMToffset;
  int time_ago = XNow - 60*60*24*Ndays;
  int Xsrtime  = 0;
  unsigned int year  = 0;
  unsigned int month = 0;
  unsigned int day   = 0;
  double Hour = 0.0;
  double Min  = 0.0;
  int hour = 0;
  int min  = 0;
  int sec  = 0;

  run     = 0;
  subrun  = 0;
  LastRun = 0;
  LastSR  = 0;

  unsigned int yearEnd  = 0;
  unsigned int monthEnd = 0;
  unsigned int dayEnd   = 0;
  double       HourEnd  = 0.0;

  // Root is stupid and you have to make it think that it is in a different time zone so that it will draw any
  // plots with time on the X-axis correctly.
  gSystem->Setenv("TZ","UTC");
  gStyle->SetTimeOffset(0);



  // Define the maximum number of points to go into a TGraph (assuming a max
  // of one subrun every minute.)
  const unsigned int Npoint = 1440 * Ndays;



  //
  // Make the plots!!!
  //



  // Open list of input files...
  char filelist_title[128];
  sprintf(filelist_title,"/home/lbnedaq/nearline/temp/35t_%.dDay_Nearline_File_List.txt",Ndays);
  std::cout << "\n\nOpening list of input files:\n" << filelist_title << "\n\n";
  inFile.open(filelist_title);


  std::vector<NearlinePlot*> NearlinePlotVec;
  std::vector<int> channelVec;
  
  channelVec.push_back(20);
  channelVec.push_back(548);
  channelVec.push_back(1297);
  channelVec.push_back(1697);
  channelVec.push_back(1838);
  channelVec.push_back(1482);
  channelVec.push_back(952);
  channelVec.push_back(454);
  for(size_t index=0;index<channelVec.size();index++){
    char hist_name[256];
    char hist_title[256];
    char hist_output_name[256];
    char graph_output_name[256];
    char graph_rms_output_name[256];
    int channel = channelVec.at(index);
    char can_metric_name[256];
    char can_metric_rms_name[256];

    if(debug) sprintf(hist_output_name,"%s/ADCSpecChan%04i_%.3u_days.png", PLOT_DIR_DEBUG.c_str(), channel, Ndays);
    else sprintf(hist_output_name,"%s/ADCSpecChan%04i_%.3u_days.png", PLOT_DIR.c_str(), channel, Ndays);

    if(debug) sprintf(graph_output_name,"%s/MeanPedChan%04i_%.3u_days", PLOT_DIR_DEBUG.c_str(), channel, Ndays);
    else sprintf(graph_output_name,"%s/MeanPedChan%04i_%.3u_days", PLOT_DIR.c_str(), channel, Ndays);

    if(debug) sprintf(graph_rms_output_name,"%s/RMSPedChan%04i_%.3u_days", PLOT_DIR_DEBUG.c_str(), channel, Ndays);
    else sprintf(graph_rms_output_name,"%s/RMSPedChan%04i_%.3u_days", PLOT_DIR.c_str(), channel, Ndays);
    
    sprintf(hist_name, "hped_per_tick_chan_%i", channel);
    sprintf(hist_title, "ADC per Tick - Channel %i", channel);
    sprintf(can_metric_name, "Mean Pedestal per Event for Channel %i", channel);
    sprintf(can_metric_rms_name, "Pedestal RMS for Channel %i", channel);

    NearlinePlot* this_plot = new NearlinePlot(Npoint, hist_name, hist_output_name, graph_output_name, graph_rms_output_name, can_metric_name, can_metric_rms_name, hist_title, 128, 0, 4096);

    NearlinePlotVec.push_back(this_plot);
    
  }




  // These variables are meant to be an expression of general DAQ health
  float *RunVSYearYear  = new float[Npoint];
  float *RunVSYearRun   = new float[Npoint];
  int    RunVSYearCount = 0;



  // initalize histos



  //
  // Looping over the list of input files...
  //
  while(inFile.good()) {
    
    // Open the Nth file...
    char filename[512];
    inFile >> filename;
    if(!inFile.good()) continue; // prevent code from running over the last file twice...
    TFile file(filename);
    file.cd("nearlineana");

    // Reset the time variables...
    Xsrtime = 0;
    year    = 0;
    month   = 0;
    day     = 0;
    Hour    = 0.0;
    Min     = 0.0;
    hour    = 0;
    min     = 0;
    sec     = 0;
    
    TTree *header  = (TTree *)file.FindObjectAny("Header");
    if(header != 0) {
      header->SetBranchAddress("StartYear",  &year);
      header->SetBranchAddress("StartMonth", &month);
      header->SetBranchAddress("StartDay",   &day);
      header->SetBranchAddress("StartHour",  &Hour);
      header->SetBranchAddress("EndYear",    &yearEnd);
      header->SetBranchAddress("EndMonth",   &monthEnd);
      header->SetBranchAddress("EndDay",     &dayEnd);
      header->SetBranchAddress("EndHour",    &HourEnd);
      header->SetBranchAddress("Run",        &run);
      header->SetBranchAddress("Subrun",     &subrun);
      header->GetEntry(0);
      
      hour = Hour;
      Min  = (Hour-hour)*60.0;
      min  = (Hour-hour)*60.0;
      sec  = (Min-min)*60.0;

      // Ana/DAQ keeps track of GMT time, so convert to local time
      SRtime->Set(year,month,day,hour,min,sec);
      Xsrtime = SRtime->Convert() - GMToffset;
      cout << "time:\t" << year << " " << month << " " << day << " " << Hour << endl;
      if(run == 0) std::cout << "\nWARNING: Run number zero for file:\t" << filename << "\n\n" << endl;
      hour = HourEnd;
      Min  = (HourEnd-hour)*60.0;
      min  = (HourEnd-hour)*60.0;
      sec  = (Min-min)*60.0;
      
      RunVSYearYear[RunVSYearCount] = year;
      RunVSYearRun [RunVSYearCount] = run;
      RunVSYearCount++;            

      if(Xsrtime < time_ago) continue; 
      
      if(run >= LastRun) {
	if(run > LastRun) {
	  LastRun = run;
	  LastSR  = subrun;
	}
	else if(subrun > LastSR) {
	  LastRun = run;
	  LastSR  = subrun;
	}
      }
    
    } // end if header != 0


    for(size_t index=0;index<NearlinePlotVec.size(); index++){
      NearlinePlot* this_plot = NearlinePlotVec.at(index);
      //  bool AddHistogram(TFile const & file, TTree* header, int Xstrtime, int Xsrtime, int XNow, int GMToffset);
      this_plot->AddHistogram(file, header, Xsrtime, XNow, GMToffset);

    }//loop over plots

    
    

    file.Close();

  } // end while loop over input files



  //
  // Make time/date stamp for each plot
  //
  TDatime *Ttemp = new TDatime;  // finish time
  int Xfin = Ttemp->Convert() - GMToffset;
  TDatime *Tfinish = new TDatime(Xfin);  // finish time



  // Make conditional standard time axis labels depending on time period.
  string taxis_labels;
  if(Ndays <= 2){
    taxis_labels = "%H:%M";
  }
  else{
    taxis_labels = "%m/%d";
  }




  //
  // Make the plots...
  //

  // Make the update text info.
  TPaveText *UpdateText = new TPaveText(0.1, 0.0, 0.5, 0.05, "NDC");
  UpdateText->SetLineColor(0);
  UpdateText->SetFillColor(0);
  UpdateText->SetBorderSize(1);
  UpdateText->SetMargin(0.0);
  UpdateText->SetTextAlign(11);
  char buff1[256];
  char buff2[256];
  sprintf(buff1, "Last updated on:    %s (central time)", Tfinish->AsString());
  UpdateText->AddText(buff1);
  sprintf(buff2, "Last run / subrun:   %i / %i", LastRun, LastSR);
  UpdateText->AddText(buff2);

    






  for(size_t index=0;index<NearlinePlotVec.size();index++){
    NearlinePlot* this_plot = NearlinePlotVec.at(index);
    //  TCanvas* makeHistoCanvas(std::string can_name, std::string can_title, int width, int height, TPaveText* updateText);

    TCanvas *canHist = this_plot->makeHistoCanvas("can_histo", "can_histo", 1200, 800, UpdateText);
    canHist->Print(this_plot->fHistOutputName.c_str());
    delete canHist;
    //  TCanvas* makeGraphMetricTimeCanvas(std::string can_name, std::string can_title, int width, int height, TPaveText* updateText, std::string taxis_labels, int time_ago, int XNow);
    

    char output_name[256];
    bool zoom;
    TCanvas *canMetricTime = this_plot->makeGraphMetricTimeCanvas(1200, 800, UpdateText, taxis_labels, time_ago, XNow);
    sprintf(output_name, "%s.png", this_plot->fGraphOutputName.c_str());
    canMetricTime->Print(output_name);
    delete canMetricTime;

    zoom=true;
    canMetricTime = this_plot->makeGraphMetricTimeCanvas(1200, 800, UpdateText, taxis_labels, time_ago, XNow, zoom);
    sprintf(output_name, "%s_zoom.png", this_plot->fGraphOutputName.c_str());
    canMetricTime->Print(output_name);      
    delete canMetricTime;

    TCanvas *canMetricRmsTime = this_plot->makeGraphMetricRmsTimeCanvas(1200, 800, UpdateText, taxis_labels, time_ago, XNow);    
    sprintf(output_name, "%s.png", this_plot->fGraphRmsOutputName.c_str());
    canMetricRmsTime->Print(output_name);
    delete canMetricRmsTime;

    zoom=true;
    canMetricRmsTime = this_plot->makeGraphMetricRmsTimeCanvas(1200, 800, UpdateText, taxis_labels, time_ago, XNow, zoom);
    sprintf(output_name, "%s_zoom.png", this_plot->fGraphRmsOutputName.c_str());
    canMetricRmsTime->Print(output_name);
    delete canMetricRmsTime;

  }//loop over plots



  //
  // Draw pretty canvases...
  //
  int maxtime = 0;
  double max  = 0.0, ave = 0.0;
  TPaveText *LastPoint = new TPaveText(0.3,0.88,0.93,0.93,"NDC");
  LastPoint->SetLineColor(1);
  LastPoint->SetFillColor(0);
  LastPoint->SetBorderSize(1);
  char lptext[128];


  char filename[128];
  char title[128];


  // A plot of general DAQ health...
  int nzero = 0;
  int firstrun = 1e9;
  int lastrun  = 0;
  for(int i = 0; i < RunVSYearCount; ++i) {
    if(RunVSYearRun[i] == 0) nzero++;
    if(RunVSYearRun[i] > lastrun  && RunVSYearRun[i] != 0) lastrun  = RunVSYearRun[i];
    if(RunVSYearRun[i] < firstrun && RunVSYearRun[i] != 0) firstrun = RunVSYearRun[i];
  }

  sprintf(lptext,"Number of points with run number = 0 : %u",nzero);
  LastPoint->Clear();
  LastPoint->AddText(lptext);


  TCanvas *cRunVSYear = new TCanvas("cRunVSYear","Year VS Run",1200,800);
  cRunVSYear->cd();
  gPad->SetGridx();
  TGraph *gRunVSYear = new TGraph(RunVSYearCount,RunVSYearRun,RunVSYearYear);
  sprintf(title,"Year VS Run");
  gRunVSYear->SetTitle(title);
  gRunVSYear->SetMarkerColor(kBlue);
  gRunVSYear->GetXaxis()->SetTitle("Run Number");
  gRunVSYear->GetYaxis()->SetTitle("Year");
  gRunVSYear->GetXaxis()->SetLimits(firstrun,lastrun);
  gRunVSYear->Draw("A*");
  UpdateText->Draw();
  LastPoint->Draw();
  if(debug)  sprintf(filename,"%s/RunVSYear_%.3u_days.png",PLOT_DIR_DEBUG.c_str(),Ndays);
  else  sprintf(filename,"%s/RunVSYear_%.3u_days.png",PLOT_DIR.c_str(),Ndays);
  cRunVSYear->Print(filename);



  // Done with everything...
  return 0;
  
} // End script.
