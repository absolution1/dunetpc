#include "NearlinePlotMaker.h"
#include "TMath.h"

const std::string PLOT_DIR_PRODUCTION = "/web/sites/lbne-dqm.fnal.gov/htdocs/NearlineMonitoring/plots";
const std::string PLOT_DIR_DEBUG = "/web/sites/lbne-dqm.fnal.gov/htdocs/NearlineMonitoring/plots_testing";
//const std::string PLOT_DIR_DEBUG = "/lbne/data2/users/lbnedaq/nearline_webpage_test/plots_testing";

std::string PLOT_DIR;

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

  if( num_entries == 0){
    old_pad->cd();
    delete can_temp;
    return;
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

Long64_t NearlinePlotMaker(int Ndays, bool debug=false);

Long64_t NearlinePlotMakerDev(int Ndays, bool debug=false){

  NearlinePlotMaker(Ndays, debug);

}

Long64_t NearlinePlotMaker(int Ndays, bool debug){

  std::cout << "\n\nMaking 35t Nearline plots for " << Ndays << " days...\n\n";
  if(debug) PLOT_DIR=PLOT_DIR_DEBUG;
  else PLOT_DIR=PLOT_DIR_PRODUCTION;


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
  //  sprintf(filelist_title,"/home/lbnedaq/nearline/temp/35t_%.dDay_Nearline_File_List_test.txt",Ndays);
  //  sprintf(filelist_title,"/lbne/app/users/jpdavies/dunetpc-nearline/srcs/dunetpc/35t_%.dDay_Nearline_File_List.txt",Ndays);
  std::cout << "\n\nOpening list of input files:\n" << filelist_title << "\n\n";
  inFile.open(filelist_title);
  

  std::vector<NearlinePlot*> NearlinePlotVec;
  
  {
    int channel;
    NearlinePlotInfo this_plot_info;
    char hist_name [256];
    NearlinePlot* this_plot;
    std::string metric_details;

    channel=20;
    metric_details = "(APA-3-U-plane)";
    this_plot_info = NearlinePlotInfo("ADC", channel, Ndays, "png", metric_details);
    sprintf(hist_name, "hped_per_tick_chan_%i", channel);
    this_plot = new NearlinePlot(hist_name, this_plot_info);
    NearlinePlotVec.push_back(this_plot);

    channel=548;
    metric_details = "(APA-2-Z-plane)";
    this_plot_info = NearlinePlotInfo("ADC", channel, Ndays, "png", metric_details);
    sprintf(hist_name, "hped_per_tick_chan_%i", channel);
    this_plot = new NearlinePlot(hist_name, this_plot_info);
    NearlinePlotVec.push_back(this_plot);

    channel=1297;
    metric_details = "(APA-1-V-plane)";
    this_plot_info = NearlinePlotInfo("ADC", channel, Ndays, "png", metric_details);
    sprintf(hist_name, "hped_per_tick_chan_%i", channel);
    this_plot = new NearlinePlot(hist_name, this_plot_info);
    NearlinePlotVec.push_back(this_plot);

    channel=1697;
    metric_details = "(APA-0-Z-plane)";
    this_plot_info = NearlinePlotInfo("ADC", channel, Ndays, "png", metric_details);
    sprintf(hist_name, "hped_per_tick_chan_%i", channel);
    this_plot = new NearlinePlot(hist_name, this_plot_info);
    NearlinePlotVec.push_back(this_plot);

    channel=1838;
    metric_details = "(APA-3-Z-plane)";
    this_plot_info = NearlinePlotInfo("ADC", channel, Ndays, "png", metric_details);
    sprintf(hist_name, "hped_per_tick_chan_%i", channel);
    this_plot = new NearlinePlot(hist_name, this_plot_info);
    NearlinePlotVec.push_back(this_plot);

    channel=1482;
    metric_details = "(APA-2-U-plane)";
    this_plot_info = NearlinePlotInfo("ADC", channel, Ndays, "png", metric_details);
    sprintf(hist_name, "hped_per_tick_chan_%i", channel);
    this_plot = new NearlinePlot(hist_name, this_plot_info);
    NearlinePlotVec.push_back(this_plot);

    channel=952;
    metric_details = "(APA-1-Z-plane)";
    this_plot_info = NearlinePlotInfo("ADC", channel, Ndays, "png", metric_details);
    sprintf(hist_name, "hped_per_tick_chan_%i", channel);
    this_plot = new NearlinePlot(hist_name, this_plot_info);
    NearlinePlotVec.push_back(this_plot);

    channel=454;
    metric_details = "(APA-0-U-plane)";
    this_plot_info = NearlinePlotInfo("ADC", channel, Ndays, "png", metric_details);
    sprintf(hist_name, "hped_per_tick_chan_%i", channel);
    this_plot = new NearlinePlot(hist_name, this_plot_info);
    NearlinePlotVec.push_back(this_plot);
  }


  {
    NearlinePlot *this_plot;
    NearlinePlotInfo this_plot_info;
    int channel = -1;
    bool normalise_histo1d=true;
    bool make_average_rms_time_plots = false;
    bool make_2d_histos = true;
    bool make_bin_by_bin_plots = true;
    NearlinePlotEnables this_plot_enables = NearlinePlotEnables(normalise_histo1d, make_average_rms_time_plots, make_2d_histos, make_bin_by_bin_plots);

    bool metric_time_graph_log=false;
    bool histo_1D_log=false;
    bool histo_2D_log=false;;
    bool bin_by_bin_log=false;

    NearlinePlotLogScale this_plot_log_scale = NearlinePlotLogScale(metric_time_graph_log, histo_1D_log, histo_2D_log, bin_by_bin_log);
    
    this_plot_info = NearlinePlotInfo("TSU Counter Rates", channel, Ndays, "png");
    this_plot = new NearlinePlot("TSUs", this_plot_info, this_plot_enables, this_plot_log_scale);
    NearlinePlotVec.push_back(this_plot);

    this_plot_info = NearlinePlotInfo("BSU Counter Rates", channel, Ndays, "png");
    this_plot = new NearlinePlot("BSUs", this_plot_info, this_plot_enables, this_plot_log_scale);
    NearlinePlotVec.push_back(this_plot);

    this_plot_info = NearlinePlotInfo("Counter Trigger Rates", channel, Ndays, "png");
    this_plot = new NearlinePlot("Triggers", this_plot_info, this_plot_enables, this_plot_log_scale);
    NearlinePlotVec.push_back(this_plot);
    
  }

  // for(int index=0;index<16;index++){
  //   int channel = index*128;
  //   NearlinePlotInfo this_plot_info("Hits", channel, Ndays, "png");
  //   char hist_name[256];
  //   sprintf(hist_name, "hhits_per_event_chan_%i", channel);
  //   NearlinePlot *this_plot = new NearlinePlot(hist_name, this_plot_info);
  //   NearlinePlotVec.push_back(this_plot);
  // }

  for(size_t index=0;index<NearlinePlotVec.size();index++){
    std::cerr << "NearlinePlot: " << (NearlinePlotVec.at(index))->fHistName << std::endl;
  }

  // These variables are meant to be an expression of general DAQ health
  float *RunVSYearYear  = new float[Npoint];
  float *RunVSYearRun   = new float[Npoint];
  int    RunVSYearCount = 0;

  NearlineProcessingTimePlot nearline_processing_time_plot;
  NearlineProcessingTimePlot nearline_evd_processing_time_plot;
  nearline_evd_processing_time_plot.fIsEVD=true;
  NearlineProcessingVersion nearline_processing_version_plot;
  NearlineProcessingPedestal nearline_processing_pedestal_plot;
  //
  // Looping over the list of input files...
  //
  while(inFile.good()) {
    
    // Open the Nth file...
    char filename[512];
    inFile >> filename;
    
    std::cerr << "INFO : " << filename << std::endl;

    if(!inFile.good()) continue; // prevent code from running over the last file twice...
    TFile file(filename);
    //    file.cd("nearlineana");

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
      cout << "INFO : time:\t" << year << " " << month << " " << day << " " << Hour << endl;
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
      this_plot->AddHistogram(file, header, Xsrtime, XNow, GMToffset, time_ago);
      //      this_plot->AddHistogram2D(file, header, Xsrtime, XNow, GMToffset, time_ago);

    }//loop over plots
    
    std::string done_file_name = NearlineProcessingTime::GetDoneFileName(filename);
    std::string evd_done_file_name = NearlineProcessingTime::GetEVDDoneFileName(filename);
    nearline_processing_time_plot.AddFile(done_file_name, run, Ndays);
    nearline_evd_processing_time_plot.AddFile(evd_done_file_name, run, Ndays);
    nearline_processing_version_plot.AddFile(filename, run);
    nearline_processing_pedestal_plot.AddFile(done_file_name, run);

    //    std::cerr << "JPD: done_file_name = " << done_file_name << std::endl;

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

    
  std::string html_string = NearlineHTML::MakePageHeader(Ndays);;
  std::string html_string_bin_by_bin = NearlineHTML::MakePageHeader(Ndays);


  for(size_t index=0;index<NearlinePlotVec.size();index++){
    NearlinePlot* this_plot = NearlinePlotVec.at(index);
    this_plot->printPlots(PLOT_DIR, UpdateText, time_ago, XNow);

    //Make HTML
    std::string plot_location;
    if(debug) plot_location = "plots_testing/";
    else plot_location = "plots/";

    //    std::cerr << NearlineHTML::MakePlotSet(plot_location, this_plot->fPlotInfo);

    html_string += NearlineHTML::MakePlotSet(plot_location, this_plot);

    if(this_plot->fPlotEnables.fMakeBinByBinPlots) html_string_bin_by_bin += NearlineHTML::MakeBinByBinGraphs(plot_location, this_plot);

  }

  //Print the Nearline Processing Time plot and add to the webpage
  nearline_processing_time_plot.PrintTimePlots(PLOT_DIR,  Ndays,  UpdateText, time_ago, XNow);
  if(debug) html_string += nearline_processing_time_plot.MakeTimePlotsHTML("plots_testing/", Ndays);
  else  html_string += nearline_processing_time_plot.MakeTimePlotsHTML("plots/", Ndays);

  //Print the Nearline Processing Event Display Time plot and add to the webpage
  nearline_evd_processing_time_plot.PrintTimePlots(PLOT_DIR,  Ndays,  UpdateText, time_ago, XNow);
  if(debug) html_string += nearline_evd_processing_time_plot.MakeTimePlotsHTML("plots_testing/", Ndays);
  else  html_string += nearline_evd_processing_time_plot.MakeTimePlotsHTML("plots/", Ndays);

  //Print the Nearline Processing Version plot and add to the webpage
  nearline_processing_version_plot.PrintVersionPlots(PLOT_DIR,  Ndays,  UpdateText, time_ago, XNow);
  if(debug) html_string += nearline_processing_version_plot.MakeVersionPlotsHTML("plots_testing/", Ndays);
  else  html_string += nearline_processing_version_plot.MakeVersionPlotsHTML("plots/", Ndays);

  //Print the Nearline Processing Pedestal plot and add to the webpage
  nearline_processing_pedestal_plot.PrintPedestalPlots(PLOT_DIR,  Ndays,  UpdateText, time_ago, XNow);
  if(debug) html_string += nearline_processing_pedestal_plot.MakePedestalPlotsHTML("plots_testing/", Ndays);
  else  html_string += nearline_processing_pedestal_plot.MakePedestalPlotsHTML("plots/", Ndays);



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

  sprintf(filename,"%s/RunVSYear_%.3u_days.png",PLOT_DIR.c_str(),Ndays);
  cRunVSYear->Print(filename);


  //Add this plot to the HTML page
  if(debug) sprintf(filename,"%s/RunVSYear_%.3u_days.png","plots_testing/",Ndays);
  else sprintf(filename,"%s/RunVSYear_%.3u_days.png","plots/",Ndays);
  html_string += NearlineHTML::MakeStartYearPlot(filename);
  
  //Now create the html pages
  std::string html_file_name;
  html_file_name = PLOT_DIR + "/../";

  if(Ndays <= 2) html_file_name += "indexDay";
  else if(Ndays==7) html_file_name += "indexWeek";
  else if(Ndays==31) html_file_name += "indexMonth";
  else html_file_name += "indexUnknown";

  if(debug) html_file_name += "_test";
  html_file_name += ".html";

  std::cerr << "INFO: Creating HTML file \"" << html_file_name << "\"\n";

  std::ofstream html_file;;
  html_file.open (html_file_name);

  html_file << html_string;
  
  html_file.close();

  std::cerr << "INFO: Done Creating HTML file \"" << html_file_name << "\"\n";


  //Create BinByBin html pages
  html_file_name = PLOT_DIR + "/../";

  if(Ndays <= 2) html_file_name += "indexDay";
  else if(Ndays==7) html_file_name += "indexWeek";
  else if(Ndays==31) html_file_name += "indexMonth";
  else html_file_name += "indexUnknown";

  html_file_name += "BinByBin";

  if(debug) html_file_name += "_test";
  html_file_name += ".html";

  std::cerr << "INFO: Creating HTML file \"" << html_file_name << "\"\n";

  html_file.open (html_file_name);

  html_file << html_string_bin_by_bin;
  
  html_file.close();

  std::cerr << "INFO: Done Creating HTML file \"" << html_file_name << "\"\n";

  // Done with everything...
  return 0;
  
} // End script.
