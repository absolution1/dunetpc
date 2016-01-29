#include "NearlinePlotMaker.h"

const std::string PLOT_DIR = "/web/sites/lbne-dqm.fnal.gov/htdocs/NearlineMonitoring/plots";
const std::string PLOT_DIR_DEBUG = "/web/sites/lbne-dqm.fnal.gov/htdocs/NearlineMonitoring/plots_testing";

void histogramZoom(TH1* hist, double n_sigma);
void graphZoom(TGraph* gr, double n_sigma);

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
  
  //Get the mean and RMS of the y-axis
  double mean = gr->GetMean(2);
  double rms = gr->GetRMS(2);

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



  // NOTE: Much of the code below is not written in a smart, flexible way. For future users, I'd advise
  //       makeing this a bit more general and less "copy and paste" style...

  // Book histos...
  TH1F *hped_per_tick_chan_20   = new TH1F("hped_per_tick_chan_20",  "ADC per Tick - Channel 20;ADC"  ,128,0,4096);
  TH1F *hped_per_tick_chan_548  = new TH1F("hped_per_tick_chan_548", "ADC per Tick - Channel 548;ADC" ,128,0,4096);
  TH1F *hped_per_tick_chan_1297 = new TH1F("hped_per_tick_chan_1297","ADC per Tick - Channel 1297;ADC",128,0,4096);
  TH1F *hped_per_tick_chan_1697 = new TH1F("hped_per_tick_chan_1697","ADC per Tick - Channel 1697;ADC",128,0,4096);

  TH1F *hped_per_tick_chan_1838 = new TH1F("hped_per_tick_chan_1838","ADC per Tick - Channel 1838;ADC",128,0,4096);
  TH1F *hped_per_tick_chan_1482 = new TH1F("hped_per_tick_chan_1482","ADC per Tick - Channel 1482;ADC",128,0,4096);
  TH1F *hped_per_tick_chan_952  = new TH1F("hped_per_tick_chan_952", "ADC per Tick - Channel 952;ADC", 128,0,4096);
  TH1F *hped_per_tick_chan_454  = new TH1F("hped_per_tick_chan_454", "ADC per Tick - Channel 454;ADC", 128,0,4096);


  
  // Declare variables used to make the TGraphs...
  float *MeanPedChan0020      = new float[Npoint];
  float *MeanPedChan0020Time  = new float[Npoint];
  int    MeanPedChan0020Count = 0;

  float *MeanPedChan0548      = new float[Npoint];
  float *MeanPedChan0548Time  = new float[Npoint];
  int    MeanPedChan0548Count = 0;

  float *MeanPedChan1297      = new float[Npoint];
  float *MeanPedChan1297Time  = new float[Npoint];
  int    MeanPedChan1297Count = 0;

  float *MeanPedChan1697      = new float[Npoint];
  float *MeanPedChan1697Time  = new float[Npoint];
  int    MeanPedChan1697Count = 0;

  float *MeanPedChan1838      = new float[Npoint];
  float *MeanPedChan1838Time  = new float[Npoint];
  int    MeanPedChan1838Count = 0;

  float *MeanPedChan1482      = new float[Npoint];
  float *MeanPedChan1482Time  = new float[Npoint];
  int    MeanPedChan1482Count = 0;

  float *MeanPedChan0952      = new float[Npoint];
  float *MeanPedChan0952Time  = new float[Npoint];
  int    MeanPedChan0952Count = 0;

  float *MeanPedChan0454      = new float[Npoint];
  float *MeanPedChan0454Time  = new float[Npoint];
  int    MeanPedChan0454Count = 0;



  float *RMSPedChan0020      = new float[Npoint];
  float *RMSPedChan0020Time  = new float[Npoint];
  int    RMSPedChan0020Count = 0;

  float *RMSPedChan0548      = new float[Npoint];
  float *RMSPedChan0548Time  = new float[Npoint];
  int    RMSPedChan0548Count = 0;

  float *RMSPedChan1297      = new float[Npoint];
  float *RMSPedChan1297Time  = new float[Npoint];
  int    RMSPedChan1297Count = 0;

  float *RMSPedChan1697      = new float[Npoint];
  float *RMSPedChan1697Time  = new float[Npoint];
  int    RMSPedChan1697Count = 0;

  float *RMSPedChan1838      = new float[Npoint];
  float *RMSPedChan1838Time  = new float[Npoint];
  int    RMSPedChan1838Count = 0;

  float *RMSPedChan1482      = new float[Npoint];
  float *RMSPedChan1482Time  = new float[Npoint];
  int    RMSPedChan1482Count = 0;

  float *RMSPedChan0952      = new float[Npoint];
  float *RMSPedChan0952Time  = new float[Npoint];
  int    RMSPedChan0952Count = 0;

  float *RMSPedChan0454      = new float[Npoint];
  float *RMSPedChan0454Time  = new float[Npoint];
  int    RMSPedChan0454Count = 0;



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

    
    
    // Lift out the histos that we want...

    //
    // NOTE: In the if statements below, there occurs this: " && Xsrtime != XNow - GMToffset"
    //       I don't know why, but somehow a bogus point was being included at the end of the
    //       TGraphs with a value of zero and a time exactly equal to (now - 6 hours). This
    //       likely comes somehow from looping over the last file twice or from somehow getting
    //       empty information and defaulting to these values, but I couldn't track it down in
    //       timely fashion. So it will remain a problem for the next person to find. For now,
    //       this if statement catches it. If this is confusing, e-mail mbaird42@FNAL.GOV with
    //       questions.)
    TH1F *hped_per_tick_chan_20_temp = (TH1F*)file.FindObjectAny("hped_per_tick_chan_20");    
    if(hped_per_tick_chan_20_temp != 0){
      hped_per_tick_chan_20->Add(hped_per_tick_chan_20_temp,1.0);
    }    
    if(hped_per_tick_chan_20_temp != 0 && header != 0 && Xsrtime != XNow - GMToffset) {
      MeanPedChan0020Time[MeanPedChan0020Count] = Xsrtime;
      MeanPedChan0020    [MeanPedChan0020Count] = hped_per_tick_chan_20_temp->GetMean(1);
      MeanPedChan0020Count++;      

      RMSPedChan0020Time[RMSPedChan0020Count] = Xsrtime;
      RMSPedChan0020    [RMSPedChan0020Count] = hped_per_tick_chan_20_temp->GetRMS(1);
      RMSPedChan0020Count++;      
    }
    
    TH1F *hped_per_tick_chan_548_temp = (TH1F*)file.FindObjectAny("hped_per_tick_chan_548");    
    if(hped_per_tick_chan_548_temp != 0){
      hped_per_tick_chan_548->Add(hped_per_tick_chan_548_temp,1.0);
    }    
    if(hped_per_tick_chan_548_temp != 0 && header != 0 && Xsrtime != XNow - GMToffset) {
      MeanPedChan0548Time[MeanPedChan0548Count] = Xsrtime;
      MeanPedChan0548    [MeanPedChan0548Count] = hped_per_tick_chan_548_temp->GetMean(1);
      MeanPedChan0548Count++;      

      RMSPedChan0548Time[RMSPedChan0548Count] = Xsrtime;
      RMSPedChan0548    [RMSPedChan0548Count] = hped_per_tick_chan_548_temp->GetRMS(1);
      RMSPedChan0548Count++;      
    }
    
    TH1F *hped_per_tick_chan_1297_temp = (TH1F*)file.FindObjectAny("hped_per_tick_chan_1297");    
    if(hped_per_tick_chan_1297_temp != 0){
      hped_per_tick_chan_1297->Add(hped_per_tick_chan_1297_temp,1.0);
    }    
    if(hped_per_tick_chan_1297_temp != 0 && header != 0 && Xsrtime != XNow - GMToffset) {
      MeanPedChan1297Time[MeanPedChan1297Count] = Xsrtime;
      MeanPedChan1297    [MeanPedChan1297Count] = hped_per_tick_chan_1297_temp->GetMean(1);
      MeanPedChan1297Count++;      

      RMSPedChan1297Time[RMSPedChan1297Count] = Xsrtime;
      RMSPedChan1297    [RMSPedChan1297Count] = hped_per_tick_chan_1297_temp->GetRMS(1);
      RMSPedChan1297Count++;      
    }
    
    TH1F *hped_per_tick_chan_1697_temp = (TH1F*)file.FindObjectAny("hped_per_tick_chan_1697");    
    if(hped_per_tick_chan_1697_temp != 0){
      hped_per_tick_chan_1697->Add(hped_per_tick_chan_1697_temp,1.0);
    }    
    if(hped_per_tick_chan_1697_temp != 0 && header != 0 && Xsrtime != XNow - GMToffset) {
      MeanPedChan1697Time[MeanPedChan1697Count] = Xsrtime;
      MeanPedChan1697    [MeanPedChan1697Count] = hped_per_tick_chan_1697_temp->GetMean(1);
      MeanPedChan1697Count++;      

      RMSPedChan1697Time[RMSPedChan1697Count] = Xsrtime;
      RMSPedChan1697    [RMSPedChan1697Count] = hped_per_tick_chan_1697_temp->GetRMS(1);
      RMSPedChan1697Count++;      
    }
        
    TH1F *hped_per_tick_chan_1838_temp = (TH1F*)file.FindObjectAny("hped_per_tick_chan_1838");    
    if(hped_per_tick_chan_1838_temp != 0){
      hped_per_tick_chan_1838->Add(hped_per_tick_chan_1838_temp,1.0);
    }    
    if(hped_per_tick_chan_1838_temp != 0 && header != 0 && Xsrtime != XNow - GMToffset) {
      MeanPedChan1838Time[MeanPedChan1838Count] = Xsrtime;
      MeanPedChan1838    [MeanPedChan1838Count] = hped_per_tick_chan_1838_temp->GetMean(1);
      MeanPedChan1838Count++;      

      RMSPedChan1838Time[RMSPedChan1838Count] = Xsrtime;
      RMSPedChan1838    [RMSPedChan1838Count] = hped_per_tick_chan_1838_temp->GetRMS(1);
      RMSPedChan1838Count++;      
    }
        
    TH1F *hped_per_tick_chan_1482_temp = (TH1F*)file.FindObjectAny("hped_per_tick_chan_1482");    
    if(hped_per_tick_chan_1482_temp != 0){
      hped_per_tick_chan_1482->Add(hped_per_tick_chan_1482_temp,1.0);
    }    
    if(hped_per_tick_chan_1482_temp != 0 && header != 0 && Xsrtime != XNow - GMToffset) {
      MeanPedChan1482Time[MeanPedChan1482Count] = Xsrtime;
      MeanPedChan1482    [MeanPedChan1482Count] = hped_per_tick_chan_1482_temp->GetMean(1);
      MeanPedChan1482Count++;      

      RMSPedChan1482Time[RMSPedChan1482Count] = Xsrtime;
      RMSPedChan1482    [RMSPedChan1482Count] = hped_per_tick_chan_1482_temp->GetRMS(1);
      RMSPedChan1482Count++;      
    }
        
    TH1F *hped_per_tick_chan_952_temp = (TH1F*)file.FindObjectAny("hped_per_tick_chan_952");    
    if(hped_per_tick_chan_952_temp != 0){
      hped_per_tick_chan_952->Add(hped_per_tick_chan_952_temp,1.0);
    }    
    if(hped_per_tick_chan_952_temp != 0 && header != 0 && Xsrtime != XNow - GMToffset) {
      MeanPedChan0952Time[MeanPedChan0952Count] = Xsrtime;
      MeanPedChan0952    [MeanPedChan0952Count] = hped_per_tick_chan_952_temp->GetMean(1);
      MeanPedChan0952Count++;      

      RMSPedChan0952Time[RMSPedChan0952Count] = Xsrtime;
      RMSPedChan0952    [RMSPedChan0952Count] = hped_per_tick_chan_952_temp->GetRMS(1);
      RMSPedChan0952Count++;      
    }
        
    TH1F *hped_per_tick_chan_454_temp = (TH1F*)file.FindObjectAny("hped_per_tick_chan_454");    
    if(hped_per_tick_chan_454_temp != 0){
      hped_per_tick_chan_454->Add(hped_per_tick_chan_454_temp,1.0);
    }    
    if(hped_per_tick_chan_454_temp != 0 && header != 0 && Xsrtime != XNow - GMToffset) {
      MeanPedChan0454Time[MeanPedChan0454Count] = Xsrtime;
      MeanPedChan0454    [MeanPedChan0454Count] = hped_per_tick_chan_454_temp->GetMean(1);
      MeanPedChan0454Count++;      

      RMSPedChan0454Time[RMSPedChan0454Count] = Xsrtime;
      RMSPedChan0454    [RMSPedChan0454Count] = hped_per_tick_chan_454_temp->GetRMS(1);
      RMSPedChan0454Count++;      
    }
        


    file.Close();

  } // end while loop over input files



  //
  // Make time/date stamp for each plot
  //
  TDatime *Ttemp = new TDatime;  // finish time
  int Xfin = Ttemp->Convert() - GMToffset;
  TDatime *Tfinish = new TDatime(Xfin);  // finish time



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

  char filename[128];
  char title[128];

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

    












  TCanvas *cADCSpecChan0020 = new TCanvas("cADCSpecChan0020","ADC Spectrum - Channel 20",1200,800);
  cADCSpecChan0020->cd();
  cADCSpecChan0020->SetLogy();
  gStyle->SetOptStat(111111);
  // hped_per_tick_chan_20->SetAxisRange(0.0,2048.0,"X");
  hped_per_tick_chan_20->SetLineWidth(2);
  hped_per_tick_chan_20->SetLineColor(kRed);
  hped_per_tick_chan_20->Draw();
  UpdateText->Draw();
  if(debug) sprintf(filename,"%s/ADCSpecChan0020_%.3u_days.png",PLOT_DIR_DEBUG.c_str(),Ndays);
  else sprintf(filename,"%s/ADCSpecChan0020_%.3u_days.png",PLOT_DIR.c_str(),Ndays);
  cADCSpecChan0020->Print(filename);

  maxtime = 0;
  max = 0.0;
  ave = 0.0;
  for(int i = 0; i < MeanPedChan0020Count; ++i) {
    ave += (double)MeanPedChan0020[i];
    if(MeanPedChan0020Time[i] > maxtime) {
      maxtime = MeanPedChan0020Time[i];
      max     = MeanPedChan0020[i];
    }
  }
  if(MeanPedChan0020Count > 0) ave = ave/(double)MeanPedChan0020Count;
  sprintf(lptext,"Last Point = %f  /  Average = %f",max,ave);
  LastPoint->Clear();
  LastPoint->AddText(lptext);
  TCanvas *cMeanPedChan0020 = new TCanvas("cMeanPedChan0020","Mean Pedestal for Channel 20",1200,800);
  cMeanPedChan0020->cd();
  gPad->SetGridx();
  TGraph *gMeanPedChan0020 = new TGraph(MeanPedChan0020Count,MeanPedChan0020Time,MeanPedChan0020);
  sprintf(title,"Mean Pedestal per Event for Channel 20");
  gMeanPedChan0020->SetTitle(title);
  gMeanPedChan0020->SetMarkerColor(kBlue);
  gMeanPedChan0020->GetXaxis()->SetTimeDisplay(1);
  gMeanPedChan0020->GetXaxis()->SetLabelSize(0.03);
  gMeanPedChan0020->GetXaxis()->SetTimeFormat(taxis_labels.c_str());
  gMeanPedChan0020->GetXaxis()->SetLimits(time_ago,XNow);
  gMeanPedChan0020->GetXaxis()->SetTitle("(central time)");
  gMeanPedChan0020->Draw("A*");
  UpdateText->Draw();
  LastPoint->Draw();
  if(debug)  sprintf(filename,"%s/MeanPedChan0020_%.3u_days.png",PLOT_DIR_DEBUG.c_str(),Ndays);
  else  sprintf(filename,"%s/MeanPedChan0020_%.3u_days.png",PLOT_DIR.c_str(),Ndays);
  cMeanPedChan0020->Print(filename);

  maxtime = 0;
  max = 0.0;
  ave = 0.0;
  for(int i = 0; i < RMSPedChan0020Count; ++i) {
    ave += (double)RMSPedChan0020[i];
    if(RMSPedChan0020Time[i] > maxtime) {
      maxtime = RMSPedChan0020Time[i];
      max     = RMSPedChan0020[i];
    }
  }
  if(RMSPedChan0020Count > 0) ave = ave/(double)RMSPedChan0020Count;
  sprintf(lptext,"Last Point = %f  /  Average = %f",max,ave);
  LastPoint->Clear();
  LastPoint->AddText(lptext);
  TCanvas *cRMSPedChan0020 = new TCanvas("cRMSPedChan0020","Pedestal RMS for Channel 20",1200,800);
  cRMSPedChan0020->cd();
  gPad->SetGridx();
  TGraph *gRMSPedChan0020 = new TGraph(RMSPedChan0020Count,RMSPedChan0020Time,RMSPedChan0020);
  sprintf(title,"Pedestal RMS for Channel 20");
  gRMSPedChan0020->SetTitle(title);
  gRMSPedChan0020->SetMarkerColor(kBlue);
  gRMSPedChan0020->GetXaxis()->SetTimeDisplay(1);
  gRMSPedChan0020->GetXaxis()->SetLabelSize(0.03);
  gRMSPedChan0020->GetXaxis()->SetTimeFormat(taxis_labels.c_str());
  gRMSPedChan0020->GetXaxis()->SetLimits(time_ago,XNow);
  gRMSPedChan0020->GetXaxis()->SetTitle("(central time)");
  gRMSPedChan0020->Draw("A*");
  UpdateText->Draw();
  LastPoint->Draw();
  if(debug)  sprintf(filename,"%s/RMSPedChan0020_%.3u_days.png",PLOT_DIR_DEBUG.c_str(),Ndays);
  else  sprintf(filename,"%s/RMSPedChan0020_%.3u_days.png",PLOT_DIR.c_str(),Ndays);
  cRMSPedChan0020->Print(filename);
  sprintf(title,"Pedestal RMS for Channel 20 - Zoom");
  gRMSPedChan0020->SetTitle(title);
  graphZoom(gRMSPedChan0020, 2.0);
  gRMSPedChan0020->Draw("A*");
  UpdateText->Draw();
  LastPoint->Draw();
  if(debug)  sprintf(filename,"%s/RMSPedChan0020_%.3u_days_zoom.png",PLOT_DIR_DEBUG.c_str(),Ndays);
  else  sprintf(filename,"%s/RMSPedChan0020_%.3u_days_zoom.png",PLOT_DIR.c_str(),Ndays);
  cRMSPedChan0020->Print(filename);




  TCanvas *cADCSpecChan0548 = new TCanvas("cADCSpecChan0548","ADC Spectrum - Channel 548",1200,800);
  cADCSpecChan0548->cd();
  cADCSpecChan0548->SetLogy();
  gStyle->SetOptStat(111111);
  // hped_per_tick_chan_548->SetAxisRange(0.0,2048.0,"X");
  hped_per_tick_chan_548->SetLineWidth(2);
  hped_per_tick_chan_548->SetLineColor(kRed);
  hped_per_tick_chan_548->Draw();
  UpdateText->Draw();
  if(debug)  sprintf(filename,"%s/ADCSpecChan0548_%.3u_days.png",PLOT_DIR_DEBUG.c_str(),Ndays);
  else  sprintf(filename,"%s/ADCSpecChan0548_%.3u_days.png",PLOT_DIR.c_str(),Ndays);
  cADCSpecChan0548->Print(filename);

  maxtime = 0;
  max = 0.0;
  ave = 0.0;
  for(int i = 0; i < MeanPedChan0548Count; ++i) {
    ave += (double)MeanPedChan0548[i];
    if(MeanPedChan0548Time[i] > maxtime) {
      maxtime = MeanPedChan0548Time[i];
      max     = MeanPedChan0548[i];
    }
  }
  if(MeanPedChan0548Count > 0) ave = ave/(double)MeanPedChan0548Count;
  sprintf(lptext,"Last Point = %f  /  Average = %f",max,ave);
  LastPoint->Clear();
  LastPoint->AddText(lptext);
  TCanvas *cMeanPedChan0548 = new TCanvas("cMeanPedChan0548","Mean Pedestal for Channel 548",1200,800);
  cMeanPedChan0548->cd();
  gPad->SetGridx();
  TGraph *gMeanPedChan0548 = new TGraph(MeanPedChan0548Count,MeanPedChan0548Time,MeanPedChan0548);
  sprintf(title,"Mean Pedestal for Channel 548");
  gMeanPedChan0548->SetTitle(title);
  gMeanPedChan0548->SetMarkerColor(kBlue);
  gMeanPedChan0548->GetXaxis()->SetTimeDisplay(1);
  gMeanPedChan0548->GetXaxis()->SetLabelSize(0.03);
  gMeanPedChan0548->GetXaxis()->SetTimeFormat(taxis_labels.c_str());
  gMeanPedChan0548->GetXaxis()->SetLimits(time_ago,XNow);
  gMeanPedChan0548->GetXaxis()->SetTitle("(central time)");
  gMeanPedChan0548->Draw("A*");
  UpdateText->Draw();
  LastPoint->Draw();
  if(debug)  sprintf(filename,"%s/MeanPedChan0548_%.3u_days.png",PLOT_DIR_DEBUG.c_str(),Ndays);
  else  sprintf(filename,"%s/MeanPedChan0548_%.3u_days.png",PLOT_DIR.c_str(),Ndays);
  cMeanPedChan0548->Print(filename);

  maxtime = 0;
  max = 0.0;
  ave = 0.0;
  for(int i = 0; i < RMSPedChan0548Count; ++i) {
    ave += (double)RMSPedChan0548[i];
    if(RMSPedChan0548Time[i] > maxtime) {
      maxtime = RMSPedChan0548Time[i];
      max     = RMSPedChan0548[i];
    }
  }
  if(RMSPedChan0548Count > 0) ave = ave/(double)RMSPedChan0548Count;
  sprintf(lptext,"Last Point = %f  /  Average = %f",max,ave);
  LastPoint->Clear();
  LastPoint->AddText(lptext);
  TCanvas *cRMSPedChan0548 = new TCanvas("cRMSPedChan0548","Pedestal RMS for Channel 548",1200,800);
  cRMSPedChan0548->cd();
  gPad->SetGridx();
  TGraph *gRMSPedChan0548 = new TGraph(RMSPedChan0548Count,RMSPedChan0548Time,RMSPedChan0548);
  sprintf(title,"Pedestal RMS for Channel 548");
  gRMSPedChan0548->SetTitle(title);
  gRMSPedChan0548->SetMarkerColor(kBlue);
  gRMSPedChan0548->GetXaxis()->SetTimeDisplay(1);
  gRMSPedChan0548->GetXaxis()->SetLabelSize(0.03);
  gRMSPedChan0548->GetXaxis()->SetTimeFormat(taxis_labels.c_str());
  gRMSPedChan0548->GetXaxis()->SetLimits(time_ago,XNow);
  gRMSPedChan0548->GetXaxis()->SetTitle("(central time)");
  gRMSPedChan0548->Draw("A*");
  UpdateText->Draw();
  LastPoint->Draw();
  if(debug)  sprintf(filename,"%s/RMSPedChan0548_%.3u_days.png",PLOT_DIR_DEBUG.c_str(),Ndays);
  else  sprintf(filename,"%s/RMSPedChan0548_%.3u_days.png",PLOT_DIR.c_str(),Ndays);
  cRMSPedChan0548->Print(filename);
  sprintf(title,"Pedestal RMS for Channel 548 - Zoom");
  gRMSPedChan0548->SetTitle(title);
  graphZoom(gRMSPedChan0548, 2.0);
  gRMSPedChan0548->Draw("A*");
  UpdateText->Draw();
  LastPoint->Draw();
  if(debug)  sprintf(filename,"%s/RMSPedChan0548_%.3u_days_zoom.png",PLOT_DIR_DEBUG.c_str(),Ndays);
  else  sprintf(filename,"%s/RMSPedChan0548_%.3u_days_zoom.png",PLOT_DIR.c_str(),Ndays);
  cRMSPedChan0548->Print(filename);




  TCanvas *cADCSpecChan1297 = new TCanvas("cADCSpecChan1297","ADC Spectrum - Channel 1297",1200,800);
  cADCSpecChan1297->cd();
  cADCSpecChan1297->SetLogy();
  gStyle->SetOptStat(111111);
  // hped_per_tick_chan_1297->SetAxisRange(0.0,2048.0,"X");
  hped_per_tick_chan_1297->SetLineWidth(2);
  hped_per_tick_chan_1297->SetLineColor(kRed);
  hped_per_tick_chan_1297->Draw();
  UpdateText->Draw();
  if(debug)  sprintf(filename,"%s/ADCSpecChan1297_%.3u_days.png",PLOT_DIR_DEBUG.c_str(),Ndays);
  else  sprintf(filename,"%s/ADCSpecChan1297_%.3u_days.png",PLOT_DIR.c_str(),Ndays);
  cADCSpecChan1297->Print(filename);

  maxtime = 0;
  max = 0.0;
  ave = 0.0;
  for(int i = 0; i < MeanPedChan1297Count; ++i) {
    ave += (double)MeanPedChan1297[i];
    if(MeanPedChan1297Time[i] > maxtime) {
      maxtime = MeanPedChan1297Time[i];
      max     = MeanPedChan1297[i];
    }
  }
  if(MeanPedChan1297Count > 0) ave = ave/(double)MeanPedChan1297Count;
  sprintf(lptext,"Last Point = %f  /  Average = %f",max,ave);
  LastPoint->Clear();
  LastPoint->AddText(lptext);
  TCanvas *cMeanPedChan1297 = new TCanvas("cMeanPedChan1297","Mean Pedestal for Channel 1297",1200,800);
  cMeanPedChan1297->cd();
  gPad->SetGridx();
  TGraph *gMeanPedChan1297 = new TGraph(MeanPedChan1297Count,MeanPedChan1297Time,MeanPedChan1297);
  sprintf(title,"Mean Pedestal for Channel 1297");
  gMeanPedChan1297->SetTitle(title);
  gMeanPedChan1297->SetMarkerColor(kBlue);
  gMeanPedChan1297->GetXaxis()->SetTimeDisplay(1);
  gMeanPedChan1297->GetXaxis()->SetLabelSize(0.03);
  gMeanPedChan1297->GetXaxis()->SetTimeFormat(taxis_labels.c_str());
  gMeanPedChan1297->GetXaxis()->SetLimits(time_ago,XNow);
  gMeanPedChan1297->GetXaxis()->SetTitle("(central time)");
  gMeanPedChan1297->Draw("A*");
  UpdateText->Draw();
  LastPoint->Draw();
  if(debug)  sprintf(filename,"%s/MeanPedChan1297_%.3u_days.png",PLOT_DIR_DEBUG.c_str(),Ndays);
  else  sprintf(filename,"%s/MeanPedChan1297_%.3u_days.png",PLOT_DIR.c_str(),Ndays);
  cMeanPedChan1297->Print(filename);

  maxtime = 0;
  max = 0.0;
  ave = 0.0;
  for(int i = 0; i < RMSPedChan1297Count; ++i) {
    ave += (double)RMSPedChan1297[i];
    if(RMSPedChan1297Time[i] > maxtime) {
      maxtime = RMSPedChan1297Time[i];
      max     = RMSPedChan1297[i];
    }
  }
  if(RMSPedChan1297Count > 0) ave = ave/(double)RMSPedChan1297Count;
  sprintf(lptext,"Last Point = %f  /  Average = %f",max,ave);
  LastPoint->Clear();
  LastPoint->AddText(lptext);
  TCanvas *cRMSPedChan1297 = new TCanvas("cRMSPedChan1297","Pedestal RMS for Channel 1297",1200,800);
  cRMSPedChan1297->cd();
  gPad->SetGridx();
  TGraph *gRMSPedChan1297 = new TGraph(RMSPedChan1297Count,RMSPedChan1297Time,RMSPedChan1297);
  sprintf(title,"Pedestal RMS for Channel 1297");
  gRMSPedChan1297->SetTitle(title);
  gRMSPedChan1297->SetMarkerColor(kBlue);
  gRMSPedChan1297->GetXaxis()->SetTimeDisplay(1);
  gRMSPedChan1297->GetXaxis()->SetLabelSize(0.03);
  gRMSPedChan1297->GetXaxis()->SetTimeFormat(taxis_labels.c_str());
  gRMSPedChan1297->GetXaxis()->SetLimits(time_ago,XNow);
  gRMSPedChan1297->GetXaxis()->SetTitle("(central time)");
  gRMSPedChan1297->Draw("A*");
  UpdateText->Draw();
  LastPoint->Draw();
  if(debug)  sprintf(filename,"%s/RMSPedChan1297_%.3u_days.png",PLOT_DIR_DEBUG.c_str(),Ndays);
  else  sprintf(filename,"%s/RMSPedChan1297_%.3u_days.png",PLOT_DIR.c_str(),Ndays);
  cRMSPedChan1297->Print(filename);
  sprintf(title,"Pedestal RMS for Channel 1297 - Zoom");
  gRMSPedChan1297->SetTitle(title);
  graphZoom(gRMSPedChan1297, 2.0);
  gRMSPedChan1297->Draw("A*");
  UpdateText->Draw();
  LastPoint->Draw();
  if(debug)  sprintf(filename,"%s/RMSPedChan1297_%.3u_days_zoom.png",PLOT_DIR_DEBUG.c_str(),Ndays);
  else  sprintf(filename,"%s/RMSPedChan1297_%.3u_days_zoom.png",PLOT_DIR.c_str(),Ndays);
  cRMSPedChan1297->Print(filename);



  TCanvas *cADCSpecChan1697 = new TCanvas("cADCSpecChan1697","ADC Spectrum - Channel 1697",1200,800);
  cADCSpecChan1697->cd();
  cADCSpecChan1697->SetLogy();
  gStyle->SetOptStat(111111);
  // hped_per_tick_chan_1697->SetAxisRange(0.0,2048.0,"X");
  hped_per_tick_chan_1697->SetLineWidth(2);
  hped_per_tick_chan_1697->SetLineColor(kRed);
  hped_per_tick_chan_1697->Draw();
  UpdateText->Draw();
  if(debug)  sprintf(filename,"%s/ADCSpecChan1697_%.3u_days.png",PLOT_DIR_DEBUG.c_str(),Ndays);
  else  sprintf(filename,"%s/ADCSpecChan1697_%.3u_days.png",PLOT_DIR.c_str(),Ndays);
  cADCSpecChan1697->Print(filename);

  maxtime = 0;
  max = 0.0;
  ave = 0.0;
  for(int i = 0; i < MeanPedChan1697Count; ++i) {
    ave += (double)MeanPedChan1697[i];
    if(MeanPedChan1697Time[i] > maxtime) {
      maxtime = MeanPedChan1697Time[i];
      max     = MeanPedChan1697[i];
    }
  }
  if(MeanPedChan1697Count > 0) ave = ave/(double)MeanPedChan1697Count;
  sprintf(lptext,"Last Point = %f  /  Average = %f",max,ave);
  LastPoint->Clear();
  LastPoint->AddText(lptext);
  TCanvas *cMeanPedChan1697 = new TCanvas("cMeanPedChan1697","Mean Pedestal for Channel 1697",1200,800);
  cMeanPedChan1697->cd();
  gPad->SetGridx();
  TGraph *gMeanPedChan1697 = new TGraph(MeanPedChan1697Count,MeanPedChan1697Time,MeanPedChan1697);
  sprintf(title,"Mean Pedestal for Channel 1697");
  gMeanPedChan1697->SetTitle(title);
  gMeanPedChan1697->SetMarkerColor(kBlue);
  gMeanPedChan1697->GetXaxis()->SetTimeDisplay(1);
  gMeanPedChan1697->GetXaxis()->SetLabelSize(0.03);
  gMeanPedChan1697->GetXaxis()->SetTimeFormat(taxis_labels.c_str());
  gMeanPedChan1697->GetXaxis()->SetLimits(time_ago,XNow);
  gMeanPedChan1697->GetXaxis()->SetTitle("(central time)");
  gMeanPedChan1697->Draw("A*");
  UpdateText->Draw();
  LastPoint->Draw();
  if(debug)  sprintf(filename,"%s/MeanPedChan1697_%.3u_days.png",PLOT_DIR_DEBUG.c_str(),Ndays);
  else  sprintf(filename,"%s/MeanPedChan1697_%.3u_days.png",PLOT_DIR.c_str(),Ndays);
  cMeanPedChan1697->Print(filename);

  maxtime = 0;
  max = 0.0;
  ave = 0.0;
  for(int i = 0; i < RMSPedChan1697Count; ++i) {
    ave += (double)RMSPedChan1697[i];
    if(RMSPedChan1697Time[i] > maxtime) {
      maxtime = RMSPedChan1697Time[i];
      max     = RMSPedChan1697[i];
    }
  }
  if(RMSPedChan1697Count > 0) ave = ave/(double)RMSPedChan1697Count;
  sprintf(lptext,"Last Point = %f  /  Average = %f",max,ave);
  LastPoint->Clear();
  LastPoint->AddText(lptext);
  TCanvas *cRMSPedChan1697 = new TCanvas("cRMSPedChan1697","Pedestal RMS for Channel 1697",1200,800);
  cRMSPedChan1697->cd();
  gPad->SetGridx();
  TGraph *gRMSPedChan1697 = new TGraph(RMSPedChan1697Count,RMSPedChan1697Time,RMSPedChan1697);
  sprintf(title,"Pedestal RMS for Channel 1697");
  gRMSPedChan1697->SetTitle(title);
  gRMSPedChan1697->SetMarkerColor(kBlue);
  gRMSPedChan1697->GetXaxis()->SetTimeDisplay(1);
  gRMSPedChan1697->GetXaxis()->SetLabelSize(0.03);
  gRMSPedChan1697->GetXaxis()->SetTimeFormat(taxis_labels.c_str());
  gRMSPedChan1697->GetXaxis()->SetLimits(time_ago,XNow);
  gRMSPedChan1697->GetXaxis()->SetTitle("(central time)");
  gRMSPedChan1697->Draw("A*");
  UpdateText->Draw();
  LastPoint->Draw();
  if(debug)  sprintf(filename,"%s/RMSPedChan1697_%.3u_days.png",PLOT_DIR_DEBUG.c_str(),Ndays);
  else  sprintf(filename,"%s/RMSPedChan1697_%.3u_days.png",PLOT_DIR.c_str(),Ndays);
  cRMSPedChan1697->Print(filename);
  sprintf(title,"Pedestal RMS for Channel 1697 - Zoom");
  gRMSPedChan1697->SetTitle(title);
  graphZoom(gRMSPedChan1697, 2.0);
  gRMSPedChan1697->Draw("A*");
  UpdateText->Draw();
  LastPoint->Draw();
  if(debug)  sprintf(filename,"%s/RMSPedChan1697_%.3u_days_zoom.png",PLOT_DIR_DEBUG.c_str(),Ndays);
  else  sprintf(filename,"%s/RMSPedChan1697_%.3u_days_zoom.png",PLOT_DIR.c_str(),Ndays);
  cRMSPedChan1697->Print(filename);



  TCanvas *cADCSpecChan1838 = new TCanvas("cADCSpecChan1838","ADC Spectrum - Channel 1838",1200,800);
  cADCSpecChan1838->cd();
  cADCSpecChan1838->SetLogy();
  gStyle->SetOptStat(111111);
  // hped_per_tick_chan_1838->SetAxisRange(0.0,2048.0,"X");
  hped_per_tick_chan_1838->SetLineWidth(2);
  hped_per_tick_chan_1838->SetLineColor(kRed);
  hped_per_tick_chan_1838->Draw();
  UpdateText->Draw();
  if(debug)  sprintf(filename,"%s/ADCSpecChan1838_%.3u_days.png",PLOT_DIR_DEBUG.c_str(),Ndays);
  else  sprintf(filename,"%s/ADCSpecChan1838_%.3u_days.png",PLOT_DIR.c_str(),Ndays);
  cADCSpecChan1838->Print(filename);

  maxtime = 0;
  max = 0.0;
  ave = 0.0;
  for(int i = 0; i < MeanPedChan1838Count; ++i) {
    ave += (double)MeanPedChan1838[i];
    if(MeanPedChan1838Time[i] > maxtime) {
      maxtime = MeanPedChan1838Time[i];
      max     = MeanPedChan1838[i];
    }
  }
  if(MeanPedChan1838Count > 0) ave = ave/(double)MeanPedChan1838Count;
  sprintf(lptext,"Last Point = %f  /  Average = %f",max,ave);
  LastPoint->Clear();
  LastPoint->AddText(lptext);
  TCanvas *cMeanPedChan1838 = new TCanvas("cMeanPedChan1838","Mean Pedestal for Channel 1838",1200,800);
  cMeanPedChan1838->cd();
  gPad->SetGridx();
  TGraph *gMeanPedChan1838 = new TGraph(MeanPedChan1838Count,MeanPedChan1838Time,MeanPedChan1838);
  sprintf(title,"Mean Pedestal per Event for Channel 1838");
  gMeanPedChan1838->SetTitle(title);
  gMeanPedChan1838->SetMarkerColor(kBlue);
  gMeanPedChan1838->GetXaxis()->SetTimeDisplay(1);
  gMeanPedChan1838->GetXaxis()->SetLabelSize(0.03);
  gMeanPedChan1838->GetXaxis()->SetTimeFormat(taxis_labels.c_str());
  gMeanPedChan1838->GetXaxis()->SetLimits(time_ago,XNow);
  gMeanPedChan1838->GetXaxis()->SetTitle("(central time)");
  gMeanPedChan1838->Draw("A*");
  UpdateText->Draw();
  LastPoint->Draw();
  if(debug)  sprintf(filename,"%s/MeanPedChan1838_%.3u_days.png",PLOT_DIR_DEBUG.c_str(),Ndays);
  else  sprintf(filename,"%s/MeanPedChan1838_%.3u_days.png",PLOT_DIR.c_str(),Ndays);
  cMeanPedChan1838->Print(filename);

  maxtime = 0;
  max = 0.0;
  ave = 0.0;
  for(int i = 0; i < RMSPedChan1838Count; ++i) {
    ave += (double)RMSPedChan1838[i];
    if(RMSPedChan1838Time[i] > maxtime) {
      maxtime = RMSPedChan1838Time[i];
      max     = RMSPedChan1838[i];
    }
  }
  if(RMSPedChan1838Count > 0) ave = ave/(double)RMSPedChan1838Count;
  sprintf(lptext,"Last Point = %f  /  Average = %f",max,ave);
  LastPoint->Clear();
  LastPoint->AddText(lptext);
  TCanvas *cRMSPedChan1838 = new TCanvas("cRMSPedChan1838","Pedestal RMS for Channel 1838",1200,800);
  cRMSPedChan1838->cd();
  gPad->SetGridx();
  TGraph *gRMSPedChan1838 = new TGraph(RMSPedChan1838Count,RMSPedChan1838Time,RMSPedChan1838);
  sprintf(title,"Pedestal RMS for Channel 1838");
  gRMSPedChan1838->SetTitle(title);
  gRMSPedChan1838->SetMarkerColor(kBlue);
  gRMSPedChan1838->GetXaxis()->SetTimeDisplay(1);
  gRMSPedChan1838->GetXaxis()->SetLabelSize(0.03);
  gRMSPedChan1838->GetXaxis()->SetTimeFormat(taxis_labels.c_str());
  gRMSPedChan1838->GetXaxis()->SetLimits(time_ago,XNow);
  gRMSPedChan1838->GetXaxis()->SetTitle("(central time)");
  gRMSPedChan1838->Draw("A*");
  UpdateText->Draw();
  LastPoint->Draw();
  if(debug)  sprintf(filename,"%s/RMSPedChan1838_%.3u_days.png",PLOT_DIR_DEBUG.c_str(),Ndays);
  else  sprintf(filename,"%s/RMSPedChan1838_%.3u_days.png",PLOT_DIR.c_str(),Ndays);
  cRMSPedChan1838->Print(filename);
  sprintf(title,"Pedestal RMS for Channel 1838 - Zoom");
  gRMSPedChan1838->SetTitle(title);
  graphZoom(gRMSPedChan1838, 2.0);
  gRMSPedChan1838->Draw("A*");
  UpdateText->Draw();
  LastPoint->Draw();
  if(debug)  sprintf(filename,"%s/RMSPedChan1838_%.3u_days_zoom.png",PLOT_DIR_DEBUG.c_str(),Ndays);
  else  sprintf(filename,"%s/RMSPedChan1838_%.3u_days_zoom.png",PLOT_DIR.c_str(),Ndays);
  cRMSPedChan1838->Print(filename);



  TCanvas *cADCSpecChan1482 = new TCanvas("cADCSpecChan1482","ADC Spectrum - Channel 1482",1200,800);
  cADCSpecChan1482->cd();
  cADCSpecChan1482->SetLogy();
  gStyle->SetOptStat(111111);
  // hped_per_tick_chan_1482->SetAxisRange(0.0,2048.0,"X");
  hped_per_tick_chan_1482->SetLineWidth(2);
  hped_per_tick_chan_1482->SetLineColor(kRed);
  hped_per_tick_chan_1482->Draw();
  UpdateText->Draw();
  if(debug)  sprintf(filename,"%s/ADCSpecChan1482_%.3u_days.png",PLOT_DIR_DEBUG.c_str(),Ndays);
  else  sprintf(filename,"%s/ADCSpecChan1482_%.3u_days.png",PLOT_DIR.c_str(),Ndays);
  cADCSpecChan1482->Print(filename);

  maxtime = 0;
  max = 0.0;
  ave = 0.0;
  for(int i = 0; i < MeanPedChan1482Count; ++i) {
    ave += (double)MeanPedChan1482[i];
    if(MeanPedChan1482Time[i] > maxtime) {
      maxtime = MeanPedChan1482Time[i];
      max     = MeanPedChan1482[i];
    }
  }
  if(MeanPedChan1482Count > 0) ave = ave/(double)MeanPedChan1482Count;
  sprintf(lptext,"Last Point = %f  /  Average = %f",max,ave);
  LastPoint->Clear();
  LastPoint->AddText(lptext);
  TCanvas *cMeanPedChan1482 = new TCanvas("cMeanPedChan1482","Mean Pedestal for Channel 1482",1200,800);
  cMeanPedChan1482->cd();
  gPad->SetGridx();
  TGraph *gMeanPedChan1482 = new TGraph(MeanPedChan1482Count,MeanPedChan1482Time,MeanPedChan1482);
  sprintf(title,"Mean Pedestal per Event for Channel 1482");
  gMeanPedChan1482->SetTitle(title);
  gMeanPedChan1482->SetMarkerColor(kBlue);
  gMeanPedChan1482->GetXaxis()->SetTimeDisplay(1);
  gMeanPedChan1482->GetXaxis()->SetLabelSize(0.03);
  gMeanPedChan1482->GetXaxis()->SetTimeFormat(taxis_labels.c_str());
  gMeanPedChan1482->GetXaxis()->SetLimits(time_ago,XNow);
  gMeanPedChan1482->GetXaxis()->SetTitle("(central time)");
  gMeanPedChan1482->Draw("A*");
  UpdateText->Draw();
  LastPoint->Draw();
  if(debug)  sprintf(filename,"%s/MeanPedChan1482_%.3u_days.png",PLOT_DIR_DEBUG.c_str(),Ndays);
  else  sprintf(filename,"%s/MeanPedChan1482_%.3u_days.png",PLOT_DIR.c_str(),Ndays);
  cMeanPedChan1482->Print(filename);

  maxtime = 0;
  max = 0.0;
  ave = 0.0;
  for(int i = 0; i < RMSPedChan1482Count; ++i) {
    ave += (double)RMSPedChan1482[i];
    if(RMSPedChan1482Time[i] > maxtime) {
      maxtime = RMSPedChan1482Time[i];
      max     = RMSPedChan1482[i];
    }
  }
  if(RMSPedChan1482Count > 0) ave = ave/(double)RMSPedChan1482Count;
  sprintf(lptext,"Last Point = %f  /  Average = %f",max,ave);
  LastPoint->Clear();
  LastPoint->AddText(lptext);
  TCanvas *cRMSPedChan1482 = new TCanvas("cRMSPedChan1482","Pedestal RMS for Channel 1482",1200,800);
  cRMSPedChan1482->cd();
  gPad->SetGridx();
  TGraph *gRMSPedChan1482 = new TGraph(RMSPedChan1482Count,RMSPedChan1482Time,RMSPedChan1482);
  sprintf(title,"Pedestal RMS for Channel 1482");
  gRMSPedChan1482->SetTitle(title);
  gRMSPedChan1482->SetMarkerColor(kBlue);
  gRMSPedChan1482->GetXaxis()->SetTimeDisplay(1);
  gRMSPedChan1482->GetXaxis()->SetLabelSize(0.03);
  gRMSPedChan1482->GetXaxis()->SetTimeFormat(taxis_labels.c_str());
  gRMSPedChan1482->GetXaxis()->SetLimits(time_ago,XNow);
  gRMSPedChan1482->GetXaxis()->SetTitle("(central time)");
  gRMSPedChan1482->Draw("A*");
  UpdateText->Draw();
  LastPoint->Draw();
  if(debug)  sprintf(filename,"%s/RMSPedChan1482_%.3u_days.png",PLOT_DIR_DEBUG.c_str(),Ndays);
  else  sprintf(filename,"%s/RMSPedChan1482_%.3u_days.png",PLOT_DIR.c_str(),Ndays);
  cRMSPedChan1482->Print(filename);
  sprintf(title,"Pedestal RMS for Channel 1482 - Zoom");
  gRMSPedChan1482->SetTitle(title);
  graphZoom(gRMSPedChan1482, 2.0);
  gRMSPedChan1482->Draw("A*");
  UpdateText->Draw();
  LastPoint->Draw();
  if(debug)  sprintf(filename,"%s/RMSPedChan1482_%.3u_days_zoom.png",PLOT_DIR_DEBUG.c_str(),Ndays);
  else  sprintf(filename,"%s/RMSPedChan1482_%.3u_days_zoom.png",PLOT_DIR.c_str(),Ndays);
  cRMSPedChan1482->Print(filename);



  TCanvas *cADCSpecChan0952 = new TCanvas("cADCSpecChan0952","ADC Spectrum - Channel 952",1200,800);
  cADCSpecChan0952->cd();
  cADCSpecChan0952->SetLogy();
  gStyle->SetOptStat(111111);
  // hped_per_tick_chan_952->SetAxisRange(0.0,2048.0,"X");
  hped_per_tick_chan_952->SetLineWidth(2);
  hped_per_tick_chan_952->SetLineColor(kRed);
  hped_per_tick_chan_952->Draw();
  UpdateText->Draw();
  if(debug)  sprintf(filename,"%s/ADCSpecChan0952_%.3u_days.png",PLOT_DIR_DEBUG.c_str(),Ndays);
  else  sprintf(filename,"%s/ADCSpecChan0952_%.3u_days.png",PLOT_DIR.c_str(),Ndays);
  cADCSpecChan0952->Print(filename);

  maxtime = 0;
  max = 0.0;
  ave = 0.0;
  for(int i = 0; i < MeanPedChan0952Count; ++i) {
    ave += (double)MeanPedChan0952[i];
    if(MeanPedChan0952Time[i] > maxtime) {
      maxtime = MeanPedChan0952Time[i];
      max     = MeanPedChan0952[i];
    }
  }
  if(MeanPedChan0952Count > 0) ave = ave/(double)MeanPedChan0952Count;
  sprintf(lptext,"Last Point = %f  /  Average = %f",max,ave);
  LastPoint->Clear();
  LastPoint->AddText(lptext);
  TCanvas *cMeanPedChan0952 = new TCanvas("cMeanPedChan0952","Mean Pedestal for Channel 952",1200,800);
  cMeanPedChan0952->cd();
  gPad->SetGridx();
  TGraph *gMeanPedChan0952 = new TGraph(MeanPedChan0952Count,MeanPedChan0952Time,MeanPedChan0952);
  sprintf(title,"Mean Pedestal per Event for Channel 952");
  gMeanPedChan0952->SetTitle(title);
  gMeanPedChan0952->SetMarkerColor(kBlue);
  gMeanPedChan0952->GetXaxis()->SetTimeDisplay(1);
  gMeanPedChan0952->GetXaxis()->SetLabelSize(0.03);
  gMeanPedChan0952->GetXaxis()->SetTimeFormat(taxis_labels.c_str());
  gMeanPedChan0952->GetXaxis()->SetLimits(time_ago,XNow);
  gMeanPedChan0952->GetXaxis()->SetTitle("(central time)");
  gMeanPedChan0952->Draw("A*");
  UpdateText->Draw();
  LastPoint->Draw();
  if(debug)  sprintf(filename,"%s/MeanPedChan0952_%.3u_days.png",PLOT_DIR_DEBUG.c_str(),Ndays);
  else  sprintf(filename,"%s/MeanPedChan0952_%.3u_days.png",PLOT_DIR.c_str(),Ndays);
  cMeanPedChan0952->Print(filename);

  maxtime = 0;
  max = 0.0;
  ave = 0.0;
  for(int i = 0; i < RMSPedChan0952Count; ++i) {
    ave += (double)RMSPedChan0952[i];
    if(RMSPedChan0952Time[i] > maxtime) {
      maxtime = RMSPedChan0952Time[i];
      max     = RMSPedChan0952[i];
    }
  }
  if(RMSPedChan0952Count > 0) ave = ave/(double)RMSPedChan0952Count;
  sprintf(lptext,"Last Point = %f  /  Average = %f",max,ave);
  LastPoint->Clear();
  LastPoint->AddText(lptext);
  TCanvas *cRMSPedChan0952 = new TCanvas("cRMSPedChan0952","Pedestal RMS for Channel 952",1200,800);
  cRMSPedChan0952->cd();
  gPad->SetGridx();
  TGraph *gRMSPedChan0952 = new TGraph(RMSPedChan0952Count,RMSPedChan0952Time,RMSPedChan0952);
  sprintf(title,"Pedestal RMS for Channel 952");
  gRMSPedChan0952->SetTitle(title);
  gRMSPedChan0952->SetMarkerColor(kBlue);
  gRMSPedChan0952->GetXaxis()->SetTimeDisplay(1);
  gRMSPedChan0952->GetXaxis()->SetLabelSize(0.03);
  gRMSPedChan0952->GetXaxis()->SetTimeFormat(taxis_labels.c_str());
  gRMSPedChan0952->GetXaxis()->SetLimits(time_ago,XNow);
  gRMSPedChan0952->GetXaxis()->SetTitle("(central time)");
  gRMSPedChan0952->Draw("A*");
  UpdateText->Draw();
  LastPoint->Draw();
  if(debug)  sprintf(filename,"%s/RMSPedChan0952_%.3u_days.png",PLOT_DIR_DEBUG.c_str(),Ndays);
  else  sprintf(filename,"%s/RMSPedChan0952_%.3u_days.png",PLOT_DIR.c_str(),Ndays);
  cRMSPedChan0952->Print(filename);
  sprintf(title,"Pedestal RMS for Channel 952 - Zoom");
  gRMSPedChan0952->SetTitle(title);
  graphZoom(gRMSPedChan0952, 2.0);
  gRMSPedChan0952->Draw("A*");
  UpdateText->Draw();
  LastPoint->Draw();
  if(debug)  sprintf(filename,"%s/RMSPedChan0952_%.3u_days_zoom.png",PLOT_DIR_DEBUG.c_str(),Ndays);
  else  sprintf(filename,"%s/RMSPedChan0952_%.3u_days_zoom.png",PLOT_DIR.c_str(),Ndays);
  cRMSPedChan0952->Print(filename);



  TCanvas *cADCSpecChan0454 = new TCanvas("cADCSpecChan0454","ADC Spectrum - Channel 454",1200,800);
  cADCSpecChan0454->cd();
  cADCSpecChan0454->SetLogy();
  gStyle->SetOptStat(111111);
  // hped_per_tick_chan_454->SetAxisRange(0.0,2048.0,"X");
  hped_per_tick_chan_454->SetLineWidth(2);
  hped_per_tick_chan_454->SetLineColor(kRed);
  hped_per_tick_chan_454->Draw();
  UpdateText->Draw();
  if(debug)  sprintf(filename,"%s/ADCSpecChan0454_%.3u_days.png",PLOT_DIR_DEBUG.c_str(),Ndays);
  else  sprintf(filename,"%s/ADCSpecChan0454_%.3u_days.png",PLOT_DIR.c_str(),Ndays);
  cADCSpecChan0454->Print(filename);

  maxtime = 0;
  max = 0.0;
  ave = 0.0;
  for(int i = 0; i < MeanPedChan0454Count; ++i) {
    ave += (double)MeanPedChan0454[i];
    if(MeanPedChan0454Time[i] > maxtime) {
      maxtime = MeanPedChan0454Time[i];
      max     = MeanPedChan0454[i];
    }
  }
  if(MeanPedChan0454Count > 0) ave = ave/(double)MeanPedChan0454Count;
  sprintf(lptext,"Last Point = %f  /  Average = %f",max,ave);
  LastPoint->Clear();
  LastPoint->AddText(lptext);
  TCanvas *cMeanPedChan0454 = new TCanvas("cMeanPedChan0454","Mean Pedestal for Channel 454",1200,800);
  cMeanPedChan0454->cd();
  gPad->SetGridx();
  TGraph *gMeanPedChan0454 = new TGraph(MeanPedChan0454Count,MeanPedChan0454Time,MeanPedChan0454);
  sprintf(title,"Mean Pedestal per Event for Channel 454");
  gMeanPedChan0454->SetTitle(title);
  gMeanPedChan0454->SetMarkerColor(kBlue);
  gMeanPedChan0454->GetXaxis()->SetTimeDisplay(1);
  gMeanPedChan0454->GetXaxis()->SetLabelSize(0.03);
  gMeanPedChan0454->GetXaxis()->SetTimeFormat(taxis_labels.c_str());
  gMeanPedChan0454->GetXaxis()->SetLimits(time_ago,XNow);
  gMeanPedChan0454->GetXaxis()->SetTitle("(central time)");
  gMeanPedChan0454->Draw("A*");
  UpdateText->Draw();
  LastPoint->Draw();
  if(debug)  sprintf(filename,"%s/MeanPedChan0454_%.3u_days.png",PLOT_DIR_DEBUG.c_str(),Ndays);
  else  sprintf(filename,"%s/MeanPedChan0454_%.3u_days.png",PLOT_DIR.c_str(),Ndays);
  cMeanPedChan0454->Print(filename);

  maxtime = 0;
  max = 0.0;
  ave = 0.0;
  for(int i = 0; i < RMSPedChan0454Count; ++i) {
    ave += (double)RMSPedChan0454[i];
    if(RMSPedChan0454Time[i] > maxtime) {
      maxtime = RMSPedChan0454Time[i];
      max     = RMSPedChan0454[i];
    }
  }
  if(RMSPedChan0454Count > 0) ave = ave/(double)RMSPedChan0454Count;
  sprintf(lptext,"Last Point = %f  /  Average = %f",max,ave);
  LastPoint->Clear();
  LastPoint->AddText(lptext);
  TCanvas *cRMSPedChan0454 = new TCanvas("cRMSPedChan0454","Pedestal RMS for Channel 454",1200,800);
  cRMSPedChan0454->cd();
  gPad->SetGridx();
  TGraph *gRMSPedChan0454 = new TGraph(RMSPedChan0454Count,RMSPedChan0454Time,RMSPedChan0454);
  sprintf(title,"Pedestal RMS for Channel 454");
  gRMSPedChan0454->SetTitle(title);
  gRMSPedChan0454->SetMarkerColor(kBlue);
  gRMSPedChan0454->GetXaxis()->SetTimeDisplay(1);
  gRMSPedChan0454->GetXaxis()->SetLabelSize(0.03);
  gRMSPedChan0454->GetXaxis()->SetTimeFormat(taxis_labels.c_str());
  gRMSPedChan0454->GetXaxis()->SetLimits(time_ago,XNow);
  gRMSPedChan0454->GetXaxis()->SetTitle("(central time)");
  gRMSPedChan0454->Draw("A*");
  UpdateText->Draw();
  LastPoint->Draw();
  if(debug)  sprintf(filename,"%s/RMSPedChan0454_%.3u_days.png",PLOT_DIR_DEBUG.c_str(),Ndays);
  else  sprintf(filename,"%s/RMSPedChan0454_%.3u_days.png",PLOT_DIR.c_str(),Ndays);
  cRMSPedChan0454->Print(filename);
  sprintf(title,"Pedestal RMS for Channel 454 - Zoom");
  gRMSPedChan0454->SetTitle(title);
  graphZoom(gRMSPedChan0454, 2.0);
  gRMSPedChan0454->Draw("A*");
  UpdateText->Draw();
  LastPoint->Draw();
  if(debug)  sprintf(filename,"%s/RMSPedChan0454_%.3u_days_zoom.png",PLOT_DIR_DEBUG.c_str(),Ndays);
  else  sprintf(filename,"%s/RMSPedChan0454_%.3u_days_zoom.png",PLOT_DIR.c_str(),Ndays);
  cRMSPedChan0454->Print(filename);




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
