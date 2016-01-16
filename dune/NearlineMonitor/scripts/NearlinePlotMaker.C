#include "NearlinePlotMaker.h"

const std::string PLOT_DIR = "/web/sites/lbne-dqm.fnal.gov/htdocs/NearlineMonitoring/plots";

//
// Make the plots for the nearline webpage from the nearline output files.
//

Long64_t NearlinePlotMaker(int Ndays){

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



  // Book histos...
  TH1F *hped_per_tick_chan_20   = new TH1F("hped_per_tick_chan_20",  "ADC per Tick - Channel 20;ADC"  ,100,0,2048);
  TH1F *hped_per_tick_chan_548  = new TH1F("hped_per_tick_chan_548", "ADC per Tick - Channel 548;ADC" ,100,0,2048);
  TH1F *hped_per_tick_chan_1297 = new TH1F("hped_per_tick_chan_1297","ADC per Tick - Channel 1297;ADC",100,0,2048);
  TH1F *hped_per_tick_chan_1697 = new TH1F("hped_per_tick_chan_1697","ADC per Tick - Channel 1697;ADC",100,0,2048);


  
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
      hour = HourEnd;
      Min  = (HourEnd-hour)*60.0;
      min  = (HourEnd-hour)*60.0;
      sec  = (Min-min)*60.0;

      // Get the end time and compute subrun duration.
      SRtime->Set(yearEnd,monthEnd,dayEnd,hour,min,sec);
      
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
    TH1F *hped_per_tick_chan_20_temp = (TH1F*)file.FindObjectAny("hped_per_tick_chan_20");    
    if(hped_per_tick_chan_20_temp != 0){
      hped_per_tick_chan_20->Add(hped_per_tick_chan_20_temp,1.0);
    }    
    if(hped_per_tick_chan_20_temp != 0 && header != 0) {
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
    if(hped_per_tick_chan_548_temp != 0 && header != 0) {
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
    if(hped_per_tick_chan_1297_temp != 0 && header != 0) {
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
    if(hped_per_tick_chan_1697_temp != 0 && header != 0) {
      MeanPedChan1697Time[MeanPedChan1697Count] = Xsrtime;
      MeanPedChan1697    [MeanPedChan1697Count] = hped_per_tick_chan_1697_temp->GetMean(1);
      MeanPedChan1697Count++;      

      RMSPedChan1697Time[RMSPedChan1697Count] = Xsrtime;
      RMSPedChan1697    [RMSPedChan1697Count] = hped_per_tick_chan_1697_temp->GetRMS(1);
      RMSPedChan1697Count++;      
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
  sprintf(filename,"%s/ADCSpecChan0020_%.3u_days.png",PLOT_DIR.c_str(),Ndays);
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
  sprintf(filename,"%s/MeanPedChan0020_%.3u_days.png",PLOT_DIR.c_str(),Ndays);
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
  sprintf(filename,"%s/RMSPedChan0020_%.3u_days.png",PLOT_DIR.c_str(),Ndays);
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
  sprintf(filename,"%s/ADCSpecChan0548_%.3u_days.png",PLOT_DIR.c_str(),Ndays);
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
  sprintf(filename,"%s/MeanPedChan0548_%.3u_days.png",PLOT_DIR.c_str(),Ndays);
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
  sprintf(filename,"%s/RMSPedChan0548_%.3u_days.png",PLOT_DIR.c_str(),Ndays);
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
  sprintf(filename,"%s/ADCSpecChan1297_%.3u_days.png",PLOT_DIR.c_str(),Ndays);
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
  sprintf(filename,"%s/MeanPedChan1297_%.3u_days.png",PLOT_DIR.c_str(),Ndays);
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
  sprintf(filename,"%s/RMSPedChan1297_%.3u_days.png",PLOT_DIR.c_str(),Ndays);
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
  sprintf(filename,"%s/ADCSpecChan1697_%.3u_days.png",PLOT_DIR.c_str(),Ndays);
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
  sprintf(filename,"%s/MeanPedChan1697_%.3u_days.png",PLOT_DIR.c_str(),Ndays);
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
  sprintf(filename,"%s/RMSPedChan1697_%.3u_days.png",PLOT_DIR.c_str(),Ndays);
  cRMSPedChan1697->Print(filename);



  // Done with everything...
  return 0;
  
} // End script.
