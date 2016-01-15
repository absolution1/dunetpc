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
  TH1F *hped_per_tick_chan_0    = new TH1F("hped_per_tick_chan_0","ADC per Tick - Channel 0;ADC",100,0,2048);
  TH1F *hped_per_tick_chan_128  = new TH1F("hped_per_tick_chan_128","ADC per Tick - Channel 128;ADC",100,0,2048);
  TH1F *hped_per_tick_chan_256  = new TH1F("hped_per_tick_chan_256","ADC per Tick - Channel 256;ADC",100,0,2048);
  TH1F *hped_per_tick_chan_384  = new TH1F("hped_per_tick_chan_384","ADC per Tick - Channel 384;ADC",100,0,2048);


  
  // Declare variables used to make the TGraphs...
  float *MeanPedChan0000      = new float[Npoint];
  float *MeanPedChan0000Time  = new float[Npoint];
  int    MeanPedChan0000Count = 0;

  float *MeanPedChan0128      = new float[Npoint];
  float *MeanPedChan0128Time  = new float[Npoint];
  int    MeanPedChan0128Count = 0;

  float *MeanPedChan0256      = new float[Npoint];
  float *MeanPedChan0256Time  = new float[Npoint];
  int    MeanPedChan0256Count = 0;

  float *MeanPedChan0384      = new float[Npoint];
  float *MeanPedChan0384Time  = new float[Npoint];
  int    MeanPedChan0384Count = 0;



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
    TH1F *hped_per_tick_chan_0_temp = (TH1F*)file.FindObjectAny("hped_per_tick_chan_0");    
    if(hped_per_tick_chan_0_temp != 0){
      hped_per_tick_chan_0->Add(hped_per_tick_chan_0_temp,1.0);
    }    
    if(hped_per_tick_chan_0_temp != 0 && header != 0) {
      MeanPedChan0000Time[MeanPedChan0000Count] = Xsrtime;
      MeanPedChan0000    [MeanPedChan0000Count] = hped_per_tick_chan_0_temp->GetMean(1);
      MeanPedChan0000Count++;      
    }
    
    TH1F *hped_per_tick_chan_128_temp = (TH1F*)file.FindObjectAny("hped_per_tick_chan_128");    
    if(hped_per_tick_chan_128_temp != 0){
      hped_per_tick_chan_128->Add(hped_per_tick_chan_128_temp,1.0);
    }    
    if(hped_per_tick_chan_128_temp != 0 && header != 0) {
      MeanPedChan0128Time[MeanPedChan0128Count] = Xsrtime;
      MeanPedChan0128    [MeanPedChan0128Count] = hped_per_tick_chan_128_temp->GetMean(1);
      MeanPedChan0128Count++;      
    }
    
    TH1F *hped_per_tick_chan_256_temp = (TH1F*)file.FindObjectAny("hped_per_tick_chan_256");    
    if(hped_per_tick_chan_256_temp != 0){
      hped_per_tick_chan_256->Add(hped_per_tick_chan_256_temp,1.0);
    }    
    if(hped_per_tick_chan_256_temp != 0 && header != 0) {
      MeanPedChan0256Time[MeanPedChan0256Count] = Xsrtime;
      MeanPedChan0256    [MeanPedChan0256Count] = hped_per_tick_chan_256_temp->GetMean(1);
      MeanPedChan0256Count++;      
    }
    
    TH1F *hped_per_tick_chan_384_temp = (TH1F*)file.FindObjectAny("hped_per_tick_chan_384");    
    if(hped_per_tick_chan_384_temp != 0){
      hped_per_tick_chan_384->Add(hped_per_tick_chan_384_temp,1.0);
    }    
    if(hped_per_tick_chan_384_temp != 0 && header != 0) {
      MeanPedChan0384Time[MeanPedChan0384Count] = Xsrtime;
      MeanPedChan0384    [MeanPedChan0384Count] = hped_per_tick_chan_384_temp->GetMean(1);
      MeanPedChan0384Count++;      
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

    












  TCanvas *cADCSpecChan0000 = new TCanvas("cADCSpecChan0000","ADC Spectrum - Channel 0",1200,800);
  cADCSpecChan0000->cd();
  cADCSpecChan0000->SetLogy();
  gStyle->SetOptStat(111111);
  // hped_per_tick_chan_0->SetAxisRange(0.0,2048.0,"X");
  hped_per_tick_chan_0->SetLineWidth(2);
  hped_per_tick_chan_0->SetLineColor(kRed);
  hped_per_tick_chan_0->Draw();
  UpdateText->Draw();
  sprintf(filename,"%s/ADCSpecChan0000_%.3u_days.png",PLOT_DIR.c_str(),Ndays);
  cADCSpecChan0000->Print(filename);

  maxtime = 0;
  max = 0.0;
  ave = 0.0;
  for(int i = 0; i < MeanPedChan0000Count; ++i) {
    ave += (double)MeanPedChan0000[i];
    if(MeanPedChan0000Time[i] > maxtime) {
      maxtime = MeanPedChan0000Time[i];
      max     = MeanPedChan0000[i];
    }
  }
  if(MeanPedChan0000Count > 0) ave = ave/(double)MeanPedChan0000Count;
  sprintf(lptext,"Last Point = %f  /  Average = %f",max,ave);
  LastPoint->Clear();
  LastPoint->AddText(lptext);
  TCanvas *cMeanPedChan0000 = new TCanvas("cMeanPedChan0000","Mean Pedestal for Channel 0",1200,800);
  cMeanPedChan0000->cd();
  gPad->SetGridx();
  TGraph *gMeanPedChan0000 = new TGraph(MeanPedChan0000Count,MeanPedChan0000Time,MeanPedChan0000);
  sprintf(title,"Mean Pedestal per Event for Channel 0");
  gMeanPedChan0000->SetTitle(title);
  gMeanPedChan0000->SetMarkerColor(kBlue);
  gMeanPedChan0000->GetXaxis()->SetTimeDisplay(1);
  gMeanPedChan0000->GetXaxis()->SetLabelSize(0.03);
  gMeanPedChan0000->GetXaxis()->SetTimeFormat(taxis_labels.c_str());
  gMeanPedChan0000->GetXaxis()->SetLimits(time_ago,XNow);
  gMeanPedChan0000->GetXaxis()->SetTitle("(central time)");
  gMeanPedChan0000->Draw("A*");
  UpdateText->Draw();
  LastPoint->Draw();
  sprintf(filename,"%s/MeanPedChan0000_%.3u_days.png",PLOT_DIR.c_str(),Ndays);
  cMeanPedChan0000->Print(filename);



  TCanvas *cADCSpecChan0128 = new TCanvas("cADCSpecChan0128","ADC Spectrum - Channel 128",1200,800);
  cADCSpecChan0128->cd();
  cADCSpecChan0128->SetLogy();
  gStyle->SetOptStat(111111);
  // hped_per_tick_chan_128->SetAxisRange(0.0,2048.0,"X");
  hped_per_tick_chan_128->SetLineWidth(2);
  hped_per_tick_chan_128->SetLineColor(kRed);
  hped_per_tick_chan_128->Draw();
  UpdateText->Draw();
  sprintf(filename,"%s/ADCSpecChan0128_%.3u_days.png",PLOT_DIR.c_str(),Ndays);
  cADCSpecChan0128->Print(filename);

  maxtime = 0;
  max = 0.0;
  ave = 0.0;
  for(int i = 0; i < MeanPedChan0128Count; ++i) {
    ave += (double)MeanPedChan0128[i];
    if(MeanPedChan0128Time[i] > maxtime) {
      maxtime = MeanPedChan0128Time[i];
      max     = MeanPedChan0128[i];
    }
  }
  if(MeanPedChan0128Count > 0) ave = ave/(double)MeanPedChan0128Count;
  sprintf(lptext,"Last Point = %f  /  Average = %f",max,ave);
  LastPoint->Clear();
  LastPoint->AddText(lptext);
  TCanvas *cMeanPedChan0128 = new TCanvas("cMeanPedChan0128","Mean Pedestal for Channel 128",1200,800);
  cMeanPedChan0128->cd();
  gPad->SetGridx();
  TGraph *gMeanPedChan0128 = new TGraph(MeanPedChan0128Count,MeanPedChan0128Time,MeanPedChan0128);
  sprintf(title,"Mean Pedestal per Event for Channel 128");
  gMeanPedChan0128->SetTitle(title);
  gMeanPedChan0128->SetMarkerColor(kBlue);
  gMeanPedChan0128->GetXaxis()->SetTimeDisplay(1);
  gMeanPedChan0128->GetXaxis()->SetLabelSize(0.03);
  gMeanPedChan0128->GetXaxis()->SetTimeFormat(taxis_labels.c_str());
  gMeanPedChan0128->GetXaxis()->SetLimits(time_ago,XNow);
  gMeanPedChan0128->GetXaxis()->SetTitle("(central time)");
  gMeanPedChan0128->Draw("A*");
  UpdateText->Draw();
  LastPoint->Draw();
  sprintf(filename,"%s/MeanPedChan0128_%.3u_days.png",PLOT_DIR.c_str(),Ndays);
  cMeanPedChan0128->Print(filename);



  TCanvas *cADCSpecChan0256 = new TCanvas("cADCSpecChan0256","ADC Spectrum - Channel 256",1200,800);
  cADCSpecChan0256->cd();
  cADCSpecChan0256->SetLogy();
  gStyle->SetOptStat(111111);
  // hped_per_tick_chan_256->SetAxisRange(0.0,2048.0,"X");
  hped_per_tick_chan_256->SetLineWidth(2);
  hped_per_tick_chan_256->SetLineColor(kRed);
  hped_per_tick_chan_256->Draw();
  UpdateText->Draw();
  sprintf(filename,"%s/ADCSpecChan0256_%.3u_days.png",PLOT_DIR.c_str(),Ndays);
  cADCSpecChan0256->Print(filename);

  maxtime = 0;
  max = 0.0;
  ave = 0.0;
  for(int i = 0; i < MeanPedChan0256Count; ++i) {
    ave += (double)MeanPedChan0256[i];
    if(MeanPedChan0256Time[i] > maxtime) {
      maxtime = MeanPedChan0256Time[i];
      max     = MeanPedChan0256[i];
    }
  }
  if(MeanPedChan0256Count > 0) ave = ave/(double)MeanPedChan0256Count;
  sprintf(lptext,"Last Point = %f  /  Average = %f",max,ave);
  LastPoint->Clear();
  LastPoint->AddText(lptext);
  TCanvas *cMeanPedChan0256 = new TCanvas("cMeanPedChan0256","Mean Pedestal for Channel 256",1200,800);
  cMeanPedChan0256->cd();
  gPad->SetGridx();
  TGraph *gMeanPedChan0256 = new TGraph(MeanPedChan0256Count,MeanPedChan0256Time,MeanPedChan0256);
  sprintf(title,"Mean Pedestal per Event for Channel 256");
  gMeanPedChan0256->SetTitle(title);
  gMeanPedChan0256->SetMarkerColor(kBlue);
  gMeanPedChan0256->GetXaxis()->SetTimeDisplay(1);
  gMeanPedChan0256->GetXaxis()->SetLabelSize(0.03);
  gMeanPedChan0256->GetXaxis()->SetTimeFormat(taxis_labels.c_str());
  gMeanPedChan0256->GetXaxis()->SetLimits(time_ago,XNow);
  gMeanPedChan0256->GetXaxis()->SetTitle("(central time)");
  gMeanPedChan0256->Draw("A*");
  UpdateText->Draw();
  LastPoint->Draw();
  sprintf(filename,"%s/MeanPedChan0256_%.3u_days.png",PLOT_DIR.c_str(),Ndays);
  cMeanPedChan0256->Print(filename);



  TCanvas *cADCSpecChan0384 = new TCanvas("cADCSpecChan0384","ADC Spectrum - Channel 384",1200,800);
  cADCSpecChan0384->cd();
  cADCSpecChan0384->SetLogy();
  gStyle->SetOptStat(111111);
  // hped_per_tick_chan_384->SetAxisRange(0.0,2048.0,"X");
  hped_per_tick_chan_384->SetLineWidth(2);
  hped_per_tick_chan_384->SetLineColor(kRed);
  hped_per_tick_chan_384->Draw();
  UpdateText->Draw();
  sprintf(filename,"%s/ADCSpecChan0384_%.3u_days.png",PLOT_DIR.c_str(),Ndays);
  cADCSpecChan0384->Print(filename);

  maxtime = 0;
  max = 0.0;
  ave = 0.0;
  for(int i = 0; i < MeanPedChan0384Count; ++i) {
    ave += (double)MeanPedChan0384[i];
    if(MeanPedChan0384Time[i] > maxtime) {
      maxtime = MeanPedChan0384Time[i];
      max     = MeanPedChan0384[i];
    }
  }
  if(MeanPedChan0384Count > 0) ave = ave/(double)MeanPedChan0384Count;
  sprintf(lptext,"Last Point = %f  /  Average = %f",max,ave);
  LastPoint->Clear();
  LastPoint->AddText(lptext);
  TCanvas *cMeanPedChan0384 = new TCanvas("cMeanPedChan0384","Mean Pedestal for Channel 384",1200,800);
  cMeanPedChan0384->cd();
  gPad->SetGridx();
  TGraph *gMeanPedChan0384 = new TGraph(MeanPedChan0384Count,MeanPedChan0384Time,MeanPedChan0384);
  sprintf(title,"Mean Pedestal per Event for Channel 384");
  gMeanPedChan0384->SetTitle(title);
  gMeanPedChan0384->SetMarkerColor(kBlue);
  gMeanPedChan0384->GetXaxis()->SetTimeDisplay(1);
  gMeanPedChan0384->GetXaxis()->SetLabelSize(0.03);
  gMeanPedChan0384->GetXaxis()->SetTimeFormat(taxis_labels.c_str());
  gMeanPedChan0384->GetXaxis()->SetLimits(time_ago,XNow);
  gMeanPedChan0384->GetXaxis()->SetTitle("(central time)");
  gMeanPedChan0384->Draw("A*");
  UpdateText->Draw();
  LastPoint->Draw();
  sprintf(filename,"%s/MeanPedChan0384_%.3u_days.png",PLOT_DIR.c_str(),Ndays);
  cMeanPedChan0384->Print(filename);



  // Done with everything...
  return 0;
  
} // End script.
