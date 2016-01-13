#include "NearlinePlotMaker.h"

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
  sprintf(filelist_title,"/dune/app/users/mbaird42/35t_nearline_releases/temp/FileList.txt");
  inFile.open(filelist_title);



  // Book histos...
  // FIX: add spectra histos...
  TH1F *hped_per_event_chan_0;


  
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

    
    
    // FIX: make the spectra plots... (book the histo above...)
    TH1F *hped_per_event_chan_0_temp = (TH1F*)file.FindObjectAny("hped_per_event_chan_0");    
    if(hped_per_event_chan_0_temp != 0){
      // hped_per_event_chan_0->Add(hped_per_event_chan_0_temp,1.0);
    }    
    if(hped_per_event_chan_0_temp != 0 && header != 0) {
      MeanPedChan0000Time[MeanPedChan0000Count] = Xsrtime;
      MeanPedChan0000    [MeanPedChan0000Count] = hped_per_event_chan_0_temp->GetMean(1);
      MeanPedChan0000Count++;      
    }
    
    // FIX: make the spectra plots... (book the histo above...)
    TH1F *hped_per_event_chan_128_temp = (TH1F*)file.FindObjectAny("hped_per_event_chan_128");    
    if(hped_per_event_chan_128_temp != 0){
      // hped_per_event_chan_128->Add(hped_per_event_chan_128_temp,1.0);
    }    
    if(hped_per_event_chan_128_temp != 0 && header != 0) {
      MeanPedChan0128Time[MeanPedChan0128Count] = Xsrtime;
      MeanPedChan0128    [MeanPedChan0128Count] = hped_per_event_chan_128_temp->GetMean(1);
      MeanPedChan0128Count++;      
    }
    
    // FIX: make the spectra plots... (book the histo above...)
    TH1F *hped_per_event_chan_256_temp = (TH1F*)file.FindObjectAny("hped_per_event_chan_256");    
    if(hped_per_event_chan_256_temp != 0){
      // hped_per_event_chan_256->Add(hped_per_event_chan_256_temp,1.0);
    }    
    if(hped_per_event_chan_256_temp != 0 && header != 0) {
      MeanPedChan0256Time[MeanPedChan0256Count] = Xsrtime;
      MeanPedChan0256    [MeanPedChan0256Count] = hped_per_event_chan_256_temp->GetMean(1);
      MeanPedChan0256Count++;      
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

    












  // FIX: example of making a spectrum plot
  /*
  TCanvas *cnSliceday = new TCanvas("cnSliceday","# Slice Day",1200,800);
  cnSliceday->cd();
  cnSliceday->SetLogy();
  gStyle->SetOptStat(111111);
  hNSliceDay->SetAxisRange(0.0,50.0,"X");
  hNSliceDay->SetLineWidth(2);
  hNSliceDay->SetLineColor(kRed);
  hNSliceDay->Draw();
  UpdateText->Draw();

  sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1uNSlice%s.png",det.c_str(),trig.c_str(),p,period.c_str());
  cnSliceday->Print(filename);
  */

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
  sprintf(filename,"/dune/app/users/mbaird42/35t_nearline_releases/MeanPedChan0000.png");
  cMeanPedChan0000->Print(filename);



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
  sprintf(filename,"/dune/app/users/mbaird42/35t_nearline_releases/MeanPedChan0128.png");
  cMeanPedChan0128->Print(filename);



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
  sprintf(filename,"/dune/app/users/mbaird42/35t_nearline_releases/MeanPedChan0256.png");
  cMeanPedChan0256->Print(filename);



  // Done with everything...
  return 0;
  
} // End script.
