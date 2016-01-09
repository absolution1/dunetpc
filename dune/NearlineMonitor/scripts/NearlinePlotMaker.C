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
  sprintf(filelist_title,"35t_%uDay_Nearline_File_List.txt",Ndays);
  inFile.open(filelist_title);



  // Book histos...
  // FIX:
  // TH1F *blah...;


  
  // Declare variables used to make the TGraphs...
  float *MeanADCChannel00      = new float;
  float *MeanADCChannel00time  = new float;
  int    MeanADCChennel00Count = 0;



  /*
  // FIX this later...
  string time;
  if(period == "Day") time = "past 24 hrs.";
  if(period == "Week") time = "past 7 days";
  if(period == "Month") time = "past 30 days";
  */

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
      header->SetBranchAddress("Nevents",    &nevent);
      header->GetEntry(0);
      
      hour = Hour;
      Min  = (Hour-hour)*60.0;
      min  = (Hour-hour)*60.0;
      sec  = (Min-min)*60.0;

      // Ana/DAQ keeps track of GMT time, so convert to local time
      SRtime->Set(year,month,day,hour,min,sec);
      Xsrtime = SRtime->Convert() - GMToffset;
      cout << "time " << year << " " << day << " " << Hour << endl;
      hour = HourEnd;
      Min  = (HourEnd-hour)*60.0;
      min  = (HourEnd-hour)*60.0;
      sec  = (Min-min)*60.0;

      // Get the end time and compute subrun duration.
      SRtime->Set(yearEnd,monthEnd,dayEnd,hour,min,sec);
      
      if(Xsrtime < time_ago) continue; 
      
      if(run >= Last_Run) {
	if(run > Last_Run) {
	  Last_Run = run;
	  Last_SR  = subrun;
	}
	else if(subrun > Last_SR) {
	  Last_Run = run;
	  Last_SR  = subrun;
	}
      }
    
    }

    

    // FIX: Example of lifting out a histo and recording the TGraph info.

    TH1F *hnSlicetemp = (TH1F*)file.FindObjectAny("fNumSlices");
    if(hnSlicetemp != 0){
            hNSliceDay->Add(hnSlicetemp,1.0);
    }
    if(hnSlicetemp != 0 && header !=0) {
      NSlicetimeDay[NSliceDayCount] = Xsrtime;
      NSliceDay    [NSliceDayCount] = hnSlicetemp->GetMean(1);
      NSliceDayCount++;      
    }




    file.Close();

  } // end while loop



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
  sprintf(buff2, "Last run / subrun:   %i / %i", Last_Run, Last_SR);
  UpdateText->AddText(buff2);

    

  // FIX: example of making a spectrum plot
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



  // FIX: example of TGraph version of the plot
  maxtime = 0;
  max = 0.0;
  ave = 0.0;
  for(int i = 0; i < NSliceDayCount; ++i) {
    ave += (double)NSliceDay[i];
    if(NSlicetimeDay[i] > maxtime) {
      maxtime = NSlicetimeDay[i];
      max     = NSliceDay[i];
    }
  }
  if(NSliceDayCount > 0) ave = ave/(double)NSliceDayCount;
  sprintf(lptext,"Last Point = %f  /  Average = %f",max,ave);
  LastPoint->Clear();
  LastPoint->AddText(lptext);
  TCanvas *ctnSliceday = new TCanvas("cnSliceday","# Non-Noise Slices day",1200,800);
  ctnSliceday->cd();
  gPad->SetGridx();
  TGraph *gnSliceday = new TGraph(NSliceDayCount,NSlicetimeDay,NSliceDay);
  sprintf(title,"Number of Slices per Subrun");
  gnSliceday->SetTitle(title);
  gnSliceday->SetMarkerColor(kBlue);
  gnSliceday->GetXaxis()->SetTimeDisplay(1);
  gnSliceday->GetXaxis()->SetLabelSize(0.03);
  gnSliceday->GetXaxis()->SetTimeFormat(taxis_labels.c_str());
  gnSliceday->GetXaxis()->SetLimits(time_ago,XNow);
  gnSliceday->GetXaxis()->SetTitle("(central time)");
  gnSliceday->Draw("A*");
  UpdateText->Draw();
  LastPoint->Draw();
  sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1unSliceperSR%s.png",det.c_str(),trig.c_str(),p,period.c_str());
  ctnSliceday->Print(filename);



  // Done with everything...
  return 0;
  
} // End script.
