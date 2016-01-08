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


  // Define the number of partitions in the data.  Make this number one more
  // than the actual number so that there can be a default unassigned value.
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
  // TH1F *blah...;


  
  // Declare variables used to make the TGraphs...
  float *MeanADCChannel00      = new float;
  float *MeanADCChannel00time  = new float;
  int    MeanADCChennel00Count = 0;



  /*
  // fix this later...
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
      
      // MICHAEL

      // Storing last run and sub run for each partition...tis a test.
      /*if(par==0){
	Last_Run[0] = LastRun;
	Last_SR[0] = LastSR;
      }
      if(par==1){
	Last_Run[1] = LastRun;
	Last_SR[1] = LastSR;
      }
      if(par==2){
	Last_Run[2] = LastRun;
	Last_SR[2] = LastSR;
      }
      if(par==3){
	Last_Run[3] = LastRun;
	Last_SR[3] = LastSR;
      }
      if(par==4){
	Last_Run[4] = LastRun;
	Last_SR[4] = LastSR;
	}*/
       
    }

    if(t_prd == "Day") std::cout << "\nDAY LOOP:   adding run " << run << ", subrun " << subrun << endl;
    if(t_prd == "Week") std::cout << "\nWEEK LOOP:   adding run " << run << ", subrun " << subrun << endl; 
    if(t_prd == "Month") std::cout << "\nMONTH LOOP:   adding run " << run << ", subrun " << subrun << endl;
    // Find DCM sync histograms from reco_hist root files      
    TH2F *hit   = (TH2F*)fileD.FindObjectAny("H3DHitOnTrack"); 
    TH2F *instr = (TH2F*)fileD.FindObjectAny("instrumented");   //MAP OF INSTRUMENTED DETECTOR (i.e. DCM WITH > 50 APDS)
    
    TH2F *trk3  = (TH2F*)fileD.FindObjectAny("H3DTrackNum");      
    TH2F *trk2  = (TH2F*)fileD.FindObjectAny("H2DTrackNum");                
    
    evoindex=0.;    
    fullyinstrdib[par] = 0;
    diblocklast[par] = 0;
    

    if((hit !=0)){      
      /*
	if(!hit){                       // TEMPORAL??
	std::cerr<<" There is no Hits from 3D Track Histogram!! Skipping...  "<<std::endl;
	continue;
	}
	else{
      */
      
      //NUEVO: xdiblock<=N_diblock
      for(unsigned int xdiblock =1; xdiblock <=14; xdiblock++){		
	sum_hits_col[xdiblock]=0.;
	
	for(unsigned int yldc =1; yldc <=N_dcms; yldc++){	  
	  sum_hits_col[xdiblock] += hit->GetBinContent(xdiblock,yldc);
	}	
	if(sum_hits_col[xdiblock] > 0){
	  ++fullyinstrdib[par];
	  diblocklast[par] = xdiblock; /// TEST	  
	}
      } //for xdiblock
      
    } //else
    

    //tr3Dfull[par][count[par]]=0.;
    //tr2Dfull[par][count[par]]=0.;
    
    // IMPORTANT: MUST GUARATNEE LIVETIME IS NEVER ZERO
    if(livetime == 0) livetime = 550000*1.0E-9*nevents;
        
    int init_dib = diblocklast[par] - fullyinstrdib[par] +1; // counts 1st db of partition

    MinDCM  = 99999.99;
    
    // NUEVO   
        
    // The negative instead...
    
    bool noInstr = false;
    
    if((!instr) || instr->GetEntries() == 0){
      std::cout<<"==============================================================================================="<<std::endl;
      std::cout<<"==============================================================================================="<<std::endl;
      std::cout<<"                                                                                               "<<std::endl;
      std::cout<<"                                                                                               "<<std::endl;
      std::cout<<" NearlineAna files have empty OR missing Instrumented Histograms!! processing will continue... "<<std::endl;
      std::cout<<"                                                                                               "<<std::endl;
      std::cout<<"                                                                                               "<<std::endl;
      std::cout<<"==============================================================================================="<<std::endl;
      std::cout<<"==============================================================================================="<<std::endl;
      
      noInstr = true;
    }

    // Dealing with funky instr histos
    if(isND && (instr->GetEntries() > 14)) noInstr = true;
    //end NUEVO
    
    
    for(int xbin= init_dib; xbin <= diblocklast[par]; xbin++){
      
	sum3dtrk_col[par][xbin]=0.;
	sum2dtrk_col[par][xbin]=0.;

	//NUEVO

	//NUEVO: ybin<=N_dcms
	for(int ybin=1; ybin<=N_dcms; ++ybin){               //counts for each dcm
	  	  
	  if(noInstr){

            if(isND){
	      
              // Last two dcms on db 4 are not instrumented in ND	      
              if((xbin==4 && ybin==3) ||(xbin ==4 && ybin ==4)) continue;

              evoindex=float(N_dcms)*(xbin-1.)+ybin-1.;

              //GETTING: PER DIBLOCK SUM OF TRACKS	      
              sum3dtrk_col[par][xbin] += trk3->GetBinContent(xbin,ybin);    //get the sum of 3D tracks for each diblock	      
              sum2dtrk_col[par][xbin] += trk2->GetBinContent(xbin,ybin);    //get the sum of 2D tracks for each diblock	      

              if((hit->GetBinContent(xbin, ybin))/(livetime) < MinDCM){
                MinDCM = (hit->GetBinContent(xbin, ybin))/(livetime);                   // rate = value/livetime   		
              }

              else{
                continue;
              }
	      
	      //TIME EVOLUTION-STRIP PLOT	      
              hNormalization[par]->Fill(Xsrtime);     //NORMALIZATION FOR STRIP PLOT                                                        	      
              tr3devo = hit->GetBinContent(hit->FindBin(xbin,ybin))/(livetime);    //Hits from 3D Tracks                                    	      
              //if(tr3devo > 2.5*MinDCM ){                                     // IS RELATED TO Min3DRate PLOT OR WHERE GAP BEGINS  	      
              if(tr3devo < lowratemin ){
                g3DRateEvo[par]->Fill(Xsrtime,evoindex,tr3devo);
              }

            }//isND
	    else if(isFD){

              evoindex=float(N_dcms)*(xbin-1.)+ybin-1.;
              //GETTING: PER DIBLOCK SUM OF TRACKS 
              sum3dtrk_col[par][xbin] += trk3->GetBinContent(xbin,ybin);    //get the sum of 3D tracks for each diblock 	      
              sum2dtrk_col[par][xbin] += trk2->GetBinContent(xbin,ybin);    //get the sum of 2D tracks for each diblock 	      
	      
              if((hit->GetBinContent(xbin, ybin))/(livetime) < MinDCM){		

                MinDCM = (hit->GetBinContent(xbin, ybin))/(livetime);                   // rate = value/livetime  		

              }
              else{
                continue;
              }
	      
              //TIME EVOLUTION-STRIP PLOT   	      
              hNormalization[par]->Fill(Xsrtime);     //NORMALIZATION FOR STRIP PLOT  	      
              tr3devo = hit->GetBinContent(hit->FindBin(xbin,ybin))/(livetime);    //Hits from 3D Tracks  	      
              //if(tr3devo > 2.5*MinDCM ){                                     // IS RELATED TO Min3DRate PLOT OR WHERE GAP BEGINS 	      

              if(tr3devo < lowratemin ){
                g3DRateEvo[par]->Fill(Xsrtime,evoindex,tr3devo);
              }

            }//end isFD
	    else{
	      continue;
	    }	    	    
	  }//end noInstr    //end NUEVO
	  else if(instr-> GetBinContent(xbin,ybin) > 0){   //Get those DCMs with +50 APDs only and neglect the rest
	    //NUEVO: float(N_dcms) 
	    evoindex=float(N_dcms)*(xbin-1.)+ybin-1.;          
	    	    	    
	    sum3dtrk_col[par][xbin] += trk3->GetBinContent(xbin,ybin);	//get the sum of 3D tracks for each diblock  
	    sum2dtrk_col[par][xbin] += trk2->GetBinContent(xbin,ybin);	//get the sum of 2D tracks for each diblock  
	    
	    //PER DCM RATE OF TRACKS
	    //if(t_prd!="Month"){
	    //tr3D[par][xbin][ybin][count[par]] = trk3->GetBinContent(xbin,ybin)/(livetime);	  
	    //tr2D[par][xbin][ybin][count[par]] = trk2->GetBinContent(xbin,ybin)/(livetime);
	    //}
	    ////////////////////////////////////////////////////////////////////////////
	    
	    //TIME EVOLUTION-STRIP PLOT
	    hNormalization[par]->Fill(Xsrtime);                                  //NORMALIZATION FOR STRIP PLOT
	    tr3devo = hit->GetBinContent(hit->FindBin(xbin,ybin))/(livetime);    //Hits from 3D Tracks
	  	  	  	  
	    //PLOTTING THE NEGATIVE: IF LOW RATE THEN PRINT IT
	    if(tr3devo < lowratemin ){	                                 // IS RELATED TO Min3DRate PLOT OR WHERE GAP BEGINS 
	      g3DRateEvo[par]->Fill(Xsrtime,evoindex,tr3devo);  
	    }   	    
	    //===========================================================================================	    
	    
	    //MINIMUM HISTOGRAM
	    if((hit->GetBinContent(xbin, ybin))/(livetime) < MinDCM){	      
	      
	      MinDCM = (hit->GetBinContent(xbin, ybin))/(livetime);                   // rate = value/livetime

	    }
	    else{
	      continue;
	    }
    	    	    
	    //===========================================================================================
	  }//end if/else
	}//end for        //counter here gives 168 = total no. dcms in FarDet       
    }//end for        //gives 14 = no. of diblocks
    //}

    minDCM[par][count[par]] = MinDCM;    
    xtime[par][count[par]] = Xsrtime;
    count[par]++;
    


    TH1F *hnSlicetemp = (TH1F*)fileD.FindObjectAny("fNumSlices");
    if(hnSlicetemp != 0){
            hNSliceDay[par]->Add(hnSlicetemp,1.0);
    }
    if(hnSlicetemp != 0 && header !=0) {
      NSlicetimeDay[par][NSliceDayCount[par]] = Xsrtime;
      NSliceDay    [par][NSliceDayCount[par]] = hnSlicetemp->GetMean(1);
      NSliceDayCount[par]++;      
    }
    //new beam stuff stuff
   TH1F *hTimeDifftemp = (TH1F*)fileD.FindObjectAny("fTimeDiffNanoSec");
    if(hTimeDifftemp != 0){
             hTimeDiffNanoSecDay[par]->Add(hTimeDifftemp,1.0);
    }
    if(hTimeDifftemp != 0 && header !=0) {
      NTimeDiffNanoSectimeDay[par][NTimeDiffDayCount[par]] = Xsrtime;
      NTimeDiffNanoSecDay    [par][NTimeDiffDayCount[par]] = hTimeDifftemp->GetMean(1);
      NTimeDiffDayCount[par]++;      
    }
   
    TH1F *hPOTSumtemp = (TH1F*)fileD.FindObjectAny("fPOTSum");
    if(hPOTSumtemp != 0){
             hPOTSumDay[par]->Add(hPOTSumtemp,1.0);
    }
    if(hPOTSumtemp != 0 && header !=0) {
      NPOTSumtimeDay[par][NPOTSumDayCount[par]] = Xsrtime;
      NPOTSumDay    [par][NPOTSumDayCount[par]] = hPOTSumtemp->GetBinContent(1);//we want the value in the first bin (this is the POT)
      NPOTSumDayCount[par]++;      
    }

    TH1F *hspillPOTtemp = (TH1F*)fileD.FindObjectAny("fspillPOT");
    if(hspillPOTtemp != 0){
             hspillPOTDay[par]->Add(hspillPOTtemp,1.0);
    }
    if(hspillPOTtemp != 0 && header !=0) {
      NspillPOTtimeDay[par][NspillPOTDayCount[par]] = Xsrtime;
      NspillPOTDay    [par][NspillPOTDayCount[par]] = hspillPOTtemp->GetMean(1);
      NspillPOTDayCount[par]++;      
    }
    
    TH1F *hXPositiontemp = (TH1F*)fileD.FindObjectAny("fXPosition");
    if(hXPositiontemp != 0){
             hXPositionDay[par]->Add(hXPositiontemp,1.0);
    }
    if(hXPositiontemp != 0 && header !=0) {
      NXPositiontimeDay[par][NXPositionDayCount[par]] = Xsrtime;
      NXPositionDay    [par][NXPositionDayCount[par]] = hXPositiontemp->GetMean(1);
      NXPositionDayCount[par]++;      
    }

    TH1F *hYPositiontemp = (TH1F*)fileD.FindObjectAny("fYPosition");
    if(hYPositiontemp != 0){
             hYPositionDay[par]->Add(hYPositiontemp,1.0);
    }
    if(hYPositiontemp != 0 && header !=0) {
      NYPositiontimeDay[par][NYPositionDayCount[par]] = Xsrtime;
      NYPositionDay    [par][NYPositionDayCount[par]] = hYPositiontemp->GetMean(1);
      NYPositionDayCount[par]++;      
    }
    
    TH1F *hXWidthtemp = (TH1F*)fileD.FindObjectAny("fXWidth");
    if(hXWidthtemp != 0){
             hXWidthDay[par]->Add(hXWidthtemp,1.0);
    }
    if(hXWidthtemp != 0 && header !=0) {
      NXWidthtimeDay[par][NXWidthDayCount[par]] = Xsrtime;
      NXWidthDay    [par][NXWidthDayCount[par]] = hXWidthtemp->GetMean(1);
      NXWidthDayCount[par]++;      
    }

    TH1F *hYWidthtemp = (TH1F*)fileD.FindObjectAny("fYWidth");
    if(hYWidthtemp != 0){
             hYWidthDay[par]->Add(hYWidthtemp,1.0);
    }
    if(hYWidthtemp != 0 && header !=0) {
      NYWidthtimeDay[par][NYWidthDayCount[par]] = Xsrtime;
      NYWidthDay    [par][NYWidthDayCount[par]] = hYWidthtemp->GetMean(1);
      NYWidthDayCount[par]++;      
    }

     
    TH1F *hGoodBeamtemp = (TH1F*)fileD.FindObjectAny("fGoodBeam");
     if(hGoodBeamtemp != 0){
              hGoodBeamDay[par]->Add(hGoodBeamtemp,1.0);
     }
     if(hGoodBeamtemp != 0 && header !=0) {
       NGoodBeamtimeDay[par][NGoodBeamDayCount[par]] = Xsrtime;
       NGoodBeamDay    [par][NGoodBeamDayCount[par]] = hGoodBeamtemp->GetMean(1);
       NGoodBeamDayCount[par]++;      
     }
     
    TH1F *hBadSpillstemp = (TH1F*)fileD.FindObjectAny("fBadSpills");
    if(hBadSpillstemp != 0){
             hBadSpillsDay[par]->Add(hBadSpillstemp,1.0);
    }
   
    TH1F *hHornCurrenttemp = (TH1F*)fileD.FindObjectAny("fHornCurrent");
    if(hHornCurrenttemp != 0){
             hHornCurrentDay[par]->Add(hHornCurrenttemp,1.0);
    }
    if(hHornCurrenttemp != 0 && header !=0) {
      NHornCurrenttimeDay[par][NHornCurrentDayCount[par]] = Xsrtime;
      NHornCurrentDay    [par][NHornCurrentDayCount[par]] = hHornCurrenttemp->GetMean(1);
      NHornCurrentDayCount[par]++;      
    }
    
    //end new beam plots
    TH1F *hnSliceHittemp = (TH1F*)fileD.FindObjectAny("fNumHits");
    if(hnSliceHittemp != 0){
       hNNonNoiseSliceHitDay[par]->Add(hnSliceHittemp,1.0);
    }
    if(hnSliceHittemp != 0 && header !=0) {
      NNonNoiseSliceHittimeDay[par][NNonNoiseSliceHitDayCount[par]] = Xsrtime;
      NNonNoiseSliceHitDay    [par][NNonNoiseSliceHitDayCount[par]] = hnSliceHittemp->GetMean(1);
      NNonNoiseSliceHitDayCount[par]++;      
    }

    TH1F *hnNoiseSliceHittemp = (TH1F*)fileD.FindObjectAny("fNumHitsNoise");
    if(hnNoiseSliceHittemp != 0){
      hNNoiseSliceHitDay[par]->Add(hnNoiseSliceHittemp,1.0);
    }
    if(hnNoiseSliceHittemp != 0 && header !=0) {
      NNoiseSliceHittimeDay[par][NNoiseSliceHitDayCount[par]] = Xsrtime;
      NNoiseSliceHitDay    [par][NNoiseSliceHitDayCount[par]] = hnNoiseSliceHittemp->GetMean(1);
      NNoiseSliceHitDayCount[par]++;      
    }


    TH1F *hNonNoiseSlicePEtemp = (TH1F*)fileD.FindObjectAny("fPENonNoiseHits");
    if(hNonNoiseSlicePEtemp != 0){
      hNonNoiseSlicePEDay[par]->Add(hNonNoiseSlicePEtemp,1.0);
    }

 
    TH1F *hNoiseSlicePEtemp = (TH1F*)fileD.FindObjectAny("fPENoiseHits");
    if(hNoiseSlicePEtemp != 0){
      hNoiseSlicePEDay[par]->Add(hNoiseSlicePEtemp,1.0);
    }

    TH1F *hTrackPEtemp = (TH1F*)fileD.FindObjectAny("fPETrackHitsAll3D");
    if(hTrackPEtemp != 0){
      hTrackPEDay[par]->Add(hTrackPEtemp,1.0);
    }
 
    TH1F *hTsdtemp = (TH1F*)fileD.FindObjectAny("fTsd");
    if(hTsdtemp != 0){
      hTsdDay[par]->Add(hTsdtemp,1.0);
    }
    if(hTsdtemp != 0 && header !=0) {
      NTsdtimeDay[par][NTsdDayCount[par]] = Xsrtime;
      NTsdDay    [par][NTsdDayCount[par]] = hTsdtemp->GetMean(1);
      NTsdDayCount[par]++;      
    }

    TH1F *hnTracktemp = (TH1F*)fileD.FindObjectAny("fNumTracksAll3D"); 
    if(hnTracktemp != 0){
      hnTrackDay[par]->Add(hnTracktemp,1.0);
    }
    if(hnTracktemp != 0 && header !=0) {
      NTracktimeDay[par][NTrackDayCount[par]] = Xsrtime;
      NTrackDay    [par][NTrackDayCount[par]] = hnTracktemp->GetMean(1);
      NTrackDayCount[par]++;      
    }

    TH1F *hTrackLentemp = (TH1F*)fileD.FindObjectAny("fTrackLengthAll3D");
    if(hTrackLentemp != 0){
      hTrackLenDay[par]->Add(hTrackLentemp,1.0);
    }
    if(hTrackLentemp != 0 && header !=0) {
      NTrackLentimeDay[par][NTrackLenDayCount[par]] = Xsrtime;
      NTrackLenDay    [par][NTrackLenDayCount[par]] = hTrackLentemp->GetMean(1);
      NTrackLenDayCount[par]++;      
    }
    
    TH1F *hCosNumitemp = (TH1F*)fileD.FindObjectAny("fTrackCosNumiAll3D");
    if(hCosNumitemp != 0){
      hCosNumiDay[par]->Add(hCosNumitemp,1.0);
    }
    
    TH2F *hPCXavetemp = (TH2F*)fileD.FindObjectAny("fPCXaveCR");
    if(hPCXavetemp != 0) {
      hPCXaveDay[par]->Add(hPCXavetemp,1.0);
    }

    TH2F *hPCYavetemp = (TH2F*)fileD.FindObjectAny("fPCYaveCR");
    if(hPCYavetemp != 0) {
      hPCYaveDay[par]->Add(hPCYavetemp,1.0);
    }
    
    TH1F *hnTrackAll3Dtemp = (TH1F*)fileD.FindObjectAny("fNumTracksAll3D");
    if(hnTrackAll3Dtemp != 0){
      hnTrackAll3DDay[par]->Add(hnTrackAll3Dtemp,1.0);
    }

    TH2F *hTrackStartXZAll3Dtemp = (TH2F*)fileD.FindObjectAny("fStartPointTrackXZCRAll3D");
    if(hTrackStartXZAll3Dtemp != 0) {
      hTrackStartXZAll3DDay[par]->Add(hTrackStartXZAll3Dtemp,1.0);
    }
    
    TH2F *hTrackStartYZAll3Dtemp = (TH2F*)fileD.FindObjectAny("fStartPointTrackYZCRAll3D");
    if(hTrackStartYZAll3Dtemp != 0) {
      hTrackStartYZAll3DDay[par]->Add(hTrackStartYZAll3Dtemp,1.0);
    }

    TH2F *hTrackStopXZAll3Dtemp = (TH2F*)fileD.FindObjectAny("fStopPointTrackXZCRAll3D");
    if(hTrackStopXZAll3Dtemp != 0) {
      hTrackStopXZAll3DDay[par]->Add(hTrackStopXZAll3Dtemp,1.0);
    }

    TH2F *hTrackStopYZAll3Dtemp = (TH2F*)fileD.FindObjectAny("fStopPointTrackYZCRAll3D");
    if(hTrackStopYZAll3Dtemp != 0) {
      hTrackStopYZAll3DDay[par]->Add(hTrackStopYZAll3Dtemp,1.0);
    }

    TH1F *hTrackLengthAll3Dtemp = (TH1F*)fileD.FindObjectAny("fTrackLengthAll3D");
    if(hTrackLengthAll3Dtemp != 0){
      hTrackLenAll3DDay[par]->Add(hTrackLengthAll3Dtemp,1.0);
    }
    
    TH1F *hnTrackCont3Dtemp = (TH1F*)fileD.FindObjectAny("fNumTracksCont3D");
    if(hnTrackCont3Dtemp != 0){
      hnTrackCont3DDay[par]->Add(hnTrackCont3Dtemp,1.0);
    }

    TH2F *hTrackStartXZCont3Dtemp = (TH2F*)fileD.FindObjectAny("fStartPointTrackXZCRCont3D");
    if(hTrackStartXZCont3Dtemp != 0) {
      hTrackStartXZCont3DDay[par]->Add(hTrackStartXZCont3Dtemp,1.0);
    }
    
    TH2F *hTrackStartYZCont3Dtemp = (TH2F*)fileD.FindObjectAny("fStartPointTrackYZCRCont3D");
    if(hTrackStartYZCont3Dtemp != 0) {
      hTrackStartYZCont3DDay[par]->Add(hTrackStartYZCont3Dtemp,1.0);
    }

    TH2F *hTrackStopXZCont3Dtemp = (TH2F*)fileD.FindObjectAny("fStopPointTrackXZCRCont3D");
    if(hTrackStopXZCont3Dtemp != 0) {
      hTrackStopXZCont3DDay[par]->Add(hTrackStopXZCont3Dtemp,1.0);
    }

    TH2F *hTrackStopYZCont3Dtemp = (TH2F*)fileD.FindObjectAny("fStopPointTrackYZCRCont3D");
    if(hTrackStopYZCont3Dtemp != 0) {
      hTrackStopYZCont3DDay[par]->Add(hTrackStopYZCont3Dtemp,1.0);
    }

    TH1F *hTrackLengthCont3Dtemp = (TH1F*)fileD.FindObjectAny("fTrackLengthCont3D");
    if(hTrackLengthCont3Dtemp != 0){
      hTrackLenCont3DDay[par]->Add(hTrackLengthCont3Dtemp,1.0);
    }

    /*if(hTrackLenAll3Dtemp != 0 && header !=0) {
      NTrackLenAll3DtimeDay[par][NTrackLenAll3DDayCount[par]] = Xsrtime;
      NTrackLenAll3DDay    [par][NTrackLenAll3DDayCount[par]] = hTrackLenAll3Dtemp->GetMean(1);
      NTrackLenAll3DDayCount[par]++;      
      }*/

    // NEW PLOTS
    TH1F *hDSlTrktemp = (TH1F*)fileD.FindObjectAny("fDeltaSliceTrack");
    if(hDSlTrktemp != 0) {
      hDeltaSliceTrackDay[par]->Add(hDSlTrktemp,1.0);
    }

    TH1F *hSlTrkRtemp = (TH1F*)fileD.FindObjectAny("fSliceTrackRatio");
    if(hSlTrkRtemp != 0) {
      hSliceTrackRatioDay[par]->Add(hSlTrkRtemp,1.0);
    }

    TH1F *hSlTrkNHRtemp = (TH1F*)fileD.FindObjectAny("fSliceTrackNHitRatio");
    if(hSlTrkNHRtemp != 0) {
      hSliceTrackNHitRatioDay[par]->Add(hSlTrkNHRtemp,1.0);
    }

    TH1F *hTrackFracAll3Dtemp  = (TH1F*)fileD.FindObjectAny("fTrackFractionAll3D");
    TH1F *hTrackFracAll2Dtemp  = (TH1F*)fileD.FindObjectAny("fTrackFractionAll2D");
    TH1F *hTrackFracCont3Dtemp = (TH1F*)fileD.FindObjectAny("fTrackFractionCont3D");
    if(hTrackFracAll3Dtemp  != 0 && hTrackFracAll2Dtemp != 0 &&
       hTrackFracCont3Dtemp != 0 && header !=0) {
      TrackFracAll3DtimeDay[par][TrackFracDayCount[par]] = Xsrtime;
      TrackFracAll3DDay    [par][TrackFracDayCount[par]] = hTrackFracAll3Dtemp->GetMean(1);
      //TrackFracAll3DDay[par][TrackFracDayCount[par]] = hTrackFracAll3Dtemp->GetMean(1);
      

      TrackFracAll2DtimeDay[par][TrackFracDayCount[par]] = Xsrtime;
      TrackFracAll2DDay    [par][TrackFracDayCount[par]] = hTrackFracAll2Dtemp->GetMean(1);

      TrackFracCont3DtimeDay[par][TrackFracDayCount[par]] = Xsrtime;
      TrackFracCont3DDay    [par][TrackFracDayCount[par]] = hTrackFracCont3Dtemp->GetMean(1);

      TrackFracDayCount[par]++;
    }    

    fileD.cd("hitefficiencyana");
    TH2F *hHitsVwcelltemp = (TH2F*)fileD.FindObjectAny("hitsVwcell");
    if(hHitsVwcelltemp != 0){
      hHitsVwcellDay[par]->Add(hHitsVwcelltemp,1.0);
    }

    TH2F *hCellsVwcelltemp = (TH2F*)fileD.FindObjectAny("cellsVwcell");
    if(hCellsVwcelltemp != 0){
      hCellsVwcellDay[par]->Add(hCellsVwcelltemp,1.0);
    }
        

    fileD.Close();
  } // end while loop

  /////////////////////////////////////////////////////////////////////////////////////////////////////
  //NORMALIZATION STRIP PLOT FOR DCM SYNC STUFF
  for(unsigned int p=0; p < Npar; ++p){
    for(int i = 0; i < g3DRateEvo[p]->GetNbinsX(); ++i){      
      double norma = hNormalization[p]->GetBinContent(hNormalization[p]->GetBin(i));    
      for(int j =0; j< g3DRateEvo[p]->GetNbinsY(); ++j){
	double content = 0.0;
	if(norma > 0.0) content = (g3DRateEvo[p]->GetBinContent(g3DRateEvo[p]->GetBin(i,j)))/(norma);
	g3DRateEvo[p]->SetBinContent(g3DRateEvo[p]->GetBin(i,j),content);
      }
    }
  }

  // Define last run and subrun unique to each partition.  Define here as to avoid overwriting same value many times.
  // Last_Run[par] = LastRun;
  // Last_SR[par] = LastSR;

  //
  // Make time/date stamp for each plot
  //
  TDatime *Ttemp = new TDatime;  // finish time
  int Xfin = Ttemp->Convert() - GMToffset;
  TDatime *Tfinish = new TDatime(Xfin);  // finish time
  // std::cout<<Tfinish->AsString()<<endl;
  /* TPaveText *UpdateText = new TPaveText(0.1, 0.0, 0.5, 0.05, "NDC");
  UpdateText->SetLineColor(0);
  UpdateText->SetFillColor(0);
  UpdateText->SetBorderSize(1);
  UpdateText->SetMargin(0.0);
  UpdateText->SetTextAlign(11);*/
  char buff1[256];
  //sprintf(buff1, "Last updated on:    %s (central time)", Tfinish->AsString());
  //UpdateText->AddText(buff1);
  char buff2[256];
  //sprintf(buff2, "Last run / subrun:   %d / %d", LastRun, LastSR);
  //UpdateText->AddText(buff2);*/

  TPaveText *StatusText = new TPaveText(0.1, 0.1, 0.9, 0.9, "NDC");
  StatusText->SetLineColor(1);
  StatusText->SetFillColor(0);
  StatusText->SetBorderSize(1);
  // StatusText->SetMargin(0.0);
  StatusText->AddText(buff1);
  StatusText->AddText(buff2);

  TCanvas *cStatus = new TCanvas("cStatus","FD Status",800,200);
  cStatus->cd();
  StatusText->Draw();
  //cStatus->Print("/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/DCMsyncTEST/FarDetAnaStatus.png");
  cStatus->Print("/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/FarDetAnaStatus.png");  

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

  /* 
  TLine *L1x = new TLine(0,0,0,384);
  TLine *L2x = new TLine(0,384,896,384);
  TLine *L3x = new TLine(0,0,896,0);
  TLine *L4x = new TLine(896,0,896,384);
  TLine *L1y = new TLine(0,0,0,384);
  TLine *L2y = new TLine(0,384,896,384);
  TLine *L3y = new TLine(0,0,896,0);
  TLine *L4y = new TLine(896,0,896,384);
  L1x->SetLineWidth(2);
  L2x->SetLineWidth(2);
  L3x->SetLineWidth(2);
  L4x->SetLineWidth(2);
  L1y->SetLineWidth(2);
  L2y->SetLineWidth(2);
  L3y->SetLineWidth(2);
  L4y->SetLineWidth(2);
  */

  // Make conditional standard time axis labels depending on period (Day, Week, Month)
  string taxis_labels;
  if(period == "Day"){
    taxis_labels = "%H:%M";
  }
  else{
    taxis_labels = "%m/%d";
  }

  TLine *L1x;
  TLine *L2x;
  TLine *L3x;
  TLine *L4x;
  TLine *L1y;
  TLine *L2y;
  TLine *L3y;
  TLine *L4y;
  

  // New lines defined in Control Room view...
  if(det_type=="NearDet"){
    L1x = new TLine(0,0,0,-95);
    L2x = new TLine(0,-95,-213,-95);
    L3x = new TLine(0,0,-213,0);
    L4x = new TLine(-213,0,-213,-95);
    L1y = new TLine(0,0,0,95);
    L2y = new TLine(0,95,-213,95);
    L3y = new TLine(0,0,-213,0);
    L4y = new TLine(-213,0,-213,95);
  }
  else{
    L1x = new TLine(0,0,0,-384);
    L2x = new TLine(0,-384,-896,-384);
    L3x = new TLine(0,0,-896,0);
    L4x = new TLine(-896,0,-896,-384);
    L1y = new TLine(0,0,0,384);
    L2y = new TLine(0,384,-896,384);
    L3y = new TLine(0,0,-896,0);
    L4y = new TLine(-896,0,-896,384);
  }
  L1x->SetLineWidth(2);
  L2x->SetLineWidth(2);
  L3x->SetLineWidth(2);
  L4x->SetLineWidth(2);
  L1y->SetLineWidth(2);
  L2y->SetLineWidth(2);
  L3y->SetLineWidth(2);
  L4y->SetLineWidth(2);

  // loop over all partitions to make 24 hour plots
  for(unsigned int p = 0; p < Npar; ++p) {
    char filename[128];
    char title[128];
    
    TPaveText *UpdateText = new TPaveText(0.1, 0.0, 0.5, 0.05, "NDC");
    UpdateText->SetLineColor(0);
    UpdateText->SetFillColor(0);
    UpdateText->SetBorderSize(1);
    UpdateText->SetMargin(0.0);
    UpdateText->SetTextAlign(11);
    char buff3[256];
    char buff4[256];
    sprintf(buff3, "Last updated on:    %s (central time)", Tfinish->AsString());
    UpdateText->AddText(buff3);

    //UInt_t lcl_lastrun = Last_Run[p];
    //UInt_t lcl_lastsr = Last_SR[p];
    // Adding last run and subrun per partition inside of par loop.
    sprintf(buff4, "Last run / subrun:   %i / %i", Last_Run[p], Last_SR[p]);
    UpdateText->AddText(buff4);
    //delete buff3;
    
    // Draw DCM sync plots
    //////////////////////////////////////////////////////
    
    maxtime = 0;
    max = 0.0;
    ave = 0.0;
    for(int s = 0; s < count[p]; ++s) {
      ave += (double)minDCM[p][s];
      if(minDCM[p][s] > maxtime) {
        maxtime = xtime[p][s];
        max     = minDCM[p][s];
      }
    }
    if(count[p] > 0) ave = ave/(double)count[p];
    sprintf(lptext,"Last Point = %.f  /  Average = %f",max,ave);
    LastPoint->Clear();
    LastPoint->AddText(lptext);
    
    //MINIMUM OF THE RATE OF HITS ON 3D TRACKS (DMC-SYNC)
    //if(t_prd!="Month"){
    TCanvas *cMin3D = new TCanvas("cMin3D","Minimum of the Rate of Hits on 3D Tracks",1000,700);
    cMin3D->cd();
    cMin3D->SetGridy();
    
    TGraph *gmin = new TGraph(count[p],xtime[p],minDCM[p]); //minDCMError was used but set to 0 for now //dim
    sprintf(title,"%s Minimum of the Rate of Hits on 3D Tracks vs Time (DCM Sync) %s  - partition %2.2d",trigname.c_str(),time.c_str(), p);     
    gmin->SetTitle(title);
    gmin->GetXaxis()->SetTimeDisplay(1);
    gmin->GetXaxis()->SetLabelSize(0.03);
    gmin->GetXaxis()->SetTimeFormat(taxis_labels.c_str());
    gmin->GetXaxis()->SetLimits(time_ago,XNow);        
    gmin->SetMarkerStyle(7);   //prev 21
    
    gmin->GetYaxis()->SetTitle("MinHit3D [Hz]");      
    gmin->GetYaxis()->SetTitleOffset(1.45);
    
    gmin->GetXaxis()->SetTitle("(central time)");
    gmin->GetXaxis()->SetLabelOffset(0.015);
    gmin->GetXaxis()->SetTitleOffset(1.25);
    gmin->SetMarkerColor(4);  
    
    gmin->Draw("AP");  //"A*"
    UpdateText->Draw();
    LastPoint->Draw();

    sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1uMin_rate_3DTrks_%s.png",det.c_str(),trig.c_str(),p,period.c_str());
    cMin3D->Print(filename);
    delete cMin3D;
    cMin3D = 0;  
    
    ///////////////////////////////////////////////////////////
    //STRIP PLOT  i.e. EVOLUTION PLOT
    
    TCanvas *cRate3Devo = new TCanvas("cRate3Devo","Rate of Hits from 3D Tracks per DCM vs Time",1100,700);
    cRate3Devo->cd();
    cRate3Devo->SetLogz();    
    
    gPad->SetGridx();   
    
    g3DRateEvo[p]->GetXaxis()->SetTimeDisplay(1);
    g3DRateEvo[p]->GetXaxis()->SetLabelSize(0.03);
    g3DRateEvo[p]->GetXaxis()->SetTimeFormat(taxis_labels.c_str());
    g3DRateEvo[p]->GetXaxis()->SetLimits(time_ago,XNow);        
    g3DRateEvo[p]->GetXaxis()->SetLabelOffset(0.015);
    
    g3DRateEvo[p]->GetXaxis()->SetTitleOffset(1.25);    
    g3DRateEvo[p]->GetYaxis()->SetTitleOffset(1.05);
        
    g3DRateEvo[p]->GetYaxis()->SetTickLength(0);
    g3DRateEvo[p]->GetYaxis()->SetLabelColor(0);
    
    g3DRateEvo[p]->SetStats(kFALSE);
    
    g3DRateEvo[p]->Draw("colz");
    UpdateText->Draw();
    
    //DIBLOCK LABELS FOR STRIP PLOT    
    for(int i=1; i<14; ++i){
      //l = new TLine(XMonthAgo,i*12.,XNow,i*12.);
      l = new TLine(time_ago,i*12.,XNow,i*12.);
      l->SetLineColor(kBlack);
      l->SetLineWidth(2);
      l->SetBit(kCannotPick);
      l->Draw("same");      
    }
    
    for (int i=0; i<14; ++i) {
      sprintf(buff,"DB%.2d",i+1);
      buff[15]=0;      
      //t = new TText(XMonthAgo-10000.,12*(i+0.5),buff);
      t = new TText(time_ago-10000.,12*(i+0.5),buff);
      t->SetTextAlign(32);
      t->SetTextSize(0.03);
      t->SetBit(kCannotPick);
      //fText.push_back(t);
      t->Draw("same");
    }

    sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1u3DHitTrackRateEvolution_%s.png",det.c_str(),trig.c_str(),p,period.c_str());
    cRate3Devo->Print(filename);
    delete cRate3Devo;
    cRate3Devo = 0;
  
    //new beam plots
    TCanvas *cTimeDiffday = new TCanvas("cTimeDiffday","Time Difference Day",1200,800);
    cTimeDiffday->cd();
    gStyle->SetOptStat(111111);
    hTimeDiffNanoSecDay[p]->SetLineWidth(2);
    hTimeDiffNanoSecDay[p]->Draw();
    cTimeDiffday->Update();
    double lmax=gPad->GetUymax();
    TLine* lTimeDiff=new TLine(0.5,0,0.5,lmax);
    lTimeDiff->SetLineColor(kRed);
    lTimeDiff->SetLineWidth(3);
    lTimeDiff->SetLineStyle(9);
    lTimeDiff->Draw();
    hTimeDiffNanoSecDay[p]->Draw("same");
    UpdateText->Draw();
    sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1uTimeDiffDBBeamSpill_%s.png",det.c_str(),trig.c_str(),p,period.c_str());
    cTimeDiffday->Print(filename);
    delete cTimeDiffday;
    cTimeDiffday = 0;

    TCanvas *cPOTSumday = new TCanvas("cPOTSumday","POT Sum Day",1200,800);
    cPOTSumday->cd();
    gStyle->SetOptStat(111111);
    hPOTSumDay[p]->SetAxisRange(0.0,0.25,"X");
    hPOTSumDay[p]->SetLineWidth(2);
    hPOTSumDay[p]->GetYaxis()->SetTitle("(POT)");
    hPOTSumDay[p]->GetXaxis()->SetTitle("");
    hPOTSumDay[p]->Draw();
    UpdateText->Draw();
    sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1uPOTSum_%s.png",det.c_str(),trig.c_str(),p,period.c_str());
    cPOTSumday->Print(filename);
    delete cPOTSumday;
    cPOTSumday = 0;
    
    TCanvas *cspillPOTday = new TCanvas("cspillPOTday","spill POT Day",1200,800);
    cspillPOTday->cd();
    gStyle->SetOptStat(111111);
    hspillPOTDay[p]->SetLineWidth(2);
    hspillPOTDay[p]->Draw();
    UpdateText->Draw();
    cspillPOTday->Update();
    lmax=gPad->GetUymax();
    TLine* lspillPOT=new TLine(2e12,0,2e12,lmax);
    lspillPOT->SetLineColor(kRed);
    lspillPOT->SetLineWidth(3);
    lspillPOT->SetLineStyle(9);
    lspillPOT->Draw();
    hspillPOTDay[p]->Draw("same");
    sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1uspillPOT_%s.png",det.c_str(),trig.c_str(),p,period.c_str());
    cspillPOTday->Print(filename);
    delete cspillPOTday;
    cspillPOTday = 0;
 
    TCanvas *cHornCurrentday = new TCanvas("cHornCurrentday","Time Difference Day",1200,800);
    cHornCurrentday->cd();
    gStyle->SetOptStat(111111);
    hHornCurrentDay[p]->SetAxisRange(-202.5,-197.5,"X");
    hHornCurrentDay[p]->SetLineWidth(2);
    hHornCurrentDay[p]->Draw();
    cHornCurrentday->Update();
    lmax=gPad->GetUymax();
    TLine* lHornCurrent=new TLine(-202,0,-202,lmax);
    TLine* l2HornCurrent=new TLine(-198,0,-198,lmax);
    lHornCurrent->SetLineColor(kRed);
    lHornCurrent->SetLineWidth(3);
    lHornCurrent->SetLineStyle(9);
    lHornCurrent->Draw();
    l2HornCurrent->SetLineColor(kRed);
    l2HornCurrent->SetLineWidth(3);
    l2HornCurrent->SetLineStyle(9);
    lHornCurrent->Draw();
    l2HornCurrent->Draw();
    hHornCurrentDay[p]->Draw("same");
    UpdateText->Draw();
    sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1uHornCurrent_%s.png",det.c_str(),trig.c_str(),p,period.c_str());
    cHornCurrentday->Print(filename);
    delete cHornCurrentday;
    cHornCurrentday = 0;
    
    TCanvas *cXPositionday = new TCanvas("cXPositionday","X Position Day",1200,800);
    cXPositionday->cd();
    gStyle->SetOptStat(111111);
    hXPositionDay[p]->SetLineWidth(2);
    hXPositionDay[p]->Draw();
    UpdateText->Draw();
    cXPositionday->Update();
    lmax=gPad->GetUymax();
    TLine* lXPosition=new TLine(2,0,2,lmax);
    TLine* l2XPosition=new TLine(0.02,0,0.02,lmax);
    lXPosition->SetLineColor(kRed);
    lXPosition->SetLineWidth(3);
    lXPosition->SetLineStyle(9);
    lXPosition->Draw();
    l2XPosition->SetLineColor(kRed);
    l2XPosition->SetLineWidth(3);
    l2XPosition->SetLineStyle(9);
    l2XPosition->Draw();
    hXPositionDay[p]->Draw("same");
    sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1uXPosition_%s.png",det.c_str(),trig.c_str(),p,period.c_str());
    cXPositionday->Print(filename);
    delete cXPositionday;
    cXPositionday = 0;

     TCanvas *cYPositionday = new TCanvas("cYPositionday","Y Position Day",1200,800);
    cYPositionday->cd();
    gStyle->SetOptStat(111111);
    hYPositionDay[p]->SetLineWidth(2);
    hYPositionDay[p]->Draw();
    UpdateText->Draw();
    cYPositionday->Update();
    lmax=gPad->GetUymax();
    TLine* lYPosition=new TLine(2,0,2,lmax);
    TLine* l2YPosition=new TLine(0.02,0,0.02,lmax);
    lYPosition->SetLineColor(kRed);
    lYPosition->SetLineWidth(3);
    lYPosition->SetLineStyle(9);
    lYPosition->Draw();
    l2YPosition->SetLineColor(kRed);
    l2YPosition->SetLineWidth(3);
    l2YPosition->SetLineStyle(9);
    l2YPosition->Draw();
    hYPositionDay[p]->Draw("same");
    sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1uYPosition_%s.png",det.c_str(),trig.c_str(),p,period.c_str());
    cYPositionday->Print(filename);
    delete cYPositionday;
    cYPositionday = 0;

    TCanvas *cXWidthday = new TCanvas("cXWidthday","X Width Day",1200,800);
    cXWidthday->cd();
    gStyle->SetOptStat(111111);
    hXWidthDay[p]->SetLineWidth(2);
    hXWidthDay[p]->Draw();
    UpdateText->Draw();
    cXWidthday->Update();
    lmax=gPad->GetUymax();
    TLine* lXWidth=new TLine(0.57,0,0.57,lmax);
    TLine* l2XWidth=new TLine(1.58,0,1.58,lmax);
    lXWidth->SetLineColor(kRed);
    lXWidth->SetLineWidth(3);
    lXWidth->SetLineStyle(9);
    lXWidth->Draw();
    l2XWidth->SetLineColor(kRed);
    l2XWidth->SetLineWidth(3);
    l2XWidth->SetLineStyle(9);
    l2XWidth->Draw();
    hXWidthDay[p]->Draw("same");
    sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1uXWidth_%s.png",det.c_str(),trig.c_str(),p,period.c_str());
    cXWidthday->Print(filename);
    delete cXWidthday;
    cXWidthday = 0;
    
    TCanvas *cYWidthday = new TCanvas("cYWidthday","Y Width Day",1200,800);
    cYWidthday->cd();
    gStyle->SetOptStat(111111);
    hYWidthDay[p]->SetLineWidth(2);
    hYWidthDay[p]->Draw();
    UpdateText->Draw();
    cYWidthday->Update();
    lmax=gPad->GetUymax();
    TLine* lYWidth=new TLine(0.57,0,0.57,lmax);
    TLine* l2YWidth=new TLine(1.58,0,1.58,lmax);
    lYWidth->SetLineColor(kRed);
    lYWidth->SetLineWidth(3);
    lYWidth->SetLineStyle(9);
    lYWidth->Draw();
    l2YWidth->SetLineColor(kRed);
    l2YWidth->SetLineWidth(3);
    l2YWidth->SetLineStyle(9);
    l2YWidth->Draw();
    hYWidthDay[p]->Draw("same");
    sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1uYWidth_%s.png",det.c_str(),trig.c_str(),p,period.c_str());
    cYWidthday->Print(filename);
    delete cYWidthday;
    cYWidthday = 0;
 
    TCanvas *cGoodBeamday = new TCanvas("cGoodBeamday","Y Good Beam Day",1200,800);
    cGoodBeamday->cd();
    gStyle->SetOptStat(111111);
    hGoodBeamDay[p]->SetLineWidth(2);
    hGoodBeamDay[p]->Draw();
    TAxis *GoodBeamaxis = hGoodBeamDay[p]->GetXaxis();
    GoodBeamaxis->SetBinLabel(1,"bad beam");
    GoodBeamaxis->SetBinLabel(2,"good beam");
    UpdateText->Draw();
    sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1uGoodBeam_%s.png",det.c_str(),trig.c_str(),p,period.c_str());
    cGoodBeamday->Print(filename);
    delete cGoodBeamday;
    cGoodBeamday = 0;

    TCanvas *cBadSpillsday = new TCanvas("cBadSpillsday","Y Bad Spills Day",1200,800);
    cBadSpillsday->cd();
    hBadSpillsDay[p]->SetLineWidth(2);
    gStyle->SetOptStat(0);
    hBadSpillsDay[p]->Draw();
    UpdateText->Draw();
    sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1uBadSpills_%s.png",det.c_str(),trig.c_str(),p,period.c_str());
    cBadSpillsday->Print(filename);
    delete cBadSpillsday;
    cBadSpillsday = 0;
    //end new

    TCanvas *cnSliceday = new TCanvas("cnSliceday","# Slice Day",1200,800);
    cnSliceday->cd();
    cnSliceday->SetLogy();
    gStyle->SetOptStat(111111);
    if(det_type=="NearDet") hNSliceDay[p]->SetAxisRange(0.0,50.0,"X");
    else hNSliceDay[p]->SetAxisRange(0.0,200.0,"X");
    hNSliceDay[p]->SetLineWidth(2);
    hNSliceDay[p]->SetLineColor(kRed);
    hNSliceDay[p]->Draw();
    UpdateText->Draw();

    sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1uNSlice%s.png",det.c_str(),trig.c_str(),p,period.c_str());
    cnSliceday->Print(filename);
    delete cnSliceday;
    cnSliceday = 0;

    TCanvas *cnSliceHitday = new TCanvas("cnSliceHitday","# Non-Noise Slice Hits Day",1200,800);
    cnSliceHitday->cd();
    cnSliceHitday->SetLogy();
    gStyle->SetOptStat(111111);
    gPad->SetLogy();
    gPad->SetLogx();
    hNNonNoiseSliceHitDay[p]->SetAxisRange(0.0,3000.0,"X");
    hNNonNoiseSliceHitDay[p]->SetLineWidth(2);
    hNNonNoiseSliceHitDay[p]->SetLineColor(kRed);
    hNNonNoiseSliceHitDay[p]->Draw();
    UpdateText->Draw();
    sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1uNSliceHit%s.png",det.c_str(),trig.c_str(),p,period.c_str());
    cnSliceHitday->Print(filename);
    delete cnSliceHitday;
    cnSliceHitday = 0;

    TCanvas *cNNoiseSliceHitday = new TCanvas("cNNoiseSliceHitday","# Noise Slice Hits Day",1200,800);
    cNNoiseSliceHitday->cd();
    cNNoiseSliceHitday->SetLogy();
    gStyle->SetOptStat(111111);
    if(det_type=="NearDet") hNNoiseSliceHitDay[p]->SetAxisRange(0.0,2000.0,"X");
    // hNNoiseSliceHitDay[p]->SetAxisRange(0.0,10000.0,"X");
    hNNoiseSliceHitDay[p]->SetLineWidth(2);
    hNNoiseSliceHitDay[p]->SetLineColor(kRed);
    hNNoiseSliceHitDay[p]->Draw();
    UpdateText->Draw();
    sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1uNNoiseSliceHit%s.png",det.c_str(),trig.c_str(),p,period.c_str());
    cNNoiseSliceHitday->Print(filename);
    delete cNNoiseSliceHitday;
    cNNoiseSliceHitday = 0;

    TCanvas *cPEday = new TCanvas("cPEday","PE distributions by Day",1200,800);
    cPEday->cd();
    gStyle->SetOptStat(0);
    gPad->SetLogy();
    gPad->SetLogx();
    sprintf(filename,"%s PE Distribution (%s) - partition %.1u",trigname.c_str(),time.c_str(),p);
    hNoiseSlicePEDay[p]->SetTitle(filename);
    hNonNoiseSlicePEDay[p]->SetTitle(filename);
    hTrackPEDay[p]->SetTitle(filename);
    hNoiseSlicePEDay[p]->SetLineWidth(2);
    hNoiseSlicePEDay[p]->SetLineColor(kBlack);
    hNoiseSlicePEDay[p]->Draw();
    hNonNoiseSlicePEDay[p]->SetLineWidth(2);
    hNonNoiseSlicePEDay[p]->SetLineColor(kRed);
    hNonNoiseSlicePEDay[p]->Draw("sames");
    hTrackPEDay[p]->SetLineWidth(2);
    hTrackPEDay[p]->SetLineColor(kBlue);
    hTrackPEDay[p]->Draw("sames");
    UpdateText->Draw();
    TLegend *legend = new TLegend(0.6,0.6,0.85,0.85);
    legend->AddEntry(hNoiseSlicePEDay[p],"Hits in Noise Slice","l");
    legend->AddEntry(hNonNoiseSlicePEDay[p],"Hits in Non-Noise Slices","l");
    legend->AddEntry(hTrackPEDay[p],"Hits on Tracks","l");
    legend->SetFillStyle(4000);
    legend->SetLineColor(0);
    legend->SetFillColor(0);
    legend->Draw();
    sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1uPE%s.png",det.c_str(),trig.c_str(),p,period.c_str());
    cPEday->Print(filename);
    delete cPEday;
    cPEday = 0;
    delete legend;
    legend = 0;

    TCanvas *cTsdday = new TCanvas("cTsdday","Slice Time Standard Dev.",1200,800);
    cTsdday->cd();
    cTsdday->SetLogy();
    gStyle->SetOptStat(111111);
    if(det_type=="NearDet") hTsdDay[p]->SetAxisRange(0.0,600.0,"X");
    else hTsdDay[p]->SetAxisRange(0.0,1000.0,"X");
    hTsdDay[p]->SetLineWidth(2);
    hTsdDay[p]->SetLineColor(kRed);
    hTsdDay[p]->Draw();
    UpdateText->Draw();
    sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1uTsd%s.png",det.c_str(),trig.c_str(),p,period.c_str());
    cTsdday->Print(filename);
    delete cTsdday;
    cTsdday = 0;

    TCanvas *cnTrackday = new TCanvas("cnTrackday"," # Cosmic Tracks",1200,800);
    cnTrackday->cd();
    cnTrackday->SetLogy();
    gStyle->SetOptStat(111111);
    hnTrackDay[p]->SetAxisRange(0.0,300.0,"X");
    hnTrackDay[p]->SetLineWidth(2);
    hnTrackDay[p]->SetLineColor(kRed);
    hnTrackDay[p]->Draw();
    UpdateText->Draw();
    sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1unTrackAll3D%s.png",det.c_str(),trig.c_str(),p,period.c_str());
    cnTrackday->Print(filename);
    delete cnTrackday;
    cnTrackday = 0;

    TCanvas *cTrackLenday = new TCanvas("cTrackLenday","Cosmic Track Length",1200,800);
    cTrackLenday->cd();
    cTrackLenday->SetLogy();
    gStyle->SetOptStat(111111);
    hTrackLenDay[p]->SetAxisRange(0.0,6000.0,"X");
    hTrackLenDay[p]->SetLineWidth(2);
    hTrackLenDay[p]->SetLineColor(kRed);
    hTrackLenDay[p]->Draw();
    UpdateText->Draw();
    sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1uTrackLenAll3D%s.png",det.c_str(),trig.c_str(),p,period.c_str());
    cTrackLenday->Print(filename);
    delete cTrackLenday;
    cTrackLenday = 0;

    TCanvas *cCosNumiday = new TCanvas("cCosNumiday","Cos Theta(numi)",1200,800);
    cCosNumiday->cd();
    gStyle->SetOptStat(111111);
    hCosNumiDay[p]->SetAxisRange(-1.0,1.0,"X");
    hCosNumiDay[p]->SetLineWidth(2);
    hCosNumiDay[p]->SetLineColor(kRed);
    hCosNumiDay[p]->Draw();
    UpdateText->Draw();
    sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1uCosNumiAll3D%s.png",det.c_str(),trig.c_str(),p,period.c_str());
    cCosNumiday->Print(filename);
    delete cCosNumiday;
    cCosNumiday = 0;

    TCanvas *cPCaveday = new TCanvas("cPCaveday","PCave day",1200,1200);
    cPCaveday->Divide(1,2);
    cPCaveday->cd(1);
    gStyle->SetOptStat(0);
    gPad->SetLogz();
    if(det_type=="NearDet"){
      hPCXaveDay[p]->SetAxisRange(-225.0,20.0,"X");
      hPCXaveDay[p]->SetAxisRange(-110.0,10.0,"Y");
    }
    else{
      hPCXaveDay[p]->SetAxisRange(-920.0,20.0,"X");
      hPCXaveDay[p]->SetAxisRange(-400.0,10.0,"Y");
    }
    hPCXaveDay[p]->Draw("colz");
    L1x->Draw();
    L2x->Draw();
    L3x->Draw();
    L4x->Draw();
    cPCaveday->cd(2);
    gStyle->SetOptStat(0);
    gPad->SetLogz();
    if(det_type=="NearDet"){
      hPCYaveDay[p]->SetAxisRange(-225.0,20.0,"X");
      hPCYaveDay[p]->SetAxisRange(-10.0,110.0,"Y");
    }
    else{
      hPCYaveDay[p]->SetAxisRange(-920.0,20.0,"X");
      hPCYaveDay[p]->SetAxisRange(-10.0,400.0,"Y");
    }
    hPCYaveDay[p]->Draw("colz");
    L1y->Draw();
    L2y->Draw();
    L3y->Draw();
    L4y->Draw();
    UpdateText->Draw();
    sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1uPCave%s.png",det.c_str(),trig.c_str(),p,period.c_str());
    cPCaveday->Print(filename);
    delete cPCaveday;
    cPCaveday = 0;

    TCanvas *cnTrackAll3Dday = new TCanvas("cnTrackday"," # 3D Cosmic Tracks",1200,800);
    cnTrackAll3Dday->cd();
    cnTrackAll3Dday->SetLogy();
    gStyle->SetOptStat(111111);
    if(det_type=="NearDet") hnTrackAll3DDay[p]->SetAxisRange(0.0,50.0,"X");
    else hnTrackAll3DDay[p]->SetAxisRange(0.0,300.0,"X");
    hnTrackAll3DDay[p]->SetLineWidth(2);
    hnTrackAll3DDay[p]->SetLineColor(kRed);
    hnTrackAll3DDay[p]->Draw();
    UpdateText->Draw();
    sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1unTrackAll3D%s.png",det.c_str(),trig.c_str(),p,period.c_str());
    cnTrackAll3Dday->Print(filename);
    delete cnTrackAll3Dday;
    cnTrackAll3Dday = 0;

    TCanvas *cTrackStartAll3Dday = new TCanvas("cTrackStartAll3Dday","All 3D Track Start day",1200,1200);
    cTrackStartAll3Dday->Divide(1,2);
    cTrackStartAll3Dday->cd(1);
    gStyle->SetOptStat(0);
    gPad->SetLogz();
    if(det_type=="NearDet"){
      hTrackStartXZAll3DDay[p]->SetAxisRange(-1750.0,0.0,"X");
      hTrackStartXZAll3DDay[p]->SetAxisRange(-250.0,250.0,"Y");
    }
    hTrackStartXZAll3DDay[p]->Draw("colz");
    cTrackStartAll3Dday->cd(2);
    gStyle->SetOptStat(0);
    gPad->SetLogz();
    if(det_type=="NearDet"){
      hTrackStartYZAll3DDay[p]->SetAxisRange(-1750.0,0.0,"X");
      hTrackStartYZAll3DDay[p]->SetAxisRange(-250.0,250.0,"Y");
    }
    hTrackStartYZAll3DDay[p]->Draw("colz");
    UpdateText->Draw();
    sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1uTrackStartAll3D%s.png",det.c_str(),trig.c_str(),p,period.c_str());
    cTrackStartAll3Dday->Print(filename);
    delete cTrackStartAll3Dday;
    cTrackStartAll3Dday = 0;

    TCanvas *cTrackStopAll3Dday = new TCanvas("cTrackStopAll3Dday","All 3D Track Stop day",1200,1200);
    cTrackStopAll3Dday->Divide(1,2);
    cTrackStopAll3Dday->cd(1);
    gStyle->SetOptStat(0);
    gPad->SetLogz();
    if(det_type=="NearDet"){
      hTrackStopXZAll3DDay[p]->SetAxisRange(-1750.0,0.0,"X");
      hTrackStopXZAll3DDay[p]->SetAxisRange(-250.0,250.0,"Y");
    }
    hTrackStopXZAll3DDay[p]->Draw("colz");
    cTrackStopAll3Dday->cd(2);
    gStyle->SetOptStat(0);
    gPad->SetLogz();
    if(det_type=="NearDet"){
      hTrackStopYZAll3DDay[p]->SetAxisRange(-1750.0,0.0,"X");
      hTrackStopYZAll3DDay[p]->SetAxisRange(-250.0,250.0,"Y");
    }
    hTrackStopYZAll3DDay[p]->Draw("colz");
    UpdateText->Draw();
    sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1uTrackStopAll3D%s.png",det.c_str(),trig.c_str(),p,period.c_str());
    cTrackStopAll3Dday->Print(filename);
    delete cTrackStopAll3Dday;
    cTrackStopAll3Dday = 0;

    TCanvas *cTrackLenAll3Dday = new TCanvas("cTrackLenAll3Dday"," 3D Cosmic Track Length",1200,800);
    cTrackLenAll3Dday->cd();
    cTrackLenAll3Dday->SetLogy();
    gStyle->SetOptStat(111111);
    if(det_type=="NearDet") hTrackLenAll3DDay[p]->SetAxisRange(0.0,2000.0,"X");
    else hTrackLenAll3DDay[p]->SetAxisRange(0.0,6000.0,"X");
    //hTrackLenAll3DDay[p]->GetXaxis()->SetRangeUser(0.,6000.);
    hTrackLenAll3DDay[p]->SetLineWidth(2);
    hTrackLenAll3DDay[p]->SetLineColor(kRed);
    hTrackLenAll3DDay[p]->Draw();
    UpdateText->Draw();
    sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1uTrackLenAll3D%s.png",det.c_str(),trig.c_str(),p,period.c_str());
    cTrackLenAll3Dday->Print(filename);
    delete cTrackLenAll3Dday;
    cTrackLenAll3Dday = 0;

    TCanvas *cnTrackCont3Dday = new TCanvas("cnTrackCont3Dday"," # Contained 3D Cosmic Tracks",1200,800);
    cnTrackCont3Dday->cd();
    cnTrackCont3Dday->SetLogy();
    gStyle->SetOptStat(111111);
    hnTrackCont3DDay[p]->SetAxisRange(0.0,300.0,"X");
    hnTrackCont3DDay[p]->SetLineWidth(2);
    hnTrackCont3DDay[p]->SetLineColor(kRed);
    hnTrackCont3DDay[p]->Draw();
    UpdateText->Draw();
    sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1unTrackCont3D%s.png",det.c_str(),trig.c_str(),p,period.c_str());
    cnTrackCont3Dday->Print(filename);
    delete cnTrackCont3Dday;
    cnTrackCont3Dday = 0;

    TCanvas *cTrackStartCont3Dday = new TCanvas("cTrackStartCont3Dday","Cont 3D Track Start day",1200,1200);
    cTrackStartCont3Dday->Divide(1,2);
    cTrackStartCont3Dday->cd(1);
    gStyle->SetOptStat(0);
    gPad->SetLogz();
    if(det_type=="NearDet"){
      hTrackStartXZCont3DDay[p]->SetAxisRange(-1750.0,0.0,"X");
      hTrackStartXZCont3DDay[p]->SetAxisRange(-250.0,250.0,"Y");
    }
    hTrackStartXZCont3DDay[p]->Draw("colz");
    cTrackStartCont3Dday->cd(2);
    gStyle->SetOptStat(0);
    gPad->SetLogz();
    if(det_type=="NearDet"){
      hTrackStartYZCont3DDay[p]->SetAxisRange(-1750.0,0.0,"X");
      hTrackStartYZCont3DDay[p]->SetAxisRange(-250.0,250.0,"Y");
    }
    hTrackStartYZCont3DDay[p]->Draw("colz");
    UpdateText->Draw();
    sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1uTrackStartCont3D%s.png",det.c_str(),trig.c_str(),p,period.c_str());
    cTrackStartCont3Dday->Print(filename);
    delete cTrackStartCont3Dday;
    cTrackStartCont3Dday = 0;

    TCanvas *cTrackStopCont3Dday = new TCanvas("cTrackStopCont3Dday","Cont 3D Track Stop day",1200,1200);
    cTrackStopCont3Dday->Divide(1,2);
    cTrackStopCont3Dday->cd(1);
    gStyle->SetOptStat(0);
    gPad->SetLogz();
    if(det_type=="NearDet"){
      hTrackStopXZCont3DDay[p]->SetAxisRange(-1750.0,0.0,"X");
      hTrackStopXZCont3DDay[p]->SetAxisRange(-250.0,250.0,"Y");
    }
    hTrackStopXZCont3DDay[p]->Draw("colz");
    cTrackStopCont3Dday->cd(2);
    gStyle->SetOptStat(0);
    gPad->SetLogz();
    if(det_type=="NearDet"){
      hTrackStopYZCont3DDay[p]->SetAxisRange(-1750.0,0.0,"X");
      hTrackStopYZCont3DDay[p]->SetAxisRange(-250.0,250.0,"Y");
    }
    hTrackStopYZCont3DDay[p]->Draw("colz");
    UpdateText->Draw();
    sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1uTrackStopCont3D%s.png",det.c_str(),trig.c_str(),p,period.c_str());
    cTrackStopCont3Dday->Print(filename);
    delete cTrackStopCont3Dday;
    cTrackStopCont3Dday = 0;

    TCanvas *cTrackLenCont3Dday = new TCanvas("cTrackLenCont3Dday"," 3D Cosmic Track Length",1200,800);
    cTrackLenCont3Dday->cd();
    cTrackLenCont3Dday->SetLogy();
    gStyle->SetOptStat(111111);
    hTrackLenCont3DDay[p]->SetAxisRange(0.0,6000.0,"X");
    //hTrackLenCont3DDay[p]->GetXaxis()->SetRangeUser(0.,6000.);
    hTrackLenCont3DDay[p]->SetLineWidth(2);
    hTrackLenCont3DDay[p]->SetLineColor(kRed);
    hTrackLenCont3DDay[p]->Draw();
    UpdateText->Draw();
    sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1uTrackLenCont3D%s.png",det.c_str(),trig.c_str(),p,period.c_str());
    cTrackLenCont3Dday->Print(filename);
    delete cTrackLenCont3Dday;
    cTrackLenCont3Dday = 0;

    // NEW PLOTS
    TCanvas *cDeltaSliceTrackday = new TCanvas("cDeltaSliceTrackday","Delta Slice/Track",1200,800);
    cDeltaSliceTrackday->cd();
    cDeltaSliceTrackday->SetLogy();
    gStyle->SetOptStat(111111);
    hDeltaSliceTrackDay[p]->SetAxisRange(-15.0,15.0,"X");
    hDeltaSliceTrackDay[p]->SetLineWidth(2);
    hDeltaSliceTrackDay[p]->SetLineColor(kRed);
    hDeltaSliceTrackDay[p]->Draw();
    UpdateText->Draw();
    sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1uDeltaSliceTrack%s.png",det.c_str(),trig.c_str(),p,period.c_str());
    cDeltaSliceTrackday->Print(filename);
    delete cDeltaSliceTrackday;
    cDeltaSliceTrackday = 0;

    TCanvas *cSliceTrackRatioday = new TCanvas("cSliceTrackRatioday","Slice/Track Ratio",1200,800);
    cSliceTrackRatioday->cd();
    cSliceTrackRatioday->SetLogy();
    gStyle->SetOptStat(111111);
    hSliceTrackRatioDay[p]->SetAxisRange(0.0,1.0,"X");  // NOTE: This range will be different for the NearDet.C_Str()
    hSliceTrackRatioDay[p]->SetLineWidth(2);
    hSliceTrackRatioDay[p]->SetLineColor(kRed);
    hSliceTrackRatioDay[p]->Draw();
    UpdateText->Draw();
    sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1uSliceTrackRatio%s.png",det.c_str(),trig.c_str(),p,period.c_str());
    cSliceTrackRatioday->Print(filename);
    delete cSliceTrackRatioday;
    cSliceTrackRatioday = 0;

    TCanvas *cSliceTrackNHitRatioday = new TCanvas("cSliceTrackNHitRatioday","Slice/Track NHit Ratio",1200,800);
    cSliceTrackNHitRatioday->cd();
    cSliceTrackNHitRatioday->SetLogy();
    gStyle->SetOptStat(111111);
    hSliceTrackNHitRatioDay[p]->SetLineWidth(2);
    hSliceTrackNHitRatioDay[p]->SetLineColor(kRed);
    hSliceTrackNHitRatioDay[p]->Draw();
    UpdateText->Draw();
    sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1uSliceTrackNHitRatio%s.png",det.c_str(),trig.c_str(),p,period.c_str());
    cSliceTrackNHitRatioday->Print(filename);
    delete cSliceTrackNHitRatioday;
    cSliceTrackNHitRatioday = 0;


    //new beam plots
    maxtime = 0;
    max = 0.0;
    ave = 0.0;
    for(int i = 0; i < NTimeDiffDayCount[p]; ++i) {
      ave += (double)NTimeDiffNanoSecDay[p][i];
      if(NTimeDiffNanoSectimeDay[p][i] > maxtime) {
	maxtime = NTimeDiffNanoSectimeDay[p][i];
	max     = NTimeDiffNanoSecDay[p][i];
      }
    }
    if(NTimeDiffDayCount[p] > 0) ave = ave/(double)NTimeDiffDayCount[p];
    sprintf(lptext,"Last Point = %f  /  Average = %f",max,ave);
    LastPoint->Clear();
    LastPoint->AddText(lptext);
    TCanvas *ctnTimeDiffday = new TCanvas("cnTimeDiffday","Time Difference (spill and IFDB) day",1200,800);
    ctnTimeDiffday->cd();
    gPad->SetGridx();
    TGraph *gnTimeDiffday = new TGraph(NTimeDiffDayCount[p],NTimeDiffNanoSectimeDay[p],NTimeDiffNanoSecDay[p]);
    sprintf(title,"%s Time Difference (spill and IFDB) - partition %.1u",trigname.c_str(),p);
    gnTimeDiffday->SetTitle(title);
    gnTimeDiffday->SetMarkerColor(kBlue);
    gnTimeDiffday->GetXaxis()->SetTimeDisplay(1);
    gnTimeDiffday->GetXaxis()->SetLabelSize(0.03);
    gnTimeDiffday->GetXaxis()->SetTimeFormat(taxis_labels.c_str());
    gnTimeDiffday->GetXaxis()->SetLimits(time_ago,XNow);
    gnTimeDiffday->GetXaxis()->SetTitle("(central time)");
    gnTimeDiffday->GetYaxis()->SetTitle("(sec)");
    gnTimeDiffday->GetYaxis()->SetRangeUser(-0.1,0.6);
    gnTimeDiffday->Draw("A*"); 
    TLine* l2TimeDiff=new TLine(time_ago,0.5,XNow,0.5);
    l2TimeDiff->SetLineColor(kRed);
    l2TimeDiff->SetLineWidth(3);
    l2TimeDiff->SetLineStyle(9);
    l2TimeDiff->Draw();
    UpdateText->Draw();
    LastPoint->Draw();
    sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1unTimeDiff%s.png",det.c_str(),trig.c_str(),p,period.c_str());
    ctnTimeDiffday->Print(filename);
    delete ctnTimeDiffday;
    ctnTimeDiffday = 0;

    // Adding "zoomed in" version of above plot
    TCanvas *ctnTimeDiff_zoom = new TCanvas("cnTimeDiff_zoom","Time Difference (sill and IFDB) day",1200,800);
    ctnTimeDiff_zoom->cd();
    gPad->SetGridx();
    TGraph *gnTimeDiff_zoom = new TGraph(NTimeDiffDayCount[p],NTimeDiffNanoSectimeDay[p],NTimeDiffNanoSecDay[p]);
    sprintf(title,"%s Time Difference (spill and IFDB) - partition %.1u",trigname.c_str(),p);
    gnTimeDiff_zoom->SetTitle(title);
    gnTimeDiff_zoom->SetMarkerColor(kBlue);
    gnTimeDiff_zoom->GetXaxis()->SetTimeDisplay(1);
    gnTimeDiff_zoom->GetXaxis()->SetLabelSize(0.03);
    gnTimeDiff_zoom->GetXaxis()->SetTimeFormat(taxis_labels.c_str());
    gnTimeDiff_zoom->GetXaxis()->SetLimits(time_ago,XNow);
    gnTimeDiff_zoom->GetXaxis()->SetTitle("(central time)");
    gnTimeDiff_zoom->GetYaxis()->SetTitle("(nsec)");
    gnTimeDiff_zoom->GetYaxis()->SetRangeUser((gnTimeDiff_zoom->GetMean(2)-3.*gnTimeDiff_zoom->GetRMS(2)),(gnTimeDiff_zoom->GetMean(2)+3.*gnTimeDiff_zoom->GetRMS(2)));
    gnTimeDiff_zoom->Draw("A*");
    UpdateText->Draw();
    LastPoint->Draw();
    sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1unTimeDiffzoom%s.png",det.c_str(),trig.c_str(),p,period.c_str());
    ctnTimeDiff_zoom->Print(filename);
    delete ctnTimeDiff_zoom;
    ctnTimeDiff_zoom = 0;

    maxtime = 0;
    max = 0.0;
    ave = 0.0;
    for(int i = 0; i < NPOTSumDayCount[p]; ++i) {
      ave += (double)NPOTSumDay[p][i];
      if(NPOTSumtimeDay[p][i] > maxtime) {
	maxtime = NPOTSumtimeDay[p][i];
	max     = NPOTSumDay[p][i];
      }
    }
    if(NPOTSumDayCount[p] > 0) ave = ave/(double)NPOTSumDayCount[p];
    sprintf(lptext,"Last Point = %f  /  Average = %f",max,ave);
    LastPoint->Clear();
    LastPoint->AddText(lptext);
    TCanvas *ctnPOTSumday = new TCanvas("cnPOTSumday","POT Sum day",1200,800);
    ctnPOTSumday->cd();
    gPad->SetGridx();
    TGraph *gnPOTSumday = new TGraph(NPOTSumDayCount[p],NPOTSumtimeDay[p],NPOTSumDay[p]);
    sprintf(title,"%s POT per Subrun - partition %.1u",trigname.c_str(),p);
    gnPOTSumday->SetTitle(title);
    gnPOTSumday->SetMarkerColor(kBlue);
    gnPOTSumday->GetXaxis()->SetTimeDisplay(1);
    gnPOTSumday->GetXaxis()->SetLabelSize(0.03);
    gnPOTSumday->GetXaxis()->SetTimeFormat(taxis_labels.c_str());
    gnPOTSumday->GetXaxis()->SetLimits(time_ago,XNow);
    gnPOTSumday->GetYaxis()->SetLabelSize(0.03);
    gnPOTSumday->GetYaxis()->SetTitle("(POT)");
    gnPOTSumday->GetXaxis()->SetTitle("(central time)");
    gnPOTSumday->GetYaxis()->SetRangeUser(25E15,100E15);
    gnPOTSumday->Draw("A*");
    UpdateText->Draw();
    LastPoint->Draw();
    sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1unPOTSumperSR%s.png",det.c_str(),trig.c_str(),p,period.c_str());
    ctnPOTSumday->Print(filename);
    delete ctnPOTSumday;
    ctnPOTSumday = 0;

    // Adding "zoomed in" version of above plot
    TCanvas *ctnPOTSumday_zoom = new TCanvas("cnPOTSumday_zoom","POT Sum day",1200,800);
    ctnPOTSumday_zoom->cd();
    gPad->SetGridx();
    TGraph *gnPOTSumday_zoom = new TGraph(NPOTSumDayCount[p],NPOTSumtimeDay[p],NPOTSumDay[p]);
    sprintf(title,"%s POT per subrun - partition %.1u",trigname.c_str(),p);
    gnPOTSumday_zoom->SetTitle(title);
    gnPOTSumday_zoom->SetMarkerColor(kBlue);
    gnPOTSumday_zoom->GetXaxis()->SetTimeDisplay(1);
    gnPOTSumday_zoom->GetXaxis()->SetLabelSize(0.03);
    gnPOTSumday_zoom->GetXaxis()->SetTimeFormat(taxis_labels.c_str());
    gnPOTSumday_zoom->GetXaxis()->SetLimits(time_ago,XNow);
    gnPOTSumday_zoom->GetYaxis()->SetLabelSize(0.03);
    gnPOTSumday_zoom->GetYaxis()->SetTitle("(POT)");
    gnPOTSumday_zoom->GetXaxis()->SetTitle("(central time)");
    gnPOTSumday_zoom->GetYaxis()->SetRangeUser((gnPOTSumday_zoom->GetMean(2)-3.*gnPOTSumday_zoom->GetRMS(2)),(gnPOTSumday_zoom->GetMean(2)+3.*gnPOTSumday_zoom->GetRMS(2)));
    gnPOTSumday_zoom->Draw("A*");
    UpdateText->Draw();
    LastPoint->Draw();
    sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1unPOTSumperSRzoom%s.png",det.c_str(),trig.c_str(),p,period.c_str());
    ctnPOTSumday_zoom->Print(filename);
    delete ctnPOTSumday_zoom;
    ctnPOTSumday_zoom = 0;
    
    maxtime = 0;
    max = 0.0;
    ave = 0.0;
    for(int i = 0; i < NspillPOTDayCount[p]; ++i) {
      ave += (double)NspillPOTDay[p][i];
      if(NspillPOTtimeDay[p][i] > maxtime) {
	maxtime = NspillPOTtimeDay[p][i];
	max     = NspillPOTDay[p][i];
      }
    }
    if(NspillPOTDayCount[p] > 0) ave = ave/(double)NspillPOTDayCount[p];
    sprintf(lptext,"Last Point = %f  /  Average = %f",max,ave);
    LastPoint->Clear();
    LastPoint->AddText(lptext);
    TCanvas *ctnspillPOTday = new TCanvas("cnspillPOTday","POT Sum day",1200,800);
    ctnspillPOTday->cd();
    gPad->SetGridx();
    TGraph *gnspillPOTday = new TGraph(NspillPOTDayCount[p],NspillPOTtimeDay[p],NspillPOTDay[p]);
    sprintf(title,"%s Spill POT per Subrun - partition %.1u",trigname.c_str(),p);
    gnspillPOTday->SetTitle(title);
    gnspillPOTday->SetMarkerColor(kBlue);
    gnspillPOTday->GetXaxis()->SetTimeDisplay(1);
    gnspillPOTday->GetXaxis()->SetLabelSize(0.03);
    gnspillPOTday->GetXaxis()->SetTimeFormat(taxis_labels.c_str());
    gnspillPOTday->GetXaxis()->SetLimits(time_ago,XNow);
    gnspillPOTday->GetYaxis()->SetLabelSize(0.03);
    gnspillPOTday->GetXaxis()->SetTitle("(central time)");
    gnspillPOTday->GetYaxis()->SetRangeUser(0,50e12);
    gnspillPOTday->Draw("A*");
    TLine* l2SpillPOT=new TLine(time_ago,2e12,XNow,2e12);
    l2SpillPOT->SetLineColor(kRed);
    l2SpillPOT->SetLineWidth(3);
    l2SpillPOT->SetLineStyle(9);
    l2SpillPOT->Draw();
    UpdateText->Draw();
    LastPoint->Draw();
    sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1unspillPOTperSR%s.png",det.c_str(),trig.c_str(),p,period.c_str());
    ctnspillPOTday->Print(filename);
    delete ctnspillPOTday;
    ctnspillPOTday = 0;

    // Adding "zoomed in" version of above plot
    TCanvas *ctnspillPOTday_zoom = new TCanvas("cnspillPOTday_zoom","POT Sum day",1200,800);
    ctnspillPOTday_zoom->cd();
    gPad->SetGridx();
    TGraph *gnspillPOTday_zoom = new TGraph(NspillPOTDayCount[p],NspillPOTtimeDay[p],NspillPOTDay[p]);
    sprintf(title,"%s Spill POT per subrun - partition %.1u",trigname.c_str(),p);
    gnspillPOTday_zoom->SetTitle(title);
    gnspillPOTday_zoom->SetMarkerColor(kBlue);
    gnspillPOTday_zoom->GetXaxis()->SetTimeDisplay(1);
    gnspillPOTday_zoom->GetXaxis()->SetLabelSize(0.03);
    gnspillPOTday_zoom->GetXaxis()->SetTimeFormat(taxis_labels.c_str());
    gnspillPOTday_zoom->GetXaxis()->SetLimits(time_ago,XNow);
    gnspillPOTday_zoom->GetYaxis()->SetLabelSize(0.03);
    gnspillPOTday_zoom->GetXaxis()->SetTitle("(central time)");
    gnspillPOTday_zoom->GetYaxis()->SetRangeUser((gnspillPOTday_zoom->GetMean(2)-3.*gnspillPOTday_zoom->GetRMS(2)),(gnspillPOTday_zoom->GetMean(2)+3.*gnspillPOTday_zoom->GetRMS(2)));
    gnspillPOTday_zoom->Draw("A*");
    UpdateText->Draw();
    LastPoint->Draw();
    sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1unspillPOTperSRzoom%s.png",det.c_str(),trig.c_str(),p,period.c_str());
    ctnspillPOTday_zoom->Print(filename);
    delete ctnspillPOTday_zoom;
    ctnspillPOTday_zoom = 0;
  
    maxtime = 0;
    max = 0; 
    ave = 0;
   
    for(int i = 0; i < NHornCurrentDayCount[p]; ++i) {
      ave += (double)NHornCurrentDay[p][i];
      if(NHornCurrenttimeDay[p][i] > maxtime) {
	maxtime = NHornCurrenttimeDay[p][i];
	max     = NHornCurrentDay[p][i];
      }
    }
    if(NHornCurrentDayCount[p] > 0) ave = ave/(double)NHornCurrentDayCount[p];
    sprintf(lptext,"Last Point = %f  /  Average = %f",max,ave);
    LastPoint->Clear();
    LastPoint->AddText(lptext);
    TCanvas *ctnHornCurrentday = new TCanvas("cnHornCurrentday","Horn Current day",1200,800);
    ctnHornCurrentday->cd();
    gPad->SetGridx();
    TGraph *gnHornCurrentday = new TGraph(NHornCurrentDayCount[p],NHornCurrenttimeDay[p],NHornCurrentDay[p]);
    sprintf(title,"%s Horn Current per Subrun - partition %.1u",trigname.c_str(),p);
    gnHornCurrentday->SetTitle(title);
    gnHornCurrentday->SetMarkerColor(kBlue);
    gnHornCurrentday->GetXaxis()->SetTimeDisplay(1);
    gnHornCurrentday->GetXaxis()->SetLabelSize(0.03);
    gnHornCurrentday->GetXaxis()->SetTimeFormat(taxis_labels.c_str());
    gnHornCurrentday->GetXaxis()->SetLimits(time_ago,XNow);
    gnHornCurrentday->GetXaxis()->SetTitle("(central time)");
    gnHornCurrentday->GetYaxis()->SetTitle("(kA)");
    gnHornCurrentday->GetYaxis()->SetRangeUser(-203,-197);
    gnHornCurrentday->Draw("A*");
    TLine* l2cHornCurrent=new TLine(time_ago,-202,XNow,-202);
    l2cHornCurrent->SetLineColor(kRed);
    l2cHornCurrent->SetLineWidth(3);
    l2cHornCurrent->SetLineStyle(9);
    l2cHornCurrent->Draw();
    TLine* l2bHornCurrent=new TLine(time_ago,-198,XNow,-198);
    l2bHornCurrent->SetLineColor(kRed);
    l2bHornCurrent->SetLineWidth(3);
    l2bHornCurrent->SetLineStyle(9);
    l2bHornCurrent->Draw();
    UpdateText->Draw();
    LastPoint->Draw();
    UpdateText->Draw();
    LastPoint->Draw();
    sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1unHornCurrentperSR%s.png",det.c_str(),trig.c_str(),p,period.c_str());
    ctnHornCurrentday->Print(filename);
    delete ctnHornCurrentday;
    ctnHornCurrentday = 0;

    // Adding "zoomed in" version of above plot
    TCanvas *ctnHornCurrentday_zoom = new TCanvas("cnHornCurrentday_zoom","Horn Current day",1200,800);
    ctnHornCurrentday_zoom->cd();
    gPad->SetGridx();
    TGraph *gnHornCurrentday_zoom = new TGraph(NHornCurrentDayCount[p],NHornCurrenttimeDay[p],NHornCurrentDay[p]);
    sprintf(title,"%s Horn Current per subrun - partition %.1u",trigname.c_str(),p);
    gnHornCurrentday_zoom->SetTitle(title);
    gnHornCurrentday_zoom->SetMarkerColor(kBlue);
    gnHornCurrentday_zoom->GetXaxis()->SetTimeDisplay(1);
    gnHornCurrentday_zoom->GetXaxis()->SetLabelSize(0.03);
    gnHornCurrentday_zoom->GetXaxis()->SetTimeFormat(taxis_labels.c_str());
    gnHornCurrentday_zoom->GetXaxis()->SetLimits(time_ago,XNow);
    gnHornCurrentday_zoom->GetXaxis()->SetTitle("(central time)");
    gnHornCurrentday_zoom->GetYaxis()->SetTitle("(kA)");
    gnHornCurrentday_zoom->GetYaxis()->SetRangeUser((gnHornCurrentday_zoom->GetMean(2)-3.*gnHornCurrentday_zoom->GetRMS(2)),(gnHornCurrentday_zoom->GetMean(2)+3.*gnHornCurrentday_zoom->GetRMS(2)));
    gnHornCurrentday_zoom->Draw("A*");
    UpdateText->Draw();
    LastPoint->Draw();
    sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1unHornCurrentperSRzoom%s.png",det.c_str(),trig.c_str(),p,period.c_str());
    ctnHornCurrentday_zoom->Print(filename);
    delete ctnHornCurrentday_zoom;
    ctnHornCurrentday_zoom = 0;

    maxtime = 0;
    max = 0.0;
    ave = 0.0;
    for(int i = 0; i < NXPositionDayCount[p]; ++i) {
      ave += (double)NXPositionDay[p][i];
      if(NXPositiontimeDay[p][i] > maxtime) {
	maxtime = NXPositiontimeDay[p][i];
	max     = NXPositionDay[p][i];
      }
    }
    if(NXPositionDayCount[p] > 0) ave = ave/(double)NXPositionDayCount[p];
    sprintf(lptext,"Last Point = %f  /  Average = %f",max,ave);
    LastPoint->Clear();
    LastPoint->AddText(lptext);
    TCanvas *ctnXPositionday = new TCanvas("cnXPositionday","X Position day",1200,800);
    ctnXPositionday->cd();
    gPad->SetGridx();
    TGraph *gnXPositionday = new TGraph(NXPositionDayCount[p],NXPositiontimeDay[p],NXPositionDay[p]);
    sprintf(title,"%s X Position per Subrun - partition %.1u",trigname.c_str(),p);
    gnXPositionday->SetTitle(title);
    gnXPositionday->SetMarkerColor(kBlue);
    gnXPositionday->GetXaxis()->SetTimeDisplay(1);
    gnXPositionday->GetXaxis()->SetLabelSize(0.03);
    gnXPositionday->GetXaxis()->SetTimeFormat(taxis_labels.c_str());
    gnXPositionday->GetXaxis()->SetLimits(time_ago,XNow);
    gnXPositionday->GetXaxis()->SetTitle("(central time)");
    gnXPositionday->GetYaxis()->SetTitle("(mm)");
    gnXPositionday->GetYaxis()->SetRangeUser(0,2.1);
    gnXPositionday->Draw("A*");
    TLine* l2cXPosition=new TLine(time_ago,0.02,XNow,0.02);
    l2cXPosition->SetLineColor(kRed);
    l2cXPosition->SetLineWidth(3);
    l2cXPosition->SetLineStyle(9);
    l2cXPosition->Draw();
    TLine* l2bXPosition=new TLine(time_ago,2.0,XNow,2.0);
    l2bXPosition->SetLineColor(kRed);
    l2bXPosition->SetLineWidth(3);
    l2bXPosition->SetLineStyle(9);
    l2bXPosition->Draw();
    UpdateText->Draw();
    LastPoint->Draw();
    sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1unXPositionperSR%s.png",det.c_str(),trig.c_str(),p,period.c_str());
    ctnXPositionday->Print(filename);
    delete ctnXPositionday;
    ctnXPositionday = 0;

    // Adding "zoomed in" version of above plot
    TCanvas *ctnXPositionday_zoom = new TCanvas("cnXPositionday_zoom","X Position day",1200,800);
    ctnXPositionday_zoom->cd();
    gPad->SetGridx();
    TGraph *gnXPositionday_zoom = new TGraph(NXPositionDayCount[p],NXPositiontimeDay[p],NXPositionDay[p]);
    sprintf(title,"%s X Position per subrun - partition %.1u",trigname.c_str(),p);
    gnXPositionday_zoom->SetTitle(title);
    gnXPositionday_zoom->SetMarkerColor(kBlue);
    gnXPositionday_zoom->GetXaxis()->SetTimeDisplay(1);
    gnXPositionday_zoom->GetXaxis()->SetLabelSize(0.03);
    gnXPositionday_zoom->GetXaxis()->SetTimeFormat(taxis_labels.c_str());
    gnXPositionday_zoom->GetXaxis()->SetLimits(time_ago,XNow);
    gnXPositionday_zoom->GetXaxis()->SetTitle("(central time)");
    gnXPositionday_zoom->GetYaxis()->SetTitle("(mm)");
    gnXPositionday_zoom->GetYaxis()->SetRangeUser((gnXPositionday_zoom->GetMean(2)-3.*gnXPositionday_zoom->GetRMS(2)),(gnXPositionday_zoom->GetMean(2)+3.*gnXPositionday_zoom->GetRMS(2)));
    gnXPositionday_zoom->Draw("A*");
    UpdateText->Draw();
    LastPoint->Draw();
    sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1unXPositionperSRzoom%s.png",det.c_str(),trig.c_str(),p,period.c_str());
    ctnXPositionday_zoom->Print(filename);
    delete ctnXPositionday_zoom;
    ctnXPositionday_zoom = 0;
    
    maxtime = 0;
    max = 0.0;
    ave = 0.0;
    for(int i = 0; i < NYPositionDayCount[p]; ++i) {
      ave += (double)NYPositionDay[p][i];
      if(NYPositiontimeDay[p][i] > maxtime) {
	maxtime = NYPositiontimeDay[p][i];
	max     = NYPositionDay[p][i];
      }
    }
    if(NYPositionDayCount[p] > 0) ave = ave/(double)NYPositionDayCount[p];
    sprintf(lptext,"Last Point = %f  /  Average = %f",max,ave);
    LastPoint->Clear();
    LastPoint->AddText(lptext);
    TCanvas *ctnYPositionday = new TCanvas("cnYPositionday","Y Position day",1200,800);
    ctnYPositionday->cd();
    gPad->SetGridx();
    TGraph *gnYPositionday = new TGraph(NYPositionDayCount[p],NYPositiontimeDay[p],NYPositionDay[p]);
    sprintf(title,"%s Y Position per Subrun - partition %.1u",trigname.c_str(),p);
    gnYPositionday->SetTitle(title);
    gnYPositionday->SetMarkerColor(kBlue);
    gnYPositionday->GetXaxis()->SetTimeDisplay(1);
    gnYPositionday->GetXaxis()->SetLabelSize(0.03);
    gnYPositionday->GetXaxis()->SetTimeFormat(taxis_labels.c_str());
    gnYPositionday->GetXaxis()->SetLimits(time_ago,XNow);
    gnYPositionday->GetXaxis()->SetTitle("(central time)");
    gnYPositionday->GetYaxis()->SetTitle("(mm)");
    gnYPositionday->GetYaxis()->SetRangeUser(0,2.1);
    gnYPositionday->Draw("A*");
    TLine* l2cYPosition=new TLine(time_ago,0.02,XNow,0.02);
    l2cYPosition->SetLineColor(kRed);
    l2cYPosition->SetLineWidth(3);
    l2cYPosition->SetLineStyle(9);
    l2cYPosition->Draw();
    TLine* l2bYPosition=new TLine(time_ago,2.0,XNow,2.0);
    l2bYPosition->SetLineColor(kRed);
    l2bYPosition->SetLineWidth(3);
    l2bYPosition->SetLineStyle(9);
    l2bYPosition->Draw();
    UpdateText->Draw();
    LastPoint->Draw();
    sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1unYPositionperSR%s.png",det.c_str(),trig.c_str(),p,period.c_str());
    ctnYPositionday->Print(filename);
    delete ctnYPositionday;
    ctnYPositionday = 0;

    // Adding "zoomed in" version of above plot
    TCanvas *ctnYPositionday_zoom = new TCanvas("cnYPositionday_zoom","Y Position day",1200,800);
    ctnYPositionday_zoom->cd();
    gPad->SetGridx();
    TGraph *gnYPositionday_zoom = new TGraph(NYPositionDayCount[p],NYPositiontimeDay[p],NYPositionDay[p]);
    sprintf(title,"%s Y Position per subrun - partition %.1u",trigname.c_str(),p);
    gnYPositionday_zoom->SetTitle(title);
    gnYPositionday_zoom->SetMarkerColor(kBlue);
    gnYPositionday_zoom->GetXaxis()->SetTimeDisplay(1);
    gnYPositionday_zoom->GetXaxis()->SetLabelSize(0.03);
    gnYPositionday_zoom->GetXaxis()->SetTimeFormat(taxis_labels.c_str());
    gnYPositionday_zoom->GetXaxis()->SetLimits(time_ago,XNow);
    gnYPositionday_zoom->GetXaxis()->SetTitle("(central time)");
    gnYPositionday_zoom->GetYaxis()->SetTitle("(mm)");
    gnYPositionday_zoom->GetYaxis()->SetRangeUser((gnYPositionday_zoom->GetMean(2)-3.*gnYPositionday_zoom->GetRMS(2)),(gnYPositionday_zoom->GetMean(2)+3.*gnYPositionday_zoom->GetRMS(2)));
    gnYPositionday_zoom->Draw("A*");
    UpdateText->Draw();
    LastPoint->Draw();
    sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1unYPositionperSRzoom%s.png",det.c_str(),trig.c_str(),p,period.c_str());
    ctnYPositionday_zoom->Print(filename);
    delete ctnYPositionday_zoom;
    ctnYPositionday_zoom = 0;

    maxtime = 0;
    max = 0.0;
    ave = 0.0;
    for(int i = 0; i < NXWidthDayCount[p]; ++i) {
      ave += (double)NXWidthDay[p][i];
      if(NXWidthtimeDay[p][i] > maxtime) {
	maxtime = NXWidthtimeDay[p][i];
	max     = NXWidthDay[p][i];
      }
    }
    if(NXWidthDayCount[p] > 0) ave = ave/(double)NXWidthDayCount[p];
    sprintf(lptext,"Last Point = %f  /  Average = %f",max,ave);
    LastPoint->Clear();
    LastPoint->AddText(lptext);
    TCanvas *ctnXWidthday = new TCanvas("cnXWidthday","X Width day",1200,800);
    ctnXWidthday->cd();
    gPad->SetGridx();
    TGraph *gnXWidthday = new TGraph(NXWidthDayCount[p],NXWidthtimeDay[p],NXWidthDay[p]);
    sprintf(title,"%s X Width per Subrun - partition %.1u",trigname.c_str(),p);
    gnXWidthday->SetTitle(title);
    gnXWidthday->SetMarkerColor(kBlue);
    gnXWidthday->GetXaxis()->SetTimeDisplay(1);
    gnXWidthday->GetXaxis()->SetLabelSize(0.03);
    gnXWidthday->GetXaxis()->SetTimeFormat(taxis_labels.c_str());
    gnXWidthday->GetXaxis()->SetLimits(time_ago,XNow);
    gnXWidthday->GetXaxis()->SetTitle("(central time)");
    gnXWidthday->GetYaxis()->SetTitle("(mm)");
    gnXWidthday->GetYaxis()->SetRangeUser(0.5,1.7);
    gnXWidthday->Draw("A*");
    TLine* l2cXWidth=new TLine(time_ago,0.57,XNow,0.57);
    l2cXWidth->SetLineColor(kRed);
    l2cXWidth->SetLineWidth(3);
    l2cXWidth->SetLineStyle(9);
    l2cXWidth->Draw();
    TLine* l2bXWidth=new TLine(time_ago,1.58,XNow,1.58);
    l2bXWidth->SetLineColor(kRed);
    l2bXWidth->SetLineWidth(3);
    l2bXWidth->SetLineStyle(9);
    l2bXWidth->Draw();
    UpdateText->Draw();
    LastPoint->Draw();
    sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1unXWidthperSR%s.png",det.c_str(),trig.c_str(),p,period.c_str());
    ctnXWidthday->Print(filename);
    delete ctnXWidthday;
    ctnXWidthday = 0;

    // Adding "zoomed in" version of above plot
    TCanvas *ctnXWidthday_zoom = new TCanvas("cnXWidthday_zoom","X Width day",1200,800);
    ctnXWidthday_zoom->cd();
    gPad->SetGridx();
    TGraph *gnXWidthday_zoom = new TGraph(NXWidthDayCount[p],NXWidthtimeDay[p],NXWidthDay[p]);
    sprintf(title,"%s X Width per subrun - partition %.1u",trigname.c_str(),p);
    gnXWidthday_zoom->SetTitle(title);
    gnXWidthday_zoom->SetMarkerColor(kBlue);
    gnXWidthday_zoom->GetXaxis()->SetTimeDisplay(1);
    gnXWidthday_zoom->GetXaxis()->SetLabelSize(0.03);
    gnXWidthday_zoom->GetXaxis()->SetTimeFormat(taxis_labels.c_str());
    gnXWidthday_zoom->GetXaxis()->SetLimits(time_ago,XNow);
    gnXWidthday_zoom->GetXaxis()->SetTitle("(central time)");
    gnXWidthday_zoom->GetYaxis()->SetTitle("(mm)");
    gnXWidthday_zoom->GetYaxis()->SetRangeUser((gnXWidthday_zoom->GetMean(2)-3.*gnXWidthday_zoom->GetRMS(2)),(gnXWidthday_zoom->GetMean(2)+3.*gnXWidthday_zoom->GetRMS(2)));
    gnXWidthday_zoom->Draw("A*");
    UpdateText->Draw();
    LastPoint->Draw();
    sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1unXWidthperSRzoom%s.png",det.c_str(),trig.c_str(),p,period.c_str());
    ctnXWidthday_zoom->Print(filename);
    delete ctnXWidthday_zoom;
    ctnXWidthday_zoom = 0;
    
    maxtime = 0;
    max = 0.0;
    ave = 0.0;
    for(int i = 0; i < NYWidthDayCount[p]; ++i) {
      ave += (double)NYWidthDay[p][i];
      if(NYWidthtimeDay[p][i] > maxtime) {
	maxtime = NYWidthtimeDay[p][i];
	max     = NYWidthDay[p][i];
      }
    }
    if(NYWidthDayCount[p] > 0) ave = ave/(double)NYWidthDayCount[p];
    sprintf(lptext,"Last Point = %f  /  Average = %f",max,ave);
    LastPoint->Clear();
    LastPoint->AddText(lptext);
    TCanvas *ctnYWidthday = new TCanvas("cnYWidthday","Y Width day",1200,800);
    ctnYWidthday->cd();
    gPad->SetGridx();
    TGraph *gnYWidthday = new TGraph(NYWidthDayCount[p],NYWidthtimeDay[p],NYWidthDay[p]);
    sprintf(title,"%s Y Width per Subrun - partition %.1u",trigname.c_str(),p);
    gnYWidthday->SetTitle(title);
    gnYWidthday->SetMarkerColor(kBlue);
    gnYWidthday->GetXaxis()->SetTimeDisplay(1);
    gnYWidthday->GetXaxis()->SetLabelSize(0.03);
    gnYWidthday->GetXaxis()->SetTimeFormat(taxis_labels.c_str());
    gnYWidthday->GetXaxis()->SetLimits(time_ago,XNow);
    gnYWidthday->GetXaxis()->SetTitle("(central time)");
    gnYWidthday->GetYaxis()->SetTitle("(mm)");
    gnYWidthday->GetYaxis()->SetRangeUser(0.5,1.7);
    gnYWidthday->Draw("A*");
    TLine* l2cYWidth=new TLine(time_ago,0.57,XNow,0.57);
    l2cYWidth->SetLineColor(kRed);
    l2cYWidth->SetLineWidth(3);
    l2cYWidth->SetLineStyle(9);
    l2cYWidth->Draw();
    TLine* l2bYWidth=new TLine(time_ago,1.58,XNow,1.58);
    l2bYWidth->SetLineColor(kRed);
    l2bYWidth->SetLineWidth(3);
    l2bYWidth->SetLineStyle(9);
    l2bYWidth->Draw();
    UpdateText->Draw();
    LastPoint->Draw();
    sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1unYWidthperSR%s.png",det.c_str(),trig.c_str(),p,period.c_str());
    ctnYWidthday->Print(filename);
    delete ctnYWidthday;
    ctnYWidthday = 0;

    // Adding "zoomed in" version of above plot
    TCanvas *ctnYWidthday_zoom = new TCanvas("cnYWidthday_zoom","Y Width day",1200,800);
    ctnYWidthday_zoom->cd();
    gPad->SetGridx();
    TGraph *gnYWidthday_zoom = new TGraph(NYWidthDayCount[p],NYWidthtimeDay[p],NYWidthDay[p]);
    sprintf(title,"%s Y Width per subrun - partition %.1u",trigname.c_str(),p);
    gnYWidthday_zoom->SetTitle(title);
    gnYWidthday_zoom->SetMarkerColor(kBlue);
    gnYWidthday_zoom->GetXaxis()->SetTimeDisplay(1);
    gnYWidthday_zoom->GetXaxis()->SetLabelSize(0.03);
    gnYWidthday_zoom->GetXaxis()->SetTimeFormat(taxis_labels.c_str());
    gnYWidthday_zoom->GetXaxis()->SetLimits(time_ago,XNow);
    gnYWidthday_zoom->GetXaxis()->SetTitle("(central time)");
    gnYWidthday_zoom->GetYaxis()->SetTitle("(mm)");
    gnYWidthday_zoom->GetYaxis()->SetRangeUser((gnYWidthday_zoom->GetMean(2)-3.*gnYWidthday_zoom->GetRMS(2)),(gnYWidthday_zoom->GetMean(2)+3.*gnYWidthday_zoom->GetRMS(2)));
    gnYWidthday_zoom->Draw("A*");
    UpdateText->Draw();
    LastPoint->Draw();
    sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1unYWidthperSRzoom%s.png",det.c_str(),trig.c_str(),p,period.c_str());
    ctnYWidthday_zoom->Print(filename);
    delete ctnYWidthday_zoom;
    ctnYWidthday_zoom = 0;
    
    maxtime = 0;
    max = 0.0;
    ave = 0.0;
    for(int i = 0; i < NGoodBeamDayCount[p]; ++i) {
      ave += (double)NGoodBeamDay[p][i];
      if(NGoodBeamtimeDay[p][i] > maxtime) {
	maxtime = NGoodBeamtimeDay[p][i];
	max     = NGoodBeamDay[p][i];
      }
    }
    if(NGoodBeamDayCount[p] > 0) ave = ave/(double)NGoodBeamDayCount[p];
    sprintf(lptext,"Last Point = %f  /  Average = %f",max,ave);
    LastPoint->Clear();
    LastPoint->AddText(lptext);
    TCanvas *ctnGoodBeamday = new TCanvas("cnGoodBeamday","Good Beam day",1200,800);
    ctnGoodBeamday->cd();
    gPad->SetGridx();
    TGraph *gnGoodBeamday = new TGraph(NGoodBeamDayCount[p],NGoodBeamtimeDay[p],NGoodBeamDay[p]);
    sprintf(title,"%s Good Beam per Subrun - partition %.1u",trigname.c_str(),p);
    gnGoodBeamday->SetTitle(title);
    gnGoodBeamday->SetMarkerColor(kBlue);
    gnGoodBeamday->GetXaxis()->SetTimeDisplay(1);
    gnGoodBeamday->GetXaxis()->SetLabelSize(0.03);
    gnGoodBeamday->GetXaxis()->SetTimeFormat(taxis_labels.c_str());
    gnGoodBeamday->GetXaxis()->SetLimits(time_ago,XNow);
    gnGoodBeamday->GetXaxis()->SetTitle("(central time)");
    gnGoodBeamday->GetYaxis()->SetRangeUser(-0.5,1.5); 
    gnGoodBeamday->Draw("A*");
    UpdateText->Draw();
    LastPoint->Draw();
    sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1unGoodBeamperSR%s.png",det.c_str(),trig.c_str(),p,period.c_str());
    ctnGoodBeamday->Print(filename);
    delete ctnGoodBeamday;
    ctnGoodBeamday = 0;

    // Adding "zoomed in" version of above plot
    TCanvas *ctnGoodBeamday_zoom = new TCanvas("cnGoodBeamday_zoom","Good Beam day",1200,800);
    ctnGoodBeamday_zoom->cd();
    gPad->SetGridx();
    TGraph *gnGoodBeamday_zoom = new TGraph(NGoodBeamDayCount[p],NGoodBeamtimeDay[p],NGoodBeamDay[p]);
    sprintf(title,"%s Good Beam per subrun - partition %.1u",trigname.c_str(),p);
    gnGoodBeamday_zoom->SetTitle(title);
    gnGoodBeamday_zoom->SetMarkerColor(kBlue);
    gnGoodBeamday_zoom->GetXaxis()->SetTimeDisplay(1);
    gnGoodBeamday_zoom->GetXaxis()->SetLabelSize(0.03);
    gnGoodBeamday_zoom->GetXaxis()->SetTimeFormat(taxis_labels.c_str());
    gnGoodBeamday_zoom->GetXaxis()->SetLimits(time_ago,XNow);
    gnGoodBeamday_zoom->GetXaxis()->SetTitle("(central time)");
    gnGoodBeamday_zoom->GetYaxis()->SetRangeUser((gnGoodBeamday_zoom->GetMean(2)-3.*gnGoodBeamday_zoom->GetRMS(2)),(gnGoodBeamday_zoom->GetMean(2)+3.*gnGoodBeamday_zoom->GetRMS(2)));
    gnGoodBeamday_zoom->Draw("A*");
    UpdateText->Draw();
    LastPoint->Draw();
    sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1unGoodBeamperSRzoom%s.png",det.c_str(),trig.c_str(),p,period.c_str());
    ctnGoodBeamday_zoom->Print(filename);
    delete ctnGoodBeamday_zoom;
    ctnGoodBeamday_zoom = 0;
    //end new beam plots
    
    maxtime = 0;
    max = 0.0;
    ave = 0.0;
    for(int i = 0; i < NSliceDayCount[p]; ++i) {
      ave += (double)NSliceDay[p][i];
      if(NSlicetimeDay[p][i] > maxtime) {
	maxtime = NSlicetimeDay[p][i];
	max     = NSliceDay[p][i];
      }
    }
    if(NSliceDayCount[p] > 0) ave = ave/(double)NSliceDayCount[p];
    sprintf(lptext,"Last Point = %f  /  Average = %f",max,ave);
    LastPoint->Clear();
    LastPoint->AddText(lptext);
    TCanvas *ctnSliceday = new TCanvas("cnSliceday","# Non-Noise Slices day",1200,800);
    ctnSliceday->cd();
    gPad->SetGridx();
    TGraph *gnSliceday = new TGraph(NSliceDayCount[p],NSlicetimeDay[p],NSliceDay[p]);
    sprintf(title,"%s Number of Slices per Subrun - partition %.1u",trigname.c_str(),p);
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
    delete ctnSliceday;
    ctnSliceday = 0;

    // Adding "zoomed in" version of above plot
    TCanvas *ctnSliceday_zoom = new TCanvas("cnSliceday_zoom","# Non-Noise Slices day",1200,800);
    ctnSliceday_zoom->cd();
    gPad->SetGridx();
    TGraph *gnSliceday_zoom = new TGraph(NSliceDayCount[p],NSlicetimeDay[p],NSliceDay[p]);
    sprintf(title,"%s Number of Slices per Subrun - partition %.1u",trigname.c_str(),p);
    gnSliceday_zoom->SetTitle(title);
    gnSliceday_zoom->SetMarkerColor(kBlue);
    gnSliceday_zoom->GetXaxis()->SetTimeDisplay(1);
    gnSliceday_zoom->GetXaxis()->SetLabelSize(0.03);
    gnSliceday_zoom->GetXaxis()->SetTimeFormat(taxis_labels.c_str());
    gnSliceday_zoom->GetXaxis()->SetLimits(time_ago,XNow);
    gnSliceday_zoom->GetXaxis()->SetTitle("(central time)");
    gnSliceday_zoom->GetYaxis()->SetRangeUser((gnSliceday_zoom->GetMean(2)-3.*gnSliceday_zoom->GetRMS(2)),(gnSliceday_zoom->GetMean(2)+3.*gnSliceday_zoom->GetRMS(2)));
    gnSliceday_zoom->Draw("A*");
    UpdateText->Draw();
    LastPoint->Draw();
    sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1unSliceperSRzoom%s.png",det.c_str(),trig.c_str(),p,period.c_str());
    ctnSliceday_zoom->Print(filename);
    delete ctnSliceday_zoom;
    ctnSliceday_zoom = 0;

    maxtime = 0;
    max = 0.0;
    ave = 0.0;
    for(int i = 0; i < NNonNoiseSliceHitDayCount[p]; ++i) {
      ave += (double)NNonNoiseSliceHitDay[p][i];
      if(NNonNoiseSliceHittimeDay[p][i] > maxtime) {
	maxtime = NNonNoiseSliceHittimeDay[p][i];
	max     = NNonNoiseSliceHitDay[p][i];
      }
    }
    if(NNonNoiseSliceHitDayCount[p] > 0) ave = ave/(double)NNonNoiseSliceHitDayCount[p];
    sprintf(lptext,"Last Point = %f  /  Average = %f",max,ave);
    LastPoint->Clear();
    LastPoint->AddText(lptext);
    TCanvas *ctnSliceHitday = new TCanvas("cnSliceHitday","# Non-Noise Slice Hits day",1200,800);
    ctnSliceHitday->cd();
    //ctnSliceHitday->SetLogy();
    gPad->SetGridx();
    TGraph *gnSliceHitday = new TGraph(NNonNoiseSliceHitDayCount[p],NNonNoiseSliceHittimeDay[p],NNonNoiseSliceHitDay[p]);
    sprintf(title,"%s Number of Slice Hits per Subrun - partition %.1u",trigname.c_str(),p);
    gnSliceHitday->SetTitle(title);
    gnSliceHitday->SetMarkerColor(kBlue);
    gnSliceHitday->GetXaxis()->SetTimeDisplay(1);
    gnSliceHitday->GetXaxis()->SetLabelSize(0.03);
    gnSliceHitday->GetXaxis()->SetTimeFormat(taxis_labels.c_str());
    gnSliceHitday->GetXaxis()->SetLimits(time_ago,XNow);
    gnSliceHitday->GetXaxis()->SetTitle("(central time)");
    gnSliceHitday->Draw("A*");
    UpdateText->Draw();
    LastPoint->Draw();
    sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1unSliceHitperSR%s.png",det.c_str(),trig.c_str(),p,period.c_str());
    ctnSliceHitday->Print(filename);
    delete ctnSliceHitday;
    ctnSliceHitday = 0;

    //Adding "zoomed in" version of the above plot
    TCanvas *ctnSliceHitday_zoom = new TCanvas("cnSliceHitday_zoom","# Non-Noise Slice Hits day",1200,800);
    ctnSliceHitday_zoom->cd();
    //ctnSliceHitday_zoom->SetLogy();
    gPad->SetGridx();
    TGraph *gnSliceHitday_zoom = new TGraph(NNonNoiseSliceHitDayCount[p],NNonNoiseSliceHittimeDay[p],NNonNoiseSliceHitDay[p]);
    sprintf(title,"%s Number of Slice Hits per Subrun - partition %.1u",trigname.c_str(),p);
    gnSliceHitday_zoom->SetTitle(title);
    gnSliceHitday_zoom->SetMarkerColor(kBlue);
    gnSliceHitday_zoom->GetXaxis()->SetTimeDisplay(1);
    gnSliceHitday_zoom->GetXaxis()->SetLabelSize(0.03);
    gnSliceHitday_zoom->GetXaxis()->SetTimeFormat(taxis_labels.c_str());
    gnSliceHitday_zoom->GetXaxis()->SetLimits(time_ago,XNow);
    gnSliceHitday_zoom->GetXaxis()->SetTitle("(central time)");
    gnSliceHitday_zoom->GetYaxis()->SetRangeUser((gnSliceHitday_zoom->GetMean(2)-3.*gnSliceHitday_zoom->GetRMS(2)),(gnSliceHitday_zoom->GetMean(2)+3.*gnSliceHitday_zoom->GetRMS(2)));
    gnSliceHitday_zoom->Draw("A*");
    UpdateText->Draw();
    LastPoint->Draw();
    sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1unSliceHitperSRzoom%s.png",det.c_str(),trig.c_str(),p,period.c_str());
    ctnSliceHitday_zoom->Print(filename);
    delete ctnSliceHitday_zoom;
    ctnSliceHitday_zoom = 0;

    maxtime = 0;
    max = 0.0;
    ave = 0.0;
    for(int i = 0; i < NNoiseSliceHitDayCount[p]; ++i) {
      ave += (double)NNoiseSliceHitDay[p][i];
      if(NNoiseSliceHittimeDay[p][i] > maxtime) {
	maxtime = NNoiseSliceHittimeDay[p][i];
	max     = NNoiseSliceHitDay[p][i];
      }
    }
    if(NNoiseSliceHitDayCount[p] > 0) ave = ave/(double)NNoiseSliceHitDayCount[p];
    sprintf(lptext,"Last Point = %f  /  Average = %f",max,ave);
    LastPoint->Clear();
    LastPoint->AddText(lptext);
    TCanvas *ctNNoiseSliceHitday = new TCanvas("cNNoiseSliceHitday","# Noise Slice Hits day",1200,800);
    ctNNoiseSliceHitday->cd();
    gPad->SetGridx();
    TGraph *gNNoiseSliceHitday = new TGraph(NNoiseSliceHitDayCount[p],NNoiseSliceHittimeDay[p],NNoiseSliceHitDay[p]);
    sprintf(title,"%s Number of Noise Slice Hits per Subrun - partition %.1u",trigname.c_str(),p);
    gNNoiseSliceHitday->SetTitle(title);
    gNNoiseSliceHitday->SetMarkerColor(kBlue);
    gNNoiseSliceHitday->GetXaxis()->SetTimeDisplay(1);
    gNNoiseSliceHitday->GetXaxis()->SetLabelSize(0.03);
    gNNoiseSliceHitday->GetXaxis()->SetTimeFormat(taxis_labels.c_str());
    gNNoiseSliceHitday->GetXaxis()->SetLimits(time_ago,XNow);
    gNNoiseSliceHitday->GetXaxis()->SetTitle("(central time)");
    gNNoiseSliceHitday->Draw("A*");
    UpdateText->Draw();
    LastPoint->Draw();
    sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1uNNoiseSliceHitperSR%s.png",det.c_str(),trig.c_str(),p,period.c_str());
    ctNNoiseSliceHitday->Print(filename);
    delete ctNNoiseSliceHitday;
    ctNNoiseSliceHitday = 0;

    // Adding "zoomed in" version of the above plot
    TCanvas *ctNNoiseSliceHitday_zoom = new TCanvas("cNNoiseSliceHitday_zoom","# Noise Slice Hits day",1200,800);
    ctNNoiseSliceHitday_zoom->cd();
    gPad->SetGridx();
    TGraph *gNNoiseSliceHitday_zoom = new TGraph(NNoiseSliceHitDayCount[p],NNoiseSliceHittimeDay[p],NNoiseSliceHitDay[p]);
    sprintf(title,"%s Number of Noise Slice Hits per Subrun - partition %.1u",trigname.c_str(),p);
    gNNoiseSliceHitday_zoom->SetTitle(title);
    gNNoiseSliceHitday_zoom->SetMarkerColor(kBlue);
    gNNoiseSliceHitday_zoom->GetXaxis()->SetTimeDisplay(1);
    gNNoiseSliceHitday_zoom->GetXaxis()->SetLabelSize(0.03);
    gNNoiseSliceHitday_zoom->GetXaxis()->SetTimeFormat(taxis_labels.c_str());
    gNNoiseSliceHitday_zoom->GetXaxis()->SetLimits(time_ago,XNow);
    gNNoiseSliceHitday_zoom->GetXaxis()->SetTitle("(central time)");
    gNNoiseSliceHitday_zoom->GetYaxis()->SetRangeUser((gNNoiseSliceHitday_zoom->GetMean(2)-3.*gNNoiseSliceHitday_zoom->GetRMS(2)),(gNNoiseSliceHitday_zoom->GetMean(2)+3.*gNNoiseSliceHitday_zoom->GetRMS(2)));
    gNNoiseSliceHitday_zoom->Draw("A*");
    UpdateText->Draw();
    LastPoint->Draw();
    sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1uNNoiseSliceHitperSRzoom%s.png",det.c_str(),trig.c_str(),p,period.c_str());
    ctNNoiseSliceHitday_zoom->Print(filename);
    delete ctNNoiseSliceHitday_zoom;
    ctNNoiseSliceHitday_zoom = 0;

    maxtime = 0;
    max = 0.0;
    ave = 0.0;
    for(int i = 0; i < NTsdDayCount[p]; ++i) {
      ave += (double)NTsdDay[p][i];
      if(NTsdtimeDay[p][i] > maxtime) {
	maxtime = NTsdtimeDay[p][i];
	max     = NTsdDay[p][i];
      }
    }
    if(NTsdDayCount[p] > 0) ave = ave/(double)NTsdDayCount[p];
    sprintf(lptext,"Last Point = %f  /  Average = %f",max,ave);
    LastPoint->Clear();
    LastPoint->AddText(lptext);
    TCanvas *ctTsdday = new TCanvas("cTsdday","Slice Time Duration Standard Dev. day",1200,800);
    ctTsdday->cd();
    gPad->SetGridx();
    TGraph *gTsdday = new TGraph(NTsdDayCount[p],NTsdtimeDay[p],NTsdDay[p]);
    sprintf(title,"%s Slice Time Duration Standard Dev. by Subrun - partition %.1u",trigname.c_str(),p);
    gTsdday->SetTitle(title);
    gTsdday->SetMarkerColor(kBlue);
    gTsdday->GetXaxis()->SetTimeDisplay(1);
    gTsdday->GetXaxis()->SetLabelSize(0.03);
    gTsdday->GetXaxis()->SetTimeFormat(taxis_labels.c_str());
    gTsdday->GetXaxis()->SetLimits(time_ago,XNow);
    gTsdday->GetXaxis()->SetTitle("(central time)");
    gTsdday->Draw("A*");
    UpdateText->Draw();
    LastPoint->Draw();
    sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1uTsdperSR%s.png",det.c_str(),trig.c_str(),p,period.c_str());
    ctTsdday->Print(filename);
    delete ctTsdday;
    ctTsdday = 0;

    // Adding "zoomed in" version of the above plot
    TCanvas *ctTsdday_zoom = new TCanvas("cTsdday_zoom","Slice Time Duration Standard Dev. day",1200,800);
    ctTsdday_zoom->cd();
    gPad->SetGridx();
    TGraph *gTsdday_zoom = new TGraph(NTsdDayCount[p],NTsdtimeDay[p],NTsdDay[p]);
    sprintf(title,"%s Slice Time Duration Standard Dev. by Subrun - partition %.1u",trigname.c_str(),p);
    gTsdday_zoom->SetTitle(title);
    gTsdday_zoom->SetMarkerColor(kBlue);
    gTsdday_zoom->GetXaxis()->SetTimeDisplay(1);
    gTsdday_zoom->GetXaxis()->SetLabelSize(0.03);
    gTsdday_zoom->GetXaxis()->SetTimeFormat(taxis_labels.c_str());
    gTsdday_zoom->GetXaxis()->SetLimits(time_ago,XNow);
    gTsdday_zoom->GetXaxis()->SetTitle("(central time)");
    gTsdday_zoom->GetYaxis()->SetRangeUser((gTsdday_zoom->GetMean(2)-3.*gTsdday_zoom->GetRMS(2)),(gTsdday_zoom->GetMean(2)+3.*gTsdday_zoom->GetRMS(2)));
    gTsdday_zoom->Draw("A*");
    UpdateText->Draw();
    LastPoint->Draw();
    sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1uTsdperSRzoom%s.png",det.c_str(),trig.c_str(),p,period.c_str());
    ctTsdday_zoom->Print(filename);
    delete ctTsdday_zoom;
    ctTsdday_zoom = 0;


    maxtime = 0;
    max = 0.0;
    ave = 0.0;
    for(int i = 0; i < NTrackDayCount[p]; ++i) {
      ave += (double)NTrackDay[p][i];
      if(NTracktimeDay[p][i] > maxtime) {
	maxtime = NTracktimeDay[p][i];
	max     = NTrackDay[p][i];
      }
    }
    if(NTrackDayCount[p] > 0) ave = ave/(double)NTrackDayCount[p];
    sprintf(lptext,"Last Point = %f  /  Average = %f",max,ave);
    LastPoint->Clear();
    LastPoint->AddText(lptext);
    TCanvas *ctnTrackday = new TCanvas("cnTrackday","# 3D Tracks day",1200,800);
    ctnTrackday->cd();
    gPad->SetGridx();
    TGraph *gnTrackday = new TGraph(NTrackDayCount[p],NTracktimeDay[p],NTrackDay[p]);
    sprintf(title,"%s # Tracks by Subrun - partition %.1u",trigname.c_str(),p);
    gnTrackday->SetTitle(title);
    gnTrackday->SetMarkerColor(kBlue);
    gnTrackday->GetXaxis()->SetTimeDisplay(1);
    gnTrackday->GetXaxis()->SetLabelSize(0.03);
    gnTrackday->GetXaxis()->SetTimeFormat(taxis_labels.c_str());
    gnTrackday->GetXaxis()->SetLimits(time_ago,XNow);
    gnTrackday->GetXaxis()->SetTitle("(central time)");
    gnTrackday->Draw("A*");
    UpdateText->Draw();
    LastPoint->Draw();
    sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1unTrackAll3DperSR%s.png",det.c_str(),trig.c_str(),p,period.c_str());
    ctnTrackday->Print(filename);
    delete ctnTrackday;
    ctnTrackday = 0;

    // Adding "zoomed in" version of the above plot
    TCanvas *ctnTrackday_zoom = new TCanvas("cnTrackday_zoom","# 3D Tracks day",1200,800);
    ctnTrackday_zoom->cd();
    gPad->SetGridx();
    TGraph *gnTrackday_zoom = new TGraph(NTrackDayCount[p],NTracktimeDay[p],NTrackDay[p]);
    sprintf(title,"%s # Tracks by Subrun - partition %.1u",trigname.c_str(),p);
    gnTrackday_zoom->SetTitle(title);
    gnTrackday_zoom->SetMarkerColor(kBlue);
    gnTrackday_zoom->GetXaxis()->SetTimeDisplay(1);
    gnTrackday_zoom->GetXaxis()->SetLabelSize(0.03);
    gnTrackday_zoom->GetXaxis()->SetTimeFormat(taxis_labels.c_str());
    gnTrackday_zoom->GetXaxis()->SetLimits(time_ago,XNow);
    gnTrackday_zoom->GetXaxis()->SetTitle("(central time)");
    gnTrackday_zoom->GetYaxis()->SetRangeUser((gnTrackday_zoom->GetMean(2)-3.*gnTrackday_zoom->GetRMS(2)),(gnTrackday_zoom->GetMean(2)+3.*gnTrackday_zoom->GetRMS(2)));
    gnTrackday_zoom->Draw("A*");
    UpdateText->Draw();
    LastPoint->Draw();
    sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1unTrackAll3DperSRzoom%s.png",det.c_str(),trig.c_str(),p,period.c_str());
    ctnTrackday_zoom->Print(filename);
    delete ctnTrackday_zoom;
    ctnTrackday_zoom = 0;

    maxtime = 0;
    max = 0.0;
    ave = 0.0;
    for(int i = 0; i < NTrackLenDayCount[p]; ++i) {
      ave += (double)NTrackLenDay[p][i];
      if(NTrackLentimeDay[p][i] > maxtime) {
	maxtime = NTrackLentimeDay[p][i];
	max     = NTrackLenDay[p][i];
      }
    }
    if(NTrackLenDayCount[p] > 0) ave = ave/(double)NTrackLenDayCount[p];
    sprintf(lptext,"Last Point = %f  /  Average = %f",max,ave);
    LastPoint->Clear();
    LastPoint->AddText(lptext);
    TCanvas *ctTrackLenday = new TCanvas("cTrackLenday","3D Track Len day",1200,800);
    ctTrackLenday->cd();
    gPad->SetGridx();
    TGraph *gTrackLenday = new TGraph(NTrackLenDayCount[p],NTrackLentimeDay[p],NTrackLenDay[p]);
    sprintf(title,"%s Ave. Track Length by Subrun - partition %.1u",trigname.c_str(),p);
    gTrackLenday->SetTitle(title);
    gTrackLenday->SetMarkerColor(kBlue);
    gTrackLenday->GetXaxis()->SetTimeDisplay(1);
    gTrackLenday->GetXaxis()->SetLabelSize(0.03);
    gTrackLenday->GetXaxis()->SetTimeFormat(taxis_labels.c_str());
    gTrackLenday->GetXaxis()->SetLimits(time_ago,XNow);
    gTrackLenday->GetXaxis()->SetTitle("(central time)");
    gTrackLenday->Draw("A*");
    UpdateText->Draw();
    LastPoint->Draw();
    sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1uTrackLenAll3DperSR%s.png",det.c_str(),trig.c_str(),p,period.c_str());
    ctTrackLenday->Print(filename);
    delete ctTrackLenday;
    ctTrackLenday = 0;

    // New "zoomed-in" version of TrackLenday above...
    TCanvas *ctTrackLenday_zoom = new TCanvas("cTrackLenday_zoom","3D Track Len day",1200,800);
    ctTrackLenday_zoom->cd();
    gPad->SetGridx();
    TGraph *gTrackLenday_zoom = new TGraph(NTrackLenDayCount[p],NTrackLentimeDay[p],NTrackLenDay[p]);
    sprintf(title,"%s Ave. Track Length by Subrun - partition %.1u",trigname.c_str(),p);
    gTrackLenday_zoom->SetTitle(title);
    gTrackLenday_zoom->SetMarkerColor(kBlue);
    gTrackLenday_zoom->GetXaxis()->SetTimeDisplay(1);
    gTrackLenday_zoom->GetXaxis()->SetLabelSize(0.03);
    gTrackLenday_zoom->GetXaxis()->SetTimeFormat(taxis_labels.c_str());
    gTrackLenday_zoom->GetXaxis()->SetLimits(time_ago,XNow);
    gTrackLenday_zoom->GetXaxis()->SetTitle("(central time)");
    gTrackLenday_zoom->GetYaxis()->SetRangeUser((gTrackLenday_zoom->GetMean(2)-3.*gTrackLenday_zoom->GetRMS(2)),(gTrackLenday_zoom->GetMean(2)+3.*gTrackLenday_zoom->GetRMS(2)));
    gTrackLenday_zoom->Draw("A*");
    UpdateText->Draw();
    LastPoint->Draw();
    sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1uTrackLenAll3DperSRzoom%s.png",det.c_str(),trig.c_str(),p,period.c_str());
    ctTrackLenday_zoom->Print(filename);
    delete ctTrackLenday_zoom;
    ctTrackLenday_zoom = 0;
    
    // NEW PLOTS
    // Comment by Hayes: Am now able to draw the three graphs on the same pad and canvas (via a TMultiGraph), but cannot convert 
    // the time on the xaxis to the readable Hour Minute format.  Tried using the same commands on the multigraph as are used 
    // for each graph and it seg faults every time...
    TCanvas *cTrackFracday = new TCanvas("cTrackFracday","Track Fractions day",1200,800);
    TPad *pad1 = new TPad("pad1","",0,0,1,1);
    TPad *pad2 = new TPad("pad2","",0,0,1,1);
    pad2->SetFillStyle(4000);  // pad2 will be transparent
    pad2->SetFillColor(0);
    pad2->SetFrameFillStyle(4000);
    pad1->Draw();
    pad1->cd();
    pad1->SetGridx();
    TMultiGraph *gTF = new TMultiGraph;
    TGraph *gTFAll3Dday  = new TGraph(TrackFracDayCount[p],TrackFracAll3DtimeDay[p],TrackFracAll3DDay[p]);
    TGraph *gTFAll2Dday  = new TGraph(TrackFracDayCount[p],TrackFracAll2DtimeDay[p],TrackFracAll2DDay[p]);
    TGraph *gTFCont3Dday = new TGraph(TrackFracDayCount[p],TrackFracCont3DtimeDay[p],TrackFracCont3DDay[p]);
    sprintf(title,"%s Track Fractions - partition %.1u",trigname.c_str(),p);
    
    pad1->SetLogy();
    
    
    gTFAll2Dday->SetTitle(title);
    gTFAll2Dday->SetMarkerColor(kBlue-7);
    gTFAll2Dday->SetFillColor(kBlue-7);
    if(det == "NearDet")
      gTFAll2Dday->SetMarkerStyle(21);
    else
      gTFAll2Dday->SetMarkerStyle(7);
    
    gTFCont3Dday->SetTitle(title);
    gTFCont3Dday->SetMarkerColor(kGray+3);
    gTFCont3Dday->SetFillColor(kGray+3);
    if(det == "NearDet")
      gTFCont3Dday->SetMarkerStyle(21);
    else
      gTFCont3Dday->SetMarkerStyle(7);
    
    //gTF->Add(gTFAll3Dday);
    gTF->Add(gTFAll2Dday);
    gTF->Add(gTFCont3Dday);
    gTF->Draw("ap");

    gTF->SetTitle(title);
    //cTrackFracday->Modified();
    //cTrackFracday->Update();
    pad1->Modified();
    pad1->Update();
    
    UpdateText->Draw();    
    
    //gTF->SetTitle(title);
    gTF->GetXaxis()->SetRangeUser(0.0,1.1);
    gTF->GetXaxis()->SetTimeDisplay(1);
    gTF->GetXaxis()->SetLabelSize(0.03);
    gTF->GetXaxis()->SetTimeFormat(taxis_labels.c_str());
    gTF->GetXaxis()->SetLimits(time_ago,XNow);
    gTF->GetXaxis()->SetTitle("(central time)");
    gTF->GetYaxis()->SetRangeUser(0.001,1.0);
    
    //cTrackFracday->cd();
    
    pad2->Draw();
    pad2->cd();
    pad2->SetGridx();
    gTFAll3Dday->SetTitle("");
    gTFAll3Dday->SetMarkerColor(kMagenta-7);
    gTFAll3Dday->SetFillColor(kMagenta-7);
    if(det == "NearDet")
      gTFAll3Dday->SetMarkerStyle(21);
    else
      gTFAll3Dday->SetMarkerStyle(7);
    gTFAll3Dday->GetXaxis()->SetRangeUser(0.0,1.1);
    gTFAll3Dday->GetXaxis()->SetTimeDisplay(1);
    gTFAll3Dday->GetXaxis()->SetLabelSize(0.03);
    gTFAll3Dday->GetXaxis()->SetTimeFormat(taxis_labels.c_str());
    gTFAll3Dday->GetXaxis()->SetLimits(time_ago,XNow);
    gTFAll3Dday->GetXaxis()->SetTitle("(central time)");
    gTFAll3Dday->GetYaxis()->SetLabelColor(kMagenta-7);
    gTFAll3Dday->GetYaxis()->SetLabelSize(0.026);
    //gTFAll3Dday->GetYaxis()->SetRangeUser(0.97,1.00);
    gTFAll3Dday->Draw("apY+");
    pad2->Modified();
    pad2->Update();
    
    // Define legend
    TLegend *leg = new TLegend(0.60,0.66,0.85,0.78);
    leg->SetFillColor(0);
    leg->AddEntry(gTFAll2Dday,"#bf{All 2D Tracks}","f");
    leg->AddEntry(gTFCont3Dday,"#bf{Contained 3D Tracks}","f");
    leg->AddEntry(gTFAll3Dday,"#bf{All 3D Tracks}","f");
    leg->SetFillStyle(4000);
    leg->SetLineColor(0);
    
    leg->Draw();

    sprintf(filename,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1uTrackFractionsperSR%s.png",det.c_str(),trig.c_str(),p,period.c_str());
    cTrackFracday->Print(filename);
    delete cTrackFracday;
    cTrackFracday = 0;
    
    
    // Now start making efficiency plots using hitsVwcell and cellsVwcell
    eff_plots(p, det, period, trig, trigname);
    
    UpdateText->Clear();
    
  } // end loop over partitions
  
  return 0;
  
} // End script.


void eff_plots(unsigned int par_eff, string det, string period, string trig, string trigname){
  
  //TFile *tfile = new TFile("/bluearc/nova/data/nearline-Ana/FarDet/S14-02-24/000143/00014332/fardet_r00014332_s01_reco_hist_t02.root");
  TFile *tfile = new TFile("/bluearc/nova/data/nearline-Ana/FarDet/S14-08-01/000168/00016898/fardet_r00016898_s13_reco_hist_t02.root");
  tfile->cd("hitefficiencyana");
  TH2F* hitsVwcell = (TH2F*)tfile->FindObjectAny("hitsVwcell");
  TH2F* cellsVwcell = (TH2F*)tfile->FindObjectAny("cellsVwcell");

  // Output file list of low efficiency modules.
  ofstream modules;

  int ncellbins_local = hitsVwcell->GetNbinsX();
  double cellmin_local = hitsVwcell->GetXaxis()->GetXmin();
  double cellmax_local = hitsVwcell->GetXaxis()->GetXmax();
  int nwbins_local = hitsVwcell->GetNbinsY();
  double wmin_local = hitsVwcell->GetYaxis()->GetXmin();
  double wmax_local = hitsVwcell->GetYaxis()->GetXmax();

  double wavglength = 200.0; // low w length to average cell/module efficiencies over
  int wbinlow = 1;
  int wbinhigh = hitsVwcell->GetYaxis()->FindFixBin(hitsVwcell->GetYaxis()->GetBinCenter(wbinlow)+wavglength);
  //int neg_wbinhigh = -800;
  //int neg_wbinlow = hitsVwcell->GetYaxis()->FindFixBin(hitsVwcell->GetYaxis()->GetBinCenter(neg_wbinhigh)+250);

  char f_name[56];
  char hmod_name[56];

  // Defing binning for histograms that will be made here
  int neffbins = 101;
  double effmin = -0.05;
  double effmax = 1.05;
  int nplbins = 897;
  double plmin = -0.5;
  double plmax = nplbins-0.5;
  int nmodbins = 12;
  double modmin = 0;
  double modmax = 12;

//   TH1F* mod_Eff[12];
//   for(int imod=0;imod<12;imod++){
//     sprintf(hmod_name,"mod%ieffperplane",imod);
//     mod_Eff[imod] = new TH1F(hmod_name,hmod_name,nplbins,0,nplbins);
//   }

  TH2F* moduleEfficiencyAverage = new TH2F("moduleEfficiencyAverage",";Plane;Module",
					   nplbins,plmin,plmax,nmodbins,modmin,modmax);
  TH2F* moduleEfficiencyRMS = new TH2F("moduleEfficiencyRMS",";Plane;Module",
				       nplbins,plmin,plmax,nmodbins,modmin,modmax);
  TH2F* moduleEfficiencyAverageXView = new TH2F("moduleEfficiencyAverageXView",";Plane;Module",
						nplbins,plmin,plmax,nmodbins,modmin,modmax);
  TH2F* moduleEfficiencyRMSXView = new TH2F("moduleEfficiencyRMSXView",";Plane;Module",
					    nplbins,plmin,plmax,nmodbins,modmin,modmax);

  TH2F* moduleEfficiencyAverageYView = new TH2F("moduleEfficiencyAverageYView",";Plane;Module",
						nplbins,plmin,plmax,nmodbins,modmin,modmax);
  TH2F* moduleEfficiencyRMSYView = new TH2F("moduleEfficiencyRMSYView",";Plane;Module",
					    nplbins,plmin,plmax,nmodbins,modmin,modmax);

  // Histogram of efficiency for a module
  TH1F* modeff = new TH1F("modeff",";Efficiency;NCells",neffbins,effmin,effmax);

  // Histogram of efficiency vs w for one cell
  TH1F* eff = new TH1F("eff","",nwbins_local,wmin_local,wmax_local);

  int currmodule = 0;
  int currplane = 0;
  int currview = 0; // View convention y = 0, x = 1

  // loop over every cell in each partition
  if(par_eff==1 && period=="Month"){
    modules.open("/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/low_eff_modules.txt");
    modules<<"Partition "<<par_eff<<endl;
    modules<<endl;
  }
  for(int icell = 1; icell <= ncellbins_local; ++icell){

    // convert this to a plane and cell number
    int cell = (icell-1)%384;
    int plane = (icell-1-cell)/384;
    int modcell = cell%32;
    int module = (cell-modcell)/32;
    int view = plane%2; // view = 0 y, view = 1 x

    if(plane != currplane || module != currmodule){
      // fill in the histogram of average efficiency and rms of average for the last module
      double effavg = modeff->GetMean();
      double effrms = modeff->GetRMS();
      //double effavg_far = modeff_far->GetMean();
      int planebin = moduleEfficiencyAverage->GetXaxis()->FindFixBin(currplane);
      int modulebin = moduleEfficiencyAverage->GetYaxis()->FindFixBin((double)currmodule + 0.5);

      if(par_eff==1 && period=="Month"){
	if(effavg<0.3){
	  modules<<"module:  "<<modulebin<<"  plane:  "<<planebin<<"  efficiency:  "<<effavg<<endl;
	}
      }
      //mod_Eff[modulebin-1]->Fill(planebin,effavg);
      moduleEfficiencyAverage->SetBinContent(planebin,modulebin,effavg);
      moduleEfficiencyRMS->SetBinContent(planebin,modulebin,effrms);
      if(currview){
	// x view
	moduleEfficiencyAverageXView->SetBinContent(planebin,modulebin,effavg);
	moduleEfficiencyRMSXView->SetBinContent(planebin,modulebin,effrms);
      }
      else{
	// y view
	moduleEfficiencyAverageYView->SetBinContent(planebin,modulebin,effavg);
	moduleEfficiencyRMSYView->SetBinContent(planebin,modulebin,effrms);
      }
      // reset things for the new module
      modeff->Reset("ICES");
      currmodule = module;
      currplane = plane;
      currview = view;
    }
    

    // need to switch the module in cell for x view
 
    if(view){ modcell = 31-modcell; }

    // loop over every w bin
    double totcell = 0.0;
    double tothit = 0.0;

    // Calculate the cell efficiency vs w
    for(int iw = 1; iw <= nwbins_local; ++iw){
      // get the efficiency at this cell/w, eff = k/n
      double k = hHitsVwcellDay[par_eff]->GetBinContent(icell,iw);
      double n = hCellsVwcellDay[par_eff]->GetBinContent(icell,iw);

      totcell+=n;
      tothit+=k;
      double ef = k/n;
      if(n == 0){ ef = 0; }

      // calculate the variance of this
      double var = (k+1)*(k+2)/((n+2)*(n+3)) - (k+1)*(k+1)/((n+2)*(n+2));
      
      double w = hHitsVwcellDay[par_eff]->GetYaxis()->GetBinCenter(iw);
      
      eff->SetBinContent(iw,ef);
      eff->SetBinError(iw,sqrt(var));
           
    }


    // Average over the last 200 cm of the cell
    double efftot = eff->Integral(wbinlow,wbinhigh);
    double navgbins = wbinhigh-wbinlow+1.0;

    modeff->Fill(efftot/navgbins);
           
  }
  if(par_eff==1 && period=="Month") modules.close();
  
  char title_buff[128];

  gStyle->SetTimeOffset(0);
  string t_title;
  if(period == "Day") t_title = "(past 24 hrs.)";
  if(period == "Week") t_title = "(past week)";
  if(period == "Month") t_title = "(past month)";

  TDatime *Ttemp = new TDatime;  // finish time
  int Xfin = Ttemp->Convert() - GMToffset;
  TDatime *Tfinish = new TDatime(Xfin);  // finish time
  TPaveText *UpdateText = new TPaveText(0.1, 0.0, 0.5, 0.05, "NDC");
  UpdateText->SetLineColor(0);
  UpdateText->SetFillColor(0);
  UpdateText->SetBorderSize(1);
  UpdateText->SetMargin(0.0);
  UpdateText->SetTextAlign(11);
  char buff1[256];
  sprintf(buff1, "Last updated on:    %s (central time)", Tfinish->AsString());
  UpdateText->AddText(buff1);
  char buff2[256];
  sprintf(buff2, "Last run / subrun:   %d / %d", Last_Run[par_eff], Last_SR[par_eff]);
  UpdateText->AddText(buff2);

  TPaveText *StatusText = new TPaveText(0.1, 0.1, 0.9, 0.9, "NDC");
  StatusText->SetLineColor(1);
  StatusText->SetFillColor(0);
  StatusText->SetBorderSize(1);
  // StatusText->SetMargin(0.0);
  StatusText->AddText(buff1);
  StatusText->AddText(buff2);

  TCanvas* c1 = new TCanvas("c1","c1",1200,800);
  gStyle->SetOptStat(0);
  sprintf(title_buff,"Average Module Efficiency %s - partition %.1u",t_title.c_str(),par_eff);
  moduleEfficiencyAverage->SetTitle(title_buff);
  moduleEfficiencyAverage->Draw("colz");
  //moduleEfficiencyAverage->GetXaxis()->SetRangeUser(0.,300.);
  gPad->Update();
//   TPaveStats *stats1 = (TPaveStats*)moduleEfficiencyAverage->FindObject("stats");
//   stats1->SetX1NDC(0.65);
//   stats1->SetY1NDC(0.65);
//   stats1->SetX2NDC(0.90);
//   stats1->SetY2NDC(0.90);  
  c1->Modified();
  c1->Update();
  UpdateText->Draw();
  sprintf(f_name,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1uAvgModuleEff%s.png",det.c_str(),trig.c_str(),par_eff,period.c_str());
  c1->Print(f_name);
  delete c1;
  c1 = 0;

  //////////////////////////////////////////////////////////////////////
  // Creating list of low efficiency modules and corresponding planes.//
  //////////////////////////////////////////////////////////////////////
  /* TFile *test_out = new TFile("testout.root","RECREATE");
  modules.open("low_eff_modules.txt");
  //int Nbinsy = moduleEfficiencyAverage->GetNbinsY();
  //int Nbinsx = moduleEfficiencyAverage->GetNbinsX();
  if(det=="FarDet" && period=="Day"){
    for(int i=1;i<257;i++){
      for(int j=1;j<13;j++){
	//Int_t bin_num = moduleEfficiencyAverage->FindFixBin(double(i)+0.5,j);
	//Int_t bin_num = moduleEfficiencyAverage->GetBin(i+1,j+1);
	Double_t eff_mod = moduleEfficiencyAverage->GetBinContent(i,j);
	//cout<<bin_num<<"      "<<eff_mod<<endl;
	//if(eff_mod<0.6){
	//modules<<" plane:  "<<i<<"   module:   "<<j<<"   efficiency:  "<<eff_mod<<endl;
	modules << moduleEfficiencyAverage->GetBinContent(i,j)<<endl;
	  //}
      }
    }
  }
  modules.close();
  test_out->cd();
  for(int ii=0;ii<12;ii++){
  mod_Eff[ii]->Write();
  }
  //moduleEfficiencyAverage->Write();
  //test_out->Close();
  */

  TCanvas* c2 = new TCanvas("c2","c2",1200,800);
  gStyle->SetOptStat(0);
  sprintf(title_buff,"%s Average Module Efficiency RMS %s - partition %.1u",trigname.c_str(),t_title.c_str(),par_eff);
  moduleEfficiencyRMS->SetTitle(title_buff);
  moduleEfficiencyRMS->Draw("colz");
  //moduleEfficiencyRMS->GetXaxis()->SetRangeUser(0.,300.);
  gPad->Update();
//   TPaveStats *stats2 = (TPaveStats*)moduleEfficiencyRMS->FindObject("stats");
//   stats2->SetX1NDC(0.65);
//   stats2->SetY1NDC(0.65);
//   stats2->SetX2NDC(0.90);
//   stats2->SetY2NDC(0.90); 
  c2->Modified();
  c2->Update();
  UpdateText->Draw();
  sprintf(f_name,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1uModuleEffRMS%s.png",det.c_str(),trig.c_str(),par_eff,period.c_str());
  c2->Print(f_name);
  delete c2;
  c2 = 0;

  TCanvas* c3 = new TCanvas("c3","c3",1200,800);
  gStyle->SetOptStat(0);
  sprintf(title_buff,"%s Average Module Efficiency In X-view %s - partition %.1u",trigname.c_str(),t_title.c_str(),par_eff);
  moduleEfficiencyAverageXView->SetTitle(title_buff);
  moduleEfficiencyAverageXView->Draw("colz");
  moduleEfficiencyAverageXView->RebinX(2);
  //moduleEfficiencyAverageXView->GetXaxis()->SetRangeUser(0.,300.);
  gPad->Update();
//   TPaveStats *stats3 = (TPaveStats*)moduleEfficiencyAverageXView->FindObject("stats");
//   stats3->SetX1NDC(0.65);
//   stats3->SetY1NDC(0.65);
//   stats3->SetX2NDC(0.90);
//   stats3->SetY2NDC(0.90); 
  c3->Modified();
  c3->Update();
  UpdateText->Draw();
  sprintf(f_name,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1uAvgModuleXEff%s.png",det.c_str(),trig.c_str(),par_eff,period.c_str());
  c3->Print(f_name);
  delete c3;
  c3 = 0;

  TCanvas* c4 = new TCanvas("c4","c4",1200,800);
  gStyle->SetOptStat(0);
  sprintf(title_buff,"%s Average Module Efficiency In X-view RMS %s - partition %.1u",trigname.c_str(),t_title.c_str(),par_eff);
  moduleEfficiencyRMSXView->SetTitle(title_buff);
  moduleEfficiencyRMSXView->Draw("colz");
  moduleEfficiencyRMSXView->RebinX(2);
  //moduleEfficiencyRMSXView->GetXaxis()->SetRangeUser(0.,300.);
  gPad->Update();
//   TPaveStats *stats4 = (TPaveStats*)moduleEfficiencyRMSXView->FindObject("stats");
//   stats4->SetX1NDC(0.65);
//   stats4->SetY1NDC(0.65);
//   stats4->SetX2NDC(0.90);
//   stats4->SetY2NDC(0.90); 
  c4->Modified();
  c4->Update();
  UpdateText->Draw();
  sprintf(f_name,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1uModuleXEffRMS%s.png",det.c_str(),trig.c_str(),par_eff,period.c_str());
  c4->Print(f_name);
  delete c4;
  c4 = 0;

  TCanvas* c5 = new TCanvas("c5","c5",1200,800);
  gStyle->SetOptStat(0);
  sprintf(title_buff,"%s Average Module Efficiency In Y-view %s - partition %.1u",trigname.c_str(),t_title.c_str(),par_eff);
  moduleEfficiencyAverageYView->SetTitle(title_buff);
  moduleEfficiencyAverageYView->Draw("colz");
  moduleEfficiencyAverageYView->RebinX(2);
  //moduleEfficiencyAverageYView->GetXaxis()->SetRangeUser(0.,300.);
  gPad->Update();
//   TPaveStats *stats5 = (TPaveStats*)moduleEfficiencyAverageYView->FindObject("stats");
//   stats5->SetX1NDC(0.65);
//   stats5->SetY1NDC(0.65);
//   stats5->SetX2NDC(0.90);
//   stats5->SetY2NDC(0.90); 
  c5->Modified();
  c5->Update();
  UpdateText->Draw();
  sprintf(f_name,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1uAvgModuleYEff%s.png",det.c_str(),trig.c_str(),par_eff,period.c_str());
  c5->Print(f_name);
  delete c5;
  c5 = 0;

  TCanvas* c6 = new TCanvas("c6","c6",1200,800);
  gStyle->SetOptStat(0);
  sprintf(title_buff,"%s Average Module Efficiency In Y-view RMS %s - partition %.1u",trigname.c_str(),t_title.c_str(),par_eff);
  moduleEfficiencyRMSYView->SetTitle(title_buff);
  moduleEfficiencyRMSYView->Draw("colz");
  moduleEfficiencyRMSYView->RebinX(2);
  //moduleEfficiencyRMSYView->GetXaxis()->SetRangeUser(0.,300.);
  gPad->Update();
//   TPaveStats *stats6 = (TPaveStats*)moduleEfficiencyRMSYView->FindObject("stats");
//   stats6->SetX1NDC(0.65);
//   stats6->SetY1NDC(0.65);
//   stats6->SetX2NDC(0.90);
//   stats6->SetY2NDC(0.90); 
  c6->Modified();
  c6->Update();
  UpdateText->Draw();
  sprintf(f_name,"/nusoft/app/web/htdoc/nova/datacheck/nearline/plots/%s-%s-P%.1uModuleYEffRMS%s.png",det.c_str(),trig.c_str(),par_eff,period.c_str());
  c6->Print(f_name);
  delete c6;
  c6 = 0;
    	    
}
