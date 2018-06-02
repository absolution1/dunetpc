#ifndef CALGAUSHFDUNE10KT_H
#define CALGAUSHFDUNE10KT_H

////////////////////////////////////////////////////////////////////////
//
// CalGausHFDUNE10kt class
//
// jti3@fnal.gov
//
// 06-21-13 Branched from CalWire_module and added GausHitFinder_module code
//
// brebel@fnal.gov
//
// 11-3-09 Pulled all FFT code out and put into Utilitiess/LArFFT
//
////////////////////////////////////////////////////////////////////////

#include <string>
#include <vector>
#include <algorithm> // std::accumulate
#include <utility> // std::move

#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h" 
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "cetlib_except/exception.h"

#include "dune/Utilities/SignalShapingServiceDUNE.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larcore/Geometry/Geometry.h"
#include "larevt/Filters/ChannelFilter.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/ArtDataHelper/HitCreator.h"
#include "lardata/Utilities/LArFFT.h"

// ROOT Includes 
#include "TH1D.h"
#include "TDecompSVD.h"
#include "TMath.h"
#include "TComplex.h"
#include "TFile.h"
#include "TH2D.h"
#include "TF1.h"
#include "TTree.h"

/* unused function
namespace{
  // Fill histogram from vector (set underflow/overflow bins to zero).

  void vector_to_hist(const std::vector<double>& v, TH1D* h)
  {
    if (!h) throw cet::exception("vector_to_hist") << "no histogram specified\n";
    int nvec = v.size();
    int nbins = h->GetNbinsX();
    int nfill = std::min(nvec, nbins);
    h->SetBinContent(0, 0.);
    int zerobin = h->GetXaxis()->FindBin(0.);
    for(int i=1; i<=nfill; ++i) {
      if(i >= zerobin)
	h->SetBinContent(i, v[i - zerobin]);
      else
	h->SetBinContent(i, v[i - zerobin + nvec]);
    }
    for(int i=nfill+1; i<=nbins+1; ++i)
      h->SetBinContent(i, 0.);
  }

}
*/

///creation of calibrated signals on wires
namespace calgaushf {

  class CalGausHFDUNE10kt : public art::EDProducer {

  public:
    
    // create calibrated signals on wires. this class runs 
    // an fft to remove the electronics shaping.     
    explicit CalGausHFDUNE10kt(fhicl::ParameterSet const& pset); 
    virtual ~CalGausHFDUNE10kt();
    
    void produce(art::Event& evt); 
    void beginJob(); 
    void endJob();                 
    void reconfigure(fhicl::ParameterSet const& p);
 
  private:
    int          fDeconvKernSize;    //< Length of truncated deconvolution kernel
    //int          fDataSize;          ///< size of raw data on one wire // unused
    int          fPostsample;        ///< number of postsample bins
    std::string  fDigitModuleLabel;  ///< module that made digits
                                                       ///< constants
    std::string  fSpillName;  ///< nominal spill is an empty string
                              ///< it is set by the DigitModuleLabel
                              ///< ex.:  "daq:preSpill" for prespill data
    std::string     fCalDataModuleLabel;
    double          fMinSigInd;     ///<Induction signal height threshold 
    double          fMinSigCol;     ///<Collection signal height threshold 
    double          fIndWidth;      ///<Initial width for induction fit
    double          fColWidth;      ///<Initial width for collection fit
    double          fIndMinWidth;   ///<Minimum induction hit width
    double          fColMinWidth;   ///<Minimum collection hit width
    int             fMaxMultiHit;   ///<maximum hits for multi fit 
    int             fAreaMethod;    ///<Type of area calculation  
    std::vector<double> fAreaNorms; ///<factors for converting area to same units as peak height 
    double	    fChi2NDF;       ///maximum Chisquared / NDF allowed for a hit to be saved
    
    //double	WireNumber[100000];
    //double 	TotalSignal[100000];
    double 	StartIime;
    double 	StartTimeError;
    double 	EndTime;
    double 	EndTimeError;
    double 	RMS;
    int	NumOfHits;
    double	MeanPosition;
    double 	MeanPosError;
    double	Amp;
    double	AmpError;
    double	Charge;
    double	ChargeError;
    double	FitGoodnes;
    int   	FitNDF;

  protected: 
    
  }; // class CalGausHFDUNE10kt

  DEFINE_ART_MODULE(CalGausHFDUNE10kt)
  
  //-------------------------------------------------
  CalGausHFDUNE10kt::CalGausHFDUNE10kt(fhicl::ParameterSet const& pset)
  {
    fSpillName="";
    this->reconfigure(pset);

    // let HitCollectionCreator declare that we are going to produce
    // hits and associations with wires and raw digits
    recob::HitCollectionCreator::declare_products
      (*this, fSpillName, false /* doWireAssns */, true /* doRawDigitAssns */);
    
  }
  
  //-------------------------------------------------
  CalGausHFDUNE10kt::~CalGausHFDUNE10kt()
  {
  }

  //////////////////////////////////////////////////////
  void CalGausHFDUNE10kt::reconfigure(fhicl::ParameterSet const& p)
  {
    fDigitModuleLabel = p.get< std::string >("DigitModuleLabel", "daq");
    fPostsample       = p.get< int >        ("PostsampleBins");
    fDeconvKernSize   = p.get< int >        ("DeconvKernSize");
    
    fSpillName="";
    
    size_t pos = fDigitModuleLabel.find(":");
    if( pos!=std::string::npos ) {
      fSpillName = fDigitModuleLabel.substr( pos+1 );
      fDigitModuleLabel = fDigitModuleLabel.substr( 0, pos );
    }
  
    // Implementation of optional member function here.
    //fCalDataModuleLabel = p.get< std::string  >("CalDataModuleLabel");
    fMinSigInd          = p.get< double       >("MinSigInd");
    fMinSigCol          = p.get< double       >("MinSigCol"); 
    fIndWidth           = p.get< double       >("IndWidth");  
    fColWidth           = p.get< double       >("ColWidth");
    fIndMinWidth        = p.get< double       >("IndMinWidth");
    fColMinWidth        = p.get< double       >("ColMinWidth"); 	  	
    fMaxMultiHit        = p.get< int          >("MaxMultiHit");
    fAreaMethod         = p.get< int          >("AreaMethod");
    fAreaNorms          = p.get< std::vector< double > >("AreaNorms");
    fChi2NDF	      = p.get< double       >("Chi2NDF");
  
  
  }

  //-------------------------------------------------
  void CalGausHFDUNE10kt::beginJob()
  {  
  }

  //////////////////////////////////////////////////////
  void CalGausHFDUNE10kt::endJob()
  {  
  }
  
  //////////////////////////////////////////////////////
  void CalGausHFDUNE10kt::produce(art::Event& evt)
  {      
    // get the geometry
    art::ServiceHandle<geo::Geometry> geom;

    // Get signal shaping service.
    art::ServiceHandle<util::SignalShapingServiceDUNE> sss;


    //Gaussian hit finder initializations etc.
    TH1::AddDirectory(kFALSE);
    
    
    // ---------------------------
    // --- Seed Fit Functions  --- (This function used to "seed" the mean peak position)
    // ---------------------------
    TF1 *hit	= new TF1("hit","gaus",0,3200); // FIXME use DetectorProperties
    
    // ###############################################
    // ### Making a ptr vector to put on the event ###
    // ###############################################
    // this contains the hit collection
    // and its associations to raw digits (not to wires)
    recob::HitCollectionCreator hcol(*this,
      evt, fSpillName, false /* doWireAssns */, true /* doRawDigitAssns */
      );
    
    
    // ##########################################
    // ### Reading in the Wire List object(s) ###
    // ##########################################
    //art::Handle< std::vector<recob::Wire> > wireVecHandle;
    //evt.getByLabel(fCalDataModuleLabel,wireVecHandle);
    //art::ServiceHandle<geo::Geometry> geom;
    
    // #########################################################
    // ### List of useful variables used throughout the code ###
    // #########################################################  
    double threshold              = 0.;               // minimum signal size for id'ing a hit
    double fitWidth               = 0.;               //hit fit width initial value
    double minWidth               = 0.;               //minimum hit width
    //channel              = 0;                // channel number
    geo::SigType_t sigType;                           // Signal Type (Collection or Induction)
    geo::View_t view;                                // view
    std::vector<int> startTimes;             	    // stores time of 1st local minimum
    std::vector<int> maxTimes;    	   	    // stores time of local maximum    
    std::vector<int> endTimes;    	     	    // stores time of 2nd local minimum
    int time             = 0;                         // current time bin
    int minTimeHolder    = 0;                         // current start time
    
    std::string eqn        = "gaus(0)";      // string for equation for gaus fit
    //std::string seed        = "gaus(0)";      // string for equation seed gaus fit
    
    
    bool maxFound        = false;            // Flag for whether a peak > threshold has been found
    std::stringstream numConv;
	  

    // Read in the digit List object(s). 
    art::Handle< std::vector<raw::RawDigit> > digitVecHandle;

 
    if(fSpillName.size()>0) evt.getByLabel(fDigitModuleLabel, fSpillName, digitVecHandle);
    else evt.getByLabel(fDigitModuleLabel, digitVecHandle);

    if (!digitVecHandle->size())  return;
    mf::LogInfo("CalGausHFDUNE10kt") << "CalGausHFDUNE10kt:: digitVecHandle size is " << digitVecHandle->size();

    // Use the handle to get a particular (0th) element of collection.
    art::Ptr<raw::RawDigit> digitVec0(digitVecHandle, 0);
            
    if (digitVec0->Compression() != raw::kZeroSuppression) {
      throw art::Exception(art::errors::UnimplementedFeature)
	<< "CalGausHFDUNE only supports zero-suppressed raw digit input!";
    } // if


    raw::ChannelID_t channel = raw::InvalidChannelID; // channel number
    unsigned int bin(0);     // time bin loop variable
    
    filter::ChannelFilter *chanFilt = new filter::ChannelFilter();  

    std::vector<float> holder;                // holds signal data
    std::vector<float> rawadc_conv;           // holds signal data before time domain deconvolution
    std::vector<short> rawadc;  // vector holding uncompressed adc values

    // loop over all wires    
    for(size_t rdIter = 0; rdIter < digitVecHandle->size(); ++rdIter){ // ++ move
      holder.clear();
      
      // get the reference to the current raw::RawDigit
      art::Ptr<raw::RawDigit> digitVec(digitVecHandle, rdIter);
      channel = digitVec->Channel();

      // skip bad channels
      if(!chanFilt->BadChannel(channel)) {


 	const int numberofblocks = digitVec->ADC(1);

 	int zerosuppressedindex = numberofblocks*2 + 2;

 	//loop over all nonzero blocks
	for(int i=0; i<numberofblocks; i++){

	  const int lengthofblock = digitVec->ADC(2+numberofblocks+i);
	  rawadc.resize(lengthofblock);
	  holder.resize(lengthofblock);
	  rawadc_conv.resize(lengthofblock);

	  for (int j=0;j<lengthofblock;j++)
	    {
	      rawadc[j] = digitVec->ADC(zerosuppressedindex);
	      zerosuppressedindex++;
	    }
    	  
	  
	  
	  // loop over all adc values and subtract the pedestal
	  for(bin = 0; bin < rawadc.size(); ++bin) 
	    rawadc_conv[bin]=(rawadc[bin]-digitVec->GetPedestal());


	  sigType = geom->SignalType(channel);
	  view = geom->View(channel);

	  if(sigType == geo::kCollection){

	    //Leave collection plane signals as is
	    for(bin = 0; bin < rawadc_conv.size(); ++bin)
	      holder[bin] = rawadc_conv[bin];
	  }
	  else if(sigType == geo::kInduction){

	    //Integrate over signal in induction planes
	    for(bin = 0;  bin < rawadc_conv.size(); ++bin) 
	      holder[bin] = (bin>0) ? rawadc_conv[bin] + holder[bin-1] : rawadc_conv[bin];
	  }


	  
	  //Beginning of Gaussian hit finder loop over signal blocks


	  //##############################
	  //### Looping over the wires ###
	  //############################## 
	  
	  //for(size_t wireIter = 0; wireIter < wirecol->size(); wireIter++){
	  
	  //art::Ptr<recob::Wire> wire(wirecol, wireIter);
	  
	  //std::vector<recob::Wire>::iterator wire;
	  //for(wire=wirecol->begin(); wire != wirecol->end(); wire++){
	  // --- Setting Channel Number and Wire Number as well as signal type ---
	  //channel = wire->RawDigit()->Channel();
	  
	  
	  // -----------------------------------------------------------
	  // -- Clearing variables at the start of looping over wires --
	  // -----------------------------------------------------------
	  startTimes.clear();
	  maxTimes.clear();
	  endTimes.clear();
	  //std::vector<float> signal(wire->Signal());
	  std::vector<float> signal(holder);
	  std::vector<float>::iterator timeIter;  	    // iterator for time bins
	  time          = 0;
	  minTimeHolder = 0;
	  maxFound      = false;
	  // -- Setting the appropriate signal widths and thresholds --
	  // --    for either the collection or induction plane      --
	  if(sigType == geo::kInduction){
	    threshold     = fMinSigInd;
	    fitWidth      = fIndWidth;
	    minWidth      = fIndMinWidth;
	  }//<-- End if Induction Plane
	  else if(sigType == geo::kCollection){
	    threshold = fMinSigCol;
	    fitWidth  = fColWidth;
	    minWidth  = fColMinWidth;
	  }//<-- End if Collection Plane
	  
	  
	  // ##################################
	  // ### Looping over Signal Vector ###
	  // ##################################
	  for(timeIter = signal.begin();timeIter+2<signal.end();timeIter++){
	    // ##########################################################
	    // ###                LOOK FOR A MINIMUM                  ###
	    // ### Testing if the point timeIter+1 is at a minimum by ###
	    // ###  checking if timeIter and timeIter+1 and if it is  ###
	    // ###         then we add this to the endTimes           ###
	    // ##########################################################
	    if(*timeIter > *(timeIter+1) && *(timeIter+1) < *(timeIter+2)){
	      //--- Note: We only keep the a minimum if we've already ---
	      //---          found a point above threshold            ---
	      if(maxFound){
		endTimes.push_back(time+1);
		maxFound = false;
		//keep these in case new hit starts right away
		minTimeHolder = time+2; 
	      }
	      else minTimeHolder = time+1; 
	      
	    }//<---End Checking if this is a minimum
	    
	    
	    // ########################################################	
	    // ### Testing if the point timeIter+1 is a maximum and ###
	    // ###  if it and is above threshold then we add it to  ###
	    // ###                  the startTime                   ###
	    // ########################################################
	    //if not a minimum, test if we are at a local maximum 
	    //if so, and the max value is above threshold, add it and proceed.
	    else if(*timeIter < *(timeIter+1) && *(timeIter+1) > *(timeIter+2) && *(timeIter+1) > threshold){ 
	      maxFound = true;
	      maxTimes.push_back(time+1);
	      startTimes.push_back(minTimeHolder);          
	    }
	    
	    time++;
	    
	  }//end loop over signal vec 
	  
	  // ###########################################################    
	  // ### If there was no minimum found before the end, but a ### 
	  // ###  maximum was found then add an end point to the end ###
	  // ###########################################################
	  while( maxTimes.size() > endTimes.size() ){ 
	    endTimes.push_back(signal.size()-1);
	    
	  }//<---End maxTimes.size > endTimes.size
	  
	  // ####################################################
	  // ### If no startTime hit was found skip this wire ###
	  // ####################################################
	  if( startTimes.size() == 0 ){
	    continue;
	  }  
	  
	  /////////////////////////////////////////////////////////////////////////////
	  /////////////////////////////////////////////////////////////////////////////
	  
	  // ##########################################################
	  // ### PERFORM THE FITTING, ADDING HITS TO THE HIT VECTOR ###
	  // ##########################################################
	  
	  //All code below does the fitting, adding of hits
	  //to the hit vector and when all wires are complete 
	  //saving them 
	  
	  double totSig(0); 						//stores the total hit signal
	  double startT(0); 						//stores the start time
	  double endT(0);  						//stores the end time
	  int numHits(0);  						//number of consecutive hits being fitted
	  int size(0);     						//size of data vector for fit
	  int hitIndex(0);						//index of current hit in sequence
	  double amplitude(0);             			        //fit parameters
	  double minPeakHeight(0);  					//lowest peak height in multi-hit fit
	  
	  
	  StartIime = 0; 							// stores the start time of the hit
	  StartTimeError = 0;						// stores the error assoc. with the start time of the hit
	  EndTime = 0;							// stores the end time of the hit
	  EndTimeError = 0;						// stores the error assoc. with the end time of the hit
	  RMS = 0;         						// stores the RMS from the gaussian fit
	  MeanPosition = 0;						// stores the peak time position of the hit
	  MeanPosError = 0;						// stores the error assoc. with thte peak time of the hit
	  Charge = 0;      						// stores the total charge assoc. with the hit
	  ChargeError = 0;              					// stores the error on the charge
	  Amp = 0;							// stores the amplitude of the hit
	  AmpError = 0;							// stores the error assoc. with the amplitude
	  NumOfHits = 0;   						// stores the multiplicity of the hit
	  FitGoodnes = 0;							// stores the Chi2/NDF of the hit
	  FitNDF = -1;							// stores the NDF of the hit
	  
	  
	  
	  //stores gaussian paramters first index is the hit number
	  //the second refers to height, position, and width respectively
	  std::vector<double>  hitSig;	
	  
	  // ########################################
	  // ### Looping over the hitIndex Vector ###
	  // ########################################
	  while( hitIndex < (signed)startTimes.size() ) {
	    // --- Zeroing Start Times and End Times ---
	    startT = 0;
	    endT   = 0;
	    // Note: numHits is hardcoded to one here (Need to fix!)
	    numHits=1;
	    
	    minPeakHeight = signal[maxTimes[hitIndex]];
	    
	    // consider adding pulse to group of consecutive hits if:
	    // 1 less than max consecutive hits
	    // 2 we are not at the last point in the signal vector
	    // 3 the height of the dip between the two is greater than threshold/2
	    // 4 and there is no gap between them
	    while(numHits < fMaxMultiHit && numHits+hitIndex < (signed)endTimes.size() && 
		  signal[endTimes[hitIndex+numHits-1]] >threshold/2.0 &&  
		  startTimes[hitIndex+numHits] - endTimes[hitIndex+numHits-1] < 2){
	      if(signal[maxTimes[hitIndex+numHits]] < minPeakHeight){ 
		minPeakHeight=signal[maxTimes[hitIndex+numHits]];
		
	      }//<---Only move the mean peakHeight in this hit
	      
	      numHits++;//<---Interate the multihit function
	      
	    }//<---End multihit while loop
	    
	    
	    
	    // ###########################################################	
	    // ### Finding the first point from the beginning (startT) ###
	    // ###     which is greater then 1/2 the smallest peak     ###
	    // ###########################################################
	    startT=startTimes[hitIndex];
	    while(signal[(int)startT] < minPeakHeight/2.0)  {startT++;}
	    
	    // #########################################################	
	    // ### Finding the first point from the end (endT) which ###
	    // ###      is greater then 1/2 the smallest peak        ###
	    // #########################################################
	    endT=endTimes[hitIndex+numHits-1];
	    while(signal[(int)endT] <minPeakHeight/2.0) {endT--;}
	    
	    
	    // ###############################################################
	    // ### Putting the current considered hit into a 1-D Histogram ###
	    // ###############################################################
	    //--- Size of hit = endT - startT ---
	    size = (int)(endT-startT);
	    
	    // ###########################################################################
	    // ###    Bug Fix: For ADC counts occuring at the end of the ticks range   ###
	    // ### the hitfinder incorrectly assigns the size of the hit as a negative ###
	    // ###      number...so we fix this to be 0 so that this hit is skipped    ###
	    // ###########################################################################
	    if(size < 0){size = 0;}
	    
	    
	    // --- TH1D HitSignal ---
	    TH1D hitSignal("hitSignal","",size,startT,endT);
	    
	    for(int i = (int)startT; i < (int)endT; i++){
	      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	      // ~~~ Filling the pulse signals into TH1D histograms ~~~
	      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	      hitSignal.Fill(i,signal[i]);
	      
	      //std::cout<<"i = "<<i<<" , signal[i] = "<<signal[i]<<std::endl;
	    }//<---End for looping over startT and endT
	    
	    
	    // ##################################################################################################
	    // ### Trying to utilize the old approach of builing my TF1 and then implementing this new method ###
	    // ##################################################################################################
	    
	    // ########################################################
	    // ### Building TFormula for seedHit and basic Gaussian ###
	    // ########################################################
	    eqn = "gaus(0)";
	    for(int i = 3; i < numHits*3; i+=3){
	      eqn.append("+gaus(");
	      numConv.str("");
	      numConv << i;
	      eqn.append(numConv.str());
	      eqn.append(")");
	    }
	    
	    
	    
	    // --- TF1 function for GausHit  ---
	    TF1 Gaus("Gaus",eqn.c_str(),0,size);
	    
	    // ##############################################################################
	    // ### For multi-pulse hits we loop over all the hits and fit N Gaus function ###
	    // ##############################################################################
	    if(numHits > 1) {
	      // --- Loop over numHits ---
	      for(int i = 0; i < numHits; ++i){
		
		// ###################################################################
		// ### Set the parameters of the Gaus hit based on the seedHit fit ###
		// ###################################################################
		//amplitude = 0.5*(threshold+signal[maxTimes[hitIndex+i]]);
		
		amplitude = signal[maxTimes[hitIndex+i]];
		Gaus.SetParameter(3*i,amplitude);
		Gaus.SetParameter(1+3*i, maxTimes[hitIndex+i]);
		Gaus.SetParameter(2+3*i, fitWidth);
		Gaus.SetParLimits(3*i, 0.0, 3.0*amplitude);
		Gaus.SetParLimits(1+3*i, startT , endT);
		Gaus.SetParLimits(2+3*i, 0.0, 10.0*fitWidth);
	      }//end loop over hits
	    }//end if numHits > 1
	    
	    // ######################################################################
	    // ### For single pulse hits we perform a fit using a single Gaussian ###
	    // ######################################################################
	    else{
	      
	      Gaus.SetParameters(signal[maxTimes[hitIndex]],maxTimes[hitIndex],fitWidth);
	      Gaus.SetParLimits(0,0.0,1.5*signal[maxTimes[hitIndex]]);
	      Gaus.SetParLimits(1, startT , endT);
	      Gaus.SetParLimits(2,0.0,10.0*fitWidth);
	    }
	    
	    // ####################################################
	    // ### PERFORMING THE TOTAL GAUSSIAN FIT OF THE HIT ###
	    // ####################################################
	    hitSignal.Sumw2();

	    hitSignal.Fit(&Gaus,"QNRW","", startT, endT);
	    //hitSignal.Fit(&Gaus,"QR0LLi","", startT, endT);

	    
	    for(int hitNumber = 0; hitNumber < numHits; ++hitNumber){
	      totSig = 0;	
	      
	      
	      // ###########################################################
	      // ### Record this hit if the amplitude is > threshold/2.0 ###
	      // ###             and the width is > minWidth             ###
	      // ###########################################################
	      if(Gaus.GetParameter(3*hitNumber) > threshold/2.0 && 
		 Gaus.GetParameter(3*hitNumber+2) > minWidth){
		StartIime			= Gaus.GetParameter(3*hitNumber+1) - Gaus.GetParameter(3*hitNumber+2); // position - width
		
		EndTime			= Gaus.GetParameter(3*hitNumber+1) + Gaus.GetParameter(3*hitNumber+2); // position + width
		
		MeanPosition		= Gaus.GetParameter(3*hitNumber+1);
		StartIime			= Gaus.GetParameter(3*hitNumber+2); // width
		
		//StartTimeError			= TMath::Sqrt( (Gaus.GetParError(3*hitNumber+1)*Gaus.GetParError(3*hitNumber+1)) + 
		//					       (Gaus.GetParError(3*hitNumber+2)*Gaus.GetParError(3*hitNumber+2)));
		
		//EndTimeError			= TMath::Sqrt( (Gaus.GetParError(3*hitNumber+1)*Gaus.GetParError(3*hitNumber+1)) + 
		//			                       (Gaus.GetParError(3*hitNumber+2)*Gaus.GetParError(3*hitNumber+2)));
		
		
		//MeanPosError			= Gaus.GetParError(3*hitNumber+1);
		
		
		
		
		hitSig.resize(size);
		for(int sigPos = 0; sigPos<size; sigPos++){ //<---Loop over the size (endT - startT)
		  
		  hitSig[sigPos] = Gaus.GetParameter(3*hitNumber)*TMath::Gaus(sigPos+startT,Gaus.GetParameter(3*hitNumber+1), Gaus.GetParameter(3*hitNumber+2));
		  totSig+=hitSig[(int)sigPos];
		  
		}//<---End Signal postion loop
		
		// --- Getting the total charge using the area method ---
		if(fAreaMethod) {
		  totSig = std::sqrt(2*TMath::Pi())*Gaus.GetParameter(3*hitNumber)*Gaus.GetParameter(3*hitNumber+2)/fAreaNorms[(size_t)sigType];
		  
		}//<---End Area Method
		Charge				= totSig;
		// ---------------------------------------------------------------------------------
		// --- chargeErr=TMath::Sqrt(TMath::Pi())*(amplitudeErr*width+widthErr*amplitude); ---
		// -----------------------------------------------------------------------------------
		ChargeError			= TMath::Sqrt(TMath::Pi())*(Gaus.GetParError(3*hitNumber+0)*Gaus.GetParameter(3*hitNumber+2)+
									    Gaus.GetParError(3*hitNumber+2)*Gaus.GetParameter(3*hitNumber+0));   //estimate from area of Gaussian
	    Amp				= Gaus.GetParameter(3*hitNumber);
	    AmpError			= Gaus.GetParError(3*hitNumber);
	    NumOfHits			= numHits;
	    
	    // #######################################################
	    // ### Using seeded values to get a better estimate of ###
	    // ###        the errors associated with the fit       ###
	    // #######################################################
	    // -----------------------------------------
	    // --- Determing the goodness of the fit ---
	    // -----------------------------------------
	    hit->SetParameters(Amp,MeanPosition, Gaus.GetParameter(3*hitNumber+2));
	    hit->FixParameter(0,Amp);
	    //hit->SetParLimits(0,Amp/2, Amp*2);
	    //hit->FixParameter(1,MeanPosition);
	    hit->SetParLimits(1,MeanPosition - 3, MeanPosition + 3);
	    //hit->FixParameter(2,Gaus.GetParameter(3*hitNumber+2));//<---Should I be setting the fitWidth initally
	    // #############################
	    // ### Perform Hit on Signal ###
	    // #############################
	    
	    
	    float TempStartTime = MeanPosition - 4;
	    float TempEndTime   = MeanPosition  + 4;
	    
	    hitSignal.Fit(hit,"QNRLLi","", TempStartTime, TempEndTime);
	    
	    
	    FitGoodnes			= hit->GetChisquare() / hit->GetNDF();
	    FitNDF    			= hit->GetNDF();
	    
	    StartTimeError			= TMath::Sqrt( (hit->GetParError(1)*hit->GetParError(1)) + 
							       (hit->GetParError(2)*hit->GetParError(2)));
	    
	    EndTimeError			= TMath::Sqrt( (hit->GetParError(1)*hit->GetParError(1)) + 
							       (hit->GetParError(2)*hit->GetParError(2)));
	    
	    
	    MeanPosError			= hit->GetParError(1);
	    
	    // ####################################################################################
	    // ### Skipping hits with a chi2/NDF larger than the value defined in the .fcl file ###
	    // ####################################################################################
	    if(FitGoodnes > fChi2NDF){continue;}
	    
	    // get the WireID for this hit
	    std::vector<geo::WireID> wids = geom->ChannelToWire(channel);
	    ///\todo need to have a disambiguation algorithm somewhere in here
	    // for now, just take the first option returned from ChannelToWire
	    geo::WireID wid = wids[0];
	    
	    recob::Hit hit(
	      channel,                 // channel
	      (raw::TDCtick_t) startT, // start_tick
	      (raw::TDCtick_t) endT,   // end_tick
	      MeanPosition,            // peak_time
	      MeanPosError,            // sigma_peak_time
	      RMS,                     // rms,
	      Amp,                     // peak_amplitude,
	      AmpError,                // sigma_peak_amplitude
	      std::accumulate          // summedADC
	        (signal.begin() + TempStartTime, signal.begin() + TempEndTime, 0.),
	      Charge,                  // hit_integral
	      ChargeError,             // hit_sigma_integral
	      NumOfHits,               // multiplicity
	      hitNumber,               // local_index
	      FitGoodnes,              // goodness_of_fit
	      FitNDF,                  // dof
	      view,                    // view
	      sigType,                 // signal_type
	      wid                      // wireID
	      );
	    hcol.emplace_back(std::move(hit), digitVec);
	    
	      }//<---End looking at hits over threshold
	      
	    }//<---End Looping over numHits
	    
	    
	    hitIndex+=numHits;
	  }//<---End while hitIndex < startTimes.size()
	  
	  
	  
	  //}//<---End looping over wireIter
	  
	  
	  //std::cout << "Test1" << std::endl;
	 
	  
	} // end loop over nonzero blocks
	
	
      } // end if not a bad channel 
      
      
    } // end loop over wires
    
    // move the hit collection and the associations into the event
    hcol.put_into(evt);

    delete chanFilt;	  
    
    
    return;
  }
  
} // end namespace calgaushf


#endif // CALGAUSHFDUNE10KT_H
