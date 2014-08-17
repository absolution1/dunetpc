///////////////////////////////////////////////
// $Id: MakeResp.cxx,v 1.0 2010/07/02 20:42:31 bpage Exp $
//
// This module produces the shape file used for the 
// argoneut CalData deconvolution.
// pagebri3@msu.edu
//
////////////////////////////////////////////////////////////////////////
#ifndef MAKERESP_H
#define MAKERESP_H

extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}

// Framework includes
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Core/EDProducer.h" // include the proper bit of the framework
#include "messagefacility/MessageLogger/MessageLogger.h" 

// LArSoft includes
#include "Utilities/LArFFT.h"
#include "Utilities/LArProperties.h"
//#include "CalData/MakeResp.h"
#include "art/Framework/Core/ModuleMacros.h" 


// ROOT includes
#include <TFile.h>
#include <TH2D.h>
#include <TF1.h>
#include <TComplex.h>
#include <iostream>

namespace caldata{

  class MakeResp : public art::EDProducer { 
  
  public:
 
    explicit MakeResp(fhicl::ParameterSet const& pset);
  
    virtual ~MakeResp();
  
    void reconfigure(fhicl::ParameterSet p);

    void produce(art::Event & evt);
 
    void beginJob();
    
  private:

  protected:
  }; //class MakeResp
}

namespace caldata{

  MakeResp::MakeResp(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset);
  }
  
  MakeResp::~MakeResp()
  {
  }
 
  void MakeResp::reconfigure(fhicl::ParameterSet p)
  {
  }

  void MakeResp::produce(art::Event & evt)
  {
  }

  void MakeResp::beginJob()
  {
    art::ServiceHandle<util::LArFFT> FFT;
    int signalSize = FFT->FFTSize();
    std::vector<double> diffsig(signalSize);
    TComplex kernBin;
    int size = signalSize/2;
    int bin=0;
    std::vector<TComplex> freqSig(size+1);
    TF1 * indFilter = new TF1("indFilter","[0]*x^[1]/([2]^[1]+x^[1])*exp(-0.5*(x/[3])^[4])",0,size+1);
    indFilter->SetParameters(15200.0,1.0,50,125,4); 
    TF1 * colFilter = new TF1("colFilter","(x>0)*[0]*exp(-0.5*(((x-[1])/[2])^2)^[3])/(1.0+(x/200)^4)",0,size+1);
    colFilter->SetParameters(.698,0.0,110.0,1.5);
    std::vector<double> bipolar(signalSize);
    double rampHeight=0.0;
    double bumpHeight=5.0;
    double asymmetry=1.2;
    for(double i = 1; i < 600; i++)  {
      if(i<20)bipolar[i]=-asymmetry*(rampHeight*2*i/pow(1.0+4*i*i/361.0,1.5)+bumpHeight*tanh(2*i)*exp(-0.5*pow(2*i/16.8,8.0)));
      bipolar[signalSize-i]=rampHeight*2*i/pow(1.0+4*i*i/361.0,1.5)+bumpHeight*tanh(2*i)*exp(-0.5*pow(2*i/16.8,8.0));
    }
    std::vector<double> gaussian(signalSize);
    gaussian[signalSize-4]=1.0;
    gaussian[signalSize-3]=8.0;
    gaussian[signalSize-2]=28.0;
    gaussian[signalSize-1]=55.0;
    gaussian[0]=70.0;
    gaussian[1]=55.0;
    gaussian[2]=28.0;
    gaussian[3]=8.0;
    gaussian[4]=1.0;
    FFT->Convolute(bipolar,gaussian);
    std::vector<double> electronics(signalSize);
    TF1 * expt = new TF1("expt","[2]+((exp([0]*x)-1.0)/(1.0+exp([0]*x)))^[1]",0,2*signalSize+4);
  
    expt->SetParameters(.113429,3.7389,.335); //.335
    for(int i = 0; i < signalSize; i++){
      electronics[i] = expt->Eval(i<3000?2*i:6000);
      //std::cout<<i<<" "<<electronics[i]<<std::endl;
    }
  
    TFile * test = new TFile("/grid/fermiapp/lbne/lar/aux/ArgoResponse1.5.root","READ");
    //TFile * shapeFile = new TFile("/home/bpage/shape-argo.root","READ");
    TH1D * decay =  new TH1D(*(TH1D*) test->Get("real/decayHist"));
    double avgDecay(0.0);
    const int iRange(2);
    const int nRanges(4);
    std::vector<double> ranges(nRanges);
    double temp[nRanges]={0,48,96,144};
    for(int i = 0; i < nRanges; i++) ranges[i]=temp[i];
    std::vector<double> decays(nRanges,0);
    int j=0;
    double decayHolder;
    for(int i = iRange; i < nRanges-1; i++) {
      for(j=0; j < ranges[i+1]-ranges[i]; j++){
	//decayHolder = decay->GetBinContent(240+j);//240 is first argoneut collection wire
	decayHolder = -3e-3;
	avgDecay+=decayHolder;
	decays[i]+=decayHolder;
      }
    }
    for(int i = 0; i < nRanges-1; i++) decays[i]/=(double)(ranges[i+1]-ranges[i]);
    avgDecay/=48.0; //was 234.0
    TFile * out = new TFile("shape.root","recreate");
  
    TDirectory *real = out->mkdir("real");
    real->AddDirectory();
    real->cd();
    //art::ServiceHandle<art::TFileService> tfs;
    TH1D * elect = new TH1D("electronics","electronics shaping",signalSize,0,signalSize);

    TH1D* decayHistNew = new TH1D("decayHist","RC decay Constants",144,0,144);
    for(int i = iRange; i < nRanges-1; i++) 
      for(j=0; j < ranges[i+1]-ranges[i]; j++)
	decayHistNew->Fill(ranges[i]+j,decays[i]);
    
    decayHistNew->Write();
    delete test;
    std::vector<double> induction(signalSize);
    std::vector<double> delta(signalSize);
    delta[signalSize-1]=1.0;
    delta[0]=1.0;
    induction[0]=0;
    for(double i = 1.0; i < signalSize; i++) {
      induction[i]=electronics[i]*exp(i*avgDecay)-electronics[i-1]*exp((i-1.0)*avgDecay);
      if(i<5 && induction[i]< 0) induction[i]=0;
    }
    
    for(int i = 0; i < signalSize; i++) elect->Fill(i,induction[i]);
    FFT->Convolute(induction,bipolar);
    TH1D * indShape = new TH1D("indShape","Induction Shape",signalSize,0,signalSize);
  
    TH1D * colShape = new TH1D("colShape","Collection Shape",signalSize,0,signalSize); 
    TH1D * bipol = new TH1D("bipol","induction input shape",signalSize,0,signalSize);
    for(int i = 0; i < signalSize; i++) { 
      indShape->Fill(i,induction[i]);
      bipol->Fill(i,bipolar[i]);
    }
    TH1D * indPow = new TH1D("indPow","Induction power spectrum",size+1,0,size+1);
    TH1D * colPow = new TH1D("colPow","Collection power spectrum",size+1,0,size+1);
    TH1D * indFil = new TH1D("indFil","Induction Filter",size+1,0,size+1);
    TH1D * colFil = new TH1D("colFil","Collection Filter",size+1,0,size+1);
    for(bin = 0; bin < size+1; bin++) {
      indFil->Fill(bin,indFilter->Eval(bin));
      colFil->Fill(bin,colFilter->Eval(bin));
    }

    TH2D * RespRe = new TH2D("RespRe",
			     " Response Functions - Real Part",
			     nRanges-1,&ranges[0],size+1, 0, size+1);
    TH2D * RespIm = new TH2D("RespIm",
			     "ResponseFunctions - Imaginary Part",   
			     nRanges-1,&ranges[0],size+1, 0, size+1);

    std::vector<double> dataTemp;
    for(int i = 0; i < nRanges-1; i++) {
      if(i < iRange) {
	dataTemp = induction;
      }
      else {
	dataTemp[0]=0;
	for(double j = 1.0; j < signalSize; j++){ 
	  dataTemp[j]=electronics[j+1]*exp(decays[i]*(j+1.0))-electronics[j]*exp(decays[i]*j);
	  if(j<5 && dataTemp[j]< 0) dataTemp[j]=0;
	}
      }
      FFT->AlignedSum(dataTemp,delta,false);
      FFT->DoFFT(dataTemp,freqSig);
      for(bin = 0; bin < signalSize/2+1; bin++) {
	if(i<iRange) {
	  kernBin=indFilter->Eval(0.5*bin)/freqSig[bin];
	  indPow->Fill(bin,freqSig[bin].Rho());
	}
	else {
	  kernBin=colFilter->Eval(0.5*bin)/freqSig[bin];
	  colPow->Fill(bin,freqSig[bin].Rho());
	}
//	switch((int)ranges[i]) {
//	case 31:
//	case 41:
//	case 108:
//	case 120:
//	case 121:
//	case 124:
//	case 139:
//	case 200:
//	case 204:
//	case 212:
//	  kernBin*=1.0-1.0/(1.0+pow((double)bin/6.0,8.0));
//	  break;
//	}
	RespRe->Fill(ranges[i],bin,kernBin.Re());     
	RespIm->Fill(ranges[i],bin,kernBin.Im());
      }
    }
    elect->Write();
    indShape->Write();
    colShape->Write();
    indPow->Write();
    colPow->Write();
    indFil->Write();
    colFil->Write();
    RespRe->Write();
    RespIm->Write();
    real->cd("../");
    TF1* indFilterS = new TF1("indFilterS","(x<2048)*(1.0-gaus(3))*[0]*exp(-0.5*(x/[1])^[2])",0,size+1);
    indFilterS->SetParameters(2.34e4,125.0,4.0,1.0,0.0,15); 
    TF1 * colFilterS = new TF1("colFilterS","[0]*exp(-0.5*(((x-[1])/[2])^2)^[3])/(1.0+(x/200)^4)",0,size+1);
    colFilterS->SetParameters(2.0e7,0.0,120.0,2.0);
    TDirectory * sim = out->mkdir("sim");
    sim->cd(); 
//    TF1 * noise = new TF1("noise","[0]*x*exp(-x/[1])",0,size+1);
//    noise->SetParameters(11.0,61.3);
//    TH1D * noise1 = new TH1D("noise","noise",size+1,0,size+1);
//    for(int i = 0; i < size+1; i++) noise1->Fill(i,noise->Eval(i));
//    TH1D * simShape= (TH1D*) shapeFile->Get("shape");
//    for(int i = 0; i < signalSize; i++) electronics[i]=simShape->GetBinContent(i+1);
//    bipolar.clear();
//    bipolar.resize(signalSize);
//    std::vector<double> ramp(signalSize);
//    ranges.clear();
//    ranges.resize(3);
//    ranges[0]=0;
//    ranges[1]=240;
//    ranges[2]=480;
//    TH2D * RespRe1 = new TH2D("RespRe","Response Functions - Real Part",2,&ranges[0],signalSize/2+1,0,signalSize/2+1);
//    TH2D* RespIm1 = new TH2D("RespIm","Response Functions - Imaginary Part",2,&ranges[0],signalSize/2+1,0,signalSize/2+1);
//    TH1D * indFil1 = new TH1D("indFil","Induction Filter",size+1,0,size+1);
//    TH1D * colFil1 = new TH1D("colFil","Collection Filter",size+1,0,size+1);
//    TH1D * decayHist1 = new TH1D("decayHist","Decay Constants",480,0,480);
//    TH1D * bipol1 = new TH1D("bipol","induction input shape",signalSize,0,signalSize);
//    for(int i = ranges[1]; i < ranges[2]; i++) 
//      decayHist1->Fill(i,-.00278871);
//    art::ServiceHandle<util::LArProperties> larp;
//    double driftvelocity=larp->DriftVelocity(larp->Efield(),larp->Temperature())/1000.;  
//    double sampleRate = 198.0;
//    double nbinc = 2.5*(.4)/(driftvelocity*sampleRate); ///number of bins //KP
//    double nbini = 1.5*(.4)/(driftvelocity*sampleRate);//KP
//    TF1 rampf("rampf","[0]*x^2/(1.0+(x/[1])^[2])",0,ramp.size());
//    rampf.SetParameters(.68/pow(nbinc,3.0),2*nbini,32.0);
//    for(size_t i = 0; i < ramp.size(); i++) ramp[i]=rampf.Eval(i);
//    
//    double crossRate = driftvelocity*sampleRate;
//    TF1 bipolf("bipolf","[0]*[1]*([2]-x)*exp(-0.5*(TMath::Abs(x-[2])/[3])^[4])",0,bipolar.size());
//    bipolf.SetParameters(.3/nbini,2.0*crossRate,bipolar.size()/2.0,nbini/2.0,2.0);
//    for(size_t i = 0; i < bipolar.size(); i++) bipolar[i]=bipolf.Eval(i);
//    FFT->ShiftData(bipolar,bipolar.size()/2.0+nbini);
//    FFT->ShiftData(ramp,2.0*nbini-nbinc);
//    
//    for(int i = 0; i < signalSize; i++) bipol1->Fill(i,bipolar[i]);
//    FFT->Convolute(bipolar,electronics);
//    TH1D * indShape1 = new TH1D("indShape","Induction Shape",signalSize,0,signalSize);
//    for(int i = 0; i < signalSize; i++) indShape1->Fill(i,bipolar[i]);
//    FFT->DoFFT(bipolar,freqSig);
//    std::vector<double> indWFilter (signalSize);
//    std::vector<double> colWFilter (signalSize);
//    for(int i = 0 ; i < size+1; i++) indWFilter[i]=118.0*(0.5+0.5*tanh((double)i/5.0))*exp(-0.5*pow((double)(i-25)/65.0,2.0))/(1.0+pow((double)i/225.0,8.0));//*freqSig[i].Rho2()/(freqSig[i].Rho2()+noise->Eval(i)*noise->Eval(i));
//    TH1D * indPow1 = new TH1D("indPow","Induction power spectrum",size+1,0,size+1);
//    TH1D * colPow1 = new TH1D("colPow","Collection power spectrum",size+1,0,size+1);
//    for(int i = 0; i <signalSize/2+1; i++) {
//      indPow1->Fill(i,freqSig[i].Rho());
//      kernBin=indWFilter[i]/freqSig[i];
//      if(i==0 || i>790) kernBin=0;
//      indFil1->Fill(i,indWFilter[i]);
//      RespRe1->Fill(0.0,i,kernBin.Re()); 
//      RespIm1->Fill(0.0,i,kernBin.Im());
//    }
//    FFT->Convolute(ramp,electronics);
//    TH1D * colShape1 = new TH1D("colShape","Collection Shape",signalSize,0,signalSize);
//    for(int i = 0; i < signalSize; i++) colShape1->Fill(i,ramp[i]);
//    FFT->DoFFT(ramp,freqSig);
//    for(int i = 0 ; i < signalSize; i++) colWFilter[i]=406.0*exp(-0.5*pow((double)(i-0)/350.0,2.0))/(1.0+pow((double)i/400.0,8.0))*freqSig[i].Rho2()/(freqSig[i].Rho2()+noise->Eval(i)*noise->Eval(i));
//    for(int i = 0; i <signalSize/2+1; i++) {
//      if(i==0) kernBin=0;
//      colPow1->Fill(i,freqSig[i].Rho());
//      kernBin=colWFilter[i]/freqSig[i];
//      colFil1->Fill(i,colWFilter[i]);
//      RespRe1->Fill(241.0,i,kernBin.Re()); 
//      RespIm1->Fill(241.0,i,kernBin.Im());
//    }
    indShape->Write();
    colShape->Write();
    indPow->Write();
    colPow->Write();
    indFil->Write();
    colFil->Write();
    RespRe->Write();
    RespIm->Write();
    decayHistNew->Write();
    //noise->Write();
    out->Write();
  }
}
namespace caldata{

  DEFINE_ART_MODULE(MakeResp);
  
} // end namespace caldata

#endif // MAKERESP_H
