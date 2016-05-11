/////////////////////////////////////////////////////////////////////////////////////////////////////
// Mike Wallbank (m.wallbank@sheffield.ac.uk), August 2015
//
// Macro to produce energy-charge conversion for the DUNE FD/35t.
//
// Useage:
//  root getEnergyConversion.C
//
//  -- runs over the output file (EMEnergyCalib.root) produced from the EMEnergyCalib_module art
//     analyser
//  -- produces an output root file which contains all of the plots used to find the conversion
//     from deposited charge to true energy
//  -- function is of the form (straight line)
//         energy = A + (B * charge) and A and B are determined from this script
//  -- linear plots showing relation are found in the output root file
//  -- the values of A and B for each plane are also printed upon completion of the script and can
//     be used in reconstruction; in DUNE the values are placed in the fhicl file:
//         larreco/larreco/RecoAlg/showeralgorithms.fcl
/////////////////////////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <iostream>
#include <utility>
#include <map>
#include <algorithm>
#include "TTree.h"
#include "TMath.h"

const int kMaxHits = 10000;

class EMEnergyConversion {
public:

  EMEnergyConversion(TTree* tree);
  ~EMEnergyConversion();
  double ConvertChargeToEnergy(double charge, int plane);
  void MakeFits();
  std::pair<double,double> MakeFit(int num, TString plane);
  void Run();
  void SaveHists();
  void SetBranchAddresses();

private:

  TTree* fTree;

  double TrueEnergy;
  double DepositU;
  double DepositV;
  double DepositZ;
  double VertexDetectorDist;
  int    NHits;
  int    Hit_TPC      [kMaxHits];
  int    Hit_Plane    [kMaxHits];
  int    Hit_Wire     [kMaxHits];
  int    Hit_Channel  [kMaxHits];
  double Hit_PeakT    [kMaxHits];
  double Hit_Charge   [kMaxHits];
  int    Hit_ClusterID[kMaxHits];

  TH2D* ChargeDepositEnergyU;
  TH2D* ChargeDepositEnergyV;
  TH2D* ChargeDepositEnergyZ;
  TH2D* EnergyDepositUDistance;
  TH2D* EnergyDepositVDistance;
  TH2D* EnergyDepositZDistance;
  TH1D* EnergyCompleteness;

  TFile *outFile;

  // Convert from charge to energy
  static const double Uintercept = -1519.33, Ugradient = 148867;
  static const double Vintercept = -1234.91, Vgradient = 149458;
  static const double Zintercept = -1089.73, Zgradient = 145372;

};

EMEnergyConversion::EMEnergyConversion(TTree* tree) {
  fTree = tree;
  ChargeDepositEnergyU = new TH2D("ChargeDepositEnergyU","Charge v Depositied Energy on U;Deposited Energy (GeV);Charge (ADC);",50,0,1,100,0,200000);
  ChargeDepositEnergyV = new TH2D("ChargeDepositEnergyV","Charge v Depositied Energy on V;Deposited Energy (GeV);Charge (ADC);",50,0,1,100,0,200000);
  ChargeDepositEnergyZ = new TH2D("ChargeDepositEnergyZ","Charge v Depositied Energy on Z;Deposited Energy (GeV);Charge (ADC);",50,0,1,100,0,200000);
  EnergyDepositUDistance = new TH2D("EnergyDepositUDistance","Deposited Energy on U vs Distance from Detector Edge;Distance (cm);Fraction of Energy Deposited;",100,0,220,50,0,1.1);
  EnergyDepositVDistance = new TH2D("EnergyDepositVDistance","Deposited Energy on V vs Distance from Detector Edge;Distance (cm);Fraction of Energy Deposited;",100,0,220,50,0,1.1);
  EnergyDepositZDistance = new TH2D("EnergyDepositZDistance","Deposited Energy on Z vs Distance from Detector Edge;Distance (cm);Fraction of Energy Deposited;",100,0,220,50,0,1.1);
  EnergyCompleteness = new TH1D("EnergyCompleteness","Fraction of Energy Reconstructed in Largest Cluster;Energy Completeness;",101,0,1.01);
  outFile = TFile::Open("EnergyFits.root","RECREATE");
}

EMEnergyConversion::~EMEnergyConversion() {
  outFile->Close();
  delete outFile;
}

double EMEnergyConversion::ConvertChargeToEnergy(double charge, int plane) {

  double energy;

  switch (plane) {
  case 0:
    energy = (double)(charge - Uintercept)/(double)Ugradient;
    break;
  case 1:
    energy = (double)(charge - Vintercept)/(double)Vgradient;
    break;
  case 2:
    energy = (double)(charge - Zintercept)/(double)Zgradient;
    break;
  }

  return energy;

}

void EMEnergyConversion::SetBranchAddresses() {
  fTree->SetBranchAddress("TrueEnergy",&TrueEnergy);
  fTree->SetBranchAddress("DepositU",&DepositU);
  fTree->SetBranchAddress("DepositV",&DepositV);
  fTree->SetBranchAddress("DepositZ",&DepositZ);
  fTree->SetBranchAddress("VertexDetectorDist",&VertexDetectorDist);
  fTree->SetBranchAddress("NHits",&NHits);
  fTree->SetBranchAddress("Hit_TPC",&Hit_TPC);
  fTree->SetBranchAddress("Hit_Plane",&Hit_Plane);
  fTree->SetBranchAddress("Hit_Wire",&Hit_Wire);
  fTree->SetBranchAddress("Hit_Channel",&Hit_Channel);
  fTree->SetBranchAddress("Hit_PeakT",&Hit_PeakT);
  fTree->SetBranchAddress("Hit_Charge",&Hit_Charge);
  fTree->SetBranchAddress("Hit_ClusterID",&Hit_ClusterID);
}

void EMEnergyConversion::MakeFits() {

  // Fit for U,V,Z plane
  std::pair<double,double> fitParamsU = this->MakeFit(0, TString("U"));
  std::pair<double,double> fitParamsV = this->MakeFit(1, TString("V"));
  std::pair<double,double> fitParamsZ = this->MakeFit(2, TString("Z"));

  std::cout << "U plane fit has intercept " << fitParamsU.first << " and gradient " << fitParamsU.second << std::endl;
  std::cout << "V plane fit has intercept " << fitParamsV.first << " and gradient " << fitParamsV.second << std::endl;
  std::cout << "Z plane fit has intercept " << fitParamsZ.first << " and gradient " << fitParamsZ.second << std::endl;

}

std::pair<double,double> EMEnergyConversion::MakeFit(int num, TString plane) {

  // Make a load of histograms!
  TH1D* ChargeDist1 = new TH1D(TString("ChargeDist1Plane")+plane,";Total ADC;",200,0,30000);
  TH1D* ChargeDist2 = new TH1D(TString("ChargeDist2Plane")+plane,";Total ADC;",200,10000,40000);
  TH1D* ChargeDist3 = new TH1D(TString("ChargeDist3Plane")+plane,";Total ADC;",200,20000,60000);
  TH1D* ChargeDist4 = new TH1D(TString("ChargeDist4Plane")+plane,";Total ADC;",200,40000,70000);
  TH1D* ChargeDist5 = new TH1D(TString("ChargeDist5Plane")+plane,";Total ADC;",200,50000,90000);
  TH1D* ChargeDist6 = new TH1D(TString("ChargeDist6Plane")+plane,";Total ADC;",200,60000,110000);
  TH1D* ChargeDist7 = new TH1D(TString("ChargeDist7Plane")+plane,";Total ADC;",200,80000,120000);
  TH1D* ChargeDist8 = new TH1D(TString("ChargeDist8Plane")+plane,";Total ADC;",200,90000,140000);
  TH1D* ChargeDist9 = new TH1D(TString("ChargeDist9Plane")+plane,";Total ADC;",200,110000,150000);
  TH1D* ChargeDist10 = new TH1D(TString("ChargeDist10Plane")+plane,";Total ADC;",200,120000,170000);
  TH1D* EnergyDist1 = new TH1D(TString("EnergyDist1Plane")+plane,";Energy (GeV);",100,0.05,0.1);
  TH1D* EnergyDist2 = new TH1D(TString("EnergyDist2Plane")+plane,";Energy (GeV);",100,0.15,0.2);
  TH1D* EnergyDist3 = new TH1D(TString("EnergyDist3Plane")+plane,";Energy (GeV);",100,0.25,0.3);
  TH1D* EnergyDist4 = new TH1D(TString("EnergyDist4Plane")+plane,";Energy (GeV);",100,0.35,0.4);
  TH1D* EnergyDist5 = new TH1D(TString("EnergyDist5Plane")+plane,";Energy (GeV);",100,0.45,0.5);
  TH1D* EnergyDist6 = new TH1D(TString("EnergyDist6Plane")+plane,";Energy (GeV);",100,0.55,0.6);
  TH1D* EnergyDist7 = new TH1D(TString("EnergyDist7Plane")+plane,";Energy (GeV);",100,0.65,0.7);
  TH1D* EnergyDist8 = new TH1D(TString("EnergyDist8Plane")+plane,";Energy (GeV);",100,0.75,0.8);
  TH1D* EnergyDist9 = new TH1D(TString("EnergyDist9Plane")+plane,";Energy (GeV);",100,0.85,0.9);
  TH1D* EnergyDist10 = new TH1D(TString("EnergyDist10Plane")+plane,";Energy (GeV);",100,0.95,1.0);

  for (unsigned int event = 0; event < fTree->GetEntriesFast(); ++event) {
    fTree->GetEntry(event);
    double deposit;
    switch (num) {
    case 0: deposit = DepositU;
    case 1: deposit = DepositV;
    case 2: deposit = DepositZ;
    }
    double planeCharge = 0;
    for (unsigned int hit = 0; hit < NHits; ++hit)
      if (Hit_Plane[hit] == num) planeCharge += (Hit_Charge[hit] * TMath::Exp((500 * Hit_PeakT[hit])/3e6));
    if (deposit >= 0.05 && deposit <= 0.1) { ChargeDist1->Fill(planeCharge); EnergyDist1->Fill(deposit); }
    if (deposit >= 0.15 && deposit <= 0.2) { ChargeDist2->Fill(planeCharge); EnergyDist2->Fill(deposit); }
    if (deposit >= 0.25 && deposit <= 0.3) { ChargeDist3->Fill(planeCharge); EnergyDist3->Fill(deposit); }
    if (deposit >= 0.35 && deposit <= 0.4) { ChargeDist4->Fill(planeCharge); EnergyDist4->Fill(deposit); }
    if (deposit >= 0.45 && deposit <= 0.5) { ChargeDist5->Fill(planeCharge); EnergyDist5->Fill(deposit); }
    if (deposit >= 0.55 && deposit <= 0.6) { ChargeDist6->Fill(planeCharge); EnergyDist6->Fill(deposit); }
    if (deposit >= 0.65 && deposit <= 0.7) { ChargeDist7->Fill(planeCharge); EnergyDist7->Fill(deposit); }
    if (deposit >= 0.75 && deposit <= 0.8) { ChargeDist8->Fill(planeCharge); EnergyDist8->Fill(deposit); }
    if (deposit >= 0.85 && deposit <= 0.9) { ChargeDist9->Fill(planeCharge); EnergyDist9->Fill(deposit); }
    if (deposit >= 0.95 && deposit <= 1.0) { ChargeDist10->Fill(planeCharge); EnergyDist10->Fill(deposit); }
  }

  TF1* fit;
  double charge[10], energy[10];
  outFile->cd();

  ChargeDist1->Fit("gaus"); fit = ChargeDist1->GetFunction("gaus"); charge[0] = fit->GetParameter(1); ChargeDist1->Write();
  ChargeDist2->Fit("gaus"); fit = ChargeDist2->GetFunction("gaus"); charge[1] = fit->GetParameter(1); ChargeDist2->Write();
  ChargeDist3->Fit("gaus"); fit = ChargeDist3->GetFunction("gaus"); charge[2] = fit->GetParameter(1); ChargeDist3->Write();
  ChargeDist4->Fit("gaus"); fit = ChargeDist4->GetFunction("gaus"); charge[3] = fit->GetParameter(1); ChargeDist4->Write();
  ChargeDist5->Fit("gaus"); fit = ChargeDist5->GetFunction("gaus"); charge[4] = fit->GetParameter(1); ChargeDist5->Write();
  ChargeDist6->Fit("gaus"); fit = ChargeDist6->GetFunction("gaus"); charge[5] = fit->GetParameter(1); ChargeDist6->Write();
  ChargeDist7->Fit("gaus"); fit = ChargeDist7->GetFunction("gaus"); charge[6] = fit->GetParameter(1); ChargeDist7->Write();
  ChargeDist8->Fit("gaus"); fit = ChargeDist8->GetFunction("gaus"); charge[7] = fit->GetParameter(1); ChargeDist8->Write();
  ChargeDist9->Fit("gaus"); fit = ChargeDist9->GetFunction("gaus"); charge[8] = fit->GetParameter(1); ChargeDist9->Write();
  ChargeDist10->Fit("gaus"); fit = ChargeDist10->GetFunction("gaus"); charge[9] = fit->GetParameter(1); ChargeDist10->Write();
  energy[0] = EnergyDist1->GetMean(); EnergyDist1->Write();
  energy[1] = EnergyDist2->GetMean(); EnergyDist2->Write();
  energy[2] = EnergyDist3->GetMean(); EnergyDist3->Write();
  energy[3] = EnergyDist4->GetMean(); EnergyDist4->Write();
  energy[4] = EnergyDist5->GetMean(); EnergyDist5->Write();
  energy[5] = EnergyDist6->GetMean(); EnergyDist6->Write();
  energy[6] = EnergyDist7->GetMean(); EnergyDist7->Write();
  energy[7] = EnergyDist8->GetMean(); EnergyDist8->Write();
  energy[8] = EnergyDist9->GetMean(); EnergyDist9->Write();
  energy[9] = EnergyDist10->GetMean(); EnergyDist10->Write();

  delete fit;

  TGraph *graph = new TGraph(10, energy, charge);
  graph->SetName(TString("FitPlane")+plane);
  graph->GetXaxis()->SetTitle("Deposited Energy (GeV)");
  graph->GetYaxis()->SetTitle("Charge");
  graph->SetTitle(plane);
  graph->SetMarkerStyle(8);
  graph->SetMarkerSize(1);
  graph->Fit("pol1");
  graph->Write();
  TF1* fit = graph->GetFunction("pol1");
  std::pair<double,double> fitParameters = std::make_pair(fit->GetParameter(0), fit->GetParameter(1));
  delete fit;
  delete graph;

  return fitParameters;

}

void EMEnergyConversion::Run() {

  this->SetBranchAddresses();

  for (unsigned int event = 0; event < fTree->GetEntriesFast(); ++event) {
    if (event % 1000 == 0) std::cout << "Processing event " << event << std::endl;
    fTree->GetEntry(event);
    this->ProcessEvent();
  }

  this->MakeFits();
  this->SaveHists();

}

void EMEnergyConversion::ProcessEvent() {

  double chargeU = 0, chargeV = 0, chargeZ = 0;
  std::map<int,double> clusterChargeU, clusterChargeV, clusterChargeZ;

  for (int hit = 0; hit < NHits; ++hit) {
    switch (Hit_Plane[hit]) {
    case 0:
      chargeU += (Hit_Charge[hit] * TMath::Exp((500 * Hit_PeakT[hit])/3e6));
      clusterChargeU[Hit_ClusterID[hit]] += (Hit_Charge[hit] * TMath::Exp((500 * Hit_PeakT[hit])/3e6));
      break;
    case 1:
      chargeV += (Hit_Charge[hit] * TMath::Exp((500 * Hit_PeakT[hit])/3e6));
      clusterChargeV[Hit_ClusterID[hit]] += (Hit_Charge[hit] * TMath::Exp((500 * Hit_PeakT[hit])/3e6));
      break;
    case 2:
      chargeZ += (Hit_Charge[hit] * TMath::Exp((500 * Hit_PeakT[hit])/3e6));
      clusterChargeZ[Hit_ClusterID[hit]] += (Hit_Charge[hit] * TMath::Exp((500 * Hit_PeakT[hit])/3e6));
      break;
    }

  } // hit loop

  // Find the highest charge cluster
  double highChargeU = 0;
  for (std::map<int,double>::iterator chargeIt = clusterChargeU.begin(); chargeIt != clusterChargeU.end(); ++chargeIt)
    if (chargeIt->second > highChargeU) highChargeU = chargeIt->second;
  double highChargeV = 0;
  for (std::map<int,double>::iterator chargeIt = clusterChargeV.begin(); chargeIt != clusterChargeV.end(); ++chargeIt)
    if (chargeIt->second > highChargeV) highChargeV = chargeIt->second;
  double highChargeZ = 0;
  for (std::map<int,double>::iterator chargeIt = clusterChargeZ.begin(); chargeIt != clusterChargeZ.end(); ++chargeIt)
    if (chargeIt->second > highChargeZ) highChargeZ = chargeIt->second;

  ChargeDepositEnergyU->Fill(DepositU, chargeU);
  ChargeDepositEnergyV->Fill(DepositV, chargeV);
  ChargeDepositEnergyZ->Fill(DepositZ, chargeZ);
  EnergyDepositUDistance->Fill(VertexDetectorDist, (double)DepositU/(double)TrueEnergy);
  EnergyDepositVDistance->Fill(VertexDetectorDist, (double)DepositV/(double)TrueEnergy);
  EnergyDepositZDistance->Fill(VertexDetectorDist, (double)DepositZ/(double)TrueEnergy);
  EnergyCompleteness->Fill((double)ConvertChargeToEnergy(highChargeU,0)/(double)DepositU);
  EnergyCompleteness->Fill((double)ConvertChargeToEnergy(highChargeV,1)/(double)DepositV);
  EnergyCompleteness->Fill((double)ConvertChargeToEnergy(highChargeZ,2)/(double)DepositZ);

}

void EMEnergyConversion::SaveHists() {

  outFile->cd();

  TCanvas cChargeDepositEnergyU = TCanvas("cChargeDepositEnergyU","",800,600);
  ChargeDepositEnergyU->Draw("colz");
  cChargeDepositEnergyU.Write("ChargeDepositEnergyU");
  TCanvas cChargeDepositEnergyV = TCanvas("cChargeDepositEnergyV","",800,600);
  ChargeDepositEnergyV->Draw("colz");
  cChargeDepositEnergyV.Write("ChargeDepositEnergyV");
  TCanvas cChargeDepositEnergyZ = TCanvas("cChargeDepositEnergyZ","",800,600);
  ChargeDepositEnergyZ->Draw("colz");
  cChargeDepositEnergyZ.Write("ChargeDepositEnergyZ");

  EnergyDepositUDistance->Write("EnergyDepositUDistance");
  EnergyDepositVDistance->Write("EnergyDepositVDistance");
  EnergyDepositZDistance->Write("EnergyDepositZDistance");

  EnergyCompleteness->Write("EnergyCompleteness");

}

void getEnergyConversion() {

  //gStyle->SetOptStat(0);
  //TFile* inFile = TFile::Open("/pnfs/lbne/persistent/users/tjyang/v04_20_00/mergeana/prod_gamma_0.1-1.0GeV_dune35t/2884577_0/EMEnergyCalib.root");
  TFile* inFile = TFile::Open("EMEnergyCalib.root");
  TTree* tree = (TTree*)inFile->Get("energyCalib/EMEnergyCalib");

  EMEnergyConversion emenergycalib = EMEnergyConversion(tree);
  emenergycalib.Run();

  inFile->Close();
  inFile->Delete();

}
