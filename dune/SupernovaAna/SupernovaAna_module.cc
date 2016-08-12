 ////////////////////////////////////////////////////////////////////////
// Class:       SupernovaAna
// Module Type: analyzer
// File:        SupernovaAna_module.cc
//
// Generated at Mon Jul 11 21:36:48 2016 by Michael Baird using the old
// copy and paste...
////////////////////////////////////////////////////////////////////////

// C++ includes

// ROOT includes
#include "TTree.h"
#include "TH1F.h"
#include "TLorentzVector.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// DUNETPC specific includes
#include "dune/DAQTriggerSim/TriggerDataProducts/TriggerTypes.h"
#include "dune/DAQTriggerSim/TriggerDataProducts/BasicTrigger.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/RecoBase/Hit.h"


const int kMaxNumParticles=500;
//const int kMaxNumHits=50000; // probably will need to be higher


class SupernovaAna;

class SupernovaAna : public art::EDAnalyzer {

public:

  explicit SupernovaAna(fhicl::ParameterSet const & p);

  // Plugins should not be copied or assigned.
  SupernovaAna(SupernovaAna const &) = delete;
  SupernovaAna(SupernovaAna &&) = delete;
  SupernovaAna & operator = (SupernovaAna const &) = delete;
  SupernovaAna & operator = (SupernovaAna &&) = delete;

  // The main guts...
  void analyze(art::Event const & e) override;

  void reconfigure(fhicl::ParameterSet const & p);

  void beginJob();



private:

  // label for the modules that made trigger objects
  std::string fTruthLabel;
  std::string fHitLabel;

  // Define an ntuple to be filled
  TTree* fNTuple;

  // Define variables to be used as branches of the nTuple
  int NHits;
  double avgRMS;
  double avgRMSloop[kMaxNumParticles];
  double avgRMSloopsum;
  double SummedADC;
  int NMCTruths;
  double MomentumYTruths[kMaxNumParticles];
  double MomentumZTruths[kMaxNumParticles];
  double NuEnergyTruths[kMaxNumParticles];
  double LeptonEnergyTruths[kMaxNumParticles];
  double ThetaTruths[kMaxNumParticles];
  int OriginTruths[kMaxNumParticles];
  int CCNCTruths[kMaxNumParticles];
  int ModeTruths[kMaxNumParticles];
  int NuPDGCodeTruths[kMaxNumParticles];  
  int LeptonPDGCodeTruths[kMaxNumParticles];  
  double xTruths[kMaxNumParticles];
  double yTruths[kMaxNumParticles];
  double zTruths[kMaxNumParticles];
  double tTruths[kMaxNumParticles];
  double EndxTruths[kMaxNumParticles];
  double EndyTruths[kMaxNumParticles];
  double EndzTruths[kMaxNumParticles];
  double EndtTruths[kMaxNumParticles];
  
  // Define 1D histograms to be filled
  TH1I *fNHits;
  TH1F *fSummedADC;
  TH1I *fNMCTruths;
  TH1F *fMomentumYTruths;
  TH1F *fMomentumZTruths;
  TH1F *fEnergyTruths;
  TH1F *fThetaTruths;
  TH1I *fOriginTruths;
  TH1I *fCCNCTruths;
  TH1I *fModeTruths;
  TH1I *fPDGCodeTruths;
  TH1F *fxTruths;

};



//......................................................
SupernovaAna::SupernovaAna(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)
{
  this->reconfigure(p);
}



//......................................................
void SupernovaAna::reconfigure(fhicl::ParameterSet const & p)
{
  fTruthLabel = p.get<std::string> ("TruthLabel");
  fHitLabel = p.get<std::string> ("HitLabel");
}



//......................................................
void SupernovaAna::beginJob()
{

  art::ServiceHandle<art::TFileService> tfs;

  // Make nTuple
  fNTuple = tfs->make<TTree>("SupernovaAna","Supernonva analysis tree");

  // Add branches to nTuple
  fNTuple->Branch("NHits",&NHits);
  fNTuple->Branch("avgRMS",&avgRMS);
  fNTuple->Branch("SummedADC",&SummedADC);
  fNTuple->Branch("NMCTruths",&NMCTruths);
  fNTuple->Branch("MomentumYTruths",MomentumYTruths,"MomentumYTruths[NMCTruths]/D");
  fNTuple->Branch("MomentumZTruths",MomentumZTruths,"MomentumZTruths[NMCTruths]/D");
  fNTuple->Branch("NuEnergyTruths",NuEnergyTruths,"NuEnergyTruths[NMCTruths]/D");
  fNTuple->Branch("LeptonEnergyTruths",LeptonEnergyTruths,"LeptonEnergyTruths[NMCTruths]/D");
  fNTuple->Branch("ThetaTruths",ThetaTruths,"ThetaTruths[NMCTruths]/D");
  fNTuple->Branch("OriginTruths",OriginTruths,"OriginTruths[NMCTruths]/I");
  fNTuple->Branch("CCNCTruths",CCNCTruths,"CCNCTruths[NMCTruths]/I");
  fNTuple->Branch("ModeTruths",ModeTruths,"ModeTruths[NMCTruths]/I");
  fNTuple->Branch("NuPDGCodeTruths",NuPDGCodeTruths,"NuPDGCodeTruths[NMCTruths]/I");
  fNTuple->Branch("LeptonPDGCodeTruths",LeptonPDGCodeTruths,"LeptonPDGCodeTruths[NMCTruths]/I");
  fNTuple->Branch("xTruths",xTruths,"xTruths[NMCTruths]/D");
  fNTuple->Branch("yTruths",yTruths,"yTruths[NMCTruths]/D");
  fNTuple->Branch("zTruths",zTruths,"zTruths[NMCTruths]/D");
  fNTuple->Branch("tTruths",tTruths,"tTruths[NMCTruths]/D");
  fNTuple->Branch("EndxTruths",EndxTruths,"EndxTruths[NMCTruths]/D");
  fNTuple->Branch("EndyTruths",EndyTruths,"EndyTruths[NMCTruths]/D");
  fNTuple->Branch("EndzTruths",EndzTruths,"EndzTruths[NMCTruths]/D");
  fNTuple->Branch("EndtTruths",EndtTruths,"EndtTruths[NMCTruths]/D");


  // Make histograms
  fNHits = tfs->make<TH1I>("fHits","Number of hits per event;Hits;Count",20001,-0.5,20000.5);

  fSummedADC = tfs->make<TH1F>("fSummedADC","Number of hits with a specific summed ADC value",1011,-10.5,1000.5);

  fNMCTruths = tfs->make<TH1I>("fNMCTruths","Number of MCTruths per event;N;count",
			       101,-0.5,100.5);

  fMomentumYTruths = tfs->make<TH1F>("fMomentumYTruths","Num of",101,-0.5,100.5);

  fMomentumZTruths = tfs->make<TH1F>("fMomentumZTruths","Num",101,-0.5,100.5);

  fEnergyTruths = tfs->make<TH1F>("fEnergyTruths","Number of neutrinos of a certain energy;Energy;count",101,-0.5,100.5);

  fThetaTruths = tfs->make<TH1F>("fThetaTruths","Number of interactions with an angle of #theta between incoming neutrino and outgoing lepton;#theta [rad];count",40,-(TMath::Pi())/40.0,41.0*(TMath::Pi())/40.0);
  
  fOriginTruths = tfs->make<TH1I>("fOriginTruths","Number of neutrinos from each source (see Origin documentation);source;count",5,-0.5,4.5);
 
  fCCNCTruths = tfs->make<TH1I>("fCCNCTruths","Number of Charged Current and Neutral Current interactions;CC or NC;count",2,-0.5,1.5);

  fModeTruths = tfs->make<TH1I>("fModeTruths","Number of interaction types;mode;count",14,-0.5,13.5);

  fPDGCodeTruths = tfs->make<TH1I>("fPDGCodeTruths","Number of particles of a type;Particle type;count",37,-18.5,18.5);
  
  fxTruths = tfs->make<TH1F>("fxTruths","Number",100,-500,500);

}



//......................................................
void SupernovaAna::analyze(art::Event const & e)
{

  
  //  art::Handle< std::vectir<recob::Vertex > > vertex_list;


  //==============================HITS================================

  // Get the Hits out of the event
  art::Handle< std::vector< recob::Hit > > hits_list;
  e.getByLabel(fHitLabel, hits_list);

  // Gets the number of hits per event and fills into a histogram
  NHits = hits_list->size();
  fNHits->Fill(hits_list->size());
  
  // 
  //for(int i=0;i<hits_list->size();i++){
  //avgRMSloop[i] = -1.0e-9;
  
  SummedADC = 0.0;
  for(unsigned int i = 0; i < hits_list->size();i++){
    recob::Hit const& hit = hits_list->at(i);  
    avgRMSloop[i] = hit.RMS();
    avgRMSloopsum += avgRMSloop[i];
    SummedADC += hit.SummedADC();
   }
  avgRMS = avgRMSloopsum/(hits_list->size());
  fSummedADC->Fill(SummedADC);

  //============================MCTRUTHS==============================

  // Get the MCTruths out of the event
  art::Handle< std::vector< simb::MCTruth > > truths_list;
  e.getByLabel(fTruthLabel, truths_list);

  

  // Gets the number of MC truths per event and fills into a histogram
  NMCTruths = truths_list->size();
  fNMCTruths->Fill(truths_list->size());
  
  // Loop through the truth arrays and sets the values of the arrays
  // to an arbitrary default (unfilled) value
  for(int i=0;i<kMaxNumParticles;i++){
    xTruths[i] = -1.0e-9;
    yTruths[i] = -1.0e-9;
    zTruths[i] = -1.0e-9;
    tTruths[i] = -1.0e-9;
    EndxTruths[i] = -1.0e-9;
    EndyTruths[i] = -1.0e-9;
    EndzTruths[i] = -1.0e-9;
    EndtTruths[i] = -1.0e-9;
    MomentumYTruths[i] = -1.0e9;
    MomentumZTruths[i] = -1.0e9;
    NuEnergyTruths[i] = -999;
    LeptonEnergyTruths[i] = -999;
    ThetaTruths[i] = -999;
    OriginTruths[i] = -999;
    CCNCTruths[i] = -999;
    ModeTruths[i] = -999;
    NuPDGCodeTruths[i] = 0;
    LeptonPDGCodeTruths[i] = 0;

  }
  
  // Detect if the number of MC truths per event is greater than 
  // or equal to the maximum number of particles per event else 
  // we'll fall off the end of the array and cause a seg fault
  if(NMCTruths >= kMaxNumParticles) { 
    std::cerr << "ERROR: NMCTruths " << NMCTruths << 
    " >= kMaxNumParticles " << kMaxNumParticles << 
    " , will cause a segmentation fault" << std::endl;}
  
  // Create nTuple data and histograms for any number of truths 
  for(unsigned int i = 0; i < truths_list->size();i++){
    simb::MCTruth const& truth = truths_list->at(i);
    

    /*fxTruths->Fill(mc_particle.Vx());
        fyTruths->Fill(mc_particle.Vy());
    fzTruths->Fill(mc_particle.Vz());
    ftTruths->Fill(mc_particle.T());    
    */
    // Get the ith particle's Origin truths and fill a histogram
    OriginTruths[i] = truth.Origin();
    fOriginTruths->Fill(truth.Origin());
           
    // Only get neutrino truths if they exist to avoid a segmentation fault
    if (truth.NeutrinoSet()){
      simb::MCNeutrino const& mc_neutrino = truth.GetNeutrino();
      // should change notation such that mc_neutrino = truth.GetNeutrino().Nu(); 
      // and then remove .Nu() when called
      simb::MCParticle const& mc_lepton = mc_neutrino.Lepton();
      
      xTruths[i] = mc_lepton.Vx();
      yTruths[i] = mc_lepton.Vy();
      zTruths[i] = mc_lepton.Vz();
      tTruths[i] = mc_lepton.T();
      EndxTruths[i] = mc_lepton.EndX();
      EndyTruths[i] = mc_lepton.EndY();
      EndzTruths[i] = mc_lepton.EndZ();
      EndtTruths[i] = mc_lepton.EndT();

      
      // Get data to be filled into nTuple
      MomentumYTruths[i] = mc_neutrino.Lepton().Py();
      MomentumZTruths[i] = mc_neutrino.Lepton().Pz();
      NuEnergyTruths[i] = mc_neutrino.Nu().E();
      LeptonEnergyTruths[i] = mc_lepton.E();
      ThetaTruths[i] = mc_neutrino.Theta();
      CCNCTruths[i] = mc_neutrino.CCNC();
      ModeTruths[i] = mc_neutrino.Mode();
      NuPDGCodeTruths[i] = mc_neutrino.Nu().PdgCode();
      LeptonPDGCodeTruths[i] = mc_lepton.PdgCode();

      // Fill histograms
      fMomentumYTruths->Fill(mc_neutrino.Lepton().Py());
      fMomentumZTruths->Fill(mc_neutrino.Lepton().Pz());
      fEnergyTruths->Fill(mc_neutrino.Nu().E());
      fThetaTruths->Fill(mc_neutrino.Theta());
      fCCNCTruths->Fill(mc_neutrino.CCNC());
      fModeTruths->Fill(mc_neutrino.Mode());
      fPDGCodeTruths->Fill(mc_neutrino.Nu().PdgCode());
    }
  }

  // Now we should have NMCTruths = truths_list->size() and 
  // PDGCodeTruths [0 - (NMCTruths - 1)] should contain PdgCodes, 
  // [ NMCTruths - 500 ] should be equal to -999, signifying that 
  // we did not fill this part
  
  // Fill nTuple
  fNTuple->Fill();
}
DEFINE_ART_MODULE(SupernovaAna)
