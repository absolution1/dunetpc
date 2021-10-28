////////////////////////////////////////////////////////////////////////
//
// \file CAFMaker_module.cc
//
// Chris Marshall's version
// Largely based on historical FDSensOpt/CAFMaker_module.cc
//
///////////////////////////////////////////////////////////////////////

#ifndef CAFMaker_H
#define CAFMaker_H

// Generic C++ includes
#include <iostream>

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h" 
#include "art/Framework/Principal/Event.h" 
#include "art/Framework/Principal/SubRun.h"
#include "fhiclcpp/ParameterSet.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "art_root_io/TFileService.h"

//#include "Utils/AppInit.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "dune/FDSensOpt/FDSensOptData/MVASelectPID.h"
#include "dune/FDSensOpt/FDSensOptData/EnergyRecoOutput.h"
#include "dune/CVN/func/InteractionType.h"
#include "dune/CVN/func/Result.h"
#include "dune/RegCNN/func/RegCNNResult.h"

// dunerw stuff
#include "systematicstools/interface/ISystProviderTool.hh"
#include "systematicstools/utility/ParameterAndProviderConfigurationUtility.hh"
#include "systematicstools/utility/exceptions.hh"
//#include "systematicstools/utility/md5.hh"

// root
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"

// pdg
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"

// genie
#include "Framework/EventGen/EventRecord.h"
#include "Framework/GHEP/GHepParticle.h"

constexpr int knShifts = 100; // number of shifts
constexpr int kmaxRwgts = 100; // Largest number of reweights in a shift

namespace dunemva {

  class CAFMaker : public art::EDAnalyzer {

    public:

      explicit CAFMaker(fhicl::ParameterSet const& pset);
      virtual ~CAFMaker();
      void beginJob() override;
      void endJob() override;
      void beginSubRun(const art::SubRun& sr) override;
      void endSubRun(const art::SubRun& sr) override;
      void reconfigure(fhicl::ParameterSet const& pset) /*override*/;
      void analyze(art::Event const & evt) override;


    private:
      std::string fMVASelectLabel;
      std::string fMVASelectNueLabel;
      std::string fMVASelectNumuLabel;

      std::string fCVNLabel;
      std::string fRegCNNLabel;

      std::string fEnergyRecoNueLabel;
      std::string fEnergyRecoNumuLabel;
      std::string fMVAMethod;

      float fOscPro;
      double fWeight;
      TTree* fTree;  
      TTree* fMetaTree;

      // Get reweight knobs from fhicl file -- no hard-coded shifts
      int fNwgt[knShifts];
      double fCvWgts[knShifts];
      double fWgts[knShifts][kmaxRwgts];

      // CAF variables
      // configuration variables
      int fIsFD, fIsFHC;
      // event accounting
      int fRun, fSubrun, fEvent;
      // Truth information
      int fIsCC, fNuPDG, fNuPDGunosc, fMode, fLepPDG; 
      double fEv, fQ2, fW, fX, fY, fNuMomX, fNuMomY, fNuMomZ, fLepMomX, fLepMomY, fLepMomZ, fLepE, fLepNuAngle;
      // True particle counts
      int nP, nN, nPip, nPim, nPi0, nKp, nKm, nK0, nEM, nOtherHad, nNucleus, nUNKNOWN;
      double eP, eN, ePip, ePim, ePi0, eOther;
      double vtx_x, vtx_y, vtx_z;

      // Reco information
      double fErecoNue;
      double fRecoLepEnNue;
      double fRecoHadEnNue;
      int fRecoMethodNue; // 1 = longest reco track + hadronic, 2 = reco shower with highest charge + hadronic, 3 = all hit charges, -1 = not set
      double fErecoNumu; 
      double fRecoLepEnNumu;
      double fRecoHadEnNumu;
      int fRecoMethodNumu; // 1 = longest reco track + hadronic, 2 = reco shower with highest charge + hadronic, 3 = all hit charges, -1 = not set
      int fLongestTrackContNumu; // 1 = contained, 0 = exiting, -1 = not set
      int fTrackMomMethodNumu; // 1 = range, 0 = MCS, -1 = not set

      double fMVAResult;
      double fMVAResultNue;
      double fMVAResultNumu;

      // CVN outputs
      double fCVNResultIsAntineutrino;
      double fCVNResultNue, fCVNResultNumu, fCVNResultNutau, fCVNResultNC; // flavour
      double fCVNResult0Protons, fCVNResult1Protons, fCVNResult2Protons, fCVNResultNProtons; // #protons
      double fCVNResult0Pions, fCVNResult1Pions, fCVNResult2Pions, fCVNResultNPions; // #pions
      double fCVNResult0Pizeros, fCVNResult1Pizeros, fCVNResult2Pizeros, fCVNResultNPizeros; // #pizeros
      double fCVNResult0Neutrons, fCVNResult1Neutrons, fCVNResult2Neutrons, fCVNResultNNeutrons; // #neutrons

      double fRegCNNNueE;

      double meta_pot;
      int meta_run, meta_subrun, meta_version;

      systtools::provider_list_t fSystProviders;

  }; // class CAFMaker


  //------------------------------------------------------------------------------
  CAFMaker::CAFMaker(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset)
  {
    this->reconfigure(pset);
  }

  dunemva::CAFMaker::~CAFMaker(){}

  //------------------------------------------------------------------------------
  void CAFMaker::reconfigure(fhicl::ParameterSet const& pset) 
  {

    fMVASelectLabel = pset.get<std::string>("MVASelectLabel");
    fMVASelectNueLabel = pset.get<std::string>("MVASelectNueLabel");
    fMVASelectNumuLabel = pset.get<std::string>("MVASelectNumuLabel");
    fCVNLabel = pset.get<std::string>("CVNLabel");
    fRegCNNLabel = pset.get<std::string>("RegCNNLabel");

    fEnergyRecoNueLabel = pset.get<std::string>("EnergyRecoNueLabel");
    fEnergyRecoNumuLabel = pset.get<std::string>("EnergyRecoNumuLabel");

    // Get DUNErw stuff from its fhicl, which should be included on the CAFMaker config file
    //if( !pset.has_key("generated_systematic_provider_configuration") ) {
    //  std::cout << "[ERROR]: Could not find producer key: "
    //               "\"generated_systematic_provider_configuration\". This should "
    //               "contain a list of configured systematic providers generated by "
    //               "GenerateSystProviderConfig." << std::endl;
    //  return;
    //}

    fhicl::ParameterSet syst_provider_config = pset.get<fhicl::ParameterSet>("generated_systematic_provider_configuration");

    fSystProviders = systtools::ConfigureISystProvidersFromParameterHeaders(syst_provider_config);
  }


  //------------------------------------------------------------------------------
  void CAFMaker::beginJob()
  {

    art::ServiceHandle<art::TFileService> tfs;
    fTree = tfs->make<TTree>("caf", "caf");
    fMetaTree = tfs->make<TTree>("meta", "meta");

    // book-keeping
    fTree->Branch("run",         &fRun,        "run/I");
    fTree->Branch("subrun",      &fSubrun,     "subrun/I");
    fTree->Branch("event",       &fEvent,      "event/I");
    fTree->Branch("isFD",        &fIsFD,       "isFD/I");
    fTree->Branch("isFHC",       &fIsFHC,      "isFHC/I");
    fTree->Branch("isCC",        &fIsCC,       "isCC/I");

    // true interaction quantities
    fTree->Branch("nuPDG",        &fNuPDG,        "nuPDG/I");
    fTree->Branch("nuPDGunosc",   &fNuPDGunosc,   "nuPDGunosc/I");
    fTree->Branch("NuMomX",       &fNuMomX,       "NuMomX/D");
    fTree->Branch("NuMomY",       &fNuMomY,       "NuMomY/D");
    fTree->Branch("NuMomZ",       &fNuMomZ,       "NuMomZ/D");
    fTree->Branch("Ev",           &fEv,           "Ev/D");
    fTree->Branch("mode",         &fMode,         "mode/I");
    fTree->Branch("LepPDG",       &fLepPDG,       "LepPDG/I");
    fTree->Branch("LepMomX",      &fLepMomX,      "LepMomX/D");
    fTree->Branch("LepMomY",      &fLepMomY,      "LepMomY/D");
    fTree->Branch("LepMomZ",      &fLepMomZ,      "LepMomZ/D");
    fTree->Branch("LepE",         &fLepE,         "LepE/D");
    fTree->Branch("LepNuAngle",   &fLepNuAngle,   "LepNuAngle/D");
    fTree->Branch("Q2",           &fQ2,           "Q2/D");
    fTree->Branch("W",            &fW,            "W/D");
    fTree->Branch("X",            &fX,            "X/D");
    fTree->Branch("Y",            &fY,            "Y/D");

    // FS particle counts
    fTree->Branch("nP",        &nP,         "nP/I");
    fTree->Branch("nN",        &nN,         "nN/I");
    fTree->Branch("nipip",     &nPip,       "nipip/I");
    fTree->Branch("nipim",     &nPim,       "nipim/I");
    fTree->Branch("nipi0",     &nPi0,       "nipi0/I");
    fTree->Branch("nikp",      &nKp,        "nikp/I");
    fTree->Branch("nikm",      &nKm,        "nikm/I");
    fTree->Branch("nik0",      &nK0,        "nik0/I");
    fTree->Branch("niem",      &nEM,        "niem/I");
    fTree->Branch("niother",   &nOtherHad,  "niother/I");
    fTree->Branch("nNucleus",  &nNucleus,   "nNucleus/I");
    fTree->Branch("nUNKNOWN",  &nUNKNOWN,   "nUNKNOWN/I");
    fTree->Branch("eP",        &eP,         "eP/D");
    fTree->Branch("eN",        &eN,         "eN/D");
    fTree->Branch("ePip",      &ePip,       "ePip/D");
    fTree->Branch("ePim",      &ePim,       "ePim/D");
    fTree->Branch("ePi0",      &ePi0,       "ePi0/D");
    fTree->Branch("eOther",    &eOther,     "eOther/D");

    // vertex position
    fTree->Branch("vtx_x",   &vtx_x,    "vtx_x/D");
    fTree->Branch("vtx_y",   &vtx_y,    "vtx_y/D");
    fTree->Branch("vtx_z",   &vtx_z,    "vtx_z/D");

    // Reco variables
    fTree->Branch("mvaresult",   &fMVAResult,  "mvaresult/D");
    fTree->Branch("mvanue",      &fMVAResultNue,  "mvanue/D");
    fTree->Branch("mvanumu",     &fMVAResultNumu, "mvanumu/D");

    fTree->Branch("cvnisantineutrino", &fCVNResultIsAntineutrino, "cvnisantineutrino/D");
    fTree->Branch("cvnnue",            &fCVNResultNue,            "cvnnue/D");
    fTree->Branch("cvnnumu",           &fCVNResultNumu,           "cvnnumu/D");
    fTree->Branch("cvnnutau",          &fCVNResultNutau,          "cvnnutau/D");
    fTree->Branch("cvnnc",             &fCVNResultNC,             "cvnnc/D");
    fTree->Branch("cvn0protons",       &fCVNResult0Protons,       "cvn0protons/D");
    fTree->Branch("cvn1protons",       &fCVNResult1Protons,       "cvn1protons/D");
    fTree->Branch("cvn2protons",       &fCVNResult2Protons,       "cvn2protons/D");
    fTree->Branch("cvnNprotons",       &fCVNResultNProtons,       "cvnNprotons/D");
    fTree->Branch("cvn0pions",         &fCVNResult0Pions,         "cvn0pions/D");
    fTree->Branch("cvn1pions",         &fCVNResult1Pions,         "cvn1pions/D");
    fTree->Branch("cvn2pions",         &fCVNResult2Pions,         "cvn2pions/D");
    fTree->Branch("cvnNpions",         &fCVNResultNPions,         "cvnNpions/D");
    fTree->Branch("cvn0pizeros",       &fCVNResult0Pizeros,       "cvn0pizeros/D");
    fTree->Branch("cvn1pizeros",       &fCVNResult1Pizeros,       "cvn1pizeros/D");
    fTree->Branch("cvn2pizeros",       &fCVNResult2Pizeros,       "cvn2pizeros/D");
    fTree->Branch("cvnNpizeros",       &fCVNResultNPizeros,       "cvnNpizeros/D");
    fTree->Branch("cvn0neutrons",      &fCVNResult0Neutrons,      "cvn0neutrons/D");
    fTree->Branch("cvn1neutrons",      &fCVNResult1Neutrons,      "cvn1neutrons/D");
    fTree->Branch("cvn2neutrons",      &fCVNResult2Neutrons,      "cvn2neutrons/D");
    fTree->Branch("cvnNneutrons",      &fCVNResultNNeutrons,      "cvnNneutrons/D");

    fTree->Branch("RegCNNNueE",  &fRegCNNNueE,   "RegCNNNueE/D");
    fTree->Branch("weight",      &fWeight,     "weight/D");
    fTree->Branch("oscpro",      &fOscPro,     "oscpro/F");

    fTree->Branch("Ev_reco_nue",      &fErecoNue,        "Ev_reco_nue/D");
    fTree->Branch("RecoLepEnNue",     &fRecoLepEnNue,    "RecoLepEnNue/D");
    fTree->Branch("RecoHadEnNue",     &fRecoHadEnNue,    "RecoHadEnNue/D");
    fTree->Branch("RecoMethodNue",    &fRecoMethodNue,   "RecoMethodNue/I");
    fTree->Branch("Ev_reco_numu",     &fErecoNumu,       "Ev_reco_numu/D");
    fTree->Branch("RecoLepEnNumu",    &fRecoLepEnNumu,   "RecoLepEnNumu/D");
    fTree->Branch("RecoHadEnNumu",    &fRecoHadEnNumu,   "RecoHadEnNumu/D");
    fTree->Branch("RecoMethodNumu",   &fRecoMethodNumu,  "RecoMethodNumu/I");
    fTree->Branch("LongestTrackContNumu",  &fLongestTrackContNumu, "LongestTrackContNumu/I");
    fTree->Branch("TrackMomMethodNumu",    &fTrackMomMethodNumu,   "TrackMomMethodNumu/I");

    fMetaTree->Branch("pot", &meta_pot, "pot/D");
    fMetaTree->Branch("run", &meta_run, "run/I");
    fMetaTree->Branch("subrun", &meta_subrun, "subrun/I");
    fMetaTree->Branch("version", &meta_version, "version/I");

    // make DUNErw variables
    for( auto &sp : fSystProviders ) {
      systtools::SystMetaData metaData = sp->GetSystMetaData();
      for( systtools::SystMetaData::iterator itMeta = metaData.begin(); itMeta != metaData.end(); ++itMeta ) {      
        systtools::SystParamHeader head = *itMeta;
        std::string name = head.prettyName;
        unsigned int parId = head.systParamId;
        std::cout << "Adding reweight branch " << parId << " for " << name << " with " << head.paramVariations.size() << " shifts" << std::endl;
        fTree->Branch( Form("%s_nshifts", name.c_str()), &fNwgt[parId], Form("%s_nshifts/I", name.c_str()) );
        fTree->Branch( Form("%s_cvwgt", name.c_str()), &fCvWgts[parId], Form("%s_cvwgt/D", name.c_str()) );
        fTree->Branch( Form("wgt_%s", name.c_str()), fWgts[parId], Form("wgt_%s[%s_nshifts]/D", name.c_str(), name.c_str()) );
      }
    }

    // initialize weight variables -- some won't ever be set
    for( int i = 0; i < knShifts; ++i ) {
      fNwgt[i] = 0;
      fCvWgts[i] = 1.;
      for( int j = 0; j < kmaxRwgts; ++j ) {
        fWgts[i][j] = 0.;
      }
    }

    meta_pot = 0.;
    meta_version = 1;
  }

  //------------------------------------------------------------------------------
  void CAFMaker::beginSubRun(const art::SubRun& sr)
  {
    auto pots = sr.getHandle< sumdata::POTSummary >("generator");
    if( pots ) meta_pot += pots->totpot;
  }

  //------------------------------------------------------------------------------
  void CAFMaker::analyze(art::Event const & evt)
  {
    auto pidin = evt.getHandle<dunemva::MVASelectPID>(fMVASelectLabel);
    auto pidinnue = evt.getHandle<dunemva::MVASelectPID>(fMVASelectNueLabel);
    auto pidinnumu = evt.getHandle<dunemva::MVASelectPID>(fMVASelectNumuLabel);
    art::InputTag itag1(fCVNLabel, "cvnresult");
    auto cvnin = evt.getHandle<std::vector<cvn::Result>>(itag1);
    art::InputTag itag2(fRegCNNLabel, "regcnnresult");
    auto regcnnin = evt.getHandle<std::vector<cnn::RegCNNResult>>(itag2);
    auto ereconuein = evt.getHandle<dune::EnergyRecoOutput>(fEnergyRecoNueLabel);
    auto ereconumuin = evt.getHandle<dune::EnergyRecoOutput>(fEnergyRecoNumuLabel);

    fRun = evt.id().run();
    fSubrun = evt.id().subRun();
    fEvent = evt.id().event();
    meta_run = fRun;
    meta_subrun = fSubrun;

    if( !pidin.failedToGet() ) {
      fMVAResult = pidin->pid;

      //Fill MVA reco stuff
      fErecoNue          = ereconuein->fNuLorentzVector.E();
      fRecoLepEnNue      = ereconuein->fLepLorentzVector.E();
      fRecoHadEnNue      = ereconuein->fHadLorentzVector.E();
      fRecoMethodNue     = ereconuein->recoMethodUsed;
      fErecoNumu         = ereconumuin->fNuLorentzVector.E();
      fRecoLepEnNumu     = ereconumuin->fLepLorentzVector.E();
      fRecoHadEnNumu     = ereconumuin->fHadLorentzVector.E();
      fRecoMethodNumu    = ereconumuin->recoMethodUsed;
      fLongestTrackContNumu  = ereconumuin->longestTrackContained;
      fTrackMomMethodNumu    = ereconumuin->trackMomMethod;
    }

    if( !pidinnue.failedToGet() ) {
      fMVAResultNue = pidinnue->pid;
    }

    if( !pidinnumu.failedToGet() ) {
      fMVAResultNumu = pidinnumu->pid;
    }

    if( !cvnin.failedToGet() ) {
      //using i = cvn::Interaction;
      //if(cvnin->empty() || (*cvnin)[0].fOutput.size() <= i::kNutauOther){
      if(cvnin->empty()){
        fCVNResultIsAntineutrino = fCVNResultNue = fCVNResultNumu = fCVNResultNutau = fCVNResultNC = \
        fCVNResult0Protons = fCVNResult1Protons = fCVNResult2Protons = fCVNResultNProtons = \
        fCVNResult0Pions = fCVNResult1Pions = fCVNResult2Pions = fCVNResultNPions = \
        fCVNResult0Pizeros = fCVNResult1Pizeros = fCVNResult2Pizeros = fCVNResultNPizeros = \
        fCVNResult0Neutrons = fCVNResult1Neutrons = fCVNResult2Neutrons = fCVNResultNNeutrons = -3;
      }
      else{
        //const std::vector<float>& v = (*cvnin)[0].fOutput;
        //fCVNResultNue = v[i::kNueQE] + v[i::kNueRes] + v[i::kNueDIS] + v[i::kNueOther];
        //fCVNResultNumu = v[i::kNumuQE] + v[i::kNumuRes] + v[i::kNumuDIS] + v[i::kNumuOther];
        //fCVNResultNutau = v[i::kNutauQE] + v[i::kNutauRes] + v[i::kNutauDIS] + v[i::kNutauOther]

        fCVNResultIsAntineutrino = (*cvnin)[0].GetIsAntineutrinoProbability();

        fCVNResultNue = (*cvnin)[0].GetNueProbability();
        fCVNResultNumu = (*cvnin)[0].GetNumuProbability();
        fCVNResultNutau = (*cvnin)[0].GetNutauProbability();
        fCVNResultNC = (*cvnin)[0].GetNCProbability();

        fCVNResult0Protons = (*cvnin)[0].Get0protonsProbability();
        fCVNResult1Protons = (*cvnin)[0].Get1protonsProbability();
        fCVNResult2Protons = (*cvnin)[0].Get2protonsProbability();
        fCVNResultNProtons = (*cvnin)[0].GetNprotonsProbability();

        fCVNResult0Pions = (*cvnin)[0].Get0pionsProbability();
        fCVNResult1Pions = (*cvnin)[0].Get1pionsProbability();
        fCVNResult2Pions = (*cvnin)[0].Get2pionsProbability();
        fCVNResultNPions = (*cvnin)[0].GetNpionsProbability();

        fCVNResult0Pizeros = (*cvnin)[0].Get0pizerosProbability();
        fCVNResult1Pizeros = (*cvnin)[0].Get1pizerosProbability();
        fCVNResult2Pizeros = (*cvnin)[0].Get2pizerosProbability();
        fCVNResultNPizeros = (*cvnin)[0].GetNpizerosProbability();

        fCVNResult0Neutrons = (*cvnin)[0].Get0neutronsProbability();
        fCVNResult1Neutrons = (*cvnin)[0].Get1neutronsProbability();
        fCVNResult2Neutrons = (*cvnin)[0].Get2neutronsProbability();
        fCVNResultNNeutrons = (*cvnin)[0].GetNneutronsProbability();
      }
    }

    fRegCNNNueE = -1.;  // initializing
    if(!regcnnin.failedToGet()){
      if (!regcnnin->empty()){
        const std::vector<float>& v = (*regcnnin)[0].fOutput;
        fRegCNNNueE = v[0];
      }
    }

    std::vector< art::Ptr<simb::MCTruth> > truth;
    auto mct = evt.getHandle< std::vector<simb::MCTruth> >("generator");
    if ( mct )
      art::fill_ptr_vector(truth, mct);
    else
      mf::LogWarning("CAFMaker") << "No MCTruth.";

    std::vector< art::Ptr<simb::MCFlux> > flux;
    auto mcf = evt.getHandle< std::vector<simb::MCFlux> >("generator");
    if ( mcf )
      art::fill_ptr_vector(flux, mcf);
    else
      mf::LogWarning("CAFMaker") << "No MCFlux.";
/*
    std::vector< art::Ptr<simb::GTruth> > gtru;
    auto gt = evt.getHandle< std::vector<simb::GTruth> >("generator");
    if ( gt )
      art::fill_ptr_vector(gtru, gt);
    else
      mf::LogWarning("CAFMaker") << "No GTruth.";
*/

    for(size_t i=0; i<truth.size(); i++){

      if(i>1){
        mf::LogWarning("CAFMaker") << "Skipping MC truth index " << i;
        continue;
      }

      fIsFD     = 1; // always FD
      fIsFHC    = 999; // don't know how to get this?
      fIsCC     = !(truth[i]->GetNeutrino().CCNC());  // ccnc is 0=CC 1=NC
      fNuPDG    = truth[i]->GetNeutrino().Nu().PdgCode();
      fNuPDGunosc = flux[i]->fntype;
      fMode     = truth[i]->GetNeutrino().Mode(); //0=QE/El, 1=RES, 2=DIS, 3=Coherent production; this is different than mode in ND
      fEv       = truth[i]->GetNeutrino().Nu().E();
      fQ2       = truth[i]->GetNeutrino().QSqr();
      fW        = truth[i]->GetNeutrino().W();
      fX        = truth[i]->GetNeutrino().X();
      fY        = truth[i]->GetNeutrino().Y();
      fNuMomX   = truth[i]->GetNeutrino().Nu().Momentum().X();
      fNuMomY   = truth[i]->GetNeutrino().Nu().Momentum().Y();
      fNuMomZ   = truth[i]->GetNeutrino().Nu().Momentum().Z();

      vtx_x     = truth[i]->GetNeutrino().Lepton().Vx();
      vtx_y     = truth[i]->GetNeutrino().Lepton().Vy();
      vtx_z     = truth[i]->GetNeutrino().Lepton().Vz();

      //Lepton stuff
      fLepPDG     = truth[i]->GetNeutrino().Lepton().PdgCode();
      fLepMomX    = truth[i]->GetNeutrino().Lepton().Momentum().X();
      fLepMomY    = truth[i]->GetNeutrino().Lepton().Momentum().Y();
      fLepMomZ    = truth[i]->GetNeutrino().Lepton().Momentum().Z();
      fLepE       = truth[i]->GetNeutrino().Lepton().Momentum().T();
      fLepNuAngle = truth[i]->GetNeutrino().Nu().Momentum().Vect().Angle(truth[i]->GetNeutrino().Lepton().Momentum().Vect());

      // TODO
      fWeight = 0.;
      fOscPro = 0.;
      //      fOscPro   = fMVAAlg.OscPro(fCCNC,fBeamPdg,fNuPDG,fEtrue);

      nP        = 0;
      nN        = 0;
      nPip      = 0;
      nPim      = 0;
      nPi0      = 0;
      nKp       = 0;
      nKm       = 0;
      nK0       = 0;
      nEM       = 0;
      nOtherHad = 0;
      nNucleus  = 0;
      nUNKNOWN  = 0;

      eP = 0.;
      eN = 0.;
      ePip = 0.;
      ePim = 0.;
      ePi0 = 0.;
      eOther = 0.;

      for( int p = 0; p < truth[i]->NParticles(); p++ ) {
        if( truth[i]->GetParticle(p).StatusCode() == genie::kIStHadronInTheNucleus ) {

          int pdg = truth[i]->GetParticle(p).PdgCode();
          double ke = truth[i]->GetParticle(p).E() - truth[i]->GetParticle(p).Mass();
          if     ( pdg == genie::kPdgProton ) {
            nP++;
            eP += ke;
          } else if( pdg == genie::kPdgNeutron ) {
            nN++;
            eN += ke;
          } else if( pdg == genie::kPdgPiP ) {
            nPip++;
            ePip += ke;
          } else if( pdg == genie::kPdgPiM ) {
            nPim++;
            ePim += ke;
          } else if( pdg == genie::kPdgPi0 ) {
            nPi0++;
            ePi0 += ke;
          } else if( pdg == genie::kPdgKP ) {
            nKp++;
            eOther += ke;
          } else if( pdg == genie::kPdgKM ) {
            nKm++;
            eOther += ke;
          } else if( pdg == genie::kPdgK0 || pdg == genie::kPdgAntiK0 || pdg == genie::kPdgK0L || pdg == genie::kPdgK0S ) {
            nK0++;
            eOther += ke;
          } else if( pdg == genie::kPdgGamma ) {
            nEM++;
            eOther += ke;
          } else if( genie::pdg::IsHadron(pdg) ) {
            nOtherHad++; // charm mesons, strange and charm baryons, antibaryons, etc.
            eOther += ke;
          } else if( genie::pdg::IsIon(pdg) ) {
            nNucleus++;
          } else {
            nUNKNOWN++;
          }

        }
      }

      // Reweighting variables
      //systtools::ScrubUnityEventResponses(er);

      // struct ParamResponses { 
      //   paramId_t pid;
      //   std::vector<double> responses;
      // }
      // typedef std::vector<ParamResponses> event_unit_response_t;
      // typedef std::vector<event_unit_response_t> EventResponse;

      for( auto &sp : fSystProviders ) {
        std::unique_ptr<systtools::EventResponse> syst_resp = sp->GetEventResponse(evt);
        if( !syst_resp ) {
          std::cout << "[ERROR]: Got nullptr systtools::EventResponse from provider "
                    << sp->GetFullyQualifiedName();
          continue;
        }

        for( systtools::EventResponse::iterator itResp = syst_resp->begin(); itResp != syst_resp->end(); ++itResp ) {
          systtools::event_unit_response_t resp = *itResp;
          for( systtools::event_unit_response_t::iterator it = resp.begin(); it != resp.end(); ++it ) {
            fNwgt[(*it).pid] = (*it).responses.size();
            //fCvWgts[(*it).pid] = (*it).CV_weight;
            for( unsigned int i = 0; i < (*it).responses.size(); ++i ) {
              fWgts[(*it).pid][i] = (*it).responses[i];
            }
          }
        }
      }
    } // loop through MC truth i

    fTree->Fill();
    return;
  }

  //------------------------------------------------------------------------------

  //------------------------------------------------------------------------------
  void CAFMaker::endSubRun(const art::SubRun& sr){
  }

  void CAFMaker::endJob()
  {
    fMetaTree->Fill();
  }

  DEFINE_ART_MODULE(CAFMaker)

} // namespace dunemva

#endif // CAFMaker_H
