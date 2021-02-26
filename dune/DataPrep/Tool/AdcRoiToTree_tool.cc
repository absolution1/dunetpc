// AdcRoiToTree_tool.cc

#include "AdcRoiToTree.h"
#include "TFile.h"
#include "TTree.h"
#include <iostream>

using std::string;
using std::cout;
using std::endl;
using fhicl::ParameterSet;

//**********************************************************************
// Class methods.
//**********************************************************************

AdcRoiToTree::AdcRoiToTree(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")),
  m_OutFile(ps.get<Name>("OutFile"))
{
  const string myname = "AdcRoiToTree::ctor: ";
  if ( m_LogLevel >=2 ) cout << myname << "Creating output file." << endl;
  TFile* pfil = TFile::Open(m_OutFile.c_str(), "RECREATE");
  if ( pfil == nullptr || ! pfil->IsOpen() ) {
    cout << myname << "ERROR: Unable to create output file " << m_OutFile << endl;
  } else {
    TreeData tdat;
    TTree* ptre = new TTree(treeName().c_str(), "ADC ROIs");
    ptre->Branch("run",     &tdat.run);
    ptre->Branch("event",   &tdat.event);
    ptre->Branch("channel", &tdat.channel);
    ptre->Branch("nroi",    &tdat.nroi);
    ptre->Branch("nsam",    &tdat.nsam[0], "nsam[nroi]/i");
    ptre->Branch("isam",    &tdat.isam[0], "isam[nroi]/i");
    ptre->Branch("qke",     &tdat.qke[0],  "qke[nroi]/F");
    ptre->ResetBranchAddresses();
    ptre->Write();
    pfil->Close();
  }
  delete pfil;
  // Display the configuration.
  if ( m_LogLevel>= 1 ) {
    cout << myname << "             LogLevel: " << m_LogLevel << endl;
    cout << myname << "              OutFile: " << m_OutFile << endl;
  }
}

//**********************************************************************

AdcRoiToTree::~AdcRoiToTree() {
  const string myname = "AdcRoiToTree::dtor: ";
  if ( m_LogLevel >= 1 ) {
    cout << myname << "Exiting." << endl;
    TFile* pfil = TFile::Open(m_OutFile.c_str(), "READ");
    if ( pfil == nullptr || ! pfil->IsOpen() ) {
      cout << myname << "ERROR: Unable to open output file " << m_OutFile << endl;
    } else {
      TTree* ptre = dynamic_cast<TTree*>(pfil->Get(treeName().c_str()));
      if ( ptre == nullptr ) {
        cout << myname << "ERROR: Unable to open tree " << treeName() << endl;
      } else {
        cout << myname << "Tree entry count is " << ptre->GetEntries() << endl;
      }
    }
    delete pfil;
  }
}

//**********************************************************************

DataMap AdcRoiToTree::viewMap(const AdcChannelDataMap& acds) const {
  const string myname = "AdcRoiToTree::viewMap: ";
  DataMap ret;
  TFile* pfil = TFile::Open(m_OutFile.c_str(), "UPDATE");
  if ( pfil == nullptr || ! pfil->IsOpen() ) {
    cout << myname << "ERROR: Unable to open output file " << m_OutFile << endl;
    return ret.setStatus(1);
  }
  TTree* ptre = dynamic_cast<TTree*>(pfil->Get(treeName().c_str()));
  if ( ptre == nullptr ) {
    cout << myname << "ERROR: Unable to open tree " << treeName() << endl;
    return ret.setStatus(2);
  }
  Index nfill = 0;
  Index maxroi = 1000;
  TreeData tdat;
  tdat.nsam.resize(maxroi);
  tdat.isam.resize(maxroi);
  tdat.qke.resize(maxroi);
  ptre->SetBranchAddress("run",     &tdat.run);
  ptre->SetBranchAddress("event",   &tdat.event);
  ptre->SetBranchAddress("channel", &tdat.channel);
  ptre->SetBranchAddress("nroi",    &tdat.nroi);
  ptre->SetBranchAddress("nsam",    &tdat.nsam[0]);
  ptre->SetBranchAddress("isam",    &tdat.isam[0]);
  ptre->SetBranchAddress("qke",     &tdat.qke[0]);
  //ptre->SetBranchAddress("nsample", pflt, "nsample[nroi]/F");
  for ( const auto& iacd : acds ) {
    const AdcChannelData& acd = iacd.second;
    tdat.run = acd.run();
    tdat.event = acd.event();
    tdat.channel = acd.channel();
    tdat.nroi = acd.rois.size();
    Index nroi = std::min(tdat.nroi, maxroi);
    for ( Index iroi=0; iroi<nroi; ++iroi ) {
      
      Index isam1 = acd.rois[iroi].first;
      Index isam2 = acd.rois[iroi].second;
      float q = 0.0;
      if ( acd.samples.size() > isam2 ) {
        for ( Index isam=isam1; isam<=isam2; ++isam ) q += acd.samples[isam];
      } else {
        cout << myname << "WARNING: Ignoring missing samples for run " << acd.run()
             << ", event " << acd.event() << ", channel " << acd.channel() << endl;
      }
      tdat.nsam[iroi] = 1 + isam2 - isam1;
      tdat.isam[iroi] = isam1;
      tdat.qke[iroi] = q;
    }
    if ( m_LogLevel >= 3 ) {
      cout << myname << "Filling run " << tdat.run << ", event " << tdat.event
           << ", channel " << tdat.channel << ", nroi " << tdat.nroi << endl;
      if ( m_LogLevel >= 4 ) {
        string spre = myname + "  nsam: [";
        for ( Index iroi=0; iroi<nroi; ++iroi ) {
          cout << spre;
          cout << tdat.nsam[iroi];
          spre = ", ";
        }
        cout << "]" << endl;
        spre = myname + "  isam: [";
        for ( Index iroi=0; iroi<nroi; ++iroi ) {
          cout << spre;
          cout << tdat.isam[iroi];
          spre = ", ";
        }
        cout << "]" << endl;
        spre = myname + "   qke: [";
        for ( Index iroi=0; iroi<nroi; ++iroi ) {
          cout << spre;
          cout << tdat.qke[iroi];
          spre = ", ";
        }
        cout << "]" << endl;
      }
    }
    ptre->Fill();
  }
  ptre->ResetBranchAddresses();
  ptre->Write();
  pfil->Close();
  delete pfil;
  ret.setInt("art_nfill", nfill);
  return ret;
}

//**********************************************************************
