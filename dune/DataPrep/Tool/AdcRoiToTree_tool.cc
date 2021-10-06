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
  m_OutFile(ps.get<Name>("OutFile")),
  m_MetadataFields(ps.get<NameVector>("MetadataFields"))
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
    ptre->Branch("status",  &tdat.status);
    tdat.mdata.resize(m_MetadataFields.size());
    float* pmd = &tdat.mdata[0];
    for ( Name mnam : m_MetadataFields ) {
      ptre->Branch(mnam.c_str(), pmd++);
    }
    ptre->Branch("nroi",    &tdat.nroi);
    ptre->Branch("nsam",    &tdat.nsam[0], "nsam[nroi]/i");
    ptre->Branch("isam",    &tdat.isam[0], "isam[nroi]/i");
    ptre->Branch("qroi",    &tdat.qroi[0],  "qroi[nroi]/F");
    ptre->Branch("hmin",    &tdat.hmin[0],  "hmin[nroi]/F");
    ptre->Branch("hmax",    &tdat.hmax[0],  "hmax[nroi]/F");
    ptre->ResetBranchAddresses();
    ptre->Write();
    pfil->Close();
  }
  delete pfil;
  // Display the configuration.
  if ( m_LogLevel>= 1 ) {
    cout << myname << "             LogLevel: " << m_LogLevel << endl;
    cout << myname << "              OutFile: " << m_OutFile << endl;
    cout << myname << "       MetadataFields: [";
    bool first = true;
    for ( Name mnam : m_MetadataFields ) {
      if ( first ) first = false;
      else cout << ", ";
      cout << mnam;
    }
    cout << "]" << endl;
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
  tdat.qroi.resize(maxroi);
  tdat.hmin.resize(maxroi);
  tdat.hmax.resize(maxroi);
  ptre->SetBranchAddress("run",     &tdat.run);
  ptre->SetBranchAddress("event",   &tdat.event);
  ptre->SetBranchAddress("channel", &tdat.channel);
  ptre->SetBranchAddress("status",  &tdat.status);
  tdat.mdata.resize(m_MetadataFields.size());
  float* pmd = &tdat.mdata[0];
  for ( Name mnam : m_MetadataFields ) {
    ptre->SetBranchAddress(mnam.c_str(), pmd++);
  }
  ptre->SetBranchAddress("nroi",    &tdat.nroi);
  ptre->SetBranchAddress("nsam",    &tdat.nsam[0]);
  ptre->SetBranchAddress("isam",    &tdat.isam[0]);
  ptre->SetBranchAddress("qroi",    &tdat.qroi[0]);
  ptre->SetBranchAddress("hmin",    &tdat.hmin[0]);
  ptre->SetBranchAddress("hmax",    &tdat.hmax[0]);
  //ptre->SetBranchAddress("nsample", pflt, "nsample[nroi]/F");
  for ( const auto& iacd : acds ) {
    const AdcChannelData& acd = iacd.second;
    tdat.run = acd.run();
    tdat.event = acd.event();
    tdat.channel = acd.channel();
    tdat.status = acd.channelStatus();
    Index idat = 0;
    for ( Name mnam : m_MetadataFields ) {
      if ( ! acd.hasMetadata(mnam) ) {
        cout << myname << "WARNING: Run/event/channel "
             << acd.run() << "/" << acd.event() << "/" << acd.channel()
             << " does not have metadata field " << mnam << endl;
      }
      tdat.mdata[idat++] = acd.getMetadata(mnam, 0.0);
    }
    tdat.nroi = acd.rois.size();
    Index nroi = std::min(tdat.nroi, maxroi);
    for ( Index iroi=0; iroi<nroi; ++iroi ) {
      
      Index isam1 = acd.rois[iroi].first;
      Index isam2 = acd.rois[iroi].second;
      float qroi = 0.0;
      float hmin = 0.0;
      float hmax = 0.0;
      bool haveSamples = false;
      if ( acd.samples.size() > isam2 ) {
        for ( Index isam=isam1; isam<=isam2; ++isam ) {
          float qsam = acd.samples[isam];
          qroi += qsam;
          if ( haveSamples ) {
            if ( qsam < hmin ) hmin = qsam;
            if ( qsam > hmax ) hmax = qsam;
          } else {
            hmin = qsam;
            hmax = qsam;
            haveSamples = true;
          }
        }
      } else {
        cout << myname << "WARNING: Ignoring missing samples for run " << acd.run()
             << ", event " << acd.event() << ", channel " << acd.channel() << endl;
      }
      tdat.nsam[iroi] = 1 + isam2 - isam1;
      tdat.isam[iroi] = isam1;
      tdat.qroi[iroi] = qroi;
      tdat.hmin[iroi] = hmin;
      tdat.hmax[iroi] = hmax;
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
        spre = myname + "  qroi: [";
        for ( Index iroi=0; iroi<nroi; ++iroi ) {
          cout << spre;
          cout << tdat.qroi[iroi];
          spre = ", ";
        }
        cout << "]" << endl;
      }
    }
    ptre->Fill();
  }
  ptre->ResetBranchAddresses();
  ptre->Write();
  gDirectory->Purge();
  pfil->Close();
  delete pfil;
  ret.setInt("art_nfill", nfill);
  return ret;
}

//**********************************************************************

DEFINE_ART_CLASS_TOOL(AdcRoiToTree)
