#ifndef PHOTOLIBFAILURESTUDY_H
#define PHOTOLIBFAILURESTUDY_H

//
// Name: PhotoLibFailureStudy.h
//
//
// Ben Jones, MIT, April 2012
//   bjpjones@mit.edu
//
#include <iostream>
#include <vector>

#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h" 

#include "CLHEP/Random/RandFlat.h"

#include "Geometry/Geometry.h"
#include "Utilities/AssociationUtil.h"
#include "Utilities/DetectorProperties.h"
#include "PhotonPropagation/PhotonVisibilityService.h"
#include "Simulation/PhotonVoxels.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"

namespace phot {

  class PhotoLibFailureStudy : public art::EDAnalyzer
  {
  public:
 
    // Constructors, destructor

    explicit PhotoLibFailureStudy(fhicl::ParameterSet const& pset);
    virtual ~PhotoLibFailureStudy();

    // Overrides.

    void reconfigure(fhicl::ParameterSet const& pset);
    void beginJob();
    void analyze(const art::Event& evt);
    
    void endJob();
  
    

  private:
      std::string fAltXAxis;
      int         fOpDet;
      float       fFailureRate;
      int         fNHardware;
      int         fGanging;
      int         fMaskNum;
      int         fMaskDenom;
  };
  
}

#endif

namespace phot {

  //----------------------------------------------------------------------------
  PhotoLibFailureStudy::PhotoLibFailureStudy(const fhicl::ParameterSet& pset)
    : EDAnalyzer(pset)
  {
    reconfigure(pset);
  }

  //----------------------------------------------------------------------------
  PhotoLibFailureStudy::~PhotoLibFailureStudy()
  {
  }

  //----------------------------------------------------------------------------
  void PhotoLibFailureStudy::reconfigure(fhicl::ParameterSet const& pset)
  {
      fAltXAxis    = pset.get<std::string>("alt_x_axis");
      fOpDet       = pset.get<int>("opdet");
      fFailureRate = pset.get<float>("failure_rate");
      fNHardware   = pset.get<int>("nhardware_chan");
      fGanging     = pset.get<int>("ganging");
      fMaskNum     = pset.get<int>("mask_num");
      fMaskDenom   = pset.get<int>("mask_denom");

  }

  //----------------------------------------------------------------------------
  void PhotoLibFailureStudy::beginJob()
  
  {
    mf::LogInfo("PhotoLibFailureStudy")<<"Analyzing photon library - begin"<<
      std::endl;

    
    art::ServiceHandle<art::TFileService> tfs;
    art::ServiceHandle<PhotonVisibilityService> pvs;
    art::ServiceHandle<geo::Geometry> geom;


    
    // Prepare structure to store failures
    int NOpDet = pvs->NOpChannels();
    std::vector<double> failrate(NOpDet, 1);
    std::vector<std::vector<int>> goodchannel;
    for (int od = 0; od < NOpDet; od++) {
        goodchannel.push_back(std::vector<int>(fNHardware, 1));
    }


    // Mask a fraction of the channels
    for (int od = 0; od < NOpDet; od++) {
        if (od % fMaskDenom < fMaskNum) {
            // Turn off this channel
            for (int hc = 0; hc < fNHardware; hc++) {
                goodchannel[od][hc] = 0;
            }
        }
    }

    
    // Choose channels with failures
    int nFailures = (int)round(fFailureRate * (float)NOpDet * (float)fNHardware);
    for (int ifail = 0; ifail < nFailures; ifail++) {
      int od = round(CLHEP::RandFlat::shootInt(NOpDet));
      int hc = round(CLHEP::RandFlat::shootInt(fNHardware));
      if ( goodchannel[od][hc] ) {
        goodchannel[od][hc] = 0;
      }
      else {
        ifail--;
      }
    }

    
    // Multiply up failures by ganging rate
    for (int od = 0; od < NOpDet; od++) {
        int groups = fNHardware/fGanging;
        for (int g = 0; g < groups; g++) {
            int AllGood = 1;
            for (int hc = g*fGanging; hc < (g+1)*fGanging; hc++){
                AllGood *= goodchannel[od][hc];
            }
            if (!AllGood) {
                for (int hc = g*fGanging; hc < (g+1)*fGanging; hc++){
                    goodchannel[od][hc] = 0;
                }
            }
        }
    }


    // Turn channel failures into failure rate
    double totalGood = 0;
    double totalAll = 0;
    for (int od = 0; od < NOpDet; od++) {
        int GoodCount = 0;
        for (auto val: goodchannel[od] ) {
            GoodCount += val;
        }
        failrate[od] = (float)GoodCount / (float)fNHardware ;
        totalGood += GoodCount;
        totalAll  += fNHardware;
    }
    
    
    

    std::cout << "============================+" << std::endl;
    std::cout << "Failed channels." <<std::endl;
    std::cout << "  Nominal failures: " << fFailureRate << std::endl;
    std::cout << "  Total failures:   " << totalGood/totalAll << std::endl;
    for (int od = 0; od < NOpDet; od++) {
        std::cout << "  " << std::setw(2) << od << ": " << std::setw(4) << failrate[od] << "  ";
        for (int hc = 0; hc < fNHardware; hc++) {
            if (goodchannel[od][hc]) std::cout << "-";
            else                     std::cout << "X";
        }
        std::cout << std::endl;
    }
    std::cout << "============================-" << std::endl;
    
    
    sim::PhotonVoxelDef TheVoxelDef = pvs->GetVoxelDef();
    TVector3 Steps = TheVoxelDef.GetSteps();
    TVector3 UpperCorner = TheVoxelDef.GetRegionUpperCorner();
    TVector3 LowerCorner = TheVoxelDef.GetRegionLowerCorner();

    mf::LogInfo("PhotoLibFailureStudy") << "UpperCorner: " << UpperCorner[0] << " " << UpperCorner[1] << " " << UpperCorner[2] << "\n"
                                         << "LowerCorner: " << LowerCorner[0] << " " << LowerCorner[1] << " " << LowerCorner[2];

    int XSteps = int(Steps[0]);
    int YSteps = int(Steps[1]);
    int ZSteps = int(Steps[2]);

    TH3D *FullVolume = tfs->make<TH3D>("FullVolume","FullVolume", 
                                       XSteps,LowerCorner[0],UpperCorner[0],
                                       YSteps,LowerCorner[1],UpperCorner[1],
                                       ZSteps,LowerCorner[2],UpperCorner[2]);
 
    
    mf::LogInfo("PhotoLibFailureStudy")<<"Analyzing photon library - making historams"<<
      std::endl;

    TH2D* XProjection;
    if (fAltXAxis == "Z") XProjection = tfs->make<TH2D>("XProjection","XProjection",ZSteps,0,ZSteps,YSteps,0,YSteps);
    else                  XProjection = tfs->make<TH2D>("XProjection","XProjection",YSteps,0,YSteps,ZSteps,0,ZSteps);
    TH2D* YProjection = tfs->make<TH2D>("YProjection","YProjection",XSteps,0,XSteps,ZSteps,0,ZSteps);
    TH2D* ZProjection = tfs->make<TH2D>("ZProjection","ZProjection",XSteps,0,XSteps,YSteps,0,YSteps);
   
    //    TH1D * PMTsNoVisibility = tfs->make<TH1D>("PMTsNoVisibility","PMTsNoVisibility", NOpDet,0,NOpDet); 

    TH1D* VisByN = tfs->make<TH1D>("VisByN","VisByN", NOpDet, 0, NOpDet);

    TH2D* XInvisibles;
    if (fAltXAxis == "Z") XInvisibles = tfs->make<TH2D>("XInvisibles","XInvisibles",ZSteps,0,ZSteps,YSteps,0,YSteps);
    else                  XInvisibles = tfs->make<TH2D>("XInvisibles","XInvisibles",YSteps,0,YSteps,ZSteps,0,ZSteps);
    TH2D* YInvisibles = tfs->make<TH2D>("YInvisibles","YInvisibles",XSteps,0,XSteps,ZSteps,0,ZSteps);
    TH2D* ZInvisibles = tfs->make<TH2D>("ZInvisibles","ZInvisibles",XSteps,0,XSteps,YSteps,0,YSteps);


  
    std::vector<TH2D*> TheXCrossSections;
    std::vector<TH2D*> TheYCrossSections;
    std::vector<TH2D*> TheZCrossSections;

    
    for(int i=0; i!=XSteps; ++i)
      {
	std::stringstream ss("");
	ss.flush();
	ss<<"projX"<<i;
        if (fAltXAxis == "Z") 	
            TheXCrossSections.push_back(tfs->make<TH2D>(ss.str().c_str(),ss.str().c_str(), ZSteps, 0,ZSteps, YSteps, 0,YSteps));
        else	
            TheXCrossSections.push_back(tfs->make<TH2D>(ss.str().c_str(),ss.str().c_str(), YSteps, 0,YSteps, ZSteps, 0,ZSteps));                  


      }

  for(int i=0; i!=YSteps; ++i)
      {
	std::stringstream ss("");
	ss.flush();
	ss<<"projY"<<i;
	TheYCrossSections.push_back(tfs->make<TH2D>(ss.str().c_str(),ss.str().c_str(), XSteps, 0,XSteps, ZSteps, 0, ZSteps));

      }

  for(int i=0; i!=ZSteps; ++i)
      {
	std::stringstream ss("");
	ss.flush();
	ss<<"projZ"<<i;
	TheZCrossSections.push_back(tfs->make<TH2D>(ss.str().c_str(),ss.str().c_str(), XSteps, 0,XSteps, YSteps, 0,YSteps));

      }
    
  mf::LogInfo("PhotoLibFailureStudy")<<"Analyzing photon library - running through voxels "<<
    std::endl;

  int reportnum=10000;

  int newX, newY;
  if (fAltXAxis == "Z") { 
      newX = 2; // Z
      newY = 1; // Y
  }
  else {
      newX = 1; // Y
      newY = 2; // Z
  }

  for(int i=0; i!=TheVoxelDef.GetNVoxels(); ++i)
    {
      if(i%reportnum==0) std::cout<<"Photon library analyzer at voxel " << i<<std::endl;
    
      std::vector<int> Coords = TheVoxelDef.GetVoxelCoords(i);
                
      const std::vector<float>* Visibilities = pvs->GetLibraryEntries(i);
      
      float TotalVis=0;
      if (fOpDet < 0) {
        for(size_t ichan=0; ichan!=Visibilities->size(); ++ichan)
        {
          TotalVis+=Visibilities->at(ichan)*failrate[ichan];	
        }
      }
      else {
        TotalVis = Visibilities->at(fOpDet)*failrate[fOpDet];
      }
      
      VisByN->Fill(Visibilities->size());
      
      if(TotalVis==0)
	{
	  XInvisibles->Fill(Coords[newX],Coords[newY]);
	  YInvisibles->Fill(Coords[0],Coords[2]);
	  ZInvisibles->Fill(Coords[0],Coords[1]);
	}

      TheXCrossSections.at(Coords.at(0))->Fill(Coords[newX],Coords[newY],TotalVis);
      TheYCrossSections.at(Coords.at(1))->Fill(Coords[0],Coords[2],TotalVis);
      TheZCrossSections.at(Coords.at(2))->Fill(Coords[0],Coords[1],TotalVis);
      XProjection->Fill(Coords[newX], Coords[newY], TotalVis);
      YProjection->Fill(Coords[0], Coords[2], TotalVis);
      ZProjection->Fill(Coords[0], Coords[1], TotalVis);
     
    }
  
  mf::LogInfo("PhotoLibFailureStudy")<<"Analyzing photon library - end"<<
    std::endl;
  }

  //----------------------------------------------------------------------------
  void PhotoLibFailureStudy::analyze(const art::Event& /*evt*/)
  {
  
  }



  
  //----------------------------------------------------------------------------
  void PhotoLibFailureStudy::endJob()
  {

  }
}


namespace phot {
  DEFINE_ART_MODULE(PhotoLibFailureStudy)
}
