////////////////////////////////////////////////////////////////////////
// Class:       HelloAuxDet
// Plugin Type: analyzer (art v2_10_03)
// File:        HelloAuxDet_module.cc
// Brief:       Demonstration of how to access AuxDetGeos and AuxDetDigits 
//              in LArSoft.  Checks for presence of AuxDetGeos in geometry 
//              service for ProtoDUNE-SP.  Absence may be due to needing 
//              a dedicated ChannelMapAlg to register AuxDets with 
//              the Geometry service.  
//
// Generated at Wed May 16 09:02:34 2018 by Andrew Olivier (aolivier@ur.rochester.edu) using cetskelgen
// from cetlib version v3_02_00.
////////////////////////////////////////////////////////////////////////

//ART includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h" 

//LArSoft includes
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/RawData/AuxDetDigit.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "nusimdata/SimulationBase/MCParticle.h"

//ROOT includes
#include "TH1D.h"

namespace ex {
  class HelloAuxDet;
}

//Helper functions
namespace
{
  /*template <class STREAM, class FLOAT>
  STREAM& PrintVec(STREAM& stream, const FLOAT& x, const FLOAT& y, const FLOAT& z)
  {
    stream << "(" << x << ", " << y << ", " << z << ")";
    return stream;
  }

  void Print(const raw::AuxDetDigit& digit)
  {
    mf::LogInfo("AuxDetDigits") << "Got a digit for a detector named " << digit.AuxDetName() << " at channel " << digit.Channel() << " at time "
                                << digit.TimeStamp() << "\n";
  }

  void Print(const sim::AuxDetSimChannel& channel)
  {
    mf::LogInfo("AuxDetSimChannels") << "Got an AuxDetSimChannel for ID " << channel.AuxDetID() << " sensitive volume "
                                     << channel.AuxDetSensitiveID() << " with " << channel.AuxDetIDEs().size() << " energy deposits:\n";
    for(const auto& ide: channel.AuxDetIDEs())
    {
      auto stream = mf::LogInfo("AuxDetIDEs");
      stream << "TrackID: " << ide.trackID
             << "energyDeposit: " << ide.energyDeposited
             << "entry point: "; 
             ::PrintVec(stream, ide.entryX, ide.entryY, ide.entryZ);
      stream << "exit point: ";
      ::PrintVec(stream, ide.exitX, ide.exitY, ide.exitZ);
    }
  }

  //TODO: Member function?
  //Print all of the data PRODUCTs of a specific type produced by each label in labels.  If labels is empty, Print all of the 
  //data PRODUCTs in the event.  
  template <class PRODUCT, class LABEL>
  void PrintAllIfEmpty(const std::vector<LABEL>& labels, const art::Event& e)
  {
    std::vector<art::Handle<std::vector<PRODUCT>>> algs;
    if(labels.empty()) //If no labels requested, print all PRODUCTs
    {
      e.getManyByType(algs);
    }
    else
    {
      for(const auto& label: labels)
      {
        algs.emplace_back();
        e.getByLabel(label, algs.back());
      }
    }
    for(const auto& prods: algs)
    {
      for(const auto& prod: *prods) ::Print(prod);
    }
  }*/  
}

class ex::HelloAuxDet : public art::EDAnalyzer {
public:
  explicit HelloAuxDet(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  HelloAuxDet(HelloAuxDet const &) = delete;
  HelloAuxDet(HelloAuxDet &&) = delete;
  HelloAuxDet & operator = (HelloAuxDet const &) = delete;
  HelloAuxDet & operator = (HelloAuxDet &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;

private:

  // Configuration data
  std::vector<art::InputTag> fDigitLabels; //Search for AuxDetDigits produced by modules with these labels
  std::vector<art::InputTag> fSimLabels; //Search for AuxDetSimChannels produced by modules with these labels
  art::InputTag fPartLabel; //Label of the module that produced simb::MCParticles

  // Observer pointers to histograms that will be written
  TH1D* fXDistToTrueTrackFront; //Distance in x from CRT AuxDetSimChannel to where true trajectory passed through a CRT module.  
  TH1D* fYDistToTrueTrackFront; //Distance in y from CRT AuxDetSimChannel to where true trajectory passed through a CRT module. 
};


ex::HelloAuxDet::HelloAuxDet(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{
  fDigitLabels = p.get<std::vector<art::InputTag>>("DigitLabels", {});
  fSimLabels = p.get<std::vector<art::InputTag>>("SimLabels", {});
  fPartLabel = p.get<art::InputTag>("PartLabel", "largeant");
}

void ex::HelloAuxDet::analyze(art::Event const & e)
{
  //Print() AuxDetDigits
  //::PrintAllIfEmpty<raw::AuxDetDigit>(fDigitLabels, e);

  //Print() AuxDetSimChannels
  //::PrintAllIfEmpty<sim::AuxDetSimChannel>(fSimLabels, e);

  //Produce a map from TrackId to MCParticle
  const auto parts = e.getValidHandle<std::vector<simb::MCParticle>>(fPartLabel);
  std::map<int, simb::MCParticle> idToPart; //Map of Geant TrackId to particle
  for(const auto& part: *parts) idToPart[part.TrackId()] = part;

  idToPart[-1] = idToPart[0]; //TODO: Remove this hack that makes sure the primary particle is recognized
  mf::LogInfo("Number of Particles") << "There are " << parts->size() << " simb::MCParticles in this event.\n";

  //Make histograms showing that channel mapping works for AuxDetSimChannels
  art::ServiceHandle<geo::Geometry> geom;
  for(const auto& label: fSimLabels)
  {
    const auto channels = e.getValidHandle<std::vector<sim::AuxDetSimChannel>>(label); 
    for(const auto& channel: *channels)
    {
      const auto& auxDet = geom->AuxDet(channel.AuxDetID());
      const auto& sens = auxDet.SensitiveVolume(channel.AuxDetSensitiveID());
      const auto center = sens.GetCenter();
      TLorentzVector centerV(center.X(), center.Y(), center.Z(), 0.); //TODO: Use GenVector when I have time to figure it out
      for(const auto& ide: channel.AuxDetIDEs())
      {
        const auto found = idToPart.find(ide.trackID);
        if(found != idToPart.end())
        {
          const auto& part = found->second;
          const auto closestPt = std::min_element(part.Trajectory().begin(), part.Trajectory().end(), 
                                                  [&centerV](const auto& first, const auto& second)
                                                  {
                                                    //SimulationBase coordinates are in mm, but Geometry service coordinates are in cm!  
                                                    //Convert to cm.
                                                    return (first.first*0.1-centerV).Vect().Mag() < (second.first*0.1-centerV).Vect().Mag();
                                                  });

          if(closestPt != part.Trajectory().end())
          {
            //Figure out the x and y positions where the CRT was hit
            //TODO: This is a hack to partially get around my lack of knowledge about the CRT channel mapping.  
            //      I'd like to figure this out by identifying the AuxDet.  
            if(sens.toWorldCoords(geo::AuxDetSensitiveGeo::LocalPoint_t(0, 0, sens.HalfLength())).Y() - center.Y() 
               > sens.toWorldCoords(geo::AuxDetSensitiveGeo::LocalPoint_t(0, 0, sens.HalfLength())).X() - center.X()) 
            //If a CRT hit that gives x position
            {
              fXDistToTrueTrackFront->Fill(centerV.X() - closestPt->first.X());
            }
            else //This is a CRT hit that gives y position
            {
              fYDistToTrueTrackFront->Fill(centerV.Y() - closestPt->first.Y());
            }
          }
          else mf::LogWarning("No Closest Point") << "No point found that is considered \"closest\" to AuxDetSensitive!\n";
        }
        else mf::LogWarning("No Such TrackID") << "No map entry for TrackID " << ide.trackID << "\n";
      }
    }
  }
}

void ex::HelloAuxDet::beginJob()
{
  // Print all of the AuxDetGeos that the Geometry service knows about.  
  const auto geom = lar::providerFrom<geo::Geometry>();

  //Get AuxDetGeos and print their names.  
  for(size_t det = 0; det < geom->NAuxDets(); ++det) 
  {
    const auto& geo = geom->AuxDet(det);
    mf::LogInfo("AuxDetGeometry") << "AuxDetGeo number " << det << ", " << geo.Name() << ", is centered at " << geo.GetCenter() 
                                  << " and has sensitive volumes:\n";
    for(size_t sens = 0; sens < geo.NSensitiveVolume(); ++sens)
    {
      const auto& sensGeo = geo.SensitiveVolume(sens);
      const auto center = sensGeo.GetCenter();
      const bool xPos = (sensGeo.toWorldCoords(geo::AuxDetSensitiveGeo::LocalPoint_t(0, 0, sensGeo.HalfLength())).Y() - center.Y() 
                         > sensGeo.toWorldCoords(geo::AuxDetSensitiveGeo::LocalPoint_t(0, 0, sensGeo.HalfLength())).X() - center.X());
      mf::LogInfo("AuxDetGeometry") << "Sensitive volume " << sens << ": " << sensGeo.TotalVolume()->GetName() << " with center " 
                                    << sensGeo.GetCenter() << " which gives information about " << std::string(xPos?"X":"Y") << " coordinates.\n";
    }
  }

  //Tell the framework that I intend to consume AuxDetDigits and AuxDetSimChannels
  for(const auto& label: fDigitLabels) consumes<std::vector<raw::AuxDetDigit>>(label);
  for(const auto& label: fSimLabels) consumes<std::vector<sim::AuxDetSimChannel>>(label);
  consumes<std::vector<simb::MCParticle>>(fPartLabel);

  //Create histograms owned by "The Framework" that I will fill
  art::ServiceHandle<art::TFileService> tfs;
  fXDistToTrueTrackFront = tfs->make<TH1D>("XDistToTrueTrackFront", "X Distance from True Trajectory to CRT Hits"
                                                                    ";X[mm];Events", 150, -30, 30);
  fYDistToTrueTrackFront = tfs->make<TH1D>("YDistToTrueTrackFront", "Y Distance from True Trajectory to CRT Hits"
                                                                    ";Y[mm];Events", 150, -30, 30);
}

DEFINE_ART_MODULE(ex::HelloAuxDet)
