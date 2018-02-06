//
// CFilter: Module to use M. Worcester's MuonCounter module as a filter
//
// M. Elnimr Dec. 2014
//
//
//
//

#include "art/Framework/Core/EDFilter.h" 
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "TMath.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "dune/Geometry/MuonCounter35Alg.h"



namespace filt{

 class CFilter : public art::EDFilter {
   public:
     explicit CFilter(fhicl::ParameterSet const& pset);
    virtual ~CFilter() { }
    virtual bool filter(art::Event& e);
    void    reconfigure(fhicl::ParameterSet const& pset);

  private:

    // the parameters we'll read from the .fcl file
    std::string fSimulationProducerLabel;
    std::string fRawDigitLabel;
   //    std::string fCounterFile;
   int fCounterType;
    int fSelectedPDG; // PDG code for primary particle 
    int fTrigger; // as per Matthew which Trigger type 
    // log file
   //    ofstream outfile; 


   int  fnumTracks;
   std::string fHitsModuleLabel;
   std::string fTrackModuleLabel;

  };

   CFilter::CFilter(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset);

  }

  void CFilter::reconfigure(fhicl::ParameterSet const& p)
  {
    
    fSimulationProducerLabel = p.get< std::string >("SimulationLabel");
    //    fCounterFile = p.get< std::string >("CounterFile");
    fSelectedPDG = p.get< int >("PDGcode");
    fCounterType=p.get< int >("CounterType");
    fTrigger=p.get< int >("Trigger");
    fnumTracks=1;

    //      fnumTracks = pset.get<int>("numTracks");
      
  }
  
  bool CFilter::filter(art::Event& event)
  {
    //muon conter geometry
    int counters_loaded=-1;
    std::vector<std::vector<double> > countergeometry;
    //load the muon counter positions from a text file

    char  counterfile[]= "/afs/fnal.gov/files/home/room1/melnimr/lr_dev_7/srcs/dunetpc/dune/Geometry/muoncounters.txt";
    counters_loaded=geo::MuonCounter35Alg::loadMuonCounterGeometry(counterfile,countergeometry);
    bool keepFlag=0;

    //    int icount=0; int trkind;
    // define the handle as an MCParticle vector and fill it with events from the simulation
    art::Handle< std::vector<simb::MCParticle> > particleHandle;
    event.getByLabel(fSimulationProducerLabel, particleHandle);
    // define a sorted map in which to put the particles
    std::map< int, const simb::MCParticle* > particleMap;
    // loop over all the particles, find the primary muon, and get the initial/final positions
    for ( auto const& particle : (*particleHandle) )
      {
        // For the methods you can call to get particle information,
        // see $NUTOOLS_INC/SimulationBase/MCParticle.h.
        int fTrackID = particle.TrackId();
        // Add the address of the MCParticle to the map, with the track ID as the key.
        particleMap[fTrackID] = &particle; 
        // PDG code of every particle in the event.
	//        int fPDG = particle.PdgCode();
        // locate the primary muon
	//        if ( particle.Process() == "primary"  &&  fPDG == fSelectedPDG )
	// {

	    // A particle has a trajectory, consisting of a set of
            // 4-positions and 4-mommenta.
            size_t numberTrajectoryPoints = particle.NumberTrajectoryPoints();
            // For trajectories, as for vectors and arrays, the
            // first point is #0, not #1.
            int last = numberTrajectoryPoints - 1;
            const TLorentzVector& positionStart = particle.Position(0);
            const TLorentzVector& positionEnd   = particle.Position(last);
            const TLorentzVector& momentumStart = particle.Momentum(0);
            const TLorentzVector& momentumEnd   = particle.Momentum(last);
            // move the initial position (for muon counter studies)
            double x_increment = 0.; // 349.666
            double z_increment = 0.; // 231.065
            TVector3 trackStart(positionStart.X()+x_increment,positionStart.Y(),positionStart.Z()+z_increment);
	    //}

	    std::vector< std::vector<double> > hitcounters;	    
	    if(counters_loaded){
	      //unsigned int counters_hit=0;
	      
	      geo::MuonCounter35Alg::testTrackInAllCounters(fTrackID,trackStart,momentumStart.Vect(),
                                                            countergeometry,hitcounters);
	    }
	    

	    // condition flags for each layer
	    bool Layer_1_2 = false;
	    bool Layer_3_4_5 = false;
	    bool Layer_E = false;
	    bool Layer_W = false;
	    bool Layer_N_U = false;
	    bool Layer_N_L = false;
	    bool Layer_S_U = false;
	    bool Layer_S_L = false;
	    int Trigger= 0;

	    for(unsigned int ii=0;ii<hitcounters.size();ii++)
	      {
		//      for(unsigned int jj=0;jj<hitcounters[ii].size();jj++)
		//	{
		std::cout<< "::Filter, counters triggered are: "<<hitcounters[ii][0] <<std::endl;
		/*	  if (hitcounters[ii][0]==fCounterType) 
		  {
		  keepFlag=1;
		  return keepFlag;
		  }*/
                // check which layer is hit
                if( 40 <= hitcounters[ii][0] && hitcounters[ii][0] <= 61){
                  Layer_1_2 = true;
                }
                if( hitcounters[ii][0] > 61){
                  Layer_3_4_5 = true;
                }
                if (14 <= hitcounters[ii][0] && hitcounters[ii][0] <=19){
                  Layer_N_U = true;
                }
                if ( 34 <= hitcounters[ii][0] && hitcounters[ii][0] <= 39){
                  Layer_S_L = true;
                }
                if (8 <= hitcounters[ii][0] && hitcounters[ii][0] <= 13){
                  Layer_N_L = true;
                }
                if (28 <= hitcounters[ii][0] && hitcounters[ii][0] <= 33){
                  Layer_S_U = true;
                }
                if (hitcounters[ii][0] <= 7){
                  Layer_E = true;
                }
                if (20 <= hitcounters[ii][0] && hitcounters[ii][0] <= 27){
                  Layer_W = true;
                }
		//	}
	      }
	    // check for a satisfied trigger condition
	    if (Layer_1_2 && Layer_3_4_5){
	      Trigger = 1;
	    }
	    if (Layer_N_U && Layer_S_L){
	      Trigger = 2;
	    }
	    if (Layer_N_L && Layer_S_U){
                Trigger = 3;
	    }
	    if (Layer_E && Layer_W){
	      Trigger = 4;
	    }
	    std::cout << "Trigger is ....." << Trigger << std::endl;
	    if(Trigger==fTrigger)
	      {
		keepFlag=1;
		return keepFlag;
	      }
      }
    return keepFlag;
  }
  // A macro required for a JobControl module.
  DEFINE_ART_MODULE(CFilter)
    
} // namespace filt
