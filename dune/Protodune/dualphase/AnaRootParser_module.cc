// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/View.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/Geometry.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/MCBase/MCStep.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/ExternalTrigger.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/BeamInfo.h"
#include "lardataobj/RawData/TriggerData.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larreco/Deprecated/BezierTrack.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "larreco/RecoAlg/TrajectoryMCSFitter.h" //uBoone fitter
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/FlashMatch.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/AnalysisBase/MVAPIDResult.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "lardata/ArtDataHelper/MVAReader.h"

#include <cstddef> // std::ptrdiff_t
#include <cstring> // std::memcpy()
#include <vector>
#include <map>
#include <iterator> // std::begin(), std::end()
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <functional> // std::mem_fun_ref
#include <typeinfo>
#include <memory> // std::unique_ptr<>

#include "TTree.h"
#include "TTimeStamp.h"

constexpr int kNplanes       = 2;     //number of wire planes
constexpr int kMaxHits       = 40000; //maximum number of hits;
constexpr int kMaxTrackHits  = 2000;  //maximum number of hits on a track
constexpr int kMaxTrackers   = 15;    //number of trackers passed into fTrackModuleLabel
constexpr int kMaxVertices   = 500;    //max number of 3D vertices
constexpr int kMaxVertexAlgos = 10;    //max number of vertex algorithms
constexpr unsigned short kMaxAuxDets = 4; ///< max number of auxiliary detector cells per MC particle
constexpr int kMaxFlashes      = 1000;   //maximum number of flashes
constexpr int kMaxExternCounts = 1000;   //maximum number of External Counters
constexpr int kMaxShowerHits   = 10000;  //maximum number of hits on a shower
constexpr int kMaxTruth        = 10;     //maximum number of neutrino truth interactions
constexpr int kMaxClusters     = 1000;   //maximum number of clusters;
constexpr int kMaxReadoutTicksInAllChannels = 76800000; //protoDUNE: 10 000 ticks * (8*960) channel
constexpr int kMaxChannels = 7680; //protoDUNE: (8*960) channel

constexpr int kMaxNDaughtersPerPFP = 10; //maximum number of daughters per PFParticle
constexpr int kMaxNClustersPerPFP  = 10; //maximum number of clusters per PFParticle
constexpr int kMaxNPFPNeutrinos    = 5;  //maximum number of reconstructed neutrino PFParticles

/// total_extent\<T\>::value has the total number of elements of an array
template <typename T>
struct total_extent {
  using value_type = size_t;
  static constexpr value_type value
    = sizeof(T) / sizeof(typename std::remove_all_extents<T>::type);
}; // total_extent<>


namespace dune {

  /// Data structure with all the tree information.
  ///
  /// Can connect to a tree, clear its fields and resize its data.
  class AnaRootParserDataStruct {
    public:

      /// A wrapper to a C array (needed to embed an array into a vector)
      template <typename Array_t>
        class BoxedArray {
          protected:
            Array_t array; // actual data

          public:
            using This_t = BoxedArray<Array_t>;
            typedef typename std::remove_all_extents<Array_t>::type Data_t;

            BoxedArray() {} // no initialization
            BoxedArray(const This_t& from)
            { std::memcpy((char*) &(data()), (char*) &(from.data()), sizeof(Array_t)); }

            Array_t& data() { return array; }
            const Array_t& data() const { return array; }

            //@{
            /// begin/end interface
            static constexpr size_t size() { return total_extent<Array_t>::value; }
            Data_t* begin() { return reinterpret_cast<Data_t*>(&array); }
            const Data_t* begin() const { return reinterpret_cast<const Data_t*>(&array); }
            Data_t* end() { return begin() + size(); }
            const Data_t* end() const { return begin() + size(); }
            //@}

            //@{
            /// Array interface
            auto operator[] (size_t index) -> decltype(*array) { return array[index]; }
            auto operator[] (size_t index) const -> decltype(*array) { return array[index]; }
            auto operator+ (ptrdiff_t index) -> decltype(&*array) { return array + index; }
            auto operator+ (ptrdiff_t index) const -> decltype(&*array) { return array + index; }
            auto operator- (ptrdiff_t index) -> decltype(&*array) { return array - index; }
            auto operator- (ptrdiff_t index) const -> decltype(&*array) { return array - index; }
            auto operator* () -> decltype(*array) { return *array; }
            auto operator* () const -> decltype(*array) { return *array; }

            operator decltype(&array[0]) () { return &array[0]; }
            operator decltype(&array[0]) () const { return &array[0]; }
            //@}

        }; // BoxedArray

      /// Tracker algorithm result
      ///
      /// Can connect to a tree, clear its fields and resize its data.
      class TrackDataStruct {
        public:
          /* Data structure size:
           *
           * TrackData_t<Short_t>                    :  2  bytes/track
           * TrackData_t<Float_t>                    :  4  bytes/track
           * TrackPlaneData_t<Float_t>, TrackPlaneData_t<Int_t>: 12  bytes/track
           * TrackHitData_t<Float_t>                      : 24k bytes/track
           * TrackHitCoordData_t<Float_t>                 : 72k bytes/track
           */
          template <typename T>
            using TrackData_t = std::vector<T>;
          template <typename T>
            using TrackPlaneData_t = std::vector<BoxedArray<T[kNplanes]>>;
          template <typename T>
            using TrackHitData_t = std::vector<BoxedArray<T[kNplanes][kMaxTrackHits]>>;
          template <typename T>
            using TrackHitCoordData_t = std::vector<BoxedArray<T[kNplanes][kMaxTrackHits][3]>>;

          size_t MaxTracks; ///< maximum number of storable tracks

          Short_t  ntracks;             //number of reconstructed tracks
          TrackPlaneData_t<Float_t>    trkke;
          TrackPlaneData_t<Float_t>    trkrange;
          TrackPlaneData_t<Int_t>      trkidtruth;  //true geant trackid
          TrackPlaneData_t<Short_t>    trkorigin;   //_ev_origin 0: unknown, 1: neutrino, 2: cosmic, 3: supernova, 4: singles
          TrackPlaneData_t<Int_t>      trkpdgtruth; //true pdg code
          TrackPlaneData_t<Float_t>    trkefftruth; //completeness
          TrackPlaneData_t<Float_t>    trkpurtruth; //purity of track
          TrackPlaneData_t<Float_t>    trkpitchc;
          TrackPlaneData_t<Short_t>    ntrkhitsperview;
          TrackHitData_t<Float_t>      trkdedx;
          TrackHitData_t<Float_t>      trkdqdx;
          TrackHitData_t<Float_t>      trkresrg;
          TrackHitData_t<Int_t>        trktpc;
          TrackHitCoordData_t<Float_t> trkxyz;

          Float_t trkdedx2[10000];
          Float_t trkdqdx2[10000];
          Float_t trktpc2[10000];
          Float_t trkx2[10000];
          Float_t trky2[10000];
          Float_t trkz2[10000];


          Float_t hittrkx[10000];
          Float_t hittrky[10000];
          Float_t hittrkz[10000];

	  Float_t hittrklocaltrackdirectionx[10000];
	  Float_t hittrklocaltrackdirectiony[10000];
	  Float_t hittrklocaltrackdirectionz[10000];
	  Float_t hittrklocaltrackdirectiontheta[10000];
	  Float_t hittrklocaltrackdirectionphi[10000];

          Float_t hittrkpitchC[10000];

          Float_t hittrkds[10000];

          Short_t hittrkchannel[10000];
          Short_t hittrktpc[10000];
          Short_t hittrkview[10000];
          Short_t hittrkwire[10000];
          Float_t hittrkpeakT[10000];
          Float_t hittrkchargeintegral[10000];
	  Float_t hittrkph[10000];
          Float_t hittrkchargesum[10000];
          Float_t hittrkstarT[10000];
          Float_t hittrkendT[10000];
          Float_t hittrkrms[10000];
          Float_t hittrkgoddnessofFit[10000];
          Short_t hittrkmultiplicity[10000];
          Int_t hittrktrueID[10000];
          Float_t hittrktrueEnergyMax[10000];
          Float_t hittrktrueEnergyFraction[10000];

          // more track info
          TrackData_t<Short_t> trkId;
          TrackData_t<Short_t> NHitsPerTrack;
          TrackData_t<Short_t> trkncosmictags_tagger;
          TrackData_t<Float_t> trkcosmicscore_tagger;
          TrackData_t<Short_t> trkcosmictype_tagger;
          TrackData_t<Short_t> trkncosmictags_containmenttagger;
          TrackData_t<Float_t> trkcosmicscore_containmenttagger;
          TrackData_t<Short_t> trkcosmictype_containmenttagger;
          TrackData_t<Short_t> trkncosmictags_flashmatch;
          TrackData_t<Float_t> trkcosmicscore_flashmatch;
          TrackData_t<Short_t> trkcosmictype_flashmatch;
          TrackData_t<Float_t> trkstartx;     // starting x position.
          TrackData_t<Float_t> trkstarty;     // starting y position.
          TrackData_t<Float_t> trkstartz;     // starting z position.
          TrackData_t<Float_t> trkstartd;     // starting distance to boundary.
          TrackData_t<Float_t> trkendx;       // ending x position.
          TrackData_t<Float_t> trkendy;       // ending y position.
          TrackData_t<Float_t> trkendz;       // ending z position.
          TrackData_t<Float_t> trkendd;       // ending distance to boundary.
          TrackData_t<Float_t> trkflashT0;   // t0 per track from matching tracks to flashes (in ns)
          TrackData_t<Float_t> trktrueT0;    // t0 per track from truth information (in ns)
          TrackData_t<Float_t> trkpurity;    // track purity based on hit information
          TrackData_t<Float_t> trkcompleteness; //track completeness based on hit information
          TrackData_t<int> trkg4id;        //true g4 track id for the reconstructed track
          TrackData_t<int> trkorig;        //origin of the track
          TrackData_t<Float_t> trkstarttheta;      // theta.
          TrackData_t<Float_t> trkstartphi;        // phi.
          TrackData_t<Float_t> trkendtheta;      // theta.
          TrackData_t<Float_t> trkendphi;        // phi.
          TrackData_t<Float_t> trkstartdirectionx;
          TrackData_t<Float_t> trkstartdirectiony;
          TrackData_t<Float_t> trkstartdirectionz;
          TrackData_t<Float_t> trkenddirectionx;
          TrackData_t<Float_t> trkenddirectiony;
          TrackData_t<Float_t> trkenddirectionz;
          TrackData_t<Float_t> trkthetaxz;    // theta_xz.
          TrackData_t<Float_t> trkthetayz;    // theta_yz.
          TrackData_t<Float_t> trkmom;        // momentum.
          TrackData_t<Float_t> trkchi2PerNDF;        // length along trajectory.
          TrackData_t<Float_t> trkNDF;        // length along trajectory.
          TrackData_t<Float_t> trklen;        // length along trajectory.
	  TrackData_t<Float_t> trklenstraightline; // shortest distance betweem start and end point of track
          TrackData_t<Float_t> trkmomrange;    // track momentum from range using CSDA tables
          TrackData_t<Float_t> trkmommschi2;   // track momentum from multiple scattering Chi2 method
          TrackData_t<Float_t> trkmommsllhd;   // track momentum from multiple scattering LLHD method
          TrackData_t<Float_t> trkmommscmic;    //TODO comment here
          TrackData_t<Float_t> trkmommscfwd;    //TODO comment here
          TrackData_t<Float_t> trkmommscbwd;    //TODO comment here
          TrackData_t<Float_t> trkmommscllfwd;  //TODO comment here
          TrackData_t<Float_t> trkmommscllbwd;  //TODO comment here
          TrackData_t<Short_t> trksvtxid;     // Vertex ID associated with the track start
          TrackData_t<Short_t> trkevtxid;     // Vertex ID associated with the track end
          TrackPlaneData_t<Int_t> trkpidpdg;       // particle PID pdg code
          TrackPlaneData_t<Float_t> trkpidchi;
          TrackPlaneData_t<Float_t> trkpidchipr;   // particle PID chisq for proton
          TrackPlaneData_t<Float_t> trkpidchika;   // particle PID chisq for kaon
          TrackPlaneData_t<Float_t> trkpidchipi;   // particle PID chisq for pion
          TrackPlaneData_t<Float_t> trkpidchimu;   // particle PID chisq for muon
          TrackPlaneData_t<Float_t> trkpidpida;    // particle PIDA
          TrackData_t<Float_t> trkpidmvamu;   // particle MVA value for muon PID
          TrackData_t<Float_t> trkpidmvae;   // particle MVA value for electron PID
          TrackData_t<Float_t> trkpidmvapich;   // particle MVA value for charged pion PID
          TrackData_t<Float_t> trkpidmvaphoton;   // particle MVA value for photon PID
          TrackData_t<Float_t> trkpidmvapr;   // particle MVA value for proton PID
          TrackData_t<Short_t> trkpidbestplane; // this is defined as the plane with most hits

          TrackData_t<Short_t> trkhasPFParticle; // whether this belongs to a PFParticle
          TrackData_t<Short_t> trkPFParticleID;  // if hasPFParticle, its ID

          /// Creates an empty tracker data structure
          TrackDataStruct(): MaxTracks(0) { Clear(); }
          /// Creates a tracker data structure allowing up to maxTracks tracks
          TrackDataStruct(size_t maxTracks): MaxTracks(maxTracks) { Clear(); }
          void Clear();
          void SetMaxTracks(size_t maxTracks)
          { MaxTracks = maxTracks; Resize(MaxTracks); }
          void Resize(size_t nTracks);
          void SetAddresses(TTree* pTree, std::string tracker, bool isCosmics);

          size_t GetMaxTracks() const { return MaxTracks; }
          size_t GetMaxPlanesPerTrack(int /* iTrack */ = 0) const
          { return (size_t) kNplanes; }
          size_t GetMaxHitsPerTrack(int /* iTrack */ = 0, int /* ipl */ = 0) const
          { return (size_t) kMaxTrackHits; }

      }; // class TrackDataStruct

      //Vertex data struct
      class VertexDataStruct {
        public:
          template <typename T>
            using VertexData_t = std::vector<T>;

          size_t MaxVertices; ///< maximum number of storable vertices

          Short_t  nvtx;             //number of reconstructed tracks
          VertexData_t<Short_t> vtxId;    // the vertex ID.
          VertexData_t<Float_t> vtxx;     // x position.
          VertexData_t<Float_t> vtxy;     // y position.
          VertexData_t<Float_t> vtxz;     // z position.

          VertexData_t<Short_t> vtxhasPFParticle; // whether this belongs to a PFParticle
          VertexData_t<Short_t> vtxPFParticleID;  // if hasPFParticle, its ID

          VertexDataStruct(): MaxVertices(0) { Clear(); }
          VertexDataStruct(size_t maxVertices): MaxVertices(maxVertices) { Clear(); }
          void Clear();
          void SetMaxVertices(size_t maxVertices)
          { MaxVertices = maxVertices; Resize(MaxVertices); }
          void Resize(size_t nVertices);
          void SetAddresses(TTree* pTree, std::string tracker, bool isCosmics);

          size_t GetMaxVertices() const { return MaxVertices; }
      }; // class VertexDataStruct


      /// Shower algorithm result
      ///
      /// Can connect to a tree, clear its fields and resize its data.
      class ShowerDataStruct {
        public:
          /* Data structure size:
           *
           * ShowerData_t<Short_t>                   :  2  bytes/shower
           * ShowerData_t<Float_t>                   :  4  bytes/shower
           * ShowerPlaneData_t<Float_t>, ShowerPlaneData_t<Int_t>: 12  bytes/shower
           * ShowerHitData_t<Float_t>                      : 24k bytes/shower
           * ShowerHitCoordData_t<Float_t>                 : 72k bytes/shower
           */
          template <typename T>
            using ShowerData_t = std::vector<T>;
          template <typename T>
            using ShowerPlaneData_t = std::vector<BoxedArray<T[kNplanes]>>;
          template <typename T>
            using ShowerHitData_t = std::vector<BoxedArray<T[kNplanes][kMaxShowerHits]>>;
          template <typename T>
            using ShowerHitCoordData_t = std::vector<BoxedArray<T[kNplanes][kMaxShowerHits][3]>>;

          std::string name; ///< name of the shower algorithm (for branch names)

          size_t MaxShowers; ///< maximum number of storable showers

          /// @{
          /// @name Branch data structures
          Short_t  nshowers;                      ///< number of showers
          ShowerData_t<Short_t>  showerID;        ///< Shower ID
          ShowerData_t<Short_t>  shwr_bestplane;  ///< Shower best plane
          ShowerData_t<Float_t>  shwr_length;     ///< Shower length
          ShowerData_t<Float_t>  shwr_startdcosx; ///< X directional cosine at start of shower
          ShowerData_t<Float_t>  shwr_startdcosy; ///< Y directional cosine at start of shower
          ShowerData_t<Float_t>  shwr_startdcosz; ///< Z directional cosine at start of shower
          ShowerData_t<Float_t>  shwr_startx;     ///< startx of shower
          ShowerData_t<Float_t>  shwr_starty;     ///< starty of shower
          ShowerData_t<Float_t>  shwr_startz;     ///< startz of shower
          ShowerPlaneData_t<Float_t>   shwr_totEng;     ///< Total energy of the shower per plane
          ShowerPlaneData_t<Float_t>   shwr_dedx;       ///< dE/dx of the shower per plane
          ShowerPlaneData_t<Float_t>   shwr_mipEng;     ///< Total MIP energy of the shower per plane

          ShowerData_t<Float_t> shwr_pidmvamu;   // particle MVA value for muon PID
          ShowerData_t<Float_t> shwr_pidmvae;   // particle MVA value for electron PID
          ShowerData_t<Float_t> shwr_pidmvapich;   // particle MVA value for charged pion PID
          ShowerData_t<Float_t> shwr_pidmvaphoton;   // particle MVA value for photon PID
          ShowerData_t<Float_t> shwr_pidmvapr;   // particle MVA value for proton PID


          ShowerData_t<Short_t>  shwr_hasPFParticle; // whether this belongs to a PFParticle
          ShowerData_t<Short_t>  shwr_PFParticleID;  // if hasPFParticle, its ID
          /// @}

          /// Creates a shower data structure allowing up to maxShowers showers
          ShowerDataStruct(std::string new_name = "", size_t maxShowers = 0):
            name(new_name), MaxShowers(maxShowers) { Clear(); }

          std::string Name() const { return name; }

          void Clear();

          /// Applies a special prescription to mark shower information as missing
          void MarkMissing(TTree* pTree);
          void SetName(std::string new_name) { name = new_name; }
          void SetMaxShowers(size_t maxShowers)
          { MaxShowers = maxShowers; Resize(MaxShowers); }
          void Resize(size_t nShowers);
          void SetAddresses(TTree* pTree);

          size_t GetMaxShowers() const { return MaxShowers; }
          size_t GetMaxPlanesPerShower(int /* iShower */ = 0) const
          { return (size_t) kNplanes; }
          size_t GetMaxHitsPerShower(int /* iShower */ = 0, int /* ipl */ = 0) const
          { return (size_t) kMaxShowerHits; }

      }; // class ShowerDataStruct

      class PFParticleDataStruct {
        public:
          /* Data structure size:
           *
           * PFParticleData_t<Short_t>   :  2  bytes/PFParticle
           * PFParticleData_t<Int_t>     :  4  bytes/PFParticle
           * DaughterData_t<Short_t>     :  20 bytes/PFParticle
           * ClusterData_t<Short_t>      :  20 bytes/PFParticle
           * Short_t [kMaxNPFPNeutrinos] :  10 bytes in total
           */
          template <typename T>
            using PFParticleData_t = std::vector<T>;
          template <typename T>
            using DaughterData_t = std::vector<BoxedArray<T[kMaxNDaughtersPerPFP]>>;
          template <typename T>
            using ClusterData_t = std::vector<BoxedArray<T[kMaxNClustersPerPFP]>>;

          size_t MaxPFParticles; ///< maximum number of storable PFParticles

          /// @{
          /// @name Branch data structures
          Short_t                   nPFParticles;     ///< the total number of PFParticles
          PFParticleData_t<Short_t> pfp_selfID;       ///< the PFParticles' own IDs
          PFParticleData_t<Short_t> pfp_isPrimary;    ///< whether the PFParticle is a primary particle

          PFParticleData_t<Short_t> pfp_numDaughters; ///< the number of daughters belonging to this PFParticle
          DaughterData_t<Short_t>   pfp_daughterIDs;  ///< the IDs of the daughter PFParticles
          PFParticleData_t<Short_t> pfp_parentID;     ///< the ID of this PFParticle's immediate parent

          PFParticleData_t<Short_t> pfp_vertexID;     ///< the ID of the vertex belonging to this PFParticle
          PFParticleData_t<Short_t> pfp_isShower;     ///< whether this PFParticle corresponds to a shower
          PFParticleData_t<Short_t> pfp_isTrack;      ///< whether this PFParticle corresponds to a track
          PFParticleData_t<Short_t> pfp_trackID;      ///< the ID of the track object corresponding to this PFParticle, if !isShower
          PFParticleData_t<Short_t> pfp_showerID;     ///< the ID of the shower object corresponding to this PFParticle, if isShower

          PFParticleData_t<Short_t> pfp_isNeutrino;   ///< whether this PFParticle is a neutrino
          PFParticleData_t<Int_t>   pfp_pdgCode;      ///< the preliminary estimate of the PFParticle type using the PDG code

          PFParticleData_t<Short_t> pfp_numClusters;  ///< the number of associated clusters
          ClusterData_t<Short_t>    pfp_clusterIDs;   ///< the IDs of any associated clusters

          Short_t                   pfp_numNeutrinos; ///< the number of reconstructed neutrinos
          Short_t pfp_neutrinoIDs[kMaxNPFPNeutrinos]; ///< the PFParticle IDs of the neutrinos
          /// @}

          /// Creates a PFParticle data structure allowing up to maxPFParticles PFParticles
          PFParticleDataStruct(size_t maxPFParticles = 0):
            MaxPFParticles(maxPFParticles) { Clear(); }

          void Clear();
          void SetMaxPFParticles(size_t maxPFParticles)
          { MaxPFParticles = maxPFParticles; Resize(MaxPFParticles); }
          void Resize(size_t numPFParticles);
          void SetAddresses(TTree* pTree);

          size_t GetMaxPFParticles() const { return MaxPFParticles; }
          size_t GetMaxDaughtersPerPFParticle(int /* iPFParticle */ = 0) const
          { return (size_t) kMaxNDaughtersPerPFP; }
          size_t GetMaxClustersPerPFParticle(int /* iPFParticle */ = 0) const
          { return (size_t) kMaxNClustersPerPFP; }

      }; // class PFParticleDataStruct

      enum DataBits_t: unsigned int {
        tdAuxDet = 0x01,
        tdCry = 0x02,
        tdGenie = 0x04,
	tdPhotons = 0x08,
        tdGenerator = 0x10,
        tdGeant = 0x20,
	tdGeantInAV = 0x40,
	tdGeantTrajectory = 0x80,
        tdHit = 0x100,
        tdTrack = 0x200,
        tdVertex = 0x400,
        tdFlash = 0x800,
        tdShower = 0x1000,
        tdMCshwr = 0x2000,
        tdMCtrk  = 0x4000,
        tdCluster = 0x8000,
        tdRawDigit = 0x10000,
	tdRecobWire = 0x20000,
        tdPandoraNuVertex = 0x40000,
        tdPFParticle = 0x80000,
        tdCount = 0x100000,
        tdProto = 0x200000,
        tdDefault = 0
      }; // DataBits_t

      /*    /// information from the run
            struct RunData_t {
            public:
            RunData_t() { Clear(); }
            void Clear() {}
            }; // struct RunData_t
            */
      /// information from the subrun
      struct SubRunData_t {
        SubRunData_t() { Clear(); }
        void Clear() {
          pot = -999.;
          potbnbETOR860 = -999.;
          potbnbETOR875 = -999.;
          potnumiETORTGT = -999.;
        }
        Double_t pot; //protons on target
        Double_t potbnbETOR860;
        Double_t potbnbETOR875;
        Double_t potnumiETORTGT;
      }; // struct SubRunData_t

      //    RunData_t    RunData; ///< run data collected at begin of run
      SubRunData_t SubRunData; ///< subrun data collected at begin of subrun

      //run information
      Int_t      run;                  //run number
      Int_t      subrun;               //subrun number
      Int_t      event;                //event number
      Int_t   evttime_seconds;              //event time in seconds
      Int_t   evttime_nanoseconds;              //event time in nanoseconcds
      Double_t   beamtime;             //beam time
      Double_t   pot;                  //protons on target moved in subrun data
      Double_t   taulife;              //electron lifetime
      Char_t     isdata;               //flag, 0=MC 1=data
      unsigned int triggernumber;      //trigger counter
      Double_t     triggertime;        //trigger time w.r.t. electronics clock T0
      Double_t     beamgatetime;       //beamgate time w.r.t. electronics clock T0
      unsigned int triggerbits;        //trigger bits
      Double_t     potbnb;             //pot per event (BNB E:TOR860)
      Double_t     potnumitgt;         //pot per event (NuMI E:TORTGT)
      Double_t     potnumi101;         //pot per event (NuMI E:TOR101)

      // hit information (non-resizeable, 45x kMaxHits = 900k bytes worth)
      Int_t    no_hits;                  //number of hits
      Int_t    NHitsInAllTracks;        //number of hits in all tracks
      Int_t    no_hits_stored;           //number of hits actually stored in the tree
      Short_t  hit_tpc[kMaxHits];        //tpc number
      Short_t  hit_view[kMaxHits];      //plane number
      Short_t  hit_wire[kMaxHits];       //wire number
      Short_t  hit_channel[kMaxHits];    //channel ID
      Float_t  hit_peakT[kMaxHits];      //peak time
      Float_t  hit_chargesum[kMaxHits];     //charge (sum)
      Float_t  hit_chargeintegral[kMaxHits];     //charge (integral)
      Float_t  hit_ph[kMaxHits];         //amplitude
      Float_t  hit_startT[kMaxHits];     //hit start time
      Float_t  hit_endT[kMaxHits];       //hit end time
      Float_t  hit_rms[kMaxHits];       //hit rms from the hit object
      Float_t  hit_goodnessOfFit[kMaxHits]; //chi2/dof goodness of fit
      Float_t  hit_fitparamampl[kMaxHits]; //dual phase hit fit
      Float_t  hit_fitparamt0[kMaxHits]; //dual phase hit fit
      Float_t  hit_fitparamtau1[kMaxHits]; //dual phase hit fit
      Float_t  hit_fitparamtau2[kMaxHits]; //dual phase hit fit
      Short_t  hit_multiplicity[kMaxHits];  //multiplicity of the given hit
      Int_t    hit_trueID[kMaxHits];  //true mctackID form backtracker
      Float_t  hit_trueEnergyMax[kMaxHits]; //energy deposited from that mctrackID
      Float_t  hit_trueEnergyFraction[kMaxHits]; //maxe/tote
      //    Float_t  hit_trueX[kMaxHits];      // hit true X (cm)
      //    Float_t  hit_nelec[kMaxHits];     //hit number of electrons
      //    Float_t  hit_energy[kMaxHits];       //hit energy
      Short_t  hit_trkid[kMaxHits];      //is this hit associated with a reco track?
      //    Short_t  hit_trkKey[kMaxHits];      //is this hit associated with a reco track,  if so associate a unique track key ID?
      Short_t  hit_clusterid[kMaxHits];  //is this hit associated with a reco cluster?
      //    Short_t  hit_clusterKey[kMaxHits];  //is this hit associated with a reco cluster, if so associate a unique cluster key ID?
/*
      Float_t rawD_ph[kMaxHits];
      Float_t rawD_peakT[kMaxHits];
      Float_t rawD_charge[kMaxHits];
      Float_t rawD_fwhh[kMaxHits];
      Double_t rawD_rms[kMaxHits];
*/
      Int_t no_channels;               //number of readout channels with raw waveform (can be different from number of max channels in simulation)
      Int_t no_ticks;                  //number of readout ticks for raw waveform
      Int_t no_ticksinallchannels;     //number of readout ticks multiplied by no_channels

      Int_t rawD_Channel[kMaxChannels];
      Short_t rawD_ADC[kMaxReadoutTicksInAllChannels];

      Int_t no_recochannels;               //number of readout channels with "reco" waveform (can be different from number of max channels in simulation)
      Int_t recoW_Channel[kMaxChannels];
      Int_t recoW_NTicks[kMaxChannels];

      Int_t no_recoticksinallchannels;     //number of readout ticks multiplied by no_channels
      Int_t recoW_Tick[kMaxReadoutTicksInAllChannels];
      Float_t recoW_ADC[kMaxReadoutTicksInAllChannels];

      //Light information
      size_t MaxPhotons = 0;
      Int_t numberofphotons;
      std::vector<Float_t> photons_time;
      std::vector<Float_t> photons_channel;

      //Pandora Nu Vertex information
      Short_t nnuvtx;
      Float_t nuvtxx[kMaxVertices];
      Float_t nuvtxy[kMaxVertices];
      Float_t nuvtxz[kMaxVertices];
      Short_t nuvtxpdg[kMaxVertices];

      //Cluster Information
      Short_t nclusters;				      //number of clusters in a given event
      Short_t clusterId[kMaxClusters];		      //ID of this cluster
      Short_t clusterView[kMaxClusters];	      //which plane this cluster belongs to
      Int_t   cluster_isValid[kMaxClusters];	      //is this cluster valid? will have a value of -1 if it is not valid
      Float_t cluster_StartCharge[kMaxClusters];	       //charge on the first wire of the cluster in ADC
      Float_t cluster_StartAngle[kMaxClusters];	      //starting angle of the cluster
      Float_t cluster_EndCharge[kMaxClusters];	      //charge on the last wire of the cluster in ADC
      Float_t cluster_EndAngle[kMaxClusters];	      //ending angle of the cluster
      Float_t cluster_Integral[kMaxClusters];	      //returns the total charge of the cluster from hit shape in ADC
      Float_t cluster_IntegralAverage[kMaxClusters];    //average charge of the cluster hits in ADC
      Float_t cluster_SummedADC[kMaxClusters];	      //total charge of the cluster from signal ADC counts
      Float_t cluster_SummedADCaverage[kMaxClusters];   //average signal ADC counts of the cluster hits.
      Float_t cluster_MultipleHitDensity[kMaxClusters]; //Density of wires in the cluster with more than one hit.
      Float_t cluster_Width[kMaxClusters];	      //cluster width in ? units
      Short_t cluster_NHits[kMaxClusters];	      //Number of hits in the cluster
      Short_t cluster_StartWire[kMaxClusters];	      //wire coordinate of the start of the cluster
      Short_t cluster_StartTick[kMaxClusters];	      //tick coordinate of the start of the cluster in time ticks
      Short_t cluster_EndWire[kMaxClusters];	      //wire coordinate of the end of the cluster
      Short_t cluster_EndTick[kMaxClusters];            //tick coordinate of the end of the cluster in time ticks
      //Cluster cosmic tagging information
      //    Short_t cluncosmictags_tagger[kMaxClusters];      //No. of cosmic tags associated to this cluster
      //    Float_t clucosmicscore_tagger[kMaxClusters];      //Cosmic score associated to this cluster. In the case of more than one tag, the first one is associated.
      //    Short_t clucosmictype_tagger[kMaxClusters];       //Cosmic tag type for this cluster.

      // flash information
      Int_t    no_flashes;                //number of flashes
      Float_t  flash_time[kMaxFlashes];   //flash time
      Float_t  flash_pe[kMaxFlashes];     //flash total PE
      Float_t  flash_ycenter[kMaxFlashes];//y center of flash
      Float_t  flash_zcenter[kMaxFlashes];//z center of flash
      Float_t  flash_ywidth[kMaxFlashes]; //y width of flash
      Float_t  flash_zwidth[kMaxFlashes]; //z width of flash
      Float_t  flash_timewidth[kMaxFlashes]; //time of flash

      // External Counter information
      Int_t   no_ExternCounts;                    // Number of External Triggers
      Float_t externcounts_time[kMaxExternCounts]; // Time of External Trigger
      Float_t externcounts_id[kMaxExternCounts];   // ID of External Trigger

      //track information
      Char_t   kNTracker;
      std::vector<TrackDataStruct> TrackData;

      //vertex information
      Char_t   kNVertexAlgos;
      std::vector<VertexDataStruct> VertexData;

      // shower information
      Char_t   kNShowerAlgos;
      std::vector<ShowerDataStruct> ShowerData;

      // PFParticle information
      PFParticleDataStruct PFParticleData;

      //mctruth information
      Int_t     mcevts_truth;    //number of neutrino Int_teractions in the spill
      Int_t     nuPDG_truth[kMaxTruth];     //neutrino PDG code
      Int_t     ccnc_truth[kMaxTruth];      //0=CC 1=NC
      Int_t     mode_truth[kMaxTruth];      //0=QE/El, 1=RES, 2=DIS, 3=Coherent production
      Float_t  enu_truth[kMaxTruth];       //true neutrino energy
      Float_t  Q2_truth[kMaxTruth];        //Momentum transfer squared
      Float_t  W_truth[kMaxTruth];         //hadronic invariant mass
      Float_t  X_truth[kMaxTruth];
      Float_t  Y_truth[kMaxTruth];
      Int_t     hitnuc_truth[kMaxTruth];    //hit nucleon
      Float_t  nuvtxx_truth[kMaxTruth];    //neutrino vertex x
      Float_t  nuvtxy_truth[kMaxTruth];    //neutrino vertex y
      Float_t  nuvtxz_truth[kMaxTruth];    //neutrino vertex z
      Float_t  nu_dcosx_truth[kMaxTruth];  //neutrino dcos x
      Float_t  nu_dcosy_truth[kMaxTruth];  //neutrino dcos y
      Float_t  nu_dcosz_truth[kMaxTruth];  //neutrino dcos z
      Float_t  lep_mom_truth[kMaxTruth];   //lepton momentum
      Float_t  lep_dcosx_truth[kMaxTruth]; //lepton dcos x
      Float_t  lep_dcosy_truth[kMaxTruth]; //lepton dcos y
      Float_t  lep_dcosz_truth[kMaxTruth]; //lepton dcos z

      //flux information
      Float_t  vx_flux[kMaxTruth];          //X position of hadron/muon decay (cm)
      Float_t  vy_flux[kMaxTruth];          //Y position of hadron/muon decay (cm)
      Float_t  vz_flux[kMaxTruth];          //Z position of hadron/muon decay (cm)
      Float_t  pdpx_flux[kMaxTruth];        //Parent X momentum at decay point (GeV)
      Float_t  pdpy_flux[kMaxTruth];        //Parent Y momentum at decay point (GeV)
      Float_t  pdpz_flux[kMaxTruth];        //Parent Z momentum at decay point (GeV)
      Float_t  ppdxdz_flux[kMaxTruth];      //Parent dxdz direction at production
      Float_t  ppdydz_flux[kMaxTruth];      //Parent dydz direction at production
      Float_t  pppz_flux[kMaxTruth];        //Parent Z momentum at production (GeV)

      Int_t    ptype_flux[kMaxTruth];        //Parent GEANT code particle ID
      Float_t  ppvx_flux[kMaxTruth];        //Parent production vertex X (cm)
      Float_t  ppvy_flux[kMaxTruth];        //Parent production vertex Y (cm)
      Float_t  ppvz_flux[kMaxTruth];        //Parent production vertex Z (cm)
      Float_t  muparpx_flux[kMaxTruth];     //Muon neutrino parent production vertex X (cm)
      Float_t  muparpy_flux[kMaxTruth];     //Muon neutrino parent production vertex Y (cm)
      Float_t  muparpz_flux[kMaxTruth];     //Muon neutrino parent production vertex Z (cm)
      Float_t  mupare_flux[kMaxTruth];      //Muon neutrino parent energy (GeV)

      Int_t    tgen_flux[kMaxTruth];        //Parent generation in cascade. (1 = primary proton,
      //2=particles produced by proton interaction, 3 = particles
      //produced by interactions of the 2's, ...
      Int_t    tgptype_flux[kMaxTruth];     //Type of particle that created a particle flying of the target
      Float_t  tgppx_flux[kMaxTruth];       //X Momentum of a particle, that created a particle that flies
      //off the target, at the interaction point. (GeV)
      Float_t  tgppy_flux[kMaxTruth];       //Y Momentum of a particle, that created a particle that flies
      //off the target, at the interaction point. (GeV)
      Float_t  tgppz_flux[kMaxTruth];       //Z Momentum of a particle, that created a particle that flies
      //off the target, at the interaction point. (GeV)
      Float_t  tprivx_flux[kMaxTruth];      //Primary particle interaction vertex X (cm)
      Float_t  tprivy_flux[kMaxTruth];      //Primary particle interaction vertex Y (cm)
      Float_t  tprivz_flux[kMaxTruth];      //Primary particle interaction vertex Z (cm)
      Float_t  dk2gen_flux[kMaxTruth];      //distance from decay to ray origin (cm)
      Float_t  gen2vtx_flux[kMaxTruth];     //distance from ray origin to event vtx (cm)

      Float_t  tpx_flux[kMaxTruth];        //Px of parent particle leaving BNB/NuMI target (GeV)
      Float_t  tpy_flux[kMaxTruth];        //Py of parent particle leaving BNB/NuMI target (GeV)
      Float_t  tpz_flux[kMaxTruth];        //Pz of parent particle leaving BNB/NuMI target (GeV)
      Int_t    tptype_flux[kMaxTruth];     //Type of parent particle leaving BNB target

      //genie information
      size_t MaxGeniePrimaries = 0;
      Int_t     genie_no_primaries;
      std::vector<Int_t>    genie_primaries_pdg;
      std::vector<Float_t>  genie_Eng;
      std::vector<Float_t>  genie_Px;
      std::vector<Float_t>  genie_Py;
      std::vector<Float_t>  genie_Pz;
      std::vector<Float_t>  genie_P;
      std::vector<Int_t>    genie_status_code;
      std::vector<Float_t>  genie_mass;
      std::vector<Int_t>    genie_trackID;
      std::vector<Int_t>    genie_ND;
      std::vector<Int_t>    genie_mother;

      //cosmic cry information
      Int_t     mcevts_truthcry;    //number of neutrino Int_teractions in the spill
      Int_t     cry_no_primaries;
      std::vector<Int_t>    cry_primaries_pdg;
      std::vector<Float_t>  cry_Eng;
      std::vector<Float_t>  cry_Px;
      std::vector<Float_t>  cry_Py;
      std::vector<Float_t>  cry_Pz;
      std::vector<Float_t>  cry_P;
      std::vector<Float_t>  cry_StartPointx;
      std::vector<Float_t>  cry_StartPointy;
      std::vector<Float_t>  cry_StartPointz;
      std::vector<Float_t>  cry_StartPointt;
      std::vector<Int_t>    cry_status_code;
      std::vector<Float_t>  cry_mass;
      std::vector<Int_t>    cry_trackID;
      std::vector<Int_t>    cry_ND;
      std::vector<Int_t>    cry_mother;

      // ProtoDUNE Beam generator information
      Int_t proto_no_primaries;
      std::vector<Int_t> proto_isGoodParticle;
      std::vector<Float_t> proto_vx;
      std::vector<Float_t> proto_vy;
      std::vector<Float_t> proto_vz;
      std::vector<Float_t> proto_t;
      std::vector<Float_t> proto_px;
      std::vector<Float_t> proto_py;
      std::vector<Float_t> proto_pz;
      std::vector<Float_t> proto_momentum;
      std::vector<Float_t> proto_energy;
      std::vector<Int_t> proto_pdg;

      //Generator and G4 MC Particle information
      size_t MaxGeneratorparticles = 0; ///! how many particles there is currently room for
      Int_t     generator_list_size;  //number of all geant particles

      size_t MaxGEANTparticles = 0; ///! how many particles there is currently room for
      Int_t     no_primaries;      //number of primary geant particles
      Int_t     geant_list_size;  //number of all geant particles
      Int_t     geant_list_size_in_tpcAV;
      std::vector<Int_t>    pdg;
      std::vector<Int_t>    status;
      std::vector<Int_t>    inTPCActive;
      std::vector<Float_t>  Eng;
      std::vector<Float_t>  EndE;
      std::vector<Float_t>  Mass;
      std::vector<Float_t>  Px;
      std::vector<Float_t>  Py;
      std::vector<Float_t>  Pz;
      std::vector<Float_t>  P;
      std::vector<Float_t>  StartPointx;
      std::vector<Float_t>  StartPointy;
      std::vector<Float_t>  StartPointz;
      std::vector<Float_t>  StartT;
      std::vector<Float_t>  EndT;
      std::vector<Float_t>  EndPointx;
      std::vector<Float_t>  EndPointy;
      std::vector<Float_t>  EndPointz;
      std::vector<Float_t>  theta;
      std::vector<Float_t>  phi;
      std::vector<Float_t>  theta_xz;
      std::vector<Float_t>  theta_yz;

      size_t MaxGEANTInAVparticles = 0;
      std::vector<Float_t>  pathlen_tpcAV;
      std::vector<Int_t>    TrackId_tpcAV;
      std::vector<Int_t>    PDGCode_tpcAV;
      std::vector<Float_t>  StartPointx_tpcAV;
      std::vector<Float_t>  StartPointy_tpcAV;
      std::vector<Float_t>  StartPointz_tpcAV;
      std::vector<Float_t>  StartT_tpcAV;
      std::vector<Float_t>  StartE_tpcAV;
      std::vector<Float_t>  StartP_tpcAV;
      std::vector<Float_t>  StartPx_tpcAV;
      std::vector<Float_t>  StartPy_tpcAV;
      std::vector<Float_t>  StartPz_tpcAV;
      std::vector<Float_t>  thetastart_tpcAV;
      std::vector<Float_t>  phistart_tpcAV;
      std::vector<Float_t>  EndPointx_tpcAV;
      std::vector<Float_t>  EndPointy_tpcAV;
      std::vector<Float_t>  EndPointz_tpcAV;
      std::vector<Float_t>  EndT_tpcAV;
      std::vector<Float_t>  EndE_tpcAV;
      std::vector<Float_t>  EndP_tpcAV;
      std::vector<Float_t>  EndPx_tpcAV;
      std::vector<Float_t>  EndPy_tpcAV;
      std::vector<Float_t>  EndPz_tpcAV;
      std::vector<Float_t>  thetaend_tpcAV;
      std::vector<Float_t>  phiend_tpcAV;

      std::vector<Float_t>  pathlen_drifted;
      std::vector<Int_t>    inTPCDrifted;
      std::vector<Float_t>  StartPointx_drifted;
      std::vector<Float_t>  StartPointy_drifted;
      std::vector<Float_t>  StartPointz_drifted;
      std::vector<Float_t>  StartT_drifted;
      std::vector<Float_t>  StartE_drifted;
      std::vector<Float_t>  StartP_drifted;
      std::vector<Float_t>  StartPx_drifted;
      std::vector<Float_t>  StartPy_drifted;
      std::vector<Float_t>  StartPz_drifted;
      std::vector<Float_t>  EndPointx_drifted;
      std::vector<Float_t>  EndPointy_drifted;
      std::vector<Float_t>  EndPointz_drifted;
      std::vector<Float_t>  EndT_drifted;
      std::vector<Float_t>  EndE_drifted;
      std::vector<Float_t>  EndP_drifted;
      std::vector<Float_t>  EndPx_drifted;
      std::vector<Float_t>  EndPy_drifted;
      std::vector<Float_t>  EndPz_drifted;
      std::vector<Int_t>    NumberDaughters;
      std::vector<Int_t>    TrackId;
      std::vector<Int_t>    Mother;
      std::vector<Int_t>    process_primary;
      std::vector<std::string> processname;
      std::vector<Int_t>    MergedId; //geant track segments, which belong to the same particle, get the same
      std::vector<Int_t>    origin;   ////0: unknown, 1: cosmic, 2: neutrino, 3: supernova, 4: singles
      std::vector<Int_t>    MCTruthIndex; //this geant particle comes from the neutrino interaction of the _truth variables with this index

      //G4 MC trajectory information
      size_t MaxGEANTtrajectorypoints = 0; ///! how many particles there is currently room for

      Int_t     geant_trajectory_size;  //number of trajectory points for all geant particles

      std::vector<Int_t> NTrajectoryPointsPerParticle;


      std::vector<Float_t>  TrajTrackId;
      std::vector<Float_t>  TrajPDGCode;

      std::vector<Float_t>  TrajX;
      std::vector<Float_t>  TrajY;
      std::vector<Float_t>  TrajZ;
      std::vector<Float_t>  TrajT;
      std::vector<Float_t>  TrajE;
      std::vector<Float_t>  TrajP;
      std::vector<Float_t>  TrajPx;
      std::vector<Float_t>  TrajPy;
      std::vector<Float_t>  TrajPz;
      std::vector<Float_t>  TrajTheta;
      std::vector<Float_t>  TrajPhi;

      //MC Shower information
      Int_t     no_mcshowers;                         //number of MC Showers in this event.
      //MC Shower particle information
      std::vector<Int_t>       mcshwr_origin;	    //MC Shower origin information.
      std::vector<Int_t>       mcshwr_pdg;	    //MC Shower particle PDG code.
      std::vector<Int_t>       mcshwr_TrackId;        //MC Shower particle G4 track ID.
      std::vector<std::string> mcshwr_Process;	    //MC Shower particle's creation process.
      std::vector<Float_t>     mcshwr_startX;	    //MC Shower particle G4 startX
      std::vector<Float_t>     mcshwr_startY;	    //MC Shower particle G4 startY
      std::vector<Float_t>     mcshwr_startZ;	    //MC Shower particle G4 startZ
      std::vector<Float_t>     mcshwr_endX;	    //MC Shower particle G4 endX
      std::vector<Float_t>     mcshwr_endY;	    //MC Shower particle G4 endY
      std::vector<Float_t>     mcshwr_endZ;	    //MC Shower particle G4 endZ
      std::vector<Float_t>    mcshwr_CombEngX;	    //MC Shower Combined energy deposition information, Start Point X Position.
      std::vector<Float_t>    mcshwr_CombEngY;	    //MC Shower Combined energy deposition information, Start Point Y Position.
      std::vector<Float_t>    mcshwr_CombEngZ;	    //MC Shower Combined energy deposition information, Start Point Z Position.
      std::vector<Float_t>     mcshwr_CombEngPx;	    //MC Shower Combined energy deposition information, Momentum X direction.
      std::vector<Float_t>     mcshwr_CombEngPy;	    //MC Shower Combined energy deposition information, Momentum X direction.
      std::vector<Float_t>     mcshwr_CombEngPz;	    //MC Shower Combined energy deposition information, Momentum X direction.
      std::vector<Float_t>     mcshwr_CombEngE;	    //MC Shower Combined energy deposition information, Energy
      std::vector<Float_t>     mcshwr_dEdx;           //MC Shower dEdx, MeV/cm
      std::vector<Float_t>     mcshwr_StartDirX;      //MC Shower Direction of begining of shower, X direction
      std::vector<Float_t>     mcshwr_StartDirY;      //MC Shower Direction of begining of shower, Y direction
      std::vector<Float_t>     mcshwr_StartDirZ;      //MC Shower Direction of begining of shower, Z direction
      std::vector<Int_t>       mcshwr_isEngDeposited;  //tells whether if this shower deposited energy in the detector or not.
      //yes = 1; no =0;
      //MC Shower mother information
      std::vector<Int_t>       mcshwr_Motherpdg;       //MC Shower's mother PDG code.
      std::vector<Int_t>       mcshwr_MotherTrkId;     //MC Shower's mother G4 track ID.
      std::vector<std::string> mcshwr_MotherProcess;   //MC Shower's mother creation process.
      std::vector<Float_t>     mcshwr_MotherstartX;    //MC Shower's mother  G4 startX .
      std::vector<Float_t>     mcshwr_MotherstartY;    //MC Shower's mother  G4 startY .
      std::vector<Float_t>     mcshwr_MotherstartZ;    //MC Shower's mother  G4 startZ .
      std::vector<Float_t>     mcshwr_MotherendX;	     //MC Shower's mother  G4 endX   .
      std::vector<Float_t>     mcshwr_MotherendY;	     //MC Shower's mother  G4 endY   .
      std::vector<Float_t>     mcshwr_MotherendZ;	     //MC Shower's mother  G4 endZ   .
      //MC Shower ancestor information
      std::vector<Int_t>       mcshwr_Ancestorpdg;       //MC Shower's ancestor PDG code.
      std::vector<Int_t>       mcshwr_AncestorTrkId;     //MC Shower's ancestor G4 track ID.
      std::vector<std::string> mcshwr_AncestorProcess;   //MC Shower's ancestor creation process.
      std::vector<Float_t>     mcshwr_AncestorstartX;    //MC Shower's ancestor  G4 startX
      std::vector<Float_t>     mcshwr_AncestorstartY;    //MC Shower's ancestor  G4 startY
      std::vector<Float_t>     mcshwr_AncestorstartZ;    //MC Shower's ancestor  G4 startZ
      std::vector<Float_t>     mcshwr_AncestorendX;      //MC Shower's ancestor  G4 endX
      std::vector<Float_t>     mcshwr_AncestorendY;      //MC Shower's ancestor  G4 endY
      std::vector<Float_t>     mcshwr_AncestorendZ;      //MC Shower's ancestor  G4 endZ

      //MC track information
      Int_t     no_mctracks;                         //number of MC tracks in this event.
      //MC track particle information
      std::vector<Int_t>       mctrk_origin;	    //MC track origin information.
      std::vector<Int_t>       mctrk_pdg;	    //MC track particle PDG code.
      std::vector<Int_t>       mctrk_TrackId;        //MC track particle G4 track ID.
      std::vector<std::string> mctrk_Process;	    //MC track particle's creation process.
      std::vector<Float_t>     mctrk_startX;	    //MC track particle G4 startX
      std::vector<Float_t>     mctrk_startY;	    //MC track particle G4 startY
      std::vector<Float_t>     mctrk_startZ;	    //MC track particle G4 startZ
      std::vector<Float_t>     mctrk_endX;		//MC track particle G4 endX
      std::vector<Float_t>     mctrk_endY;		//MC track particle G4 endY
      std::vector<Float_t>     mctrk_endZ;		//MC track particle G4 endZ
      std::vector<Float_t>     mctrk_startX_drifted;        //MC track particle first step in TPC x
      std::vector<Float_t>     mctrk_startY_drifted;        //MC track particle first step in TPC y
      std::vector<Float_t>     mctrk_startZ_drifted;        //MC track particle first step in TPC z
      std::vector<Float_t>     mctrk_endX_drifted;          //MC track particle last step in TPC x
      std::vector<Float_t>     mctrk_endY_drifted;          //MC track particle last step in TPC y
      std::vector<Float_t>     mctrk_endZ_drifted;          //MC track particle last step in TPC z
      std::vector<Float_t>     mctrk_len_drifted;           //MC track length within TPC
      std::vector<Float_t>     mctrk_p_drifted;             //MC track momentum at start point in TPC
      std::vector<Float_t>     mctrk_px_drifted;            //MC track x momentum at start point in TPC
      std::vector<Float_t>     mctrk_py_drifted;            //MC track y momentum at start point in TPC
      std::vector<Float_t>     mctrk_pz_drifted;            //MC track z momentum at start point in TPC
      //MC Track mother information
      std::vector<Int_t>       mctrk_Motherpdg;       //MC Track's mother PDG code.
      std::vector<Int_t>       mctrk_MotherTrkId;     //MC Track's mother G4 track ID.
      std::vector<std::string> mctrk_MotherProcess;   //MC Track's mother creation process.
      std::vector<Float_t>     mctrk_MotherstartX;    //MC Track's mother  G4 startX .
      std::vector<Float_t>     mctrk_MotherstartY;    //MC Track's mother  G4 startY .
      std::vector<Float_t>     mctrk_MotherstartZ;    //MC Track's mother  G4 startZ .
      std::vector<Float_t>     mctrk_MotherendX;	     //MC Track's mother  G4 endX   .
      std::vector<Float_t>     mctrk_MotherendY;	     //MC Track's mother  G4 endY   .
      std::vector<Float_t>     mctrk_MotherendZ;	     //MC Track's mother  G4 endZ   .
      //MC Track ancestor information
      std::vector<Int_t>       mctrk_Ancestorpdg;       //MC Track's ancestor PDG code.
      std::vector<Int_t>       mctrk_AncestorTrkId;     //MC Track's ancestor G4 track ID.
      std::vector<std::string> mctrk_AncestorProcess;   //MC Track's ancestor creation process.
      std::vector<Float_t>     mctrk_AncestorstartX;    //MC Track's ancestor  G4 startX
      std::vector<Float_t>     mctrk_AncestorstartY;    //MC Track's ancestor  G4 startY
      std::vector<Float_t>     mctrk_AncestorstartZ;    //MC Track's ancestor  G4 startZ
      std::vector<Float_t>     mctrk_AncestorendX;      //MC Track's ancestor  G4 endX
      std::vector<Float_t>     mctrk_AncestorendY;      //MC Track's ancestor  G4 endY
      std::vector<Float_t>     mctrk_AncestorendZ;      //MC Track's ancestor  G4 endZ

      // Auxiliary detector variables saved for each geant track
      // This data is saved as a vector (one item per GEANT particle) of C arrays
      // (wrapped in a BoxedArray for technical reasons), one item for each
      // affected detector cell (which one is saved in AuxDetID
      template <typename T>
        using AuxDetMCData_t = std::vector<BoxedArray<T[kMaxAuxDets]>>;

      std::vector<UShort_t> NAuxDets;         ///< Number of AuxDets crossed by this particle
      AuxDetMCData_t<Short_t> AuxDetID;       ///< Which AuxDet this particle went through
      AuxDetMCData_t<Float_t> entryX;         ///< Entry X position of particle into AuxDet
      AuxDetMCData_t<Float_t> entryY;         ///< Entry Y position of particle into AuxDet
      AuxDetMCData_t<Float_t> entryZ;         ///< Entry Z position of particle into AuxDet
      AuxDetMCData_t<Float_t> entryT;         ///< Entry T position of particle into AuxDet
      AuxDetMCData_t<Float_t> exitX;          ///< Exit X position of particle out of AuxDet
      AuxDetMCData_t<Float_t> exitY;          ///< Exit Y position of particle out of AuxDet
      AuxDetMCData_t<Float_t> exitZ;          ///< Exit Z position of particle out of AuxDet
      AuxDetMCData_t<Float_t> exitT;          ///< Exit T position of particle out of AuxDet
      AuxDetMCData_t<Float_t> exitPx;         ///< Exit x momentum of particle out of AuxDet
      AuxDetMCData_t<Float_t> exitPy;         ///< Exit y momentum of particle out of AuxDet
      AuxDetMCData_t<Float_t> exitPz;         ///< Exit z momentum of particle out of AuxDet
      AuxDetMCData_t<Float_t> CombinedEnergyDep; ///< Sum energy of all particles with this trackID (+ID or -ID) in AuxDet

      unsigned int bits; ///< complementary information

      /// Returns whether we have auxiliary detector data
      bool hasAuxDetector() const { return bits & tdAuxDet; }

      /// Returns whether we have Cry data
      bool hasCryInfo() const { return bits & tdCry; }

      /// Returns whether we have Genie data
      bool hasGenieInfo() const { return bits & tdGenie; }

      /// Returns whether we have MCShower data
      bool hasMCShowerInfo() const { return bits & tdMCshwr; }

      /// Returns whether we have MCTrack data
      bool hasMCTrackInfo() const { return bits & tdMCtrk; }

      /// Returns whether we have Hit data
      bool hasHitInfo() const { return bits & tdHit; }

      /// Returns whether we have RawDigit data
      bool hasRawDigitInfo() const { return bits & tdRawDigit; }

      /// Returns whether we have RecobWire data
      bool hasRecobWireInfo() const { return bits & tdRecobWire; }

      /// Returns whether we have Track data
      bool hasTrackInfo() const { return bits & tdTrack; }

      /// Returns whether we have Shower data
      bool hasShowerInfo() const { return bits & tdShower; }

      /// Returns whether we have Vertex data
      bool hasVertexInfo() const { return bits & tdVertex; }

      /// Returns whether we have PFParticle data
      bool hasPFParticleInfo() const { return bits & tdPFParticle; }

      /// Returns whether we have Cluster data
      bool hasClusterInfo() const { return bits & tdCluster; }

      /// Returns whether we have Pandora Nu Vertex data
      bool hasPandoraNuVertexInfo() const { return bits & tdPandoraNuVertex; }

      /// Returns whether we have photon data
      bool hasPhotonInfo() const { return bits & tdPhotons; }

      /// Returns whether we have Generator data
      bool hasGeneratorInfo() const { return bits & tdGenerator; }

      /// Returns whether we have Geant data
      bool hasGeantInfo() const { return bits & tdGeant; }

      /// Returns whether we have Geant data
      bool hasGeantInAVInfo() const { return bits & tdGeantInAV; }

      /// Returns whether we have Geant trajectory data
      bool hasGeantTrajectoryInfo() const { return bits & tdGeantTrajectory; }

      /// Returns whether we have Flash data
      bool hasFlashInfo() const { return bits & tdFlash; }

      /// Returns whether we have External Counter data
      bool hasExternCountInfo() const { return bits & tdCount; }

      /// Returns whether we have protoDUNE beam primaries
      bool hasProtoInfo() const {return bits & tdProto; }

      /// Sets the specified bits
      void SetBits(unsigned int setbits, bool unset = false)
      { if (unset) bits &= ~setbits; else bits |= setbits; }

      /// Constructor; clears all fields
      AnaRootParserDataStruct(size_t nTrackers = 0, size_t nVertexAlgos = 0,
          std::vector<std::string> const& ShowerAlgos = {}):
        bits(tdDefault)
        { SetTrackers(nTrackers); SetVertexAlgos(nVertexAlgos); SetShowerAlgos(ShowerAlgos); Clear(); }

      TrackDataStruct& GetTrackerData(size_t iTracker)
      { return TrackData.at(iTracker); }
      const TrackDataStruct& GetTrackerData(size_t iTracker) const
      { return TrackData.at(iTracker); }

      ShowerDataStruct& GetShowerData(size_t iShower)
      { return ShowerData.at(iShower); }
      ShowerDataStruct const& GetShowerData(size_t iShower) const
      { return ShowerData.at(iShower); }

      VertexDataStruct& GetVertexData(size_t iVertex)
      { return VertexData.at(iVertex); }
      const VertexDataStruct& GetVertexData(size_t iVertex) const
      { return VertexData.at(iVertex); }

      PFParticleDataStruct& GetPFParticleData()
      { return PFParticleData; }
      const PFParticleDataStruct& GetPFParticleData() const
      { return PFParticleData; }

      /// Clear all fields if this object (not the tracker algorithm data)
      void ClearLocalData();

      /// Clear all fields
      void Clear();


      /// Allocates data structures for the given number of trackers (no Clear())
      void SetTrackers(size_t nTrackers) { TrackData.resize(nTrackers); }

      /// Allocates data structures for the given number of vertex algos (no Clear())
      void SetVertexAlgos(size_t nVertexAlgos) { VertexData.resize(nVertexAlgos); }

      /// Allocates data structures for the given number of trackers (no Clear())
      void SetShowerAlgos(std::vector<std::string> const& ShowerAlgos);

      /// Resize the data strutcure for Generator particles
      void ResizePhotons(int nPhotons);

      /// Resize the data strutcure for Generator particles
      void ResizeGenerator(int nParticles);

      /// Resize the data strutcure for GEANT particles
      void ResizeGEANT(int nParticles);

      /// Resize the data strutcure for GEANT particles in active volume
      void ResizeGEANTInAV(int nParticlesInAV);

      /// Resize the data strutcure for GEANT trajectory info
      void ResizeGEANTTrajectory(int nTrajectoryPoints);

      /// Resize the data strutcure for Genie primaries
      void ResizeGenie(int nPrimaries);

      /// Resize the data strutcure for Cry primaries
      void ResizeCry(int nPrimaries);

      /// Resize the data structure for ProtoDUNE primaries
      void ResizeProto(int nPrimaries);

      /// Resize the data strutcure for  MC Showers
      void ResizeMCShower(int nMCShowers);

      /// Resize the data strutcure for  MC Tracks
      void ResizeMCTrack(int nMCTracks);

      /// Connect this object with a tree
      void SetAddresses(
          TTree* pTree,
          std::vector<std::string> const& trackers,
          std::vector<std::string> const& vertexalgos,
          std::vector<std::string> const& showeralgos,
          bool isCosmics
          );


      /// Returns the number of trackers for which data structures are allocated
      size_t GetNTrackers() const { return TrackData.size(); }

      /// Returns the number of Vertex algos for which data structures are allocated
      size_t GetNVertexAlgos() const { return VertexData.size(); }

      /// Returns the number of trackers for which data structures are allocated
      size_t GetNShowerAlgos() const { return ShowerData.size(); }

      /// Returns the number of hits for which memory is allocated
      size_t GetMaxHits() const { return kMaxHits; }

      /// Returns the number of trackers for which memory is allocated
      size_t GetMaxTrackers() const { return TrackData.capacity(); }

      /// Returns the number of trackers for which memory is allocated
      size_t GetMaxVertexAlgos() const { return VertexData.capacity(); }

      /// Returns the number of trackers for which memory is allocated
      size_t GetMaxShowers() const { return ShowerData.capacity(); }

      /// Returns the number of GEANT particles for which memory is allocated
      size_t GetMaxPhotons() const { return MaxPhotons; }

      /// Returns the number of GEANT particles for which memory is allocated
      size_t GetMaxGeneratorparticles() const { return MaxGeneratorparticles; }

      /// Returns the number of GEANT particles for which memory is allocated
      size_t GetMaxGEANTparticles() const { return MaxGEANTparticles; }

      /// Returns the number of GEANT particles in AV for which memory is allocated
      size_t GetMaxGEANTInAVparticles() const { return MaxGEANTInAVparticles; }

      /// Returns the number of GEANT particles for which memory is allocated
      size_t GetMaxGEANTtrajectorypoints() const { return MaxGEANTtrajectorypoints; }

      /// Returns the number of GENIE primaries for which memory is allocated
      size_t GetMaxGeniePrimaries() const { return MaxGeniePrimaries; }


    private:
      /// Little helper functor class to create or reset branches in a tree
      class BranchCreator {
        public:
          TTree* pTree; ///< the tree to be worked on
          BranchCreator(TTree* tree): pTree(tree) {}

          //@{
          /// Create a branch if it does not exist, and set its address
          void operator()
            (std::string name, void* address, std::string leaflist /*, int bufsize = 32000 */)
            {

              if (!pTree) return;
              TBranch* pBranch = pTree->GetBranch(name.c_str());

              if (!pBranch) {
                pTree->Branch(name.c_str(), address, leaflist.c_str() /*, bufsize */);

                LOG_DEBUG("AnaRootParserStructure")
                  << "Creating branch '" << name << " with leaf '" << leaflist << "'";

              }
              else if (pBranch->GetAddress() != address) {
                pBranch->SetAddress(address);
                LOG_DEBUG("AnaRootParserStructure")
                  << "Reassigning address to branch '" << name << "'";
              }
              else {
                LOG_DEBUG("AnaRootParserStructure")
                  << "Branch '" << name << "' is fine";
              }
            } // operator()
          void operator()
            (std::string name, void* address, const std::stringstream& leaflist /*, int bufsize = 32000 */)
            { return this->operator() (name, address, leaflist.str() /*, int bufsize = 32000 */); }
          template <typename T>
            void operator()
            (std::string name, std::vector<T>& data, std::string leaflist /*, int bufsize = 32000 */)
            { return this->operator() (name, (void*) data.data(), leaflist /*, int bufsize = 32000 */); }

          template <typename T>
            void operator() (std::string name, std::vector<T>& data)
            {
              // overload for a generic object expressed directly by reference
              // (as opposed to a generic object expressed by a pointer or
              // to a simple leaf sequence specification);
              // TTree::Branch(name, T* obj, Int_t bufsize, splitlevel) and
              // TTree::SetObject() are used.
              if (!pTree) return;
              TBranch* pBranch = pTree->GetBranch(name.c_str());
              if (!pBranch) {
                pTree->Branch(name.c_str(), &data);
                // ROOT needs a TClass definition for T in order to create a branch,
                // se we are sure that at this point the TClass exists
                LOG_DEBUG("AnaRootParserStructure")
                  << "Creating object branch '" << name
                  << " with " << TClass::GetClass(typeid(T))->ClassName();
              }
              else if
                (*(reinterpret_cast<std::vector<T>**>(pBranch->GetAddress())) != &data)
                {
                  // when an object is provided directly, the address of the object
                  // is assigned in TBranchElement::fObject (via TObject::SetObject())
                  // and the address itself is set to the address of the fObject
                  // member. Here we check that the address of the object in fObject
                  // is the same as the address of our current data type
                  pBranch->SetObject(&data);
                  LOG_DEBUG("AnaRootParserStructure")
                    << "Reassigning object to branch '" << name << "'";
                }
              else {
                LOG_DEBUG("AnaRootParserStructure")
                  << "Branch '" << name << "' is fine";
              }
            } // operator()
          //@}
      }; // class BranchCreator

  }; // class AnaRootParserDataStruct


  /// Contains ROOTTreeCode<>::code, ROOT tree character for branch of type T
  template <typename T> struct ROOTTreeCode; // generally undefined

  template<> struct ROOTTreeCode<Short_t>  { static constexpr char code = 'S'; };
  template<> struct ROOTTreeCode<Int_t>    { static constexpr char code = 'I'; };
  template<> struct ROOTTreeCode<Double_t> { static constexpr char code = 'D'; };


  /// Class whose "type" contains the base data type of the container
  template <typename C> struct ContainerValueType; // generally undefined

  template <typename A>
    struct ContainerValueType<std::vector<AnaRootParserDataStruct::BoxedArray<A>>>
    { using type = typename AnaRootParserDataStruct::BoxedArray<A>::value_type; };

  template <typename T>
    struct ContainerValueType<std::vector<T>>
    { using type = typename std::vector<T>::value_type; };


  /**
   * @brief Creates a simple ROOT tree with tracking and calorimetry information
   *
   * <h2>Configuration parameters</h2>
   * - <b>UseBuffers</b> (default: false): if enabled, memory is allocated for
   *   tree data for all the run; otherwise, it's allocated on each event, used
   *   and freed; use "true" for speed, "false" to save memory
   * - <b>SaveAuxDetInfo</b> (default: false): if enabled, auxiliary detector
   *   data will be extracted and included in the tree
   */
  class AnaRootParser : public art::EDAnalyzer {

    public:

      explicit AnaRootParser(fhicl::ParameterSet const& pset);
      virtual ~AnaRootParser();

      /// read access to event
      void analyze(const art::Event& evt);
      //  void beginJob() {}
      void beginSubRun(const art::SubRun& sr);
      void endSubRun(const art::SubRun& sr);

    private:

      void   HitsBackTrack( art::Ptr<recob::Hit> const& hit, int & trackid, float & maxe, float & purity );
      double length(const recob::Track& track);
      double driftedLength(const simb::MCParticle& part, TLorentzVector& start, TLorentzVector& end, unsigned int &starti, unsigned int &endi);
      double driftedLength(const sim::MCTrack& mctrack, TLorentzVector& tpcstart, TLorentzVector& tpcend, TLorentzVector& tpcmom);
      double length(const simb::MCParticle& part, TLorentzVector& start, TLorentzVector& end, unsigned int &starti, unsigned int &endi);
      double bdist(const TVector3& pos);
      int    CountHits(const art::Event& evt, const art::InputTag& which, unsigned int cryostat, unsigned int tpc, unsigned int plane);

      TTree* fTree;
      //    TTree* fPOT;
      // event information is huge and dynamic;
      // run information is much smaller and we still store it statically
      // in the event
      std::unique_ptr<AnaRootParserDataStruct> fData;
      //    AnaRootParserDataStruct::RunData_t RunData;
      AnaRootParserDataStruct::SubRunData_t SubRunData;

      int fLogLevel;
      short fEventsPerSubrun;

      std::string fRawDigitModuleLabel;
      std::string fHitsModuleLabel;
      std::string fLArG4ModuleLabel;
      std::string fCalDataModuleLabel;
      std::string fGenieGenModuleLabel;
      std::string fCryGenModuleLabel;
      std::string fProtoGenModuleLabel;
      std::string fG4ModuleLabel;
      std::string fClusterModuleLabel;
      std::string fPandoraNuVertexModuleLabel;
      std::string fOpFlashModuleLabel;
      std::string fExternalCounterModuleLabel;
      std::string fMCShowerModuleLabel;
      std::string fMCTrackModuleLabel;
      std::vector<std::string> fTrackModuleLabel;
      std::string fPFParticleModuleLabel;
      std::vector<std::string> fVertexModuleLabel;
      std::vector<std::string> fShowerModuleLabel;
      std::vector<std::string> fCalorimetryModuleLabel;
      std::vector<std::string> fParticleIDModuleLabel;
      std::vector<std::string> fMVAPIDShowerModuleLabel;
      std::vector<std::string> fMVAPIDTrackModuleLabel;
      std::vector<std::string> fFlashT0FinderLabel;
      std::vector<std::string> fMCT0FinderLabel;
      std::string fPOTModuleLabel;
      std::string fCosmicClusterTaggerAssocLabel;
      bool fIsMC; ///< whether to use a permanent buffer (faster, huge memory)
      bool fUseBuffer; ///< whether to use a permanent buffer (faster, huge memory)

//      bool fSaveRecobWireInfo; //whether to extract and save recob::wire info
      bool fSaveAuxDetInfo; ///< whether to extract and save auxiliary detector data
      bool fSaveCryInfo; ///whether to extract and save CRY particle data
      bool fSaveGenieInfo; ///whether to extract and save Genie information
      bool fSaveProtoInfo; ///whether to extract and save ProtDUNE beam simulation information
      bool fSavePhotonInfo; ///whether to extract and save collected photons
      bool fSaveGeneratorInfo; ///whether to extract and save Geant information
      bool fSaveGeantInfo; ///whether to extract and save Geant information
      bool fSaveGeantInAVInfo; ///whether to extract and save Geant information
      bool fSaveGeantTrajectoryInfo; ///whether to extract and save Geant information
      bool fSaveMCShowerInfo; ///whether to extract and save MC Shower information
      bool fSaveMCTrackInfo; ///whether to extract and save MC Track information
      bool fSaveHitInfo; ///whether to extract and save Hit information
      bool fSaveRawDigitInfo; ///whether to extract and save raw digit information
      bool fSaveRecobWireInfo; ///whether to extract and save recob wire information
      bool fSaveTrackInfo; ///whether to extract and save Track information
      bool fSaveVertexInfo; ///whether to extract and save Vertex information
      bool fSaveClusterInfo;  ///whether to extract and save Cluster information
      bool fSavePandoraNuVertexInfo; ///whether to extract and save nu vertex information from Pandora
      bool fSaveFlashInfo;  ///whether to extract and save Flash information
      bool fSaveExternCounterInfo;  ///whether to extract and save External Counter information
      bool fSaveShowerInfo;  ///whether to extract and save Shower information
      bool fSavePFParticleInfo; ///whether to extract and save PFParticle information

      std::vector<std::string> fCosmicTaggerAssocLabel;
      std::vector<std::string> fContainmentTaggerAssocLabel;
      std::vector<std::string> fFlashMatchAssocLabel;

      bool bIgnoreMissingShowers; ///< whether to ignore missing shower information

      bool isCosmics;      ///< if it contains cosmics
      bool fSaveCaloCosmics; ///< save calorimetry information for cosmics
      float fG4minE;         ///< Energy threshold to save g4 particle info

      trkf::TrajectoryMCSFitter fMCSFitter; /// objet holding the configuration of the uBoone MCS fitting alg

      double ActiveBounds[6]; // Cryostat boundaries ( neg x, pos x, neg y, pos y, neg z, pos z )

      /// Returns the number of trackers configured
      size_t GetNTrackers() const { return fTrackModuleLabel.size(); }

      size_t GetNVertexAlgos() const { return fVertexModuleLabel.size(); }

      /// Returns the number of shower algorithms configured
      size_t GetNShowerAlgos() const { return fShowerModuleLabel.size(); }

      /// Returns the name of configured shower algorithms (converted to string)
      std::vector<std::string> GetShowerAlgos() const
      { return { fShowerModuleLabel.begin(), fShowerModuleLabel.end() }; }

      /// Creates the structure for the tree data; optionally initializes it
      void CreateData(bool bClearData = false)
      {
        if (!fData) {
          fData.reset
            (new AnaRootParserDataStruct(GetNTrackers(), GetNVertexAlgos(), GetShowerAlgos()));
          fData->SetBits(AnaRootParserDataStruct::tdCry,    !fSaveCryInfo);
          fData->SetBits(AnaRootParserDataStruct::tdGenie,  !fSaveGenieInfo);
          fData->SetBits(AnaRootParserDataStruct::tdProto,  !fSaveProtoInfo);
          fData->SetBits(AnaRootParserDataStruct::tdPhotons,  !fSavePhotonInfo);
          fData->SetBits(AnaRootParserDataStruct::tdGenerator,  !fSaveGeneratorInfo);
          fData->SetBits(AnaRootParserDataStruct::tdGeant,  !fSaveGeantInfo);
          fData->SetBits(AnaRootParserDataStruct::tdGeantInAV,  !fSaveGeantInAVInfo);
          fData->SetBits(AnaRootParserDataStruct::tdGeantTrajectory,  !fSaveGeantTrajectoryInfo);
          fData->SetBits(AnaRootParserDataStruct::tdMCshwr, !fSaveMCShowerInfo);
          fData->SetBits(AnaRootParserDataStruct::tdMCtrk,  !fSaveMCTrackInfo);
          fData->SetBits(AnaRootParserDataStruct::tdHit,    !fSaveHitInfo);
          fData->SetBits(AnaRootParserDataStruct::tdRawDigit,    !fSaveRawDigitInfo);
          fData->SetBits(AnaRootParserDataStruct::tdRecobWire,    !fSaveRecobWireInfo);
          fData->SetBits(AnaRootParserDataStruct::tdFlash,  !fSaveFlashInfo);
          fData->SetBits(AnaRootParserDataStruct::tdCount,  !fSaveExternCounterInfo);
          fData->SetBits(AnaRootParserDataStruct::tdShower, !fSaveShowerInfo);
          fData->SetBits(AnaRootParserDataStruct::tdCluster,!fSaveClusterInfo);
          fData->SetBits(AnaRootParserDataStruct::tdPandoraNuVertex,!fSavePandoraNuVertexInfo);
          fData->SetBits(AnaRootParserDataStruct::tdTrack,  !fSaveTrackInfo);
          fData->SetBits(AnaRootParserDataStruct::tdVertex, !fSaveVertexInfo);
          fData->SetBits(AnaRootParserDataStruct::tdAuxDet, !fSaveAuxDetInfo);
          fData->SetBits(AnaRootParserDataStruct::tdPFParticle, !fSavePFParticleInfo);
        }
        else {
          fData->SetTrackers(GetNTrackers());
          fData->SetVertexAlgos(GetNVertexAlgos());
          fData->SetShowerAlgos(GetShowerAlgos());
          if (bClearData) fData->Clear();
        }
      } // CreateData()

      /// Sets the addresses of all the tree branches, creating the missing ones
      void SetAddresses()
      {
        CheckData(__func__); CheckTree(__func__);
        fData->SetAddresses
          (fTree, fTrackModuleLabel, fVertexModuleLabel, fShowerModuleLabel, isCosmics);
      } // SetAddresses()

      /// Sets the addresses of all the tree branches of the specified tracking algo,
      /// creating the missing ones
      void SetTrackerAddresses(size_t iTracker)
      {
        CheckData(__func__); CheckTree(__func__);
        if (iTracker >= fData->GetNTrackers()) {
          throw art::Exception(art::errors::LogicError)
            << "AnaRootParser::SetTrackerAddresses(): no tracker #" << iTracker
            << " (" << fData->GetNTrackers() << " available)";
        }
        fData->GetTrackerData(iTracker)
          .SetAddresses(fTree, fTrackModuleLabel[iTracker], isCosmics);
      } // SetTrackerAddresses()


      void SetVertexAddresses(size_t iVertexAlg)
      {
        CheckData(__func__); CheckTree(__func__);
        if (iVertexAlg >= fData->GetNVertexAlgos()) {
          throw art::Exception(art::errors::LogicError)
            << "AnaRootParser::SetVertexAddresses(): no vertex alg #" << iVertexAlg
            << " (" << fData->GetNVertexAlgos() << " available)";
        }
        fData->GetVertexData(iVertexAlg)
          .SetAddresses(fTree, fVertexModuleLabel[iVertexAlg], isCosmics);
      } // SetVertexAddresses()

      /// Sets the addresses of all the tree branches of the specified shower algo,
      /// creating the missing ones
      void SetShowerAddresses(size_t iShower)
      {
        CheckData(__func__); CheckTree(__func__);
        if (iShower >= fData->GetNShowerAlgos()) {
          throw art::Exception(art::errors::LogicError)
            << "AnaRootParser::SetShowerAddresses(): no shower algo #" << iShower
            << " (" << fData->GetNShowerAlgos() << " available)";
        }
        fData->GetShowerData(iShower).SetAddresses(fTree);
      } // SetShowerAddresses()

      /// Sets the addresses of the tree branch of the PFParticle,
      /// creating it if missing
      void SetPFParticleAddress()
      {
        CheckData(__func__); CheckTree(__func__);
        fData->GetPFParticleData().SetAddresses(fTree);
      } // SetPFParticleAddress()

      /// Create the output tree and the data structures, if needed
      void CreateTree(bool bClearData = false);

      /// Destroy the local buffers (existing branches will point to invalid address!)
      void DestroyData() { fData.reset(); }

      /// Helper function: throws if no data structure is available
      void CheckData(std::string caller) const
      {
        if (fData) return;
        throw art::Exception(art::errors::LogicError)
          << "AnaRootParser::" << caller << ": no data";
      } // CheckData()
      /// Helper function: throws if no tree is available
      void CheckTree(std::string caller) const
      {
        if (fTree) return;
        throw art::Exception(art::errors::LogicError)
          << "AnaRootParser::" << caller << ": no tree";
      } // CheckData()

      /// Stores the information of shower in slot iShower of showerData
      void FillShower(
          AnaRootParserDataStruct::ShowerDataStruct& showerData,
          size_t iShower, recob::Shower const& showers, const bool fSavePFParticleInfo,
          const std::map<Short_t, Short_t> &showerIDtoPFParticleIDMap
          ) const;

      /// Stores the information of all showers into showerData
      void FillShowers(
          AnaRootParserDataStruct::ShowerDataStruct& showerData,
          std::vector<recob::Shower> const& showers, const bool fSavePFParticleInfo,
          const std::map<Short_t, Short_t> &showerIDtoPFParticleIDMap
          ) const;

  }; // class dune::AnaRootParser
} // namespace dune


namespace { // local namespace
  /// Simple stringstream which empties its buffer on operator() call
  class AutoResettingStringSteam: public std::ostringstream {
    public:
      AutoResettingStringSteam& operator() () { str(""); return *this; }
  }; // class AutoResettingStringSteam

  /// Fills a sequence of TYPE elements
  template <typename ITER, typename TYPE>
    inline void FillWith(ITER from, ITER to, TYPE value)
    { std::fill(from, to, value); }

  /// Fills a sequence of TYPE elements
  template <typename ITER, typename TYPE>
    inline void FillWith(ITER from, size_t n, TYPE value)
    { std::fill(from, from + n, value); }

  /// Fills a container with begin()/end() interface
  template <typename CONT, typename V>
    inline void FillWith(CONT& data, const V& value)
    { FillWith(std::begin(data), std::end(data), value); }

} // local namespace


//------------------------------------------------------------------------------
//---  AnaRootParserDataStruct::TrackDataStruct
//---

void dune::AnaRootParserDataStruct::TrackDataStruct::Resize(size_t nTracks)
{
  MaxTracks = nTracks;

  trkId.resize(MaxTracks);
  NHitsPerTrack.resize(MaxTracks);
  trkncosmictags_tagger.resize(MaxTracks);
  trkcosmicscore_tagger.resize(MaxTracks);
  trkcosmictype_tagger.resize(MaxTracks);
  trkncosmictags_containmenttagger.resize(MaxTracks);
  trkcosmicscore_containmenttagger.resize(MaxTracks);
  trkcosmictype_containmenttagger.resize(MaxTracks);
  trkncosmictags_flashmatch.resize(MaxTracks);
  trkcosmicscore_flashmatch.resize(MaxTracks);
  trkcosmictype_flashmatch.resize(MaxTracks);
  trkstartx.resize(MaxTracks);
  trkstarty.resize(MaxTracks);
  trkstartz.resize(MaxTracks);
  trkstartd.resize(MaxTracks);
  trkendx.resize(MaxTracks);
  trkendy.resize(MaxTracks);
  trkendz.resize(MaxTracks);
  trkendd.resize(MaxTracks);
  trkflashT0.resize(MaxTracks);
  trktrueT0.resize(MaxTracks);
  trkstarttheta.resize(MaxTracks);
  trkstartphi.resize(MaxTracks);
  trkendtheta.resize(MaxTracks);
  trkendphi.resize(MaxTracks);
  trkstartdirectionx.resize(MaxTracks);
  trkstartdirectiony.resize(MaxTracks);
  trkstartdirectionz.resize(MaxTracks);
  trkenddirectionx.resize(MaxTracks);
  trkenddirectiony.resize(MaxTracks);
  trkenddirectionz.resize(MaxTracks);
  trkthetaxz.resize(MaxTracks);
  trkthetayz.resize(MaxTracks);
  trkmom.resize(MaxTracks);
  trkmomrange.resize(MaxTracks);
  trkmommschi2.resize(MaxTracks);
  trkmommsllhd.resize(MaxTracks);
  trkmommscmic.resize(MaxTracks);
  trkmommscfwd.resize(MaxTracks);
  trkmommscbwd.resize(MaxTracks);
  trkmommscllfwd.resize(MaxTracks);
  trkmommscllbwd.resize(MaxTracks);
  trkchi2PerNDF.resize(MaxTracks);
  trkNDF.resize(MaxTracks);
  trklen.resize(MaxTracks);
  trklenstraightline.resize(MaxTracks);
  trksvtxid.resize(MaxTracks);
  trkevtxid.resize(MaxTracks);
  // PID variables
  trkpidpdg.resize(MaxTracks);
  trkpidchi.resize(MaxTracks);
  trkpidchipr.resize(MaxTracks);
  trkpidchika.resize(MaxTracks);
  trkpidchipi.resize(MaxTracks);
  trkpidchimu.resize(MaxTracks);
  trkpidpida.resize(MaxTracks);
  trkpidbestplane.resize(MaxTracks);
  trkpidmvamu.resize(MaxTracks);
  trkpidmvae.resize(MaxTracks);
  trkpidmvapich.resize(MaxTracks);
  trkpidmvaphoton.resize(MaxTracks);
  trkpidmvapr.resize(MaxTracks);

  trkke.resize(MaxTracks);
  trkrange.resize(MaxTracks);
  trkidtruth.resize(MaxTracks);
  trkorigin.resize(MaxTracks);
  trkpdgtruth.resize(MaxTracks);
  trkefftruth.resize(MaxTracks);
  trkpurtruth.resize(MaxTracks);
  trkpurity.resize(MaxTracks);
  trkcompleteness.resize(MaxTracks);
  trkg4id.resize(MaxTracks);
  trkorig.resize(MaxTracks);
  trkpitchc.resize(MaxTracks);
  ntrkhitsperview.resize(MaxTracks);

  trkdedx.resize(MaxTracks);
  trkdqdx.resize(MaxTracks);
  trkresrg.resize(MaxTracks);
  trktpc.resize(MaxTracks);
  trkxyz.resize(MaxTracks);

  trkhasPFParticle.resize(MaxTracks);
  trkPFParticleID.resize(MaxTracks);

} // dune::AnaRootParserDataStruct::TrackDataStruct::Resize()

void dune::AnaRootParserDataStruct::TrackDataStruct::Clear() {
  Resize(MaxTracks);
  ntracks = 0;

  std::fill(trkdedx2, trkdedx2 + sizeof(trkdedx2)/sizeof(trkdedx2[0]), -999.);
  std::fill(trkdqdx2, trkdqdx2 + sizeof(trkdqdx2)/sizeof(trkdqdx2[0]), -999.);
  std::fill(trktpc2, trktpc2 + sizeof(trktpc2)/sizeof(trktpc2[0]), -999.);
  std::fill(trkx2, trkx2 + sizeof(trkx2)/sizeof(trkx2[0]), -999.);
  std::fill(trky2, trky2 + sizeof(trky2)/sizeof(trky2[0]), -999.);
  std::fill(trkz2, trkz2 + sizeof(trkz2)/sizeof(trkz2[0]), -999.);

  std::fill(hittrkx, hittrkx + sizeof(hittrkx)/sizeof(hittrkx[0]), -999.);
  std::fill(hittrky, hittrky + sizeof(hittrky)/sizeof(hittrky[0]), -999.);
  std::fill(hittrkz, hittrkz + sizeof(hittrkz)/sizeof(hittrkz[0]), -999.);

  std::fill(hittrklocaltrackdirectionx, hittrklocaltrackdirectionx + sizeof(hittrklocaltrackdirectionx)/sizeof(hittrklocaltrackdirectionx[0]), -999.);
  std::fill(hittrklocaltrackdirectiony, hittrklocaltrackdirectiony + sizeof(hittrklocaltrackdirectiony)/sizeof(hittrklocaltrackdirectiony[0]), -999.);
  std::fill(hittrklocaltrackdirectionz, hittrklocaltrackdirectionz + sizeof(hittrklocaltrackdirectionz)/sizeof(hittrklocaltrackdirectionz[0]), -999.);
  std::fill(hittrklocaltrackdirectiontheta, hittrklocaltrackdirectiontheta + sizeof(hittrklocaltrackdirectiontheta)/sizeof(hittrklocaltrackdirectiontheta[0]), -999.);
  std::fill(hittrklocaltrackdirectionphi, hittrklocaltrackdirectionphi + sizeof(hittrklocaltrackdirectionphi)/sizeof(hittrklocaltrackdirectionphi[0]), -999.);


  std::fill(hittrkpitchC, hittrkpitchC + sizeof(hittrkpitchC)/sizeof(hittrkpitchC[0]), -999.);

  std::fill(hittrkds, hittrkds + sizeof(hittrkds)/sizeof(hittrkds[0]), -999.);

  std::fill(hittrkchannel, hittrkchannel + sizeof(hittrkchannel)/sizeof(hittrkchannel[0]), -999);
  std::fill(hittrktpc, hittrktpc + sizeof(hittrktpc)/sizeof(hittrktpc[0]), -999);
  std::fill(hittrkview, hittrkview + sizeof(hittrkview)/sizeof(hittrkview[0]), -999);
  std::fill(hittrkwire, hittrkwire + sizeof(hittrkwire)/sizeof(hittrkwire[0]), -999);
  std::fill(hittrkpeakT, hittrkpeakT + sizeof(hittrkpeakT)/sizeof(hittrkpeakT[0]), -999.);
  std::fill(hittrkchargeintegral, hittrkchargeintegral + sizeof(hittrkchargeintegral)/sizeof(hittrkchargeintegral[0]), -999.);
  std::fill(hittrkph, hittrkph + sizeof(hittrkph)/sizeof(hittrkph[0]), -999.);
  std::fill(hittrkchargesum, hittrkchargesum + sizeof(hittrkchargesum)/sizeof(hittrkchargesum[0]), -999.);
  std::fill(hittrkstarT, hittrkstarT + sizeof(hittrkstarT)/sizeof(hittrkstarT[0]), -999.);
  std::fill(hittrkendT, hittrkendT + sizeof(hittrkendT)/sizeof(hittrkendT[0]), -999.);
  std::fill(hittrkrms, hittrkrms + sizeof(hittrkrms)/sizeof(hittrkrms[0]), -999.);
  std::fill(hittrkgoddnessofFit, hittrkgoddnessofFit + sizeof(hittrkgoddnessofFit)/sizeof(hittrkgoddnessofFit[0]), -999.);
  std::fill(hittrkmultiplicity, hittrkmultiplicity + sizeof(hittrkmultiplicity)/sizeof(hittrkmultiplicity[0]), -999);
  std::fill(hittrktrueID, hittrktrueID + sizeof(hittrktrueID)/sizeof(hittrktrueID[0]), -999);
  std::fill(hittrktrueEnergyMax, hittrktrueEnergyMax + sizeof(hittrktrueEnergyMax)/sizeof(hittrktrueEnergyMax[0]), -999);
  std::fill(hittrktrueEnergyFraction, hittrktrueEnergyFraction + sizeof(hittrktrueEnergyFraction)/sizeof(hittrktrueEnergyFraction[0]), -999);

  FillWith(trkId        , -999  );
  FillWith(NHitsPerTrack        , -999  );
  FillWith(trkncosmictags_tagger, -999  );
  FillWith(trkcosmicscore_tagger, -999.);
  FillWith(trkcosmictype_tagger, -999  );
  FillWith(trkncosmictags_containmenttagger, -999  );
  FillWith(trkcosmicscore_containmenttagger, -999.);
  FillWith(trkcosmictype_containmenttagger, -999  );
  FillWith(trkncosmictags_flashmatch, -999  );
  FillWith(trkcosmicscore_flashmatch, -999.);
  FillWith(trkcosmictype_flashmatch, -999  );
  FillWith(trkstartx    , -999.);
  FillWith(trkstarty    , -999.);
  FillWith(trkstartz    , -999.);
  FillWith(trkstartd    , -999.);
  FillWith(trkendx      , -999.);
  FillWith(trkendy      , -999.);
  FillWith(trkendz      , -999.);
  FillWith(trkendd      , -999.);
  FillWith(trkflashT0   , -999.);
  FillWith(trktrueT0    , -999.);
  FillWith(trkg4id      , -9999 );
  FillWith(trkpurity    , -999.);
  FillWith(trkcompleteness, -999.);
  FillWith(trkorig      , -9999 );
  FillWith(trkstarttheta     , -999.);
  FillWith(trkstartphi       , -999.);
  FillWith(trkendtheta     , -999.);
  FillWith(trkendphi       , -999.);
  FillWith(trkstartdirectionx, -999.);
  FillWith(trkstartdirectiony, -999.);
  FillWith(trkstartdirectionz, -999.);
  FillWith(trkenddirectionx  , -999.);
  FillWith(trkenddirectiony  , -999.);
  FillWith(trkenddirectionz  , -999.);
  FillWith(trkthetaxz   , -999.);
  FillWith(trkthetayz   , -999.);
  FillWith(trkmom       , -999.);
  FillWith(trkmomrange  , -999.);
  FillWith(trkmommschi2 , -999.);
  FillWith(trkmommsllhd , -999.);
  FillWith(trkmommscmic , -999.);
  FillWith(trkmommscfwd , -999.);
  FillWith(trkmommscbwd , -999.);
  FillWith(trkmommscllfwd , -999.);
  FillWith(trkmommscllbwd , -999.);
  FillWith(trkchi2PerNDF, -999.);
  FillWith(trkNDF, -999.);
  FillWith(trklen       , -999.);
  FillWith(trklenstraightline, -999.);
  FillWith(trksvtxid    , -1);
  FillWith(trkevtxid    , -1);
  FillWith(trkpidbestplane, -1);
  FillWith(trkpidmvamu  , -999.);
  FillWith(trkpidmvae   , -999.);
  FillWith(trkpidmvapich, -999.);
  FillWith(trkpidmvaphoton , -999.);
  FillWith(trkpidmvapr  , -999.);

  FillWith(trkhasPFParticle, -1);
  FillWith(trkPFParticleID , -1);

  for (size_t iTrk = 0; iTrk < MaxTracks; ++iTrk){

    // the following are BoxedArray's;
    // their iterators traverse all the array dimensions
    FillWith(trkke[iTrk]      , -999.);
    FillWith(trkrange[iTrk]   , -999.);
    FillWith(trkidtruth[iTrk] , -9999 );
    FillWith(trkorigin[iTrk]  , -1 );
    FillWith(trkpdgtruth[iTrk], -9999 );
    FillWith(trkefftruth[iTrk], -999.);
    FillWith(trkpurtruth[iTrk], -999.);
    FillWith(trkpitchc[iTrk]  , -999.);
    FillWith(ntrkhitsperview[iTrk]   ,  -999 );

    FillWith(trkdedx[iTrk], 0.);
    FillWith(trkdqdx[iTrk], 0.);
    FillWith(trkresrg[iTrk], 0.);
    FillWith(trktpc[iTrk], -1);
    FillWith(trkxyz[iTrk], 0.);

    FillWith(trkpidpdg[iTrk]    , -1);
    FillWith(trkpidchi[iTrk]    , -999.);
    FillWith(trkpidchipr[iTrk]  , -999.);
    FillWith(trkpidchika[iTrk]  , -999.);
    FillWith(trkpidchipi[iTrk]  , -999.);
    FillWith(trkpidchimu[iTrk]  , -999.);
    FillWith(trkpidpida[iTrk]   , -999.);
  } // for track

} // dune::AnaRootParserDataStruct::TrackDataStruct::Clear()


void dune::AnaRootParserDataStruct::TrackDataStruct::SetAddresses(
    TTree* pTree, std::string tracker, bool isCosmics
    ) {
  if (MaxTracks == 0) return; // no tracks, no tree!

  dune::AnaRootParserDataStruct::BranchCreator CreateBranch(pTree);

  AutoResettingStringSteam sstr;
  sstr() << kMaxTrackHits;
  std::string MaxTrackHitsIndexStr("[" + sstr.str() + "]");

  std::string TrackLabel = tracker;
  std::string BranchName;

  BranchName = "NumberOfTracks";
  CreateBranch(BranchName, &ntracks, BranchName + "/S");
  std::string NTracksIndexStr = "[" + BranchName + "]";

  BranchName = "TrackID";
  CreateBranch(BranchName, trkId, BranchName + NTracksIndexStr + "/S");

  BranchName = "Track_NumberOfHits";
  CreateBranch(BranchName, NHitsPerTrack, BranchName + NTracksIndexStr + "/S");

  BranchName = "Track_Length_Trajectory";
  CreateBranch(BranchName, trklen, BranchName + NTracksIndexStr + "/F");

  BranchName = "Track_Length_StraightLine";
  CreateBranch(BranchName, trklenstraightline, BranchName + NTracksIndexStr + "/F");

  BranchName = "Track_Chi2PerNDF";
  CreateBranch(BranchName, trkchi2PerNDF, BranchName + NTracksIndexStr + "/F");

  BranchName = "Track_NDF";
  CreateBranch(BranchName, trkNDF, BranchName + NTracksIndexStr + "/F");

  BranchName = "Track_StartPoint_X";
  CreateBranch(BranchName, trkstartx, BranchName + NTracksIndexStr + "/F");

  BranchName = "Track_StartPoint_Y";
  CreateBranch(BranchName, trkstarty, BranchName + NTracksIndexStr + "/F");

  BranchName = "Track_StartPoint_Z";
  CreateBranch(BranchName, trkstartz, BranchName + NTracksIndexStr + "/F");

  BranchName = "Track_StartPoint_DistanceToBoundary";
  CreateBranch(BranchName, trkstartd, BranchName + NTracksIndexStr + "/F");

  BranchName = "Track_EndPoint_X";
  CreateBranch(BranchName, trkendx, BranchName + NTracksIndexStr + "/F");

  BranchName = "Track_EndPoint_Y";
  CreateBranch(BranchName, trkendy, BranchName + NTracksIndexStr + "/F");

  BranchName = "Track_EndPoint_Z";
  CreateBranch(BranchName, trkendz, BranchName + NTracksIndexStr + "/F");

  BranchName = "Track_EndPoint_DistanceToBoundary";
  CreateBranch(BranchName, trkendd, BranchName + NTracksIndexStr + "/F");

  BranchName = "Track_StartDirection_Theta";
  CreateBranch(BranchName, trkstarttheta, BranchName + NTracksIndexStr + "/F");

  BranchName = "Track_StartDirection_Phi";
  CreateBranch(BranchName, trkstartphi, BranchName + NTracksIndexStr + "/F");

  BranchName = "Track_StartDirection_X";
  CreateBranch(BranchName, trkstartdirectionx, BranchName + NTracksIndexStr + "/F");

  BranchName = "Track_StartDirection_Y";
  CreateBranch(BranchName, trkstartdirectiony, BranchName + NTracksIndexStr + "/F");

  BranchName = "Track_StartDirection_Z";
  CreateBranch(BranchName, trkstartdirectionz, BranchName + NTracksIndexStr + "/F");

  BranchName = "Track_EndDirection_Theta";
  CreateBranch(BranchName, trkendtheta, BranchName + NTracksIndexStr + "/F");

  BranchName = "Track_EndDirection_Phi";
  CreateBranch(BranchName, trkendphi, BranchName + NTracksIndexStr + "/F");

  BranchName = "Track_EndDirection_X";
  CreateBranch(BranchName, trkenddirectionx, BranchName + NTracksIndexStr + "/F");

  BranchName = "Track_EndDirection_Y";
  CreateBranch(BranchName, trkenddirectiony, BranchName + NTracksIndexStr + "/F");

  BranchName = "Track_EndDirection_Z";
  CreateBranch(BranchName, trkenddirectionz, BranchName + NTracksIndexStr + "/F");



  //Branches for calculated momentum
  BranchName = "Track_Momentum";
  CreateBranch(BranchName, trkmom, BranchName + NTracksIndexStr + "/F");

  BranchName = "Track_Momentum_Range";
  CreateBranch(BranchName, trkmomrange, BranchName + NTracksIndexStr + "/F");

  BranchName = "Track_Momentum_mschi";
  CreateBranch(BranchName, trkmommschi2, BranchName + NTracksIndexStr + "/F");

  BranchName = "Track_Momentum_mscmic";
  CreateBranch(BranchName, trkmommscmic, BranchName + NTracksIndexStr + "/F");

  BranchName = "Track_Momentum_mscfwd";
  CreateBranch(BranchName, trkmommscfwd, BranchName + NTracksIndexStr + "/F");

  BranchName = "Track_Momentum_mscbwd";
  CreateBranch(BranchName, trkmommscbwd, BranchName + NTracksIndexStr + "/F");

  BranchName = "Track_Momentum_mscllfwd";
  CreateBranch(BranchName, trkmommscllfwd, BranchName + NTracksIndexStr + "/F");

  BranchName = "Track_Momentum_mscllbwd";
  CreateBranch(BranchName, trkmommscllbwd, BranchName + NTracksIndexStr + "/F");


  /*
     BranchName = "trkncosmictags_tagger";
     CreateBranch(BranchName, trkncosmictags_tagger, BranchName + NTracksIndexStr + "/S");

     BranchName = "trkcosmicscore_tagger";
     CreateBranch(BranchName, trkcosmicscore_tagger, BranchName + NTracksIndexStr + "/F");

     BranchName = "trkcosmictype_tagger";
     CreateBranch(BranchName, trkcosmictype_tagger, BranchName + NTracksIndexStr + "/S");

     BranchName = "trkncosmictags_containmenttagger";
     CreateBranch(BranchName, trkncosmictags_containmenttagger, BranchName + NTracksIndexStr + "/S");

     BranchName = "trkcosmicscore_containmenttagger";
     CreateBranch(BranchName, trkcosmicscore_containmenttagger, BranchName + NTracksIndexStr + "/F");

     BranchName = "trkcosmictype_containmenttagger";
     CreateBranch(BranchName, trkcosmictype_containmenttagger, BranchName + NTracksIndexStr + "/S");

     BranchName = "trkncosmictags_flashmatch";
     CreateBranch(BranchName, trkncosmictags_flashmatch, BranchName + NTracksIndexStr + "/S");

     BranchName = "trkcosmicscore_flashmatch";
     CreateBranch(BranchName, trkcosmicscore_flashmatch, BranchName + NTracksIndexStr + "/F");

     BranchName = "trkcosmictype_flashmatch";
     CreateBranch(BranchName, trkcosmictype_flashmatch, BranchName + NTracksIndexStr + "/S");
     */
  /*
     BranchName = "Track_KineticEnergy";
     CreateBranch(BranchName, trkke, BranchName + NTracksIndexStr + "[2]/F");

     BranchName = "Track_Range";
     CreateBranch(BranchName, trkrange, BranchName + NTracksIndexStr + "[2]/F");
     */
  /*
     BranchName = "trkidtruth";
     CreateBranch(BranchName, trkidtruth, BranchName + NTracksIndexStr + "[2]/I");

     BranchName = "trkorigin";
     CreateBranch(BranchName, trkorigin, BranchName + NTracksIndexStr + "[2]/S");

     BranchName = "trkpdgtruth";
     CreateBranch(BranchName, trkpdgtruth, BranchName + NTracksIndexStr + "[2]/I");

     BranchName = "trkefftruth";
     CreateBranch(BranchName, trkefftruth, BranchName + NTracksIndexStr + "[2]/F");

     BranchName = "trkpurtruth";
     CreateBranch(BranchName, trkpurtruth, BranchName + NTracksIndexStr + "[2]/F");
     */
  BranchName = "Track_PitchInViews";
  CreateBranch(BranchName, trkpitchc, BranchName + NTracksIndexStr + "[2]/F");

  BranchName = "Track_NumberOfHitsPerView";
  CreateBranch(BranchName, ntrkhitsperview, BranchName + NTracksIndexStr + "[2]/S");


  /*
     BranchName = "trkresrg";
     CreateBranch(BranchName, trkresrg, BranchName + "[no_hits]/F");
     */
/*
  BranchName = "Track_TPC";
  CreateBranch(BranchName, trktpc2, BranchName + "[no_hits]/F");

  BranchName = "Track_dEdx";
  CreateBranch(BranchName, trkdedx2, BranchName + "[no_hits]/F");

  BranchName = "Track_dQdx";
  CreateBranch(BranchName, trkdqdx2, BranchName + "[no_hits]/F");

  BranchName = "Track_HitX";
  CreateBranch(BranchName, trkx2, BranchName + "[no_hits]/F");

  BranchName = "Track_HitY";
  CreateBranch(BranchName, trky2, BranchName + "[no_hits]/F");

  BranchName = "Track_HitZ";
  CreateBranch(BranchName, trkz2, BranchName + "[no_hits]/F");
*/
//////////// new arrays ///////////////////////////

  BranchName = "Track_Hit_X";
  CreateBranch(BranchName, hittrkx, BranchName + "[no_hits]/F");

  BranchName = "Track_Hit_Y";
  CreateBranch(BranchName, hittrky, BranchName + "[no_hits]/F");

  BranchName = "Track_Hit_Z";
  CreateBranch(BranchName, hittrkz, BranchName + "[no_hits]/F");

  BranchName = "Track_Hit_LocalTrackDirection_X";
  CreateBranch(BranchName, hittrklocaltrackdirectionx, BranchName + "[no_hits]/F");

  BranchName = "Track_Hit_LocalTrackDirection_Y";
  CreateBranch(BranchName, hittrklocaltrackdirectiony, BranchName + "[no_hits]/F");

  BranchName = "Track_Hit_LocalTrackDirection_Z";
  CreateBranch(BranchName, hittrklocaltrackdirectionz, BranchName + "[no_hits]/F");

  BranchName = "Track_Hit_LocalTrackDirection_Theta";
  CreateBranch(BranchName, hittrklocaltrackdirectiontheta, BranchName + "[no_hits]/F");

  BranchName = "Track_Hit_LocalTrackDirection_Phi";
  CreateBranch(BranchName, hittrklocaltrackdirectionphi, BranchName + "[no_hits]/F");

  BranchName = "Track_Hit_ds_LocalTrackDirection";
  CreateBranch(BranchName, hittrkpitchC, BranchName + "[no_hits]/F");

  BranchName = "Track_Hit_ds_3DPosition";
  CreateBranch(BranchName, hittrkds, BranchName + "[no_hits]/F");

  BranchName = "Track_Hit_TPC";
  CreateBranch(BranchName, hittrktpc, BranchName + "[no_hits]/S");

  BranchName = "Track_Hit_View";
  CreateBranch(BranchName, hittrkview, BranchName + "[no_hits]/S");

  BranchName = "Track_Hit_Channel";
  CreateBranch(BranchName, hittrkchannel, BranchName + "[no_hits]/S");

//  BranchName = "Track_Hit_Wire";
//  CreateBranch(BranchName, hittrkwire, BranchName + "[no_hits]/S");

  BranchName = "Track_Hit_PeakTime";
  CreateBranch(BranchName, hittrkpeakT, BranchName + "[no_hits]/F");

  BranchName = "Track_Hit_ChargeSummedADC";
  CreateBranch(BranchName, hittrkchargesum, BranchName + "[no_hits]/F");

  BranchName = "Track_Hit_ChargeIntegral";
  CreateBranch(BranchName, hittrkchargeintegral, BranchName + "[no_hits]/F");

  BranchName = "Track_Hit_Amplitude";
  CreateBranch(BranchName, hittrkph, BranchName + "[no_hits]/F");

  BranchName = "Track_Hit_StartTime";
  CreateBranch(BranchName, hittrkstarT, BranchName + "[no_hits]/F");

  BranchName = "Track_Hit_EndTime";
  CreateBranch(BranchName, hittrkendT, BranchName + "[no_hits]/F");

  BranchName = "Track_Hit_Width";
  CreateBranch(BranchName, hittrkrms, BranchName + "[no_hits]/F");

  BranchName = "Track_Hit_GoodnessOfFit";
  CreateBranch(BranchName, hittrkgoddnessofFit, BranchName + "[no_hits]/F");

  BranchName = "Track_Hit_Multiplicity";
  CreateBranch(BranchName, hittrkmultiplicity, BranchName + "[no_hits]/S");

  BranchName = "Track_Hit_trueID";
  CreateBranch(BranchName, hittrktrueID, BranchName + "[no_hits]/I");

  BranchName = "Track_Hit_trueEnergyMax";
  CreateBranch(BranchName, hittrktrueEnergyMax, BranchName + "[no_hits]/F");

  BranchName = "Track_Hit_trueEnergyFraction";
  CreateBranch(BranchName, hittrktrueEnergyFraction, BranchName + "[no_hits]/F");

//////////// end new arrays ///////////////////////////



  /*
     if (!isCosmics){
     BranchName = "Track_dEdx";
     CreateBranch(BranchName, trkdedx, BranchName + NTracksIndexStr + "[2]" + MaxTrackHitsIndexStr + "/F");

     BranchName = "Track_dQdx";
     CreateBranch(BranchName, trkdqdx, BranchName + NTracksIndexStr + "[2]" + MaxTrackHitsIndexStr + "/F");

     BranchName = "trkresrg";
     CreateBranch(BranchName, trkresrg, BranchName + NTracksIndexStr + "[2]" + MaxTrackHitsIndexStr + "/F");

     BranchName = "Track_TPC";
     CreateBranch(BranchName, trktpc, BranchName + NTracksIndexStr + "[2]" + MaxTrackHitsIndexStr + "/I");

     BranchName = "Track_XYZ";
     CreateBranch(BranchName, trkxyz, BranchName + NTracksIndexStr + "[2]" + MaxTrackHitsIndexStr + "[3]" + "/F");
     }
     */

  /*
     BranchName = "trkflashT0";
     CreateBranch(BranchName, trkflashT0, BranchName + NTracksIndexStr + "/F");

     BranchName = "trktrueT0";
     CreateBranch(BranchName, trktrueT0, BranchName + NTracksIndexStr + "/F");

     BranchName = "trkg4id";
     CreateBranch(BranchName, trkg4id, BranchName + NTracksIndexStr + "/I");

     BranchName = "trkorig";
     CreateBranch(BranchName, trkorig, BranchName + NTracksIndexStr + "/I");

     BranchName = "trkpurity";
     CreateBranch(BranchName, trkpurity, BranchName + NTracksIndexStr + "/F");

     BranchName = "trkcompleteness";
     CreateBranch(BranchName, trkcompleteness, BranchName + NTracksIndexStr + "/F");
     */

  /*


     BranchName = "trkthetaxz";
     CreateBranch(BranchName, trkthetaxz, BranchName + NTracksIndexStr + "/F");

     BranchName = "trkthetayz";
     CreateBranch(BranchName, trkthetayz, BranchName + NTracksIndexStr + "/F");

     BranchName = "trkmom";
     CreateBranch(BranchName, trkmom, BranchName + NTracksIndexStr + "/F");

     BranchName = "trkmomrange";
     CreateBranch(BranchName, trkmomrange, BranchName + NTracksIndexStr + "/F");
     */
  /*
     BranchName = "trkmommschi2";
     CreateBranch(BranchName, trkmommschi2, BranchName + NTracksIndexStr + "/F");

     BranchName = "trkmommsllhd";
     CreateBranch(BranchName, trkmommsllhd, BranchName + NTracksIndexStr + "/F");
     */

  /*
     BranchName = "trksvtxid";
     CreateBranch(BranchName, trksvtxid, BranchName + NTracksIndexStr + "/S");

     BranchName = "trkevtxid";
     CreateBranch(BranchName, trkevtxid, BranchName + NTracksIndexStr + "/S");

     BranchName = "trkpidmvamu";
     CreateBranch(BranchName, trkpidmvamu, BranchName + NTracksIndexStr + "/F");

     BranchName = "trkpidmvae";
     CreateBranch(BranchName, trkpidmvae, BranchName + NTracksIndexStr + "/F");

     BranchName = "trkpidmvapich";
     CreateBranch(BranchName, trkpidmvapich, BranchName + NTracksIndexStr + "/F");

     BranchName = "trkpidmvaphoton";
     CreateBranch(BranchName, trkpidmvaphoton, BranchName + NTracksIndexStr + "/F");

     BranchName = "trkpidmvapr";
     CreateBranch(BranchName, trkpidmvapr, BranchName + NTracksIndexStr + "/F");

     BranchName = "trkpidpdg";
     CreateBranch(BranchName, trkpidpdg, BranchName + NTracksIndexStr + "[2]/I");

     BranchName = "trkpidchi";
     CreateBranch(BranchName, trkpidchi, BranchName + NTracksIndexStr + "[2]/F");

     BranchName = "trkpidchipr";
     CreateBranch(BranchName, trkpidchipr, BranchName + NTracksIndexStr + "[2]/F");

     BranchName = "trkpidchika";
     CreateBranch(BranchName, trkpidchika, BranchName + NTracksIndexStr + "[2]/F");

     BranchName = "trkpidchipi";
     CreateBranch(BranchName, trkpidchipi, BranchName + NTracksIndexStr + "[2]/F");

     BranchName = "trkpidchimu";
     CreateBranch(BranchName, trkpidchimu, BranchName + NTracksIndexStr + "[2]/F");

     BranchName = "trkpidpida";
     CreateBranch(BranchName, trkpidpida, BranchName + NTracksIndexStr + "[2]/F");

     BranchName = "trkpidbestview";
     CreateBranch(BranchName, trkpidbestplane, BranchName + NTracksIndexStr + "/S");

     BranchName = "trkhasPFParticle";
     CreateBranch(BranchName, trkhasPFParticle, BranchName + NTracksIndexStr + "/S");

     BranchName = "trkPFParticleID";
     CreateBranch(BranchName, trkPFParticleID, BranchName + NTracksIndexStr + "/S");
     */
} // dune::AnaRootParserDataStruct::TrackDataStruct::SetAddresses()



void dune::AnaRootParserDataStruct::VertexDataStruct::Resize(size_t nVertices)
{
  MaxVertices = nVertices;
  vtxId.resize(MaxVertices);
  vtxx.resize(MaxVertices);
  vtxy.resize(MaxVertices);
  vtxz.resize(MaxVertices);

  vtxhasPFParticle.resize(MaxVertices);
  vtxPFParticleID.resize(MaxVertices);
}

void dune::AnaRootParserDataStruct::VertexDataStruct::Clear() {
  Resize(MaxVertices);
  nvtx = -999;

  FillWith(vtxId       , -999  );
  FillWith(vtxx        , -999  );
  FillWith(vtxy        , -999  );
  FillWith(vtxz        , -999  );
  FillWith(vtxhasPFParticle, -1  );
  FillWith(vtxPFParticleID , -1  );
}

void dune::AnaRootParserDataStruct::VertexDataStruct::SetAddresses(
    TTree* pTree, std::string alg, bool isCosmics
    ) {
  if (MaxVertices == 0) return; // no tracks, no tree!

  dune::AnaRootParserDataStruct::BranchCreator CreateBranch(pTree);

  AutoResettingStringSteam sstr;

  std::string VertexLabel = alg;
  std::string BranchName;
  /*
     BranchName = "nvtx_" + VertexLabel;
     CreateBranch(BranchName, &nvtx, BranchName + "/S");
     std::string NVertexIndexStr = "[" + BranchName + "]";

     BranchName = "vtxId_" + VertexLabel;
     CreateBranch(BranchName, vtxId, BranchName + NVertexIndexStr + "/S");

     BranchName = "vtxx_" + VertexLabel;
     CreateBranch(BranchName, vtxx, BranchName + NVertexIndexStr + "/F");

     BranchName = "vtxy_" + VertexLabel;
     CreateBranch(BranchName, vtxy, BranchName + NVertexIndexStr + "/F");

     BranchName = "vtxz_" + VertexLabel;
     CreateBranch(BranchName, vtxz, BranchName + NVertexIndexStr + "/F");

     BranchName = "vtxhasPFParticle_" + VertexLabel;
     CreateBranch(BranchName, vtxhasPFParticle, BranchName + NVertexIndexStr + "/S");

     BranchName = "vtxPFParticleID_" + VertexLabel;
     CreateBranch(BranchName, vtxPFParticleID, BranchName + NVertexIndexStr + "/S");
     */
}

//------------------------------------------------------------------------------
//---  AnaRootParserDataStruct::PFParticleDataStruct
//---

void dune::AnaRootParserDataStruct::PFParticleDataStruct::Resize(size_t nPFParticles)
{
  MaxPFParticles = nPFParticles;

  pfp_selfID.resize(MaxPFParticles);
  pfp_isPrimary.resize(MaxPFParticles);
  pfp_numDaughters.resize(MaxPFParticles);
  pfp_daughterIDs.resize(MaxPFParticles);
  pfp_parentID.resize(MaxPFParticles);
  pfp_vertexID.resize(MaxPFParticles);
  pfp_isShower.resize(MaxPFParticles);
  pfp_isTrack.resize(MaxPFParticles);
  pfp_trackID.resize(MaxPFParticles);
  pfp_showerID.resize(MaxPFParticles);
  pfp_pdgCode.resize(MaxPFParticles);
  pfp_numClusters.resize(MaxPFParticles);
  pfp_clusterIDs.resize(MaxPFParticles);
  pfp_isNeutrino.resize(MaxPFParticles);
}

void dune::AnaRootParserDataStruct::PFParticleDataStruct::Clear() {
  Resize(MaxPFParticles);

  nPFParticles = -999;
  FillWith(pfp_selfID, -999);
  FillWith(pfp_isPrimary, -999);
  FillWith(pfp_numDaughters, -999);
  FillWith(pfp_parentID, -999);
  FillWith(pfp_vertexID, -999);
  FillWith(pfp_isShower, -999);
  FillWith(pfp_isTrack, -999);
  FillWith(pfp_trackID, -999);
  FillWith(pfp_showerID, -999);
  FillWith(pfp_pdgCode, -999);
  FillWith(pfp_isNeutrino, -999);
  pfp_numNeutrinos = -999;
  FillWith(pfp_neutrinoIDs, -999);

  for (size_t iPFParticle = 0; iPFParticle < MaxPFParticles; ++iPFParticle){
    // the following are BoxedArrays;
    // their iterators traverse all the array dimensions
    FillWith(pfp_daughterIDs[iPFParticle], -999);
    FillWith(pfp_clusterIDs[iPFParticle], -999);
  }
}

void dune::AnaRootParserDataStruct::PFParticleDataStruct::SetAddresses(
    TTree* pTree
    ) {

  if (MaxPFParticles == 0) { return; } // no PFParticles, no tree

  dune::AnaRootParserDataStruct::BranchCreator CreateBranch(pTree);

  AutoResettingStringSteam sstr;
  sstr() << kMaxNDaughtersPerPFP;
  std::string MaxNDaughtersIndexStr("[" + sstr.str() + "]");

  sstr.str("");
  sstr() << kMaxNClustersPerPFP;
  std::string MaxNClustersIndexStr("[" + sstr.str() + "]");

  sstr.str("");
  sstr() << kMaxNPFPNeutrinos;
  std::string MaxNNeutrinosIndexStr("[" + sstr.str() + "]");

  std::string BranchName;

  BranchName = "nPFParticles";
  CreateBranch(BranchName, &nPFParticles, BranchName + "/S");
  std::string NPFParticleIndexStr = "[" + BranchName + "]";

  BranchName = "pfp_selfID";
  CreateBranch(BranchName, pfp_selfID, BranchName + NPFParticleIndexStr + "/S");

  BranchName = "pfp_isPrimary";
  CreateBranch(BranchName, pfp_isPrimary, BranchName + NPFParticleIndexStr + "/S");

  BranchName = "pfp_numDaughters";
  CreateBranch(BranchName, pfp_numDaughters, BranchName + NPFParticleIndexStr + "/S");

  BranchName = "pfp_daughterIDs";
  CreateBranch(BranchName, pfp_daughterIDs, BranchName + NPFParticleIndexStr + MaxNDaughtersIndexStr + "/S");

  BranchName = "pfp_parentID";
  CreateBranch(BranchName, pfp_parentID, BranchName + NPFParticleIndexStr + "/S");

  BranchName = "pfp_vertexID";
  CreateBranch(BranchName, pfp_vertexID, BranchName + NPFParticleIndexStr + "/S");

  BranchName = "pfp_isShower";
  CreateBranch(BranchName, pfp_isShower, BranchName + NPFParticleIndexStr + "/S");

  BranchName = "pfp_isTrack";
  CreateBranch(BranchName, pfp_isTrack, BranchName + NPFParticleIndexStr + "/S");

  BranchName = "pfp_trackID";
  CreateBranch(BranchName, pfp_trackID, BranchName + NPFParticleIndexStr + "/S");

  BranchName = "pfp_showerID";
  CreateBranch(BranchName, pfp_showerID, BranchName + NPFParticleIndexStr + "/S");

  BranchName = "pfp_pdgCode";
  CreateBranch(BranchName, pfp_pdgCode, BranchName + NPFParticleIndexStr + "/I");

  BranchName = "pfp_numClusters";
  CreateBranch(BranchName, pfp_numClusters, BranchName + NPFParticleIndexStr + "/S");

  BranchName = "pfp_clusterIDs";
  CreateBranch(BranchName, pfp_clusterIDs, BranchName + NPFParticleIndexStr + MaxNClustersIndexStr + "/S");

  BranchName = "pfp_isNeutrino";
  CreateBranch(BranchName, pfp_isNeutrino, BranchName + NPFParticleIndexStr + "/S");

  BranchName = "pfp_numNeutrinos";
  CreateBranch(BranchName, &pfp_numNeutrinos, BranchName + "/S");

  BranchName = "pfp_neutrinoIDs";
  CreateBranch(BranchName, pfp_neutrinoIDs, BranchName + MaxNNeutrinosIndexStr + "/S");
}

//------------------------------------------------------------------------------
//---  AnaRootParserDataStruct::ShowerDataStruct
//---

  void dune::AnaRootParserDataStruct::ShowerDataStruct::Resize
(size_t nShowers)
{
  MaxShowers = nShowers;

  showerID.resize(MaxShowers);
  shwr_bestplane.resize(MaxShowers);
  shwr_length.resize(MaxShowers);
  shwr_startdcosx.resize(MaxShowers);
  shwr_startdcosy.resize(MaxShowers);
  shwr_startdcosz.resize(MaxShowers);
  shwr_startx.resize(MaxShowers);
  shwr_starty.resize(MaxShowers);
  shwr_startz.resize(MaxShowers);
  shwr_totEng.resize(MaxShowers);
  shwr_dedx.resize(MaxShowers);
  shwr_mipEng.resize(MaxShowers);
  shwr_pidmvamu.resize(MaxShowers);
  shwr_pidmvae.resize(MaxShowers);
  shwr_pidmvapich.resize(MaxShowers);
  shwr_pidmvaphoton.resize(MaxShowers);
  shwr_pidmvapr.resize(MaxShowers);

  shwr_hasPFParticle.resize(MaxShowers);
  shwr_PFParticleID.resize(MaxShowers);

} // dune::AnaRootParserDataStruct::ShowerDataStruct::Resize()

void dune::AnaRootParserDataStruct::ShowerDataStruct::Clear() {
  Resize(MaxShowers);
  nshowers = 0;

  FillWith(showerID,         -999 );
  FillWith(shwr_bestplane,   -999 );
  FillWith(shwr_length,     -999.);
  FillWith(shwr_startdcosx, -999.);
  FillWith(shwr_startdcosy, -999.);
  FillWith(shwr_startdcosz, -999.);
  FillWith(shwr_startx,     -999.);
  FillWith(shwr_starty,     -999.);
  FillWith(shwr_startz,     -999.);
  FillWith(shwr_pidmvamu,   -999.);
  FillWith(shwr_pidmvae,    -999.);
  FillWith(shwr_pidmvapich, -999.);
  FillWith(shwr_pidmvaphoton,  -999.);
  FillWith(shwr_pidmvapr,   -999.);

  FillWith(shwr_hasPFParticle, -1);
  FillWith(shwr_PFParticleID,  -1);

  for (size_t iShw = 0; iShw < MaxShowers; ++iShw){
    // the following are BoxedArray's;
    // their iterators traverse all the array dimensions
    FillWith(shwr_totEng[iShw], -999.);
    FillWith(shwr_dedx[iShw],   -999.);
    FillWith(shwr_mipEng[iShw], -999.);
  } // for shower

} // dune::AnaRootParserDataStruct::ShowerDataStruct::Clear()


  void dune::AnaRootParserDataStruct::ShowerDataStruct::MarkMissing
(TTree* pTree)
{
  // here we implement the policy prescription for a missing set of showers;
  // this means that no shower data product was found in the event,
  // yet the user has accepted to go on.
  // We now need to mark this product in a unmistakably clear way, so that it
  // is not confused with a valid collection of an event where no showers
  // were reconstructed, not as a list of valid showers.
  // The prescription currently implemented is:
  // - have only one shower in the list;
  // - set the ID of that shower as -999
  //

  // first set the data structures to contain one invalid shower:
  SetMaxShowers(1); // includes resize to a set of one shower
  Clear(); // initializes all the showers in the set (one) as invalid
  // now set the tree addresses to the newly allocated memory;
  // this creates the tree branches in case they are not there yet
  SetAddresses(pTree);

  // then, set the variables so that ROOT tree knows there is one shower only
  nshowers = 1;

} // dune::AnaRootParserDataStruct::ShowerDataStruct::MarkMissing()


  void dune::AnaRootParserDataStruct::ShowerDataStruct::SetAddresses
(TTree* pTree)
{
  if (MaxShowers == 0) return; // no showers, no tree!

  dune::AnaRootParserDataStruct::BranchCreator CreateBranch(pTree);

  AutoResettingStringSteam sstr;
  sstr() << kMaxShowerHits;
  std::string MaxShowerHitsIndexStr("[" + sstr.str() + "]");

  std::string ShowerLabel = Name();
  std::string BranchName;

  BranchName = "nshowers_" + ShowerLabel;
  CreateBranch(BranchName, &nshowers, BranchName + "/S");
  std::string NShowerIndexStr = "[" + BranchName + "]";

  BranchName = "showerID_" + ShowerLabel;
  CreateBranch(BranchName, showerID, BranchName + NShowerIndexStr + "/S");

  BranchName = "shwr_bestplane_" + ShowerLabel;
  CreateBranch(BranchName, shwr_bestplane, BranchName + NShowerIndexStr + "/S");

  BranchName = "shwr_length_" + ShowerLabel;
  CreateBranch(BranchName, shwr_length, BranchName + NShowerIndexStr + "/F");

  BranchName = "shwr_startdcosx_" + ShowerLabel;
  CreateBranch(BranchName, shwr_startdcosx, BranchName + NShowerIndexStr + "/F");

  BranchName = "shwr_startdcosy_" + ShowerLabel;
  CreateBranch(BranchName, shwr_startdcosy, BranchName + NShowerIndexStr + "/F");

  BranchName = "shwr_startdcosz_" + ShowerLabel;
  CreateBranch(BranchName, shwr_startdcosz, BranchName + NShowerIndexStr + "/F");

  BranchName = "shwr_startx_" + ShowerLabel;
  CreateBranch(BranchName, shwr_startx, BranchName + NShowerIndexStr + "/F");

  BranchName = "shwr_starty_" + ShowerLabel;
  CreateBranch(BranchName, shwr_starty, BranchName + NShowerIndexStr + "/F");

  BranchName = "shwr_startz_" + ShowerLabel;
  CreateBranch(BranchName, shwr_startz, BranchName + NShowerIndexStr + "/F");

  BranchName = "shwr_totEng_" + ShowerLabel;
  CreateBranch(BranchName, shwr_totEng, BranchName + NShowerIndexStr + "[3]/F");

  BranchName = "shwr_dedx_" + ShowerLabel;
  CreateBranch(BranchName, shwr_dedx, BranchName + NShowerIndexStr + "[3]/F");

  BranchName = "shwr_mipEng_" + ShowerLabel;
  CreateBranch(BranchName, shwr_mipEng, BranchName + NShowerIndexStr + "[3]/F");

  BranchName = "shwr_hasPFParticle_" + ShowerLabel;
  CreateBranch(BranchName, shwr_hasPFParticle, BranchName + NShowerIndexStr + "/S");

  BranchName = "shwr_PFParticleID_" + ShowerLabel;
  CreateBranch(BranchName, shwr_PFParticleID, BranchName + NShowerIndexStr + "/S");

  BranchName = "shwr_pidmvamu_" + ShowerLabel;
  CreateBranch(BranchName, shwr_pidmvamu, BranchName + NShowerIndexStr + "/F");

  BranchName = "shwr_pidmvae_" + ShowerLabel;
  CreateBranch(BranchName, shwr_pidmvae, BranchName + NShowerIndexStr + "/F");

  BranchName = "shwr_pidmvapich_" + ShowerLabel;
  CreateBranch(BranchName, shwr_pidmvapich, BranchName + NShowerIndexStr + "/F");

  BranchName = "shwr_pidmvaphoton_" + ShowerLabel;
  CreateBranch(BranchName, shwr_pidmvaphoton, BranchName + NShowerIndexStr + "/F");

  BranchName = "shwr_pidmvapr_" + ShowerLabel;
  CreateBranch(BranchName, shwr_pidmvapr, BranchName + NShowerIndexStr + "/F");

} // dune::AnaRootParserDataStruct::ShowerDataStruct::SetAddresses()

//------------------------------------------------------------------------------
//---  AnaRootParserDataStruct
//---

void dune::AnaRootParserDataStruct::ClearLocalData() {

  //  RunData.Clear();
  SubRunData.Clear();

  run = -9999;
  subrun = -9999;
  event = -9999;
  evttime_seconds = -9999;
  evttime_nanoseconds = -9999;
  beamtime = -9999;
  isdata = -99;
  taulife = -9999;
  triggernumber = 0;
  triggertime = -9999;
  beamgatetime = -9999;
  triggerbits = 0;
  potbnb = 0;
  potnumitgt = 0;
  potnumi101 = 0;

  no_hits = 0;
  no_hits_stored = 0;
  NHitsInAllTracks = 0;

  std::fill(hit_tpc, hit_tpc + sizeof(hit_tpc)/sizeof(hit_tpc[0]), -999);
  std::fill(hit_view, hit_view + sizeof(hit_view)/sizeof(hit_view[0]), -999);
  std::fill(hit_wire, hit_wire + sizeof(hit_wire)/sizeof(hit_wire[0]), -999);
  std::fill(hit_channel, hit_channel + sizeof(hit_channel)/sizeof(hit_channel[0]), -999);
  std::fill(hit_peakT, hit_peakT + sizeof(hit_peakT)/sizeof(hit_peakT[0]), -999.);
  std::fill(hit_chargesum, hit_chargesum + sizeof(hit_chargesum)/sizeof(hit_chargesum[0]), -999.);
  std::fill(hit_chargeintegral, hit_chargeintegral + sizeof(hit_chargeintegral)/sizeof(hit_chargeintegral[0]), -999.);
  std::fill(hit_ph, hit_ph + sizeof(hit_ph)/sizeof(hit_ph[0]), -999.);
  std::fill(hit_startT, hit_startT + sizeof(hit_startT)/sizeof(hit_startT[0]), -999.);
  std::fill(hit_endT, hit_endT + sizeof(hit_endT)/sizeof(hit_endT[0]), -999.);
  std::fill(hit_rms, hit_rms + sizeof(hit_rms)/sizeof(hit_rms[0]), -999.);
  //  std::fill(hit_trueX, hit_trueX + sizeof(hit_trueX)/sizeof(hit_trueX[0]), -999.);
  std::fill(hit_goodnessOfFit, hit_goodnessOfFit + sizeof(hit_goodnessOfFit)/sizeof(hit_goodnessOfFit[0]), -999.);
  std::fill(hit_fitparamampl, hit_fitparamampl + sizeof(hit_fitparamampl)/sizeof(hit_fitparamampl[0]), -999.);
  std::fill(hit_fitparamt0, hit_fitparamt0 + sizeof(hit_fitparamt0)/sizeof(hit_fitparamt0[0]), -999.);
  std::fill(hit_fitparamtau1, hit_fitparamtau1 + sizeof(hit_fitparamtau1)/sizeof(hit_fitparamtau1[0]), -999.);
  std::fill(hit_fitparamtau2, hit_fitparamtau2 + sizeof(hit_fitparamtau2)/sizeof(hit_fitparamtau2[0]), -999.);
  std::fill(hit_multiplicity, hit_multiplicity + sizeof(hit_multiplicity)/sizeof(hit_multiplicity[0]), -999.);
  std::fill(hit_trueID, hit_trueID + sizeof(hit_trueID)/sizeof(hit_trueID[0]), -999.);
  std::fill(hit_trueEnergyMax, hit_trueEnergyMax + sizeof(hit_trueEnergyMax)/sizeof(hit_trueEnergyMax[0]), -999.);
  std::fill(hit_trueEnergyFraction, hit_trueEnergyFraction + sizeof(hit_trueEnergyFraction)/sizeof(hit_trueEnergyFraction[0]), -999.);
  std::fill(hit_trkid, hit_trkid + sizeof(hit_trkid)/sizeof(hit_trkid[0]), -999);
  //  std::fill(hit_trkKey, hit_trkKey + sizeof(hit_trkKey)/sizeof(hit_trkKey[0]), -999);
  std::fill(hit_clusterid, hit_clusterid + sizeof(hit_clusterid)/sizeof(hit_clusterid[0]), -9999);
  //  std::fill(hit_clusterKey, hit_clusterKey + sizeof(hit_clusterKey)/sizeof(hit_clusterKey[0]), -999);
  //  std::fill(hit_nelec, hit_nelec + sizeof(hit_nelec)/sizeof(hit_nelec[0]), -999.);
  //  std::fill(hit_energy, hit_energy + sizeof(hit_energy)/sizeof(hit_energy[0]), -999.);
  //raw digit information
/*
  std::fill(rawD_ph, rawD_ph + sizeof(rawD_ph)/sizeof(rawD_ph[0]), -999.);
  std::fill(rawD_peakT, rawD_peakT + sizeof(rawD_peakT)/sizeof(rawD_peakT[0]), -999.);
  std::fill(rawD_charge, rawD_charge + sizeof(rawD_charge)/sizeof(rawD_charge[0]), -999.);
  std::fill(rawD_fwhh, rawD_fwhh + sizeof(rawD_fwhh)/sizeof(rawD_fwhh[0]), -999.);
  std::fill(rawD_rms, rawD_rms + sizeof(rawD_rms)/sizeof(rawD_rms[0]), -999.);
*/
  no_channels = 0;
  no_ticks = 0;
  no_ticksinallchannels = 0;
  std::fill(rawD_ADC, rawD_ADC + sizeof(rawD_ADC)/sizeof(rawD_ADC[0]), -999);
  std::fill(rawD_Channel, rawD_Channel + sizeof(rawD_Channel)/sizeof(rawD_Channel[0]), -999);

  no_recochannels=0;
  std::fill(recoW_Channel, recoW_Channel + sizeof(recoW_Channel)/sizeof(recoW_Channel[0]), -999);
  std::fill(recoW_NTicks, recoW_NTicks + sizeof(recoW_NTicks)/sizeof(recoW_NTicks[0]), -999);

  no_recoticksinallchannels=0;
  std::fill(recoW_Tick, recoW_Tick + sizeof(recoW_Tick)/sizeof(recoW_Tick[0]), -999);
  std::fill(recoW_ADC, recoW_ADC + sizeof(recoW_ADC)/sizeof(recoW_ADC[0]), -999);

  numberofphotons=0;
  FillWith(photons_time,-999);
  FillWith(photons_channel,-999);

  no_flashes = 0;
  std::fill(flash_time, flash_time + sizeof(flash_time)/sizeof(flash_time[0]), -999);
  std::fill(flash_pe, flash_pe + sizeof(flash_pe)/sizeof(flash_pe[0]), -999);
  std::fill(flash_ycenter, flash_ycenter + sizeof(flash_ycenter)/sizeof(flash_ycenter[0]), -999);
  std::fill(flash_zcenter, flash_zcenter + sizeof(flash_zcenter)/sizeof(flash_zcenter[0]), -999);
  std::fill(flash_ywidth, flash_ywidth + sizeof(flash_ywidth)/sizeof(flash_ywidth[0]), -999);
  std::fill(flash_zwidth, flash_zwidth + sizeof(flash_zwidth)/sizeof(flash_zwidth[0]), -999);
  std::fill(flash_timewidth, flash_timewidth + sizeof(flash_timewidth)/sizeof(flash_timewidth[0]), -999);

  no_ExternCounts = 0;
  std::fill(externcounts_time, externcounts_time + sizeof(externcounts_time)/sizeof(externcounts_time[0]), -999);
  std::fill(externcounts_id, externcounts_id + sizeof(externcounts_id)/sizeof(externcounts_id[0]), -999);

  nclusters = 0;
  std::fill(clusterId, clusterId + sizeof(clusterId)/sizeof(clusterId[0]), -999);
  std::fill(clusterView, clusterView + sizeof(clusterView)/sizeof(clusterView[0]), -999);
  std::fill(cluster_isValid, cluster_isValid + sizeof(cluster_isValid)/sizeof(cluster_isValid[0]), -1);
  std::fill(cluster_StartCharge, cluster_StartCharge +  sizeof(cluster_StartCharge)/sizeof(cluster_StartCharge[0]), -999.);
  std::fill(cluster_StartAngle, cluster_StartAngle + sizeof(cluster_StartAngle)/sizeof(cluster_StartAngle[0]), -999.);
  std::fill(cluster_EndCharge, cluster_EndCharge + sizeof(cluster_EndCharge)/sizeof(cluster_EndCharge[0]), -999.);
  std::fill(cluster_EndAngle , cluster_EndAngle + sizeof(cluster_EndAngle)/sizeof(cluster_EndAngle[0]), -999.);
  std::fill(cluster_Integral , cluster_Integral + sizeof(cluster_Integral)/sizeof(cluster_Integral[0]), -999.);
  std::fill(cluster_IntegralAverage, cluster_IntegralAverage + sizeof(cluster_IntegralAverage)/sizeof(cluster_IntegralAverage[0]), -999.);
  std::fill(cluster_SummedADC, cluster_SummedADC + sizeof(cluster_SummedADC)/sizeof(cluster_SummedADC[0]), -999.);
  std::fill(cluster_SummedADCaverage, cluster_SummedADCaverage + sizeof(cluster_SummedADCaverage)/sizeof(cluster_SummedADCaverage[0]), -999.);
  std::fill(cluster_MultipleHitDensity, cluster_MultipleHitDensity + sizeof(cluster_MultipleHitDensity)/sizeof(cluster_MultipleHitDensity[0]), -999.);
  std::fill(cluster_Width, cluster_Width + sizeof(cluster_Width)/sizeof(cluster_Width[0]), -999.);
  std::fill(cluster_NHits, cluster_NHits + sizeof(cluster_NHits)/sizeof(cluster_NHits[0]), -999);
  std::fill(cluster_StartWire, cluster_StartWire + sizeof(cluster_StartWire)/sizeof(cluster_StartWire[0]), -999);
  std::fill(cluster_StartTick, cluster_StartTick + sizeof(cluster_StartTick)/sizeof(cluster_StartTick[0]), -999);
  std::fill(cluster_EndWire, cluster_EndWire + sizeof(cluster_EndWire)/sizeof(cluster_EndWire[0]), -999);
  std::fill(cluster_EndTick, cluster_EndTick + sizeof(cluster_EndTick)/sizeof(cluster_EndTick[0]), -999);
  //  std::fill(cluncosmictags_tagger, cluncosmictags_tagger + sizeof(cluncosmictags_tagger)/sizeof(cluncosmictags_tagger[0]), -999);
  //  std::fill(clucosmicscore_tagger, clucosmicscore_tagger + sizeof(clucosmicscore_tagger)/sizeof(clucosmicscore_tagger[0]), -999.);
  //  std::fill(clucosmictype_tagger , clucosmictype_tagger  + sizeof(clucosmictype_tagger )/sizeof(clucosmictype_tagger [0]), -999);

  nnuvtx = 0;
  std::fill(nuvtxx, nuvtxx + sizeof(nuvtxx)/sizeof(nuvtxx[0]), -999.);
  std::fill(nuvtxy, nuvtxy + sizeof(nuvtxy)/sizeof(nuvtxy[0]), -999.);
  std::fill(nuvtxz, nuvtxz + sizeof(nuvtxz)/sizeof(nuvtxz[0]), -999.);
  std::fill(nuvtxpdg, nuvtxpdg + sizeof(nuvtxpdg)/sizeof(nuvtxpdg[0]), -9999);

  mcevts_truth = -9999;
  mcevts_truthcry = -9999;
  std::fill(nuPDG_truth, nuPDG_truth + sizeof(nuPDG_truth)/sizeof(nuPDG_truth[0]), -999.);
  std::fill(ccnc_truth, ccnc_truth + sizeof(ccnc_truth)/sizeof(ccnc_truth[0]), -999.);
  std::fill(mode_truth, mode_truth + sizeof(mode_truth)/sizeof(mode_truth[0]), -999.);
  std::fill(enu_truth, enu_truth + sizeof(enu_truth)/sizeof(enu_truth[0]), -999.);
  std::fill(Q2_truth, Q2_truth + sizeof(Q2_truth)/sizeof(Q2_truth[0]), -999.);
  std::fill(W_truth, W_truth + sizeof(W_truth)/sizeof(W_truth[0]), -999.);
  std::fill(X_truth, X_truth + sizeof(X_truth)/sizeof(X_truth[0]), -999.);
  std::fill(Y_truth, Y_truth + sizeof(Y_truth)/sizeof(Y_truth[0]), -999.);
  std::fill(hitnuc_truth, hitnuc_truth + sizeof(hitnuc_truth)/sizeof(hitnuc_truth[0]), -999.);
  std::fill(nuvtxx_truth, nuvtxx_truth + sizeof(nuvtxx_truth)/sizeof(nuvtxx_truth[0]), -999.);
  std::fill(nuvtxy_truth, nuvtxy_truth + sizeof(nuvtxy_truth)/sizeof(nuvtxy_truth[0]), -999.);
  std::fill(nuvtxz_truth, nuvtxz_truth + sizeof(nuvtxz_truth)/sizeof(nuvtxz_truth[0]), -999.);
  std::fill(nu_dcosx_truth, nu_dcosx_truth + sizeof(nu_dcosx_truth)/sizeof(nu_dcosx_truth[0]), -999.);
  std::fill(nu_dcosy_truth, nu_dcosy_truth + sizeof(nu_dcosy_truth)/sizeof(nu_dcosy_truth[0]), -999.);
  std::fill(nu_dcosz_truth, nu_dcosz_truth + sizeof(nu_dcosz_truth)/sizeof(nu_dcosz_truth[0]), -999.);
  std::fill(lep_mom_truth, lep_mom_truth + sizeof(lep_mom_truth)/sizeof(lep_mom_truth[0]), -999.);
  std::fill(lep_dcosx_truth, lep_dcosx_truth + sizeof(lep_dcosx_truth)/sizeof(lep_dcosx_truth[0]), -999.);
  std::fill(lep_dcosy_truth, lep_dcosy_truth + sizeof(lep_dcosy_truth)/sizeof(lep_dcosy_truth[0]), -999.);
  std::fill(lep_dcosz_truth, lep_dcosz_truth + sizeof(lep_dcosz_truth)/sizeof(lep_dcosz_truth[0]), -999.);

  //Flux information
  std::fill(vx_flux, vx_flux + sizeof(vx_flux)/sizeof(vx_flux[0]), -999.);
  std::fill(vy_flux, vy_flux + sizeof(vy_flux)/sizeof(vy_flux[0]), -999.);
  std::fill(vz_flux, vz_flux + sizeof(vz_flux)/sizeof(vz_flux[0]), -999.);
  std::fill(pdpx_flux, pdpx_flux + sizeof(pdpx_flux)/sizeof(pdpx_flux[0]), -999.);
  std::fill(pdpy_flux, pdpy_flux + sizeof(pdpy_flux)/sizeof(pdpy_flux[0]), -999.);
  std::fill(pdpz_flux, pdpz_flux + sizeof(pdpz_flux)/sizeof(pdpz_flux[0]), -999.);
  std::fill(ppdxdz_flux, ppdxdz_flux + sizeof(ppdxdz_flux)/sizeof(ppdxdz_flux[0]), -999.);
  std::fill(ppdydz_flux, ppdydz_flux + sizeof(ppdydz_flux)/sizeof(ppdydz_flux[0]), -999.);
  std::fill(pppz_flux, pppz_flux + sizeof(pppz_flux)/sizeof(pppz_flux[0]), -999.);
  std::fill(ptype_flux, ptype_flux + sizeof(ptype_flux)/sizeof(ptype_flux[0]), -999);
  std::fill(ppvx_flux, ppvx_flux + sizeof(ppvx_flux)/sizeof(ppvx_flux[0]), -999.);
  std::fill(ppvy_flux, ppvy_flux + sizeof(ppvy_flux)/sizeof(ppvy_flux[0]), -999.);
  std::fill(ppvz_flux, ppvz_flux + sizeof(ppvz_flux)/sizeof(ppvz_flux[0]), -999.);
  std::fill(muparpx_flux, muparpx_flux + sizeof(muparpx_flux)/sizeof(muparpx_flux[0]), -999.);
  std::fill(muparpy_flux, muparpy_flux + sizeof(muparpy_flux)/sizeof(muparpy_flux[0]), -999.);
  std::fill(muparpz_flux, muparpz_flux + sizeof(muparpz_flux)/sizeof(muparpz_flux[0]), -999.);
  std::fill(mupare_flux, mupare_flux + sizeof(mupare_flux)/sizeof(mupare_flux[0]), -999.);
  std::fill(tgen_flux, tgen_flux + sizeof(tgen_flux)/sizeof(tgen_flux[0]), -999);
  std::fill(tgptype_flux, tgptype_flux + sizeof(tgptype_flux)/sizeof(tgptype_flux[0]), -999);
  std::fill(tgppx_flux, tgppx_flux + sizeof(tgppx_flux)/sizeof(tgppx_flux[0]), -999.);
  std::fill(tgppy_flux, tgppy_flux + sizeof(tgppy_flux)/sizeof(tgppy_flux[0]), -999.);
  std::fill(tgppz_flux, tgppz_flux + sizeof(tgppz_flux)/sizeof(tgppz_flux[0]), -999.);
  std::fill(tprivx_flux, tprivx_flux + sizeof(tprivx_flux)/sizeof(tprivx_flux[0]), -999.);
  std::fill(tprivy_flux, tprivy_flux + sizeof(tprivy_flux)/sizeof(tprivy_flux[0]), -999.);
  std::fill(tprivz_flux, tprivz_flux + sizeof(tprivz_flux)/sizeof(tprivz_flux[0]), -999.);
  std::fill(dk2gen_flux, dk2gen_flux + sizeof(dk2gen_flux)/sizeof(dk2gen_flux[0]), -999.);
  std::fill(gen2vtx_flux, gen2vtx_flux + sizeof(gen2vtx_flux)/sizeof(gen2vtx_flux[0]), -999.);
  std::fill(tpx_flux, tpx_flux + sizeof(tpx_flux)/sizeof(tpx_flux[0]), -999.);
  std::fill(tpy_flux, tpy_flux + sizeof(tpy_flux)/sizeof(tpy_flux[0]), -999.);
  std::fill(tpz_flux, tpz_flux + sizeof(tpz_flux)/sizeof(tpz_flux[0]), -999.);
  std::fill(tptype_flux, tptype_flux + sizeof(tptype_flux)/sizeof(tptype_flux[0]), -999.);

  genie_no_primaries = 0;
  cry_no_primaries = 0;
  proto_no_primaries = 0;
  generator_list_size = 0;
  no_primaries = 0;
  geant_list_size=0;
  geant_list_size_in_tpcAV = 0;
  no_mcshowers = 0;
  no_mctracks = 0;

  FillWith(pdg, -9999);
  FillWith(status, -9999);
  FillWith(Mass, -999.);
  FillWith(Eng, -999.);
  FillWith(EndE, -999.);
  FillWith(Px, -999.);
  FillWith(Py, -999.);
  FillWith(Pz, -999.);
  FillWith(P, -999.);
  FillWith(StartPointx, -999.);
  FillWith(StartPointy, -999.);
  FillWith(StartPointz, -999.);
  FillWith(StartT, -99e7);
  FillWith(EndT, -999.);
  FillWith(EndPointx, -999.);
  FillWith(EndPointy, -999.);
  FillWith(EndPointz, -999.);
  FillWith(EndT, -99e7);
  FillWith(theta, -999.);
  FillWith(phi, -999.);
  FillWith(theta_xz, -999.);
  FillWith(theta_yz, -999.);

  FillWith(pathlen_tpcAV, -999.);
  FillWith(inTPCActive, -9999);
  FillWith(TrackId_tpcAV, -9999);
  FillWith(PDGCode_tpcAV, -9999);
  FillWith(StartPointx_tpcAV, -999.);
  FillWith(StartPointy_tpcAV, -999.);
  FillWith(StartPointz_tpcAV, -999.);
  FillWith(StartT_tpcAV, -99e7);
  FillWith(StartE_tpcAV, -999.);
  FillWith(StartP_tpcAV, -999.);
  FillWith(StartPx_tpcAV, -999.);
  FillWith(StartPy_tpcAV, -999.);
  FillWith(StartPz_tpcAV, -999.);
  FillWith(thetastart_tpcAV, -999.);
  FillWith(phistart_tpcAV, -999.);
  FillWith(EndPointx_tpcAV, -999.);
  FillWith(EndPointy_tpcAV, -999.);
  FillWith(EndPointz_tpcAV, -999.);
  FillWith(EndT_tpcAV, -99e7);
  FillWith(EndE_tpcAV, -999.);
  FillWith(EndP_tpcAV, -999.);
  FillWith(EndPx_tpcAV, -999.);
  FillWith(EndPy_tpcAV, -999.);
  FillWith(EndPz_tpcAV, -999.);
  FillWith(thetaend_tpcAV, -999.);
  FillWith(phiend_tpcAV, -999.);

  FillWith(pathlen_drifted, -999.);
  FillWith(inTPCDrifted, -9999);
  FillWith(StartPointx_drifted, -999.);
  FillWith(StartPointy_drifted, -999.);
  FillWith(StartPointz_drifted, -999.);
  FillWith(StartT_drifted, -99e7);
  FillWith(StartE_drifted, -999.);
  FillWith(StartP_drifted, -999.);
  FillWith(StartPx_drifted, -999.);
  FillWith(StartPy_drifted, -999.);
  FillWith(StartPz_drifted, -999.);
  FillWith(EndPointx_drifted, -999.);
  FillWith(EndPointy_drifted, -999.);
  FillWith(EndPointz_drifted, -999.);
  FillWith(EndT_drifted, -99e7);
  FillWith(EndE_drifted, -999.);
  FillWith(EndP_drifted, -999.);
  FillWith(EndPx_drifted, -999.);
  FillWith(EndPy_drifted, -999.);
  FillWith(EndPz_drifted, -999.);
  FillWith(NumberDaughters, -9999);
  FillWith(Mother, -9999);
  FillWith(TrackId, -9999);
  FillWith(process_primary, -9999);
  FillWith(processname, "noname");
  FillWith(MergedId, -9999);
  FillWith(origin, -9999);
  FillWith(MCTruthIndex, -9999);
  FillWith(genie_primaries_pdg, -9999);
  FillWith(genie_Eng, -999.);
  FillWith(genie_Px, -999.);
  FillWith(genie_Py, -999.);
  FillWith(genie_Pz, -999.);
  FillWith(genie_P, -999.);
  FillWith(genie_status_code, -9999);
  FillWith(genie_mass, -999.);
  FillWith(genie_trackID, -9999);
  FillWith(genie_ND, -9999);
  FillWith(genie_mother, -9999);
  FillWith(cry_primaries_pdg, -9999);
  FillWith(cry_Eng, -999.);
  FillWith(cry_Px, -999.);
  FillWith(cry_Py, -999.);
  FillWith(cry_Pz, -999.);
  FillWith(cry_P, -999.);
  FillWith(cry_StartPointx, -999.);
  FillWith(cry_StartPointy, -999.);
  FillWith(cry_StartPointz, -999.);
  FillWith(cry_StartPointt, -999.);
  FillWith(cry_status_code, -9999);
  FillWith(cry_mass, -999.);
  FillWith(cry_trackID, -9999);
  FillWith(cry_ND, -9999);
  FillWith(cry_mother, -9999);
  // Start of ProtoDUNE Beam generator section
  FillWith(proto_isGoodParticle,-9999);
  FillWith(proto_vx,-999.);
  FillWith(proto_vy,-999.);
  FillWith(proto_vz,-999.);
  FillWith(proto_t,-999.);
  FillWith(proto_px,-999.);
  FillWith(proto_py,-999.);
  FillWith(proto_pz,-999.);
  FillWith(proto_momentum,-999.);
  FillWith(proto_energy,-999.);
  FillWith(proto_pdg,-9999);
  // End of ProtoDUNE Beam generator section
  FillWith(mcshwr_origin, -1);
  FillWith(mcshwr_pdg, -9999);
  FillWith(mcshwr_TrackId, -9999);
  FillWith(mcshwr_Process, "noname");
  FillWith(mcshwr_startX, -999.);
  FillWith(mcshwr_startY, -999.);
  FillWith(mcshwr_startZ, -999.);
  FillWith(mcshwr_endX, -999.);
  FillWith(mcshwr_endY, -999.);
  FillWith(mcshwr_endZ, -999.);
  FillWith(mcshwr_CombEngX, -999.);
  FillWith(mcshwr_CombEngY, -999.);
  FillWith(mcshwr_CombEngZ, -999.);
  FillWith(mcshwr_CombEngPx, -999.);
  FillWith(mcshwr_CombEngPy, -999.);
  FillWith(mcshwr_CombEngPz, -999.);
  FillWith(mcshwr_CombEngE, -999.);
  FillWith(mcshwr_dEdx, -999.);
  FillWith(mcshwr_StartDirX, -999.);
  FillWith(mcshwr_StartDirY, -999.);
  FillWith(mcshwr_StartDirZ, -999.);
  FillWith(mcshwr_isEngDeposited, -999);
  FillWith(mcshwr_Motherpdg, -9999);
  FillWith(mcshwr_MotherTrkId, -9999);
  FillWith(mcshwr_MotherProcess, "noname");
  FillWith(mcshwr_MotherstartX, -999.);
  FillWith(mcshwr_MotherstartY, -999.);
  FillWith(mcshwr_MotherstartZ, -999.);
  FillWith(mcshwr_MotherendX, -999.);
  FillWith(mcshwr_MotherendY, -999.);
  FillWith(mcshwr_MotherendZ, -999.);
  FillWith(mcshwr_Ancestorpdg, -9999);
  FillWith(mcshwr_AncestorTrkId, -9999);
  FillWith(mcshwr_AncestorProcess, "noname");
  FillWith(mcshwr_AncestorstartX, -999.);
  FillWith(mcshwr_AncestorstartY, -999.);
  FillWith(mcshwr_AncestorstartZ, -999.);
  FillWith(mcshwr_AncestorendX, -999.);
  FillWith(mcshwr_AncestorendY, -999.);
  FillWith(mcshwr_AncestorendZ, -999.);

  geant_trajectory_size = 0;

  FillWith(NTrajectoryPointsPerParticle, -999);

  FillWith(TrajTrackId, -999.);
  FillWith(TrajPDGCode, -999.);

  FillWith(TrajX, -999.);
  FillWith(TrajY, -999.);
  FillWith(TrajZ, -999.);
  FillWith(TrajT, -999.);
  FillWith(TrajE, -999.);
  FillWith(TrajP, -999.);
  FillWith(TrajPx, -999.);
  FillWith(TrajPy, -999.);
  FillWith(TrajPz, -999.);
  FillWith(TrajTheta, -999.);
  FillWith(TrajPhi, -999.);


  // auxiliary detector information;
  FillWith(NAuxDets, 0);
  // - set to -999 all the values of each of the arrays in AuxDetID;
  //   this auto is BoxedArray<Short_t>
  for (auto& partInfo: AuxDetID) FillWith(partInfo, -999);
  // - pythonish C++: as the previous line, for each one in a list of containers
  //   of the same type (C++ is not python yet), using pointers to avoid copy;
  for (AuxDetMCData_t<Float_t>* cont: {
      &entryX, &entryY, &entryZ, &entryT,
      &exitX , &exitY , &exitZ, &exitT, &exitPx, &exitPy, &exitPz,
      &CombinedEnergyDep
      })
  {
    // this auto is BoxedArray<Float_t>
    for (auto& partInfo: *cont) FillWith(partInfo, -999.);
  } // for container

} // dune::AnaRootParserDataStruct::ClearLocalData()


void dune::AnaRootParserDataStruct::Clear() {
  ClearLocalData();
  std::for_each
    (TrackData.begin(), TrackData.end(), std::mem_fn(&TrackDataStruct::Clear));
  std::for_each
    (VertexData.begin(), VertexData.end(), std::mem_fn(&VertexDataStruct::Clear));
  std::for_each
    (ShowerData.begin(), ShowerData.end(), std::mem_fn(&ShowerDataStruct::Clear));
} // dune::AnaRootParserDataStruct::Clear()


  void dune::AnaRootParserDataStruct::SetShowerAlgos
(std::vector<std::string> const& ShowerAlgos)
{

  size_t const nShowerAlgos = ShowerAlgos.size();
  ShowerData.resize(nShowerAlgos);
  for (size_t iAlgo = 0; iAlgo < nShowerAlgos; ++iAlgo)
    ShowerData[iAlgo].SetName(ShowerAlgos[iAlgo]);

} // dune::AnaRootParserDataStruct::SetShowerAlgos()


void dune::AnaRootParserDataStruct::ResizePhotons(int nPhotons) {

  // minimum size is 1, so that we always have an address
  MaxPhotons = (size_t) std::max(nPhotons, 1);

  photons_time.resize(MaxPhotons);
  photons_channel.resize(MaxPhotons);
} // dune::AnaRootParserDataStruct::ResizePhotons



void dune::AnaRootParserDataStruct::ResizeGenerator(int nParticles) {

  // minimum size is 1, so that we always have an address
  MaxGeneratorparticles = (size_t) std::max(nParticles, 1);

  TrackId.resize(MaxGeneratorparticles);
  pdg.resize(MaxGeneratorparticles);
  status.resize(MaxGeneratorparticles);
  Mass.resize(MaxGeneratorparticles);
  Eng.resize(MaxGeneratorparticles);
  Px.resize(MaxGeneratorparticles);
  Py.resize(MaxGeneratorparticles);
  Pz.resize(MaxGeneratorparticles);
  P.resize(MaxGeneratorparticles);
  StartPointx.resize(MaxGeneratorparticles);
  StartPointy.resize(MaxGeneratorparticles);
  StartPointz.resize(MaxGeneratorparticles);
  StartT.resize(MaxGeneratorparticles);
  theta.resize(MaxGeneratorparticles);
  phi.resize(MaxGeneratorparticles);
} //dune::AnaRootParserDataStruct::ResizeGenerator

void dune::AnaRootParserDataStruct::ResizeGEANT(int nParticles) {

  // minimum size is 1, so that we always have an address
  MaxGEANTparticles = (size_t) std::max(nParticles, 1);

  pdg.resize(MaxGEANTparticles);
  status.resize(MaxGEANTparticles);
  inTPCActive.resize(MaxGEANTparticles);
  Mass.resize(MaxGEANTparticles);
  Eng.resize(MaxGEANTparticles);
  EndE.resize(MaxGEANTparticles);
  Px.resize(MaxGEANTparticles);
  Py.resize(MaxGEANTparticles);
  Pz.resize(MaxGEANTparticles);
  P.resize(MaxGEANTparticles);
  StartPointx.resize(MaxGEANTparticles);
  StartPointy.resize(MaxGEANTparticles);
  StartPointz.resize(MaxGEANTparticles);
  StartT.resize(MaxGEANTparticles);
  EndT.resize(MaxGEANTparticles);
  EndPointx.resize(MaxGEANTparticles);
  EndPointy.resize(MaxGEANTparticles);
  EndPointz.resize(MaxGEANTparticles);
  EndT.resize(MaxGEANTparticles);
  theta.resize(MaxGEANTparticles);
  phi.resize(MaxGEANTparticles);
  theta_xz.resize(MaxGEANTparticles);
  theta_yz.resize(MaxGEANTparticles);
  NumberDaughters.resize(MaxGEANTparticles);
  Mother.resize(MaxGEANTparticles);
  TrackId.resize(MaxGEANTparticles);
  process_primary.resize(MaxGEANTparticles);
  processname.resize(MaxGEANTparticles);
  MergedId.resize(MaxGEANTparticles);
  origin.resize(MaxGEANTparticles);
  MCTruthIndex.resize(MaxGEANTparticles);

  pathlen_drifted.resize(MaxGEANTparticles);
  inTPCDrifted.resize(MaxGEANTparticles);
  StartPointx_drifted.resize(MaxGEANTparticles);
  StartPointy_drifted.resize(MaxGEANTparticles);
  StartPointz_drifted.resize(MaxGEANTparticles);
  StartT_drifted.resize(MaxGEANTparticles);
  StartE_drifted.resize(MaxGEANTparticles);
  StartP_drifted.resize(MaxGEANTparticles);
  StartPx_drifted.resize(MaxGEANTparticles);
  StartPy_drifted.resize(MaxGEANTparticles);
  StartPz_drifted.resize(MaxGEANTparticles);
  EndPointx_drifted.resize(MaxGEANTparticles);
  EndPointy_drifted.resize(MaxGEANTparticles);
  EndPointz_drifted.resize(MaxGEANTparticles);
  EndT_drifted.resize(MaxGEANTparticles);
  EndE_drifted.resize(MaxGEANTparticles);
  EndP_drifted.resize(MaxGEANTparticles);
  EndPx_drifted.resize(MaxGEANTparticles);
  EndPy_drifted.resize(MaxGEANTparticles);
  EndPz_drifted.resize(MaxGEANTparticles);

  // auxiliary detector structure
  NAuxDets.resize(MaxGEANTparticles);
  AuxDetID.resize(MaxGEANTparticles);
  entryX.resize(MaxGEANTparticles);
  entryY.resize(MaxGEANTparticles);
  entryZ.resize(MaxGEANTparticles);
  entryT.resize(MaxGEANTparticles);
  exitX.resize(MaxGEANTparticles);
  exitY.resize(MaxGEANTparticles);
  exitZ.resize(MaxGEANTparticles);
  exitT.resize(MaxGEANTparticles);
  exitPx.resize(MaxGEANTparticles);
  exitPy.resize(MaxGEANTparticles);
  exitPz.resize(MaxGEANTparticles);
  CombinedEnergyDep.resize(MaxGEANTparticles);

  NTrajectoryPointsPerParticle.resize(MaxGEANTparticles);
} // dune::AnaRootParserDataStruct::ResizeGEANT()

void dune::AnaRootParserDataStruct::ResizeGEANTInAV(int nParticlesInAV){

  MaxGEANTInAVparticles = (size_t) std::max(nParticlesInAV, 1);

  pathlen_tpcAV.resize(MaxGEANTInAVparticles);
  TrackId_tpcAV.resize(MaxGEANTInAVparticles);
  PDGCode_tpcAV.resize(MaxGEANTInAVparticles);
  StartPointx_tpcAV.resize(MaxGEANTInAVparticles);
  StartPointy_tpcAV.resize(MaxGEANTInAVparticles);
  StartPointz_tpcAV.resize(MaxGEANTInAVparticles);
  StartT_tpcAV.resize(MaxGEANTInAVparticles);
  StartE_tpcAV.resize(MaxGEANTInAVparticles);
  StartP_tpcAV.resize(MaxGEANTInAVparticles);
  StartPx_tpcAV.resize(MaxGEANTInAVparticles);
  StartPy_tpcAV.resize(MaxGEANTInAVparticles);
  StartPz_tpcAV.resize(MaxGEANTInAVparticles);
  thetastart_tpcAV.resize(MaxGEANTInAVparticles);
  phistart_tpcAV.resize(MaxGEANTInAVparticles);
  EndPointx_tpcAV.resize(MaxGEANTInAVparticles);
  EndPointy_tpcAV.resize(MaxGEANTInAVparticles);
  EndPointz_tpcAV.resize(MaxGEANTInAVparticles);
  EndT_tpcAV.resize(MaxGEANTInAVparticles);
  EndE_tpcAV.resize(MaxGEANTInAVparticles);
  EndP_tpcAV.resize(MaxGEANTInAVparticles);
  EndPx_tpcAV.resize(MaxGEANTInAVparticles);
  EndPy_tpcAV.resize(MaxGEANTInAVparticles);
  EndPz_tpcAV.resize(MaxGEANTInAVparticles);
  thetaend_tpcAV.resize(MaxGEANTInAVparticles);
  phiend_tpcAV.resize(MaxGEANTInAVparticles);
} // dune::AnaRootParserDataStruct::ResizeGEANTInAV()

void dune::AnaRootParserDataStruct::ResizeGEANTTrajectory(int nTrajectoryPoints) {

  // minimum size is 1, so that we always have an address
  MaxGEANTtrajectorypoints = (size_t) std::max(nTrajectoryPoints, 1);

  TrajTrackId.resize(MaxGEANTtrajectorypoints);
  TrajPDGCode.resize(MaxGEANTtrajectorypoints);

  TrajX.resize(MaxGEANTtrajectorypoints);
  TrajY.resize(MaxGEANTtrajectorypoints);
  TrajZ.resize(MaxGEANTtrajectorypoints);
  TrajT.resize(MaxGEANTtrajectorypoints);
  TrajE.resize(MaxGEANTtrajectorypoints);
  TrajP.resize(MaxGEANTtrajectorypoints);
  TrajPx.resize(MaxGEANTtrajectorypoints);
  TrajPy.resize(MaxGEANTtrajectorypoints);
  TrajPz.resize(MaxGEANTtrajectorypoints);
  TrajTheta.resize(MaxGEANTtrajectorypoints);
  TrajPhi.resize(MaxGEANTtrajectorypoints);

} // dune::AnaRootParserDataStruct::ResizeGEANTTrajectory()


void dune::AnaRootParserDataStruct::ResizeGenie(int nPrimaries) {

  // minimum size is 1, so that we always have an address
  MaxGeniePrimaries = (size_t) std::max(nPrimaries, 1);
  genie_primaries_pdg.resize(MaxGeniePrimaries);
  genie_Eng.resize(MaxGeniePrimaries);
  genie_Px.resize(MaxGeniePrimaries);
  genie_Py.resize(MaxGeniePrimaries);
  genie_Pz.resize(MaxGeniePrimaries);
  genie_P.resize(MaxGeniePrimaries);
  genie_status_code.resize(MaxGeniePrimaries);
  genie_mass.resize(MaxGeniePrimaries);
  genie_trackID.resize(MaxGeniePrimaries);
  genie_ND.resize(MaxGeniePrimaries);
  genie_mother.resize(MaxGeniePrimaries);
} // dune::AnaRootParserDataStruct::ResizeGenie()

void dune::AnaRootParserDataStruct::ResizeCry(int nPrimaries) {

  cry_primaries_pdg.resize(nPrimaries);
  cry_Eng.resize(nPrimaries);
  cry_Px.resize(nPrimaries);
  cry_Py.resize(nPrimaries);
  cry_Pz.resize(nPrimaries);
  cry_P.resize(nPrimaries);
  cry_StartPointx.resize(nPrimaries);
  cry_StartPointy.resize(nPrimaries);
  cry_StartPointz.resize(nPrimaries);
  cry_StartPointt.resize(nPrimaries);
  cry_status_code.resize(nPrimaries);
  cry_mass.resize(nPrimaries);
  cry_trackID.resize(nPrimaries);
  cry_ND.resize(nPrimaries);
  cry_mother.resize(nPrimaries);

} // dune::AnaRootParserDataStruct::ResizeCry()

void dune::AnaRootParserDataStruct::ResizeProto(int nPrimaries){

  proto_isGoodParticle.resize(nPrimaries);
  proto_vx.resize(nPrimaries);
  proto_vy.resize(nPrimaries);
  proto_vz.resize(nPrimaries);
  proto_t.resize(nPrimaries);
  proto_px.resize(nPrimaries);
  proto_py.resize(nPrimaries);
  proto_pz.resize(nPrimaries);
  proto_momentum.resize(nPrimaries);
  proto_energy.resize(nPrimaries);
  proto_pdg.resize(nPrimaries);

} // dune::AnaRootParserDataStruct::ResizeProto()

void dune::AnaRootParserDataStruct::ResizeMCShower(int nMCShowers) {
  mcshwr_origin.resize(nMCShowers);
  mcshwr_pdg.resize(nMCShowers);
  mcshwr_TrackId.resize(nMCShowers);
  mcshwr_Process.resize(nMCShowers);
  mcshwr_startX.resize(nMCShowers);
  mcshwr_startY.resize(nMCShowers);
  mcshwr_startZ.resize(nMCShowers);
  mcshwr_endX.resize(nMCShowers);
  mcshwr_endY.resize(nMCShowers);
  mcshwr_endZ.resize(nMCShowers);
  mcshwr_CombEngX.resize(nMCShowers);
  mcshwr_CombEngY.resize(nMCShowers);
  mcshwr_CombEngZ.resize(nMCShowers);
  mcshwr_CombEngPx.resize(nMCShowers);
  mcshwr_CombEngPy.resize(nMCShowers);
  mcshwr_CombEngPz.resize(nMCShowers);
  mcshwr_CombEngE.resize(nMCShowers);
  mcshwr_dEdx.resize(nMCShowers);
  mcshwr_StartDirX.resize(nMCShowers);
  mcshwr_StartDirY.resize(nMCShowers);
  mcshwr_StartDirZ.resize(nMCShowers);
  mcshwr_isEngDeposited.resize(nMCShowers);
  mcshwr_Motherpdg.resize(nMCShowers);
  mcshwr_MotherTrkId.resize(nMCShowers);
  mcshwr_MotherProcess.resize(nMCShowers);
  mcshwr_MotherstartX.resize(nMCShowers);
  mcshwr_MotherstartY.resize(nMCShowers);
  mcshwr_MotherstartZ.resize(nMCShowers);
  mcshwr_MotherendX.resize(nMCShowers);
  mcshwr_MotherendY.resize(nMCShowers);
  mcshwr_MotherendZ.resize(nMCShowers);
  mcshwr_Ancestorpdg.resize(nMCShowers);
  mcshwr_AncestorTrkId.resize(nMCShowers);
  mcshwr_AncestorProcess.resize(nMCShowers);
  mcshwr_AncestorstartX.resize(nMCShowers);
  mcshwr_AncestorstartY.resize(nMCShowers);
  mcshwr_AncestorstartZ.resize(nMCShowers);
  mcshwr_AncestorendX.resize(nMCShowers);
  mcshwr_AncestorendY.resize(nMCShowers);
  mcshwr_AncestorendZ.resize(nMCShowers);

} // dune::AnaRootParserDataStruct::ResizeMCShower()

void dune::AnaRootParserDataStruct::ResizeMCTrack(int nMCTracks) {
  mctrk_origin.resize(nMCTracks);
  mctrk_pdg.resize(nMCTracks);
  mctrk_TrackId.resize(nMCTracks);
  mctrk_Process.resize(nMCTracks);
  mctrk_startX.resize(nMCTracks);
  mctrk_startY.resize(nMCTracks);
  mctrk_startZ.resize(nMCTracks);
  mctrk_endX.resize(nMCTracks);
  mctrk_endY.resize(nMCTracks);
  mctrk_endZ.resize(nMCTracks);
  mctrk_startX_drifted.resize(nMCTracks);
  mctrk_startY_drifted.resize(nMCTracks);
  mctrk_startZ_drifted.resize(nMCTracks);
  mctrk_endX_drifted.resize(nMCTracks);
  mctrk_endY_drifted.resize(nMCTracks);
  mctrk_endZ_drifted.resize(nMCTracks);
  mctrk_len_drifted.resize(nMCTracks);
  mctrk_p_drifted.resize(nMCTracks);
  mctrk_px_drifted.resize(nMCTracks);
  mctrk_py_drifted.resize(nMCTracks);
  mctrk_pz_drifted.resize(nMCTracks);
  mctrk_Motherpdg.resize(nMCTracks);
  mctrk_MotherTrkId.resize(nMCTracks);
  mctrk_MotherProcess.resize(nMCTracks);
  mctrk_MotherstartX.resize(nMCTracks);
  mctrk_MotherstartY.resize(nMCTracks);
  mctrk_MotherstartZ.resize(nMCTracks);
  mctrk_MotherendX.resize(nMCTracks);
  mctrk_MotherendY.resize(nMCTracks);
  mctrk_MotherendZ.resize(nMCTracks);
  mctrk_Ancestorpdg.resize(nMCTracks);
  mctrk_AncestorTrkId.resize(nMCTracks);
  mctrk_AncestorProcess.resize(nMCTracks);
  mctrk_AncestorstartX.resize(nMCTracks);
  mctrk_AncestorstartY.resize(nMCTracks);
  mctrk_AncestorstartZ.resize(nMCTracks);
  mctrk_AncestorendX.resize(nMCTracks);
  mctrk_AncestorendY.resize(nMCTracks);
  mctrk_AncestorendZ.resize(nMCTracks);

} // dune::AnaRootParserDataStruct::ResizeMCTrack()



void dune::AnaRootParserDataStruct::SetAddresses(
    TTree* pTree,
    const std::vector<std::string>& trackers,
    const std::vector<std::string>& vertexalgos,
    const std::vector<std::string>& showeralgos,
    bool isCosmics
    ) {
  BranchCreator CreateBranch(pTree);

  CreateBranch("Run",&run,"run/I");
  CreateBranch("Subrun",&subrun,"subrun/I");
//  CreateBranch("EventNumberInSubrun",&event,"event/I");
  CreateBranch("EventNumberInRun",&event,"event/I");
  CreateBranch("EventTimeSeconds",&evttime_seconds,"evttime_seconds/I");
  CreateBranch("EventTimeNanoseconds",&evttime_nanoseconds,"evttime_nanoseconds/I");
  //  CreateBranch("beamtime",&beamtime,"beamtime/D");
  //  CreateBranch("pot",&SubRunData.pot,"pot/D");
  CreateBranch("IsData",&isdata,"isdata/B");
  //  CreateBranch("taulife",&taulife,"taulife/D");
  //  CreateBranch("triggernumber",&triggernumber,"triggernumber/i");
  //  CreateBranch("triggertime",&triggertime,"triggertime/D");
  //  CreateBranch("beamgatetime",&beamgatetime,"beamgatetime/D");
  //  CreateBranch("triggerbits",&triggerbits,"triggerbits/i");
  //  CreateBranch("potbnb",&potbnb,"potbnb/D");
  //  CreateBranch("potnumitgt",&potnumitgt,"potnumitgt/D");
  //  CreateBranch("potnumi101",&potnumi101,"potnumi101/D");

if (hasRawDigitInfo()){
    CreateBranch("RawWaveform_NumberOfChannels",&no_channels,"no_channels/I");
    CreateBranch("RawWaveform_NumberOfTicks",&no_ticks,"no_ticks/I");

    CreateBranch("RawWaveform_Channel",rawD_Channel,"rawD_Channel[no_channels]/I");

    CreateBranch("RawWaveform_NumberOfTicksInAllChannels",&no_ticksinallchannels,"no_ticksinallchannels/I");
    CreateBranch("RawWaveform_ADC",rawD_ADC,"rawD_ADC[no_ticksinallchannels]/S");
}

if (hasRecobWireInfo()){
    CreateBranch("RecoWaveforms_NumberOfChannels",&no_recochannels,"no_recochannels/I");
    CreateBranch("RecoWaveform_Channel",recoW_Channel,"recoW_Channel[no_recochannels]/I");
    CreateBranch("RecoWaveform_NTicks",recoW_NTicks,"recoW_NTicks[no_recochannels]/I");

    CreateBranch("RecoWaveform_NumberOfTicksInAllChannels",&no_recoticksinallchannels,"no_recoticksinallchannels/I");
    CreateBranch("RecoWaveform_Tick",recoW_Tick,"recoW_Tick[no_recoticksinallchannels]/I");
    CreateBranch("RecoWaveform_ADC",recoW_ADC,"recoW_ADC[no_recoticksinallchannels]/F");
}

  if (hasHitInfo()){
    CreateBranch("NumberOfHits",&no_hits,"no_hits/I");
    //CreateBranch("NumberOfHits_Stored,&no_hits_stored,"no_hits_stored/I");
    CreateBranch("Hit_TPC",hit_tpc,"hit_tpc[no_hits]/S");
    CreateBranch("Hit_View",hit_view,"hit_view[no_hits]/S");
    //CreateBranch("hit_wire",hit_wire,"hit_wire[no_hits]/S");
    CreateBranch("Hit_Channel",hit_channel,"hit_channel[no_hits]/S");
    CreateBranch("Hit_PeakTime",hit_peakT,"hit_peakT[no_hits]/F");
    CreateBranch("Hit_ChargeSummedADC",hit_chargesum,"hit_chargesum[no_hits]/F");
    CreateBranch("Hit_ChargeIntegral",hit_chargeintegral,"hit_chargeintegral[no_hits]/F");
    CreateBranch("Hit_Amplitude",hit_ph,"hit_ph[no_hits]/F");
    CreateBranch("Hit_StartTime",hit_startT,"hit_startT[no_hits]/F");
    CreateBranch("Hit_EndTime",hit_endT,"hit_endT[no_hits]/F");
    CreateBranch("Hit_Width",hit_rms,"hit_rms[no_hits]/F");
    //CreateBranch("hit_trueX",hit_trueX,"hit_trueX[no_hits]/F");
    CreateBranch("Hit_GoodnessOfFit",hit_goodnessOfFit,"hit_goodnessOfFit[no_hits]/F");

    CreateBranch("Hit_FitParameter_Amplitude", hit_fitparamampl, "hit_fitparamamp[no_hits]/F");
    CreateBranch("Hit_FitParameter_Offset", hit_fitparamt0, "hit_fitparamt0[no_hits]/F");
    CreateBranch("Hit_FitParameter_Tau1", hit_fitparamtau1, "hit_fitparamtau1[no_hits]/F");
    CreateBranch("Hit_FitParameter_Tau2", hit_fitparamtau2, "hit_fitparamtau2[no_hits]/F");

    CreateBranch("Hit_Multiplicity",hit_multiplicity,"hit_multiplicity[no_hits]/S");
    CreateBranch("Hit_trueID", hit_trueID,  "Hit_trueID[no_hits]/I");
    CreateBranch("Hit_trueEnergyMax", hit_trueEnergyMax,  "hit_trueEnergyMax[no_hits]/F");
    CreateBranch("Hit_trueEnergyFraction", hit_trueEnergyFraction,  "hit_trueEnergyFraction[no_hits]/F");
    CreateBranch("Hit_TrackID",hit_trkid,"hit_trkid[no_hits]/S");
    //CreateBranch("hit_trkKey",hit_trkKey,"hit_trkKey[no_hits]/S");
    CreateBranch("Hit_ClusterID",hit_clusterid,"hit_clusterid[no_hits]/S");
    //CreateBranch("hit_clusterKey",hit_clusterKey,"hit_clusterKey[no_hits]/S");
    /*    if (!isCosmics){
          CreateBranch("hit_nelec",hit_nelec,"hit_nelec[no_hits]/F");
          CreateBranch("hit_energy",hit_energy,"hit_energy[no_hits]/F");
          }
          */  /*  if (hasRawDigitInfo()){
            CreateBranch("rawD_ph",rawD_ph,"rawD_ph[no_hits]/F");
            CreateBranch("rawD_peakT",rawD_peakT,"rawD_peakT[no_hits]/F");
            CreateBranch("rawD_charge",rawD_charge,"rawD_charge[no_hits]/F");
            CreateBranch("rawD_fwhh",rawD_fwhh,"rawD_fwhh[no_hits]/F");
            CreateBranch("rawD_rms",rawD_rms,"rawD_rms[no_hits]/D");
          }*/
  }
  /*
     if (hasPandoraNuVertexInfo()){
     CreateBranch("nnuvtx", &nnuvtx, "nnuvtx/S");
     CreateBranch("nuvtxx", nuvtxx, "nuvtxx[nnuvtx]/F");
     CreateBranch("nuvtxy", nuvtxy, "nuvtxy[nnuvtx]/F");
     CreateBranch("nuvtxz", nuvtxz, "nuvtxz[nnuvtx]/F");
     CreateBranch("nuvtxpdg", nuvtxpdg, "nuvtxpdg[nnuvtx]/S");
     }
     */
  if (hasClusterInfo()){
    CreateBranch("NumberOfClusters",&nclusters,"nclusters/S");
    CreateBranch("ClusterID", clusterId, "clusterId[nclusters]/S");
    CreateBranch("Cluster_NumberOfHits", cluster_NHits, "cluster_NHits[nclusters]/S");
    CreateBranch("Cluster_View", clusterView, "clusterView[nclusters]/S");
    CreateBranch("Cluster_ChargeIntegral", cluster_Integral, "cluster_Integral[nclusters]/F");
    CreateBranch("Cluster_ChargeIntegralAveragePerHit", cluster_IntegralAverage, "cluster_IntegralAverage[nclusters]/F");
    //    CreateBranch("Cluster_SummedADC", cluster_SummedADC, "cluster_SummedADC[nclusters]/F");
    //    CreateBranch("Cluster_SummedADCAveragePerHit", cluster_SummedADCaverage, "cluster_SummedADCaverage[nclusters]/F");
    //    CreateBranch("Cluster_MultipleHitDensity", cluster_MultipleHitDensity, "cluster_MultipleHitDensity[nclusters]/F");
    //    CreateBranch("Cluster_Width", cluster_Width, "cluster_Width[nclusters]/F");
    CreateBranch("Cluster_StartChannel", cluster_StartWire, "cluster_StartWire[nclusters]/S");
    CreateBranch("Cluster_StartTick", cluster_StartTick, "cluster_StartTick[nclusters]/S");
    CreateBranch("Cluster_EndChannel", cluster_EndWire, "cluster_EndWire[nclusters]/S");
    CreateBranch("Cluster_EndTick", cluster_EndTick, "cluster_EndTick[nclusters]/S");
    CreateBranch("Cluster_StartCharge", cluster_StartCharge, "cluster_StartCharge[nclusters]/F");
    CreateBranch("Cluster_StartAngle", cluster_StartAngle, "cluster_StartAngle[nclusters]/F");
    CreateBranch("Cluster_EndCharge", cluster_EndCharge, "cluster_EndCharge[nclusters]/F");
    CreateBranch("Cluster_EndAngle", cluster_EndAngle, "cluster_EndAngle[nclusters]/F");

    //    CreateBranch("cluncosmictags_tagger", cluncosmictags_tagger, "cluncosmictags_tagger[nclusters]/S");
    //    CreateBranch("clucosmicscore_tagger", clucosmicscore_tagger, "clucosmicscore_tagger[nclusters]/F");
    //    CreateBranch("clucosmictype_tagger", clucosmictype_tagger, "clucosmictype_tagger[nclusters]/S");
  }

  /*  if (hasFlashInfo()){
      CreateBranch("no_flashes",&no_flashes,"no_flashes/I");
      CreateBranch("flash_time",flash_time,"flash_time[no_flashes]/F");
      CreateBranch("flash_pe",flash_pe,"flash_pe[no_flashes]/F");
      CreateBranch("flash_ycenter",flash_ycenter,"flash_ycenter[no_flashes]/F");
      CreateBranch("flash_zcenter",flash_zcenter,"flash_zcenter[no_flashes]/F");
      CreateBranch("flash_ywidth",flash_ywidth,"flash_ywidth[no_flashes]/F");
      CreateBranch("flash_zwidth",flash_zwidth,"flash_zwidth[no_flashes]/F");
      CreateBranch("flash_timewidth",flash_timewidth,"flash_timewidth[no_flashes]/F");
      }
      */
  /*  if (hasExternCountInfo()){
      CreateBranch("no_ExternCounts",&no_ExternCounts,"no_ExternCounts/I");
      CreateBranch("externcounts_time",externcounts_time,"externcounts_time[no_ExternCounts]/F");
      CreateBranch("externcounts_id",externcounts_id,"externcounts_id[no_ExternCounts]/F");
      }
      */
  /*  if (hasTrackInfo()){
      kNTracker = trackers.size();
      CreateBranch("kNTracker",&kNTracker,"kNTracker/B");
      for(int i=0; i<kNTracker; i++){
      std::string TrackLabel = trackers[i];

  // note that if the tracker data has maximum number of tracks 0,
  // nothing is initialized (branches are not even created)
  TrackData[i].SetAddresses(pTree, TrackLabel, isCosmics);
  } // for trackers
  }
  */
  /*  if (hasVertexInfo()){
      kNVertexAlgos = vertexalgos.size();
      CreateBranch("kNVertexAlgos",&kNVertexAlgos,"kNVertexAlgos/B");
      for(int i=0; i<kNVertexAlgos; i++){
      std::string VertexLabel = vertexalgos[i];

// note that if the tracker data has maximum number of tracks 0,
// nothing is initialized (branches are not even created)
VertexData[i].SetAddresses(pTree, VertexLabel, isCosmics);
} // for trackers
}
*/
if (hasShowerInfo()){
  kNShowerAlgos = showeralgos.size();
  CreateBranch("kNShowerAlgos",&kNShowerAlgos,"kNShowerAlgos/B");
  for(int i=0; i<kNShowerAlgos; i++){
    // note that if the shower data has maximum number of showers 0,
    // nothing is initialized (branches are not even created)
    ShowerData[i].SetAddresses(pTree);
  } // for showers
} // if we have shower algos

if (hasPFParticleInfo()){
  //CreateBranch("kNVertexAlgos",&kNVertexAlgos,"kNVertexAlgos/B"); // What would be the PFParticle equivalent of this? There's only 1 algo!
  PFParticleData.SetAddresses(pTree);
}
/*
   if (hasGenieInfo()){
   CreateBranch("mcevts_truth",&mcevts_truth,"mcevts_truth/I");
   CreateBranch("nuPDG_truth",nuPDG_truth,"nuPDG_truth[mcevts_truth]/I");
   CreateBranch("ccnc_truth",ccnc_truth,"ccnc_truth[mcevts_truth]/I");
   CreateBranch("mode_truth",mode_truth,"mode_truth[mcevts_truth]/I");
   CreateBranch("enu_truth",enu_truth,"enu_truth[mcevts_truth]/F");
   CreateBranch("Q2_truth",Q2_truth,"Q2_truth[mcevts_truth]/F");
   CreateBranch("W_truth",W_truth,"W_truth[mcevts_truth]/F");
   CreateBranch("X_truth",X_truth,"X_truth[mcevts_truth]/F");
   CreateBranch("Y_truth",Y_truth,"Y_truth[mcevts_truth]/F");
   CreateBranch("hitnuc_truth",hitnuc_truth,"hitnuc_truth[mcevts_truth]/I");
   CreateBranch("nuvtxx_truth",nuvtxx_truth,"nuvtxx_truth[mcevts_truth]/F");
   CreateBranch("nuvtxy_truth",nuvtxy_truth,"nuvtxy_truth[mcevts_truth]/F");
   CreateBranch("nuvtxz_truth",nuvtxz_truth,"nuvtxz_truth[mcevts_truth]/F");
   CreateBranch("nu_dcosx_truth",nu_dcosx_truth,"nu_dcosx_truth[mcevts_truth]/F");
   CreateBranch("nu_dcosy_truth",nu_dcosy_truth,"nu_dcosy_truth[mcevts_truth]/F");
   CreateBranch("nu_dcosz_truth",nu_dcosz_truth,"nu_dcosz_truth[mcevts_truth]/F");
   CreateBranch("lep_mom_truth",lep_mom_truth,"lep_mom_truth[mcevts_truth]/F");
   CreateBranch("lep_dcosx_truth",lep_dcosx_truth,"lep_dcosx_truth[mcevts_truth]/F");
   CreateBranch("lep_dcosy_truth",lep_dcosy_truth,"lep_dcosy_truth[mcevts_truth]/F");
   CreateBranch("lep_dcosz_truth",lep_dcosz_truth,"lep_dcosz_truth[mcevts_truth]/F");

   CreateBranch("vx_flux",vx_flux,"vx_flux[mcevts_truth]/F");
   CreateBranch("vy_flux",vy_flux,"vy_flux[mcevts_truth]/F");
   CreateBranch("vz_flux",vz_flux,"vz_flux[mcevts_truth]/F");
   CreateBranch("pdpx_flux",pdpx_flux,"pdpx_flux[mcevts_truth]/F");
   CreateBranch("pdpy_flux",pdpy_flux,"pdpy_flux[mcevts_truth]/F");
   CreateBranch("pdpz_flux",pdpz_flux,"pdpz_flux[mcevts_truth]/F");
   CreateBranch("ppdxdz_flux",ppdxdz_flux,"ppdxdz_flux[mcevts_truth]/F");
   CreateBranch("ppdydz_flux",ppdydz_flux,"ppdydz_flux[mcevts_truth]/F");
   CreateBranch("pppz_flux",pppz_flux,"pppz_flux[mcevts_truth]/F");
   CreateBranch("ptype_flux",ptype_flux,"ptype_flux[mcevts_truth]/I");
   CreateBranch("ppvx_flux",ppvx_flux,"ppvx_flux[mcevts_truth]/F");
   CreateBranch("ppvy_flux",ppvy_flux,"ppvy_flux[mcevts_truth]/F");
   CreateBranch("ppvz_flux",ppvz_flux,"ppvz_flux[mcevts_truth]/F");
   CreateBranch("muparpx_flux",muparpx_flux,"muparpx_flux[mcevts_truth]/F");
   CreateBranch("muparpy_flux",muparpy_flux,"muparpy_flux[mcevts_truth]/F");
   CreateBranch("muparpz_flux",muparpz_flux,"muparpz_flux[mcevts_truth]/F");
   CreateBranch("mupare_flux",mupare_flux,"mupare_flux[mcevts_truth]/F");
   CreateBranch("tgen_flux",tgen_flux,"tgen_flux[mcevts_truth]/I");
   CreateBranch("tgptype_flux",tgptype_flux,"tgptype_flux[mcevts_truth]/I");
   CreateBranch("tgppx_flux",tgppx_flux,"tgppx_flux[mcevts_truth]/F");
   CreateBranch("tgppy_flux",tgppy_flux,"tgppy_flux[mcevts_truth]/F");
   CreateBranch("tgppz_flux",tgppz_flux,"tgppz_flux[mcevts_truth]/F");
   CreateBranch("tprivx_flux",tprivx_flux,"tprivx_flux[mcevts_truth]/F");
   CreateBranch("tprivy_flux",tprivy_flux,"tprivy_flux[mcevts_truth]/F");
   CreateBranch("tprivz_flux",tprivz_flux,"tprivz_flux[mcevts_truth]/F");
   CreateBranch("dk2gen_flux",dk2gen_flux,"dk2gen_flux[mcevts_truth]/F");
   CreateBranch("gen2vtx_flux",gen2vtx_flux,"gen2vtx_flux[mcevts_truth]/F");
   CreateBranch("tpx_flux",tpx_flux,"tpx_flux[mcevts_truth]/F");
   CreateBranch("tpy_flux",tpy_flux,"tpy_flux[mcevts_truth]/F");
   CreateBranch("tpz_flux",tpz_flux,"tpz_flux[mcevts_truth]/F");
   CreateBranch("tptype_flux",tptype_flux,"tptype_flux[mcevts_truth]/I");

   CreateBranch("genie_no_primaries",&genie_no_primaries,"genie_no_primaries/I");
   CreateBranch("genie_primaries_pdg",genie_primaries_pdg,"genie_primaries_pdg[genie_no_primaries]/I");
   CreateBranch("genie_Eng",genie_Eng,"genie_Eng[genie_no_primaries]/F");
   CreateBranch("genie_Px",genie_Px,"genie_Px[genie_no_primaries]/F");
   CreateBranch("genie_Py",genie_Py,"genie_Py[genie_no_primaries]/F");
   CreateBranch("genie_Pz",genie_Pz,"genie_Pz[genie_no_primaries]/F");
   CreateBranch("genie_P",genie_P,"genie_P[genie_no_primaries]/F");
   CreateBranch("genie_status_code",genie_status_code,"genie_status_code[genie_no_primaries]/I");
   CreateBranch("genie_mass",genie_mass,"genie_mass[genie_no_primaries]/F");
   CreateBranch("genie_trackID",genie_trackID,"genie_trackID[genie_no_primaries]/I");
   CreateBranch("genie_ND",genie_ND,"genie_ND[genie_no_primaries]/I");
   CreateBranch("genie_mother",genie_mother,"genie_mother[genie_no_primaries]/I");
   }
   */
/*
   if (hasCryInfo()){
   CreateBranch("mcevts_truthcry",&mcevts_truthcry,"mcevts_truthcry/I");
   CreateBranch("cry_no_primaries",&cry_no_primaries,"cry_no_primaries/I");
   CreateBranch("cry_primaries_pdg",cry_primaries_pdg,"cry_primaries_pdg[cry_no_primaries]/I");
   CreateBranch("cry_Eng",cry_Eng,"cry_Eng[cry_no_primaries]/F");
   CreateBranch("cry_Px",cry_Px,"cry_Px[cry_no_primaries]/F");
   CreateBranch("cry_Py",cry_Py,"cry_Py[cry_no_primaries]/F");
   CreateBranch("cry_Pz",cry_Pz,"cry_Pz[cry_no_primaries]/F");
   CreateBranch("cry_P",cry_P,"cry_P[cry_no_primaries]/F");
   CreateBranch("cry_StartPointx",cry_StartPointx,"cry_StartPointx[cry_no_primaries]/F");
   CreateBranch("cry_StartPointy",cry_StartPointy,"cry_StartPointy[cry_no_primaries]/F");
   CreateBranch("cry_StartPointz",cry_StartPointz,"cry_StartPointz[cry_no_primaries]/F");
   CreateBranch("cry_StartPointt",cry_StartPointt,"cry_StartPointt[cry_no_primaries]/F");
   CreateBranch("cry_status_code",cry_status_code,"cry_status_code[cry_no_primaries]/I");
   CreateBranch("cry_mass",cry_mass,"cry_mass[cry_no_primaries]/F");
   CreateBranch("cry_trackID",cry_trackID,"cry_trackID[cry_no_primaries]/I");
   CreateBranch("cry_ND",cry_ND,"cry_ND[cry_no_primaries]/I");
   CreateBranch("cry_mother",cry_mother,"cry_mother[cry_no_primaries]/I");
   }
   */
/*
   if (hasProtoInfo()) {
   CreateBranch("proto_no_primaries",&proto_no_primaries,"proto_no_primaries/I");
   CreateBranch("proto_isGoodParticle",proto_isGoodParticle,"proto_isGoodParticle[proto_no_primaries]/I");
   CreateBranch("proto_vx",proto_vx,"proto_vx[proto_no_primaries]/F");
   CreateBranch("proto_vy",proto_vy,"proto_vy[proto_no_primaries]/F");
   CreateBranch("proto_vz",proto_vz,"proto_vz[proto_no_primaries]/F");
   CreateBranch("proto_t",proto_t,"proto_t[proto_no_primaries]/F");
   CreateBranch("proto_px",proto_px,"proto_px[proto_no_primaries]/F");
   CreateBranch("proto_py",proto_py,"proto_py[proto_no_primaries]/F");
   CreateBranch("proto_pz",proto_pz,"proto_pz[proto_no_primaries]/F");
   CreateBranch("proto_momentum",proto_momentum,"proto_momentum[proto_no_primaries]/F");
   CreateBranch("proto_energy",proto_energy,"proto_energy[proto_no_primaries]/F");
   CreateBranch("proto_pdg",proto_pdg,"proto_pdg[proto_no_primaries]/I");
   }
   */


   if (hasGeneratorInfo()){
   CreateBranch("MCTruth_Generator_NumberOfParticles",&generator_list_size,"generator_list_size/I");
   CreateBranch("MCTruth_Generator_ParticleID",TrackId,"TrackId[generator_list_size]/I");
   CreateBranch("MCTruth_Generator_PDGCode",pdg,"pdg[generator_list_size]/I");
   CreateBranch("MCTruth_Generator_Status",status,"status[generator_list_size]/I");
   CreateBranch("MCTruth_Generator_Mass",Mass,"Mass[generator_list_size]/F");
   CreateBranch("MCTruth_Generator_StartTime",StartT,"StartT[generator_list_size]/F");
   CreateBranch("MCTruth_Generator_StartEnergy",Eng,"Eng[generator_list_size]/F");
   CreateBranch("MCTruth_Generator_StartMomentum",P,"P[generator_list_size]/F");
   CreateBranch("MCTruth_Generator_StartPoint_X",StartPointx,"StartPointx[generator_list_size]/F");
   CreateBranch("MCTruth_Generator_StartPoint_Y",StartPointy,"StartPointy[generator_list_size]/F");
   CreateBranch("MCTruth_Generator_StartPoint_Z",StartPointz,"StartPointz[generator_list_size]/F");
   CreateBranch("MCTruth_Generator_StartMomentum_X",Px,"Px[generator_list_size]/F");
   CreateBranch("MCTruth_Generator_StartMomentum_Y",Py,"Py[generator_list_size]/F");
   CreateBranch("MCTruth_Generator_StartMomentum_Z",Pz,"Pz[generator_list_size]/F");
   CreateBranch("MCTruth_Generator_StartDirection_Theta",theta,"theta[generator_list_size]/F");
   CreateBranch("MCTruth_Generator_StartDirection_Phi",phi,"phi[generator_list_size]/F");
   }

   if (hasPhotonInfo()){
   CreateBranch("MCTruth_GEANT4_NumberOfDetectedPhotons",&numberofphotons,"numberofphotons/I");
   CreateBranch("MCTruth_GEANT4_DetectedPhoton_Channel",photons_channel,"photons_channel[numberofphotons]/F");
   CreateBranch("MCTruth_GEANT4_DetectedPhoton_Time",photons_time,"photons_time[numberofphotons]/F");
   }

   if (hasGeantInfo()){
   CreateBranch("MCTruth_GEANT4_NumberOfParticles",&geant_list_size,"geant_list_size/I");
   CreateBranch("MCTruth_GEANT4_NumberOfPrimaries",&no_primaries,"no_primaries/I");
   CreateBranch("MCTruth_GEANT4_ParticleID",TrackId,"TrackId[geant_list_size]/I");
   CreateBranch("MCTruth_GEANT4_PDGCode",pdg,"pdg[geant_list_size]/I");

   CreateBranch("MCTruth_GEANT4_Status",status,"status[geant_list_size]/I");
   CreateBranch("MCTruth_GEANT4_IsInTPCAV",inTPCActive,"inTPCActive[geant_list_size]/I");
   CreateBranch("MCTruth_GEANT4_NumberOfDaughterParticles",NumberDaughters,"NumberDaughters[geant_list_size]/I");
   CreateBranch("MCTruth_GEANT4_MotherParticle",Mother,"Mother[geant_list_size]/I");
   CreateBranch("MCTruth_GEANT4_Mass",Mass,"Mass[geant_list_size]/F");
   CreateBranch("MCTruth_GEANT4_StartPoint_X",StartPointx,"StartPointx[geant_list_size]/F");
   CreateBranch("MCTruth_GEANT4_StartPoint_Y",StartPointy,"StartPointy[geant_list_size]/F");
   CreateBranch("MCTruth_GEANT4_StartPoint_Z",StartPointz,"StartPointz[geant_list_size]/F");
   CreateBranch("MCTruth_GEANT4_StartTime",StartT,"StartT[geant_list_size]/F");
   CreateBranch("MCTruth_GEANT4_StartEnergy",Eng,"Eng[geant_list_size]/F");
   CreateBranch("MCTruth_GEANT4_StartMomentum",P,"P[geant_list_size]/F");
   CreateBranch("MCTruth_GEANT4_StartMomentum_X",Px,"Px[geant_list_size]/F");
   CreateBranch("MCTruth_GEANT4_StartMomentum_Y",Py,"Py[geant_list_size]/F");
   CreateBranch("MCTruth_GEANT4_StartMomentum_Z",Pz,"Pz[geant_list_size]/F");
   CreateBranch("MCTruth_GEANT4_StartDirection_Theta",theta,"theta[geant_list_size]/F");
   CreateBranch("MCTruth_GEANT4_StartDirection_Phi",phi,"phi[geant_list_size]/F");
   CreateBranch("MCTruth_GEANT4_InTPCAV_NumberOfParticles",&geant_list_size_in_tpcAV,"geant_list_size_in_tpcAV/I");
//   CreateBranch("MCTruth_EndEnergy",EndE,"EndE[geant_list_size]/F");
//   CreateBranch("MCTruth_EndPoint_X",EndPointx,"EndPointx[geant_list_size]/F");
//   CreateBranch("MCTruth_EndPoint_Y",EndPointy,"EndPointy[geant_list_size]/F");
//   CreateBranch("MCTruth_EndPoint_Z",EndPointz,"EndPointz[geant_list_size]/F");
//   CreateBranch("MCTruth_EndTime",EndT,"EndT[geant_list_size]/F");

//   CreateBranch("theta_xz",theta_xz,"theta_xz[geant_list_size]/F");
//   CreateBranch("theta_yz",theta_yz,"theta_yz[geant_list_size]/F");

   if (hasGeantInAVInfo()){
   CreateBranch("MCTruth_GEANT4_InTPCAV_ParticleID",TrackId_tpcAV,"TrackId_tpcAV[geant_list_size_in_tpcAV]/I");
   CreateBranch("MCTruth_GEANT4_InTPCAV_PDGCode",PDGCode_tpcAV,"PDGCode_tpcAV[geant_list_size_in_tpcAV]/I");

   CreateBranch("MCTruth_GEANT4_InTPCAV_Pathlength",pathlen_tpcAV,"pathlen_tpcAV[geant_list_size_in_tpcAV]/F");
   CreateBranch("MCTruth_GEANT4_InTPCAV_StartPoint_X",StartPointx_tpcAV,"StartPointx_tpcAV[geant_list_size_in_tpcAV]/F");
   CreateBranch("MCTruth_GEANT4_InTPCAV_StartPoint_Y",StartPointy_tpcAV,"StartPointy_tpcAV[geant_list_size_in_tpcAV]/F");
   CreateBranch("MCTruth_GEANT4_InTPCAV_StartPoint_Z",StartPointz_tpcAV,"StartPointz_tpcAV[geant_list_size_in_tpcAV]/F");
   CreateBranch("MCTruth_GEANT4_InTPCAV_StartTime",StartT_tpcAV,"StartT_tpcAV[geant_list_size_in_tpcAV]/F");
   CreateBranch("MCTruth_GEANT4_InTPCAV_StartEnergy",StartE_tpcAV,"StartE_tpcAV[geant_list_size_in_tpcAV]/F");
   CreateBranch("MCTruth_GEANT4_InTPCAV_StartMomentum",StartP_tpcAV,"StartP_tpcAV[geant_list_size_in_tpcAV]/F");
   CreateBranch("MCTruth_GEANT4_InTPCAV_StartMomentum_X",StartPx_tpcAV,"StartPx_tpcAV[geant_list_size_in_tpcAV]/F");
   CreateBranch("MCTruth_GEANT4_InTPCAV_StartMomentum_Y",StartPy_tpcAV,"StartPy_tpcAV[geant_list_size_in_tpcAV]/F");
   CreateBranch("MCTruth_GEANT4_InTPCAV_StartMomentum_Z",StartPz_tpcAV,"StartPz_tpcAV[geant_list_size_in_tpcAV]/F");
   CreateBranch("MCTruth_GEANT4_InTPCAV_StartDirection_Theta",thetastart_tpcAV,"thetastart_tpcAV[geant_list_size_in_tpcAV]/F");
   CreateBranch("MCTruth_GEANT4_InTPCAV_StartDirection_Phi",phistart_tpcAV,"phistart_tpcAV[geant_list_size_in_tpcAV]/F");
   CreateBranch("MCTruth_GEANT4_InTPCAV_EndPoint_X",EndPointx_tpcAV,"EndPointx_tpcAV[geant_list_size_in_tpcAV]/F");
   CreateBranch("MCTruth_GEANT4_InTPCAV_EndPoint_Y",EndPointy_tpcAV,"EndPointy_tpcAV[geant_list_size_in_tpcAV]/F");
   CreateBranch("MCTruth_GEANT4_InTPCAV_EndPoint_Z",EndPointz_tpcAV,"EndPointz_tpcAV[geant_list_size_in_tpcAV]/F");
   CreateBranch("MCTruth_GEANT4_InTPCAV_EndTime",EndT_tpcAV,"EndT_tpcAV[geant_list_size_in_tpcAV]/F");
   CreateBranch("MCTruth_GEANT4_InTPCAV_EndEnergy",EndE_tpcAV,"EndE_tpcAV[geant_list_size_in_tpcAV]/F");
   CreateBranch("MCTruth_GEANT4_InTPCAV_EndMomentum",EndP_tpcAV,"EndP_tpcAV[geant_list_size_in_tpcAV]/F");
   CreateBranch("MCTruth_GEANT4_InTPCAV_EndMomentum_X",EndPx_tpcAV,"EndPx_tpcAV[geant_list_size_in_tpcAV]/F");
   CreateBranch("MCTruth_GEANT4_InTPCAV_EndMomentum_Y",EndPy_tpcAV,"EndPy_tpcAV[geant_list_size_in_tpcAV]/F");
   CreateBranch("MCTruth_GEANT4_InTPCAV_EndMomentum_Z",EndPz_tpcAV,"EndPz_tpcAV[geant_list_size_in_tpcAV]/F");
   CreateBranch("MCTruth_GEANT4_InTPCAV_EndDirection_Theta",thetaend_tpcAV,"thetaend_tpcAV[geant_list_size_in_tpcAV]/F");
   CreateBranch("MCTruth_GEANT4_InTPCAV_EndDirection_Phi",phiend_tpcAV,"phiend_tpcAV[geant_list_size_in_tpcAV]/F");
}


/*
   CreateBranch("pathlen_drifted",pathlen_drifted,"pathlen_drifted[geant_list_size]/F");
   CreateBranch("inTPCDrifted",inTPCDrifted,"inTPCDrifted[geant_list_size]/I");
   CreateBranch("StartPointx_drifted",StartPointx_drifted,"StartPointx_drifted[geant_list_size]/F");
   CreateBranch("StartPointy_drifted",StartPointy_drifted,"StartPointy_drifted[geant_list_size]/F");
   CreateBranch("StartPointz_drifted",StartPointz_drifted,"StartPointz_drifted[geant_list_size]/F");
   CreateBranch("StartT_drifted",StartT_drifted,"StartT_drifted[geant_list_size]/F");
   CreateBranch("StartE_drifted",StartE_drifted,"StartE_drifted[geant_list_size]/F");
   CreateBranch("StartP_drifted",StartP_drifted,"StartP_drifted[geant_list_size]/F");
   CreateBranch("StartPx_drifted",StartPx_drifted,"StartPx_drifted[geant_list_size]/F");
   CreateBranch("StartPy_drifted",StartPy_drifted,"StartPy_drifted[geant_list_size]/F");
   CreateBranch("StartPz_drifted",StartPz_drifted,"StartPz_drifted[geant_list_size]/F");
   CreateBranch("EndPointx_drifted",EndPointx_drifted,"EndPointx_drifted[geant_list_size]/F");
   CreateBranch("EndPointy_drifted",EndPointy_drifted,"EndPointy_drifted[geant_list_size]/F");
   CreateBranch("EndPointz_drifted",EndPointz_drifted,"EndPointz_drifted[geant_list_size]/F");
   CreateBranch("EndT_drifted",EndT_drifted,"EndT_drifted[geant_list_size]/F");
   CreateBranch("EndE_drifted",EndE_drifted,"EndE_drifted[geant_list_size]/F");
   CreateBranch("EndP_drifted",EndP_drifted,"EndP_drifted[geant_list_size]/F");
   CreateBranch("EndPx_drifted",EndPx_drifted,"EndPx_drifted[geant_list_size]/F");
   CreateBranch("EndPy_drifted",EndPy_drifted,"EndPy_drifted[geant_list_size]/F");
   CreateBranch("EndPz_drifted",EndPz_drifted,"EndPz_drifted[geant_list_size]/F");
*/

/*
   CreateBranch("TrackId",TrackId,"TrackId[geant_list_size]/I");
   CreateBranch("MergedId", MergedId, "MergedId[geant_list_size]/I");
   CreateBranch("origin", origin, "origin[geant_list_size]/I");
   CreateBranch("MCTruthIndex", MCTruthIndex, "MCTruthIndex[geant_list_size]/I");
   CreateBranch("process_primary",process_primary,"process_primary[geant_list_size]/I");
   CreateBranch("processname", processname);
*/
}

if(hasGeantTrajectoryInfo()){
  CreateBranch("MCTruth_GEANT4_NumberOfTrajectoryStepsPerParticle",NTrajectoryPointsPerParticle,"NTrajectoryPointsPerParticle[geant_list_size]/I");
  CreateBranch("MCTruth_GEANT4_TotalNumberOfTrajectoryStepsForAllParticles",&geant_trajectory_size,"geant_trajectory_size/I");

  CreateBranch("MCTruth_GEANT4_TrajectoryStep_ParticleID",TrajTrackId,"TrajTrackId[geant_trajectory_size]/F");
  CreateBranch("MCTruth_GEANT4_TrajectoryStep_PDGCode",TrajPDGCode,"TrajPDGCode[geant_trajectory_size]/F");

  CreateBranch("MCTruth_GEANT4_TrajectoryStep_Point_X",TrajX,"TrajX[geant_trajectory_size]/F");
  CreateBranch("MCTruth_GEANT4_TrajectoryStep_Point_Y",TrajY,"TrajY[geant_trajectory_size]/F");
  CreateBranch("MCTruth_GEANT4_TrajectoryStep_Point_Z",TrajZ,"TrajZ[geant_trajectory_size]/F");
  CreateBranch("MCTruth_GEANT4_TrajectoryStep_Time",TrajT,"TrajT[geant_trajectory_size]/F");
  CreateBranch("MCTruth_GEANT4_TrajectoryStep_Energy",TrajE,"TrajE[geant_trajectory_size]/F");
  CreateBranch("MCTruth_GEANT4_TrajectoryStep_Momentum",TrajP,"TrajP[geant_trajectory_size]/F");
  CreateBranch("MCTruth_GEANT4_TrajectoryStep_Momentum_X",TrajPx,"TrajPx[geant_trajectory_size]/F");
  CreateBranch("MCTruth_GEANT4_TrajectoryStep_Momentum_Y",TrajPy,"TrajPy[geant_trajectory_size]/F");
  CreateBranch("MCTruth_GEANT4_TrajectoryStep_Momentum_Z",TrajPz,"TrajPz[geant_trajectory_size]/F");
  CreateBranch("MCTruth_GEANT4_TrajectoryStep_Direction_Theta",TrajTheta,"TrajTheta[geant_trajectory_size]/F");
  CreateBranch("MCTruth_GEANT4_TrajectoryStep_Direction_Phi",TrajPhi,"TrajPhi[geant_trajectory_size]/F");
}

/*
   if (hasMCShowerInfo()){
   CreateBranch("no_mcshowers",&no_mcshowers,"no_mcshowers/I");
   CreateBranch("mcshwr_origin",mcshwr_origin,"mcshwr_origin[no_mcshowers]/I");
   CreateBranch("mcshwr_pdg",mcshwr_pdg,"mcshwr_pdg[no_mcshowers]/I");
   CreateBranch("mcshwr_TrackId",mcshwr_TrackId,"mcshwr_TrackId[no_mcshowers]/I");
   CreateBranch("mcshwr_Process",mcshwr_Process);
   CreateBranch("mcshwr_startX",mcshwr_startX,"mcshwr_startX[no_mcshowers]/F");
   CreateBranch("mcshwr_startY",mcshwr_startY,"mcshwr_startY[no_mcshowers]/F");
   CreateBranch("mcshwr_startZ",mcshwr_startZ,"mcshwr_startZ[no_mcshowers]/F");
   CreateBranch("mcshwr_endX",mcshwr_endX,"mcshwr_endX[no_mcshowers]/F");
   CreateBranch("mcshwr_endY",mcshwr_endY,"mcshwr_endY[no_mcshowers]/F");
   CreateBranch("mcshwr_endZ",mcshwr_endZ,"mcshwr_endZ[no_mcshowers]/F");
   CreateBranch("mcshwr_CombEngX",mcshwr_CombEngX,"mcshwr_CombEngX[no_mcshowers]/F");
   CreateBranch("mcshwr_CombEngY",mcshwr_CombEngY,"mcshwr_CombEngY[no_mcshowers]/F");
   CreateBranch("mcshwr_CombEngZ",mcshwr_CombEngZ,"mcshwr_CombEngZ[no_mcshowers]/F");
   CreateBranch("mcshwr_CombEngPx",mcshwr_CombEngPx,"mcshwr_CombEngPx[no_mcshowers]/F");
   CreateBranch("mcshwr_CombEngPy",mcshwr_CombEngPy,"mcshwr_CombEngPy[no_mcshowers]/F");
   CreateBranch("mcshwr_CombEngPz",mcshwr_CombEngPz,"mcshwr_CombEngPz[no_mcshowers]/F");
   CreateBranch("mcshwr_CombEngE",mcshwr_CombEngE,"mcshwr_CombEngE[no_mcshowers]/F");
   CreateBranch("mcshwr_dEdx",mcshwr_dEdx,"mcshwr_dEdx[no_mcshowers]/F");
   CreateBranch("mcshwr_StartDirX",mcshwr_StartDirX,"mcshwr_StartDirX[no_mcshowers]/F");
   CreateBranch("mcshwr_StartDirY",mcshwr_StartDirY,"mcshwr_StartDirY[no_mcshowers]/F");
   CreateBranch("mcshwr_StartDirZ",mcshwr_StartDirZ,"mcshwr_StartDirZ[no_mcshowers]/F");
   CreateBranch("mcshwr_isEngDeposited",mcshwr_isEngDeposited,"mcshwr_isEngDeposited[no_mcshowers]/I");
   CreateBranch("mcshwr_Motherpdg",mcshwr_Motherpdg,"mcshwr_Motherpdg[no_mcshowers]/I");
   CreateBranch("mcshwr_MotherTrkId",mcshwr_MotherTrkId,"mcshwr_MotherTrkId[no_mcshowers]/I");
   CreateBranch("mcshwr_MotherProcess",mcshwr_MotherProcess);
   CreateBranch("mcshwr_MotherstartX",mcshwr_MotherstartX,"mcshwr_MotherstartX[no_mcshowers]/F");
   CreateBranch("mcshwr_MotherstartY",mcshwr_MotherstartY,"mcshwr_MotherstartY[no_mcshowers]/F");
   CreateBranch("mcshwr_MotherstartZ",mcshwr_MotherstartZ,"mcshwr_MotherstartZ[no_mcshowers]/F");
   CreateBranch("mcshwr_MotherendX",mcshwr_MotherendX,"mcshwr_MotherendX[no_mcshowers]/F");
   CreateBranch("mcshwr_MotherendY",mcshwr_MotherendY,"mcshwr_MotherendY[no_mcshowers]/F");
   CreateBranch("mcshwr_MotherendZ",mcshwr_MotherendZ,"mcshwr_MotherendZ[no_mcshowers]/F");
   CreateBranch("mcshwr_Ancestorpdg",mcshwr_Ancestorpdg,"mcshwr_Ancestorpdg[no_mcshowers]/I");
   CreateBranch("mcshwr_AncesotorTrkId",mcshwr_AncestorTrkId,"mcshwr_AncestorTrkId[no_mcshowers]/I");
   CreateBranch("mcshwr_AncesotorProcess",mcshwr_AncestorProcess);
   CreateBranch("mcshwr_AncestorstartX",mcshwr_AncestorstartX,"mcshwr_AncestorstartX[no_mcshowers]/F");
   CreateBranch("mcshwr_AncestorstartY",mcshwr_AncestorstartY,"mcshwr_AncestorstartY[no_mcshowers]/F");
   CreateBranch("mcshwr_AncestorstartZ",mcshwr_AncestorstartZ,"mcshwr_AncestorstartZ[no_mcshowers]/F");
   CreateBranch("mcshwr_AncestorendX",mcshwr_AncestorendX,"mcshwr_AncestorendX[no_mcshowers]/F");
   CreateBranch("mcshwr_AncestorendY",mcshwr_AncestorendY,"mcshwr_AncestorendY[no_mcshowers]/F");
   CreateBranch("mcshwr_AncestorendZ",mcshwr_AncestorendZ,"mcshwr_AncestorendZ[no_mcshowers]/F");
   }
   */
/*
   if (hasMCTrackInfo()){
   CreateBranch("no_mctracks",&no_mctracks,"no_mctracks/I");
   CreateBranch("mctrk_origin",mctrk_origin,"mctrk_origin[no_mctracks]/I");
   CreateBranch("mctrk_pdg",mctrk_pdg,"mctrk_pdg[no_mctracks]/I");
   CreateBranch("mctrk_TrackId",mctrk_TrackId,"mctrk_TrackId[no_mctracks]/I");
   CreateBranch("mctrk_Process",mctrk_Process);
   CreateBranch("mctrk_startX",mctrk_startX,"mctrk_startX[no_mctracks]/F");
   CreateBranch("mctrk_startY",mctrk_startY,"mctrk_startY[no_mctracks]/F");
   CreateBranch("mctrk_startZ",mctrk_startZ,"mctrk_startZ[no_mctracks]/F");
   CreateBranch("mctrk_endX",mctrk_endX,"mctrk_endX[no_mctracks]/F");
   CreateBranch("mctrk_endY",mctrk_endY,"mctrk_endY[no_mctracks]/F");
   CreateBranch("mctrk_endZ",mctrk_endZ,"mctrk_endZ[no_mctracks]/F");
   CreateBranch("mctrk_startX_drifted",mctrk_startX_drifted,"mctrk_startX_drifted[no_mctracks]/F");
   CreateBranch("mctrk_startY_drifted",mctrk_startY_drifted,"mctrk_startY_drifted[no_mctracks]/F");
   CreateBranch("mctrk_startZ_drifted",mctrk_startZ_drifted,"mctrk_startZ_drifted[no_mctracks]/F");
   CreateBranch("mctrk_endX_drifted",mctrk_endX_drifted,"mctrk_endX_drifted[no_mctracks]/F");
   CreateBranch("mctrk_endY_drifted",mctrk_endY_drifted,"mctrk_endY_drifted[no_mctracks]/F");
   CreateBranch("mctrk_endZ_drifted",mctrk_endZ_drifted,"mctrk_endZ_drifted[no_mctracks]/F");
   CreateBranch("mctrk_len_drifted",mctrk_len_drifted,"mctrk_len_drifted[no_mctracks]/F");
   CreateBranch("mctrk_p_drifted",mctrk_p_drifted,"mctrk_p_drifted[no_mctracks]/F");
   CreateBranch("mctrk_px_drifted",mctrk_px_drifted,"mctrk_px_drifted[no_mctracks]/F");
   CreateBranch("mctrk_py_drifted",mctrk_py_drifted,"mctrk_py_drifted[no_mctracks]/F");
   CreateBranch("mctrk_pz_drifted",mctrk_pz_drifted,"mctrk_pz_drifted[no_mctracks]/F");
   CreateBranch("mctrk_Motherpdg",mctrk_Motherpdg,"mctrk_Motherpdg[no_mctracks]/I");
   CreateBranch("mctrk_MotherTrkId",mctrk_MotherTrkId,"mctrk_MotherTrkId[no_mctracks]/I");
   CreateBranch("mctrk_MotherProcess",mctrk_MotherProcess);
   CreateBranch("mctrk_MotherstartX",mctrk_MotherstartX,"mctrk_MotherstartX[no_mctracks]/F");
   CreateBranch("mctrk_MotherstartY",mctrk_MotherstartY,"mctrk_MotherstartY[no_mctracks]/F");
   CreateBranch("mctrk_MotherstartZ",mctrk_MotherstartZ,"mctrk_MotherstartZ[no_mctracks]/F");
   CreateBranch("mctrk_MotherendX",mctrk_MotherendX,"mctrk_MotherendX[no_mctracks]/F");
   CreateBranch("mctrk_MotherendY",mctrk_MotherendY,"mctrk_MotherendY[no_mctracks]/F");
   CreateBranch("mctrk_MotherendZ",mctrk_MotherendZ,"mctrk_MotherendZ[no_mctracks]/F");
   CreateBranch("mctrk_Ancestorpdg",mctrk_Ancestorpdg,"mctrk_Ancestorpdg[no_mctracks]/I");
   CreateBranch("mctrk_AncesotorTrkId",mctrk_AncestorTrkId,"mctrk_AncestorTrkId[no_mctracks]/I");
   CreateBranch("mctrk_AncesotorProcess",mctrk_AncestorProcess);
   CreateBranch("mctrk_AncestorstartX",mctrk_AncestorstartX,"mctrk_AncestorstartX[no_mctracks]/F");
   CreateBranch("mctrk_AncestorstartY",mctrk_AncestorstartY,"mctrk_AncestorstartY[no_mctracks]/F");
   CreateBranch("mctrk_AncestorstartZ",mctrk_AncestorstartZ,"mctrk_AncestorstartZ[no_mctracks]/F");
   CreateBranch("mctrk_AncestorendX",mctrk_AncestorendX,"mctrk_AncestorendX[no_mctracks]/F");
   CreateBranch("mctrk_AncestorendY",mctrk_AncestorendY,"mctrk_AncestorendY[no_mctracks]/F");
   CreateBranch("mctrk_AncestorendZ",mctrk_AncestorendZ,"mctrk_AncestorendZ[no_mctracks]/F");
   }
   */
/*
   if (hasAuxDetector()) {
// Geant information is required to fill aux detector information.
// if fSaveGeantInfo is not set to true, show an error message and quit!
if (!hasGeantInfo()){
throw art::Exception(art::errors::Configuration)
<< "Saving Auxiliary detector information requies saving GEANT information, "
<<"please set fSaveGeantInfo flag to true in your fhicl file and rerun.\n";
}
std::ostringstream sstr;
sstr << "[" << kMaxAuxDets << "]";
std::string MaxAuxDetIndexStr = sstr.str();
CreateBranch("NAuxDets",     NAuxDets, "NAuxDets[geant_list_size]/s");
CreateBranch("AuxDetID",     AuxDetID, "AuxDetID[geant_list_size]" + MaxAuxDetIndexStr + "/S");
CreateBranch("AuxDetEntryX", entryX,   "AuxDetEntryX[geant_list_size]" + MaxAuxDetIndexStr + "/F");
CreateBranch("AuxDetEntryY", entryY,   "AuxDetEntryY[geant_list_size]" + MaxAuxDetIndexStr + "/F");
CreateBranch("AuxDetEntryZ", entryZ,   "AuxDetEntryZ[geant_list_size]" + MaxAuxDetIndexStr + "/F");
CreateBranch("AuxDetEntryT", entryT,   "AuxDetEntryT[geant_list_size]" + MaxAuxDetIndexStr + "/F");
CreateBranch("AuxDetExitX",  exitX,    "AuxDetExitX[geant_list_size]"  + MaxAuxDetIndexStr + "/F");
CreateBranch("AuxDetExitY",  exitY,    "AuxDetExitY[geant_list_size]"  + MaxAuxDetIndexStr + "/F");
CreateBranch("AuxDetExitZ",  exitZ,    "AuxDetExitZ[geant_list_size]"  + MaxAuxDetIndexStr + "/F");
CreateBranch("AuxDetExitT",  exitT,    "AuxDetExitT[geant_list_size]"  + MaxAuxDetIndexStr + "/F");
CreateBranch("AuxDetExitPx", exitPx,   "AuxDetExitPx[geant_list_size]" + MaxAuxDetIndexStr + "/F");
CreateBranch("AuxDetExitPy", exitPy,   "AuxDetExitPy[geant_list_size]" + MaxAuxDetIndexStr + "/F");
CreateBranch("AuxDetExitPz", exitPz,   "AuxDetExitPz[geant_list_size]" + MaxAuxDetIndexStr + "/F");
CreateBranch("CombinedEnergyDep", CombinedEnergyDep,
"CombinedEnergyDep[geant_list_size]" + MaxAuxDetIndexStr + "/F");
} // if hasAuxDetector
*/
} // dune::AnaRootParserDataStruct::SetAddresses()


//------------------------------------------------------------------------------
//---  AnaRootParser
//---

dune::AnaRootParser::AnaRootParser(fhicl::ParameterSet const& pset) :
  EDAnalyzer(pset),
  fTree(nullptr),
  //  fPOT(nullptr),

  fLogLevel	            (pset.get< short >("LogLevel")        ),
  fEventsPerSubrun          (pset.get< short >("EventsPerSubrun")        ),
  fRawDigitModuleLabel         (pset.get< std::string >("RawDigitModuleLabel")        ),
  fHitsModuleLabel          (pset.get< std::string >("HitsModuleLabel")         ),
  fLArG4ModuleLabel         (pset.get< std::string >("LArGeantModuleLabel")     ),
  fCalDataModuleLabel       (pset.get< std::string >("CalDataModuleLabel")      ),
  fGenieGenModuleLabel      (pset.get< std::string >("GenieGenModuleLabel")     ),
  fCryGenModuleLabel        (pset.get< std::string >("CryGenModuleLabel")       ),
  fProtoGenModuleLabel      (pset.get< std::string >("ProtoGenModuleLabel")     ),
  fG4ModuleLabel            (pset.get< std::string >("G4ModuleLabel")           ),
  fClusterModuleLabel       (pset.get< std::string >("ClusterModuleLabel")     ),
  fPandoraNuVertexModuleLabel (pset.get< std::string >("PandoraNuVertexModuleLabel")     ),
  fOpFlashModuleLabel       (pset.get< std::string >("OpFlashModuleLabel")      ),
  fExternalCounterModuleLabel (pset.get< std::string >("ExternalCounterModuleLabel")      ),
  fMCShowerModuleLabel      (pset.get< std::string >("MCShowerModuleLabel")     ),
  fMCTrackModuleLabel      (pset.get< std::string >("MCTrackModuleLabel")     ),
  fTrackModuleLabel         (pset.get< std::vector<std::string> >("TrackModuleLabel")),
  fVertexModuleLabel        (pset.get< std::vector<std::string> >("VertexModuleLabel")),
  fShowerModuleLabel        (pset.get< std::vector<std::string> >("ShowerModuleLabel")),
  fCalorimetryModuleLabel   (pset.get< std::vector<std::string> >("CalorimetryModuleLabel")),
  fParticleIDModuleLabel    (pset.get< std::vector<std::string> >("ParticleIDModuleLabel")   ),
  fMVAPIDShowerModuleLabel  (pset.get< std::vector<std::string> >("MVAPIDShowerModuleLabel")   ),
  fMVAPIDTrackModuleLabel   (pset.get< std::vector<std::string> >("MVAPIDTrackModuleLabel")   ),
  fFlashT0FinderLabel       (pset.get< std::vector<std::string> >("FlashT0FinderLabel")   ),
  fMCT0FinderLabel          (pset.get< std::vector<std::string> >("MCT0FinderLabel")   ),
  fPOTModuleLabel           (pset.get< std::string >("POTModuleLabel")),
  fCosmicClusterTaggerAssocLabel (pset.get< std::string >("CosmicClusterTaggerAssocLabel")),
  fIsMC                	    (pset.get< bool >("IsMC", false)),
  fUseBuffer                (pset.get< bool >("UseBuffers", false)),
  fSaveAuxDetInfo           (pset.get< bool >("SaveAuxDetInfo", false)),
  fSaveCryInfo              (pset.get< bool >("SaveCryInfo", false)),
  fSaveGenieInfo	    (pset.get< bool >("SaveGenieInfo", false)),
  fSaveProtoInfo	    (pset.get< bool >("SaveProtoInfo", false)),
  fSavePhotonInfo	    (pset.get< bool >("SavePhotonInfo", false)),
  fSaveGeneratorInfo	    (pset.get< bool >("SaveGeneratorInfo", false)),
  fSaveGeantInfo	    (pset.get< bool >("SaveGeantInfo", false)),
  fSaveGeantInAVInfo	    (pset.get< bool >("SaveGeantInAVInfo", false)),
  fSaveGeantTrajectoryInfo  (pset.get< bool >("SaveGeantTrajectoryInfo", false)),
  fSaveMCShowerInfo	    (pset.get< bool >("SaveMCShowerInfo", false)),
  fSaveMCTrackInfo	    (pset.get< bool >("SaveMCTrackInfo", false)),
  fSaveHitInfo	            (pset.get< bool >("SaveHitInfo", false)),
  fSaveRawDigitInfo	    (pset.get< bool >("SaveRawDigitInfo", false)),
  fSaveRecobWireInfo	    (pset.get< bool >("SaveRecobWireInfo", false)),
  fSaveTrackInfo	    (pset.get< bool >("SaveTrackInfo", false)),
  fSaveVertexInfo	    (pset.get< bool >("SaveVertexInfo", false)),
  fSaveClusterInfo	    (pset.get< bool >("SaveClusterInfo", false)),
  fSavePandoraNuVertexInfo  (pset.get< bool >("SavePandoraNuVertexInfo", false)),
  fSaveFlashInfo            (pset.get< bool >("SaveFlashInfo", false)),
  fSaveExternCounterInfo    (pset.get< bool >("SaveExternCounterInfo", false)),
  fSaveShowerInfo           (pset.get< bool >("SaveShowerInfo", false)),
  fSavePFParticleInfo	    (pset.get< bool >("SavePFParticleInfo", false)),
  fCosmicTaggerAssocLabel   (pset.get<std::vector< std::string > >("CosmicTaggerAssocLabel") ),
  fContainmentTaggerAssocLabel  (pset.get<std::vector< std::string > >("ContainmentTaggerAssocLabel") ),
  fFlashMatchAssocLabel     (pset.get<std::vector< std::string > >("FlashMatchAssocLabel") ),
  bIgnoreMissingShowers     (pset.get< bool >("IgnoreMissingShowers", false)),
  isCosmics(false),
  fSaveCaloCosmics          (pset.get< bool >("SaveCaloCosmics",false)),
  fG4minE                   (pset.get< float>("G4minE",0.01)),
  fMCSFitter( pset.get<fhicl::ParameterSet>("TrajMCSFitter") )

{

  if (fSavePFParticleInfo) fPFParticleModuleLabel = pset.get<std::string>("PFParticleModuleLabel");

  if (fSaveAuxDetInfo == true) fSaveGeantInfo = true;
//  if (fSaveRawDigitInfo == true) fSaveHitInfo = true;
  mf::LogInfo("AnaRootParser") << "Configuration:"
    << "\n  UseBuffers: " << std::boolalpha << fUseBuffer
    ;
  if (GetNTrackers() > kMaxTrackers) {
    throw art::Exception(art::errors::Configuration)
      << "AnaRootParser currently supports only up to " << kMaxTrackers
      << " tracking algorithms, but " << GetNTrackers() << " are specified."
      << "\nYou can increase kMaxTrackers and recompile.";
  } // if too many trackers
  if (fTrackModuleLabel.size() != fCalorimetryModuleLabel.size()){
    throw art::Exception(art::errors::Configuration)
      << "fTrackModuleLabel.size() = "<<fTrackModuleLabel.size()<<" does not match "
      << "fCalorimetryModuleLabel.size() = "<<fCalorimetryModuleLabel.size();
  }
  if (fTrackModuleLabel.size() != fParticleIDModuleLabel.size()){
    throw art::Exception(art::errors::Configuration)
      << "fTrackModuleLabel.size() = "<<fTrackModuleLabel.size()<<" does not match "
      << "fParticleIDModuleLabel.size() = "<<fParticleIDModuleLabel.size();
  }
  if (fTrackModuleLabel.size() != fFlashT0FinderLabel.size()){
    throw art::Exception(art::errors::Configuration)
      << "fTrackModuleLabel.size() = "<<fTrackModuleLabel.size()<<" does not match "
      << "fFlashT0FinderLabel.size() = "<<fFlashT0FinderLabel.size();
  }
  if (fTrackModuleLabel.size() != fMCT0FinderLabel.size()){
    throw art::Exception(art::errors::Configuration)
      << "fTrackModuleLabel.size() = "<<fTrackModuleLabel.size()<<" does not match "
      << "fMCT0FinderLabel.size() = "<<fMCT0FinderLabel.size();
  }
  if (fTrackModuleLabel.size() != fCosmicTaggerAssocLabel.size()) {
    throw art::Exception(art::errors::Configuration)
      << "fTrackModuleLabel.size() = "<<fTrackModuleLabel.size()<<" does not match "
      << "fCosmicTaggerAssocLabel.size() = "<<fCosmicTaggerAssocLabel.size();
  }
  if (fTrackModuleLabel.size() != fContainmentTaggerAssocLabel.size()) {
    throw art::Exception(art::errors::Configuration)
      << "fTrackModuleLabel.size() = "<<fTrackModuleLabel.size()<<" does not match "
      << "fCosmicTaggerAssocLabel.size() = "<<fContainmentTaggerAssocLabel.size();
  }
  if (fTrackModuleLabel.size() != fFlashMatchAssocLabel.size()) {
    throw art::Exception(art::errors::Configuration)
      << "fTrackModuleLabel.size() = "<<fTrackModuleLabel.size()<<" does not match "
      << "fCosmicTaggerAssocLabel.size() = "<<fFlashMatchAssocLabel.size();
  }
  if (GetNVertexAlgos() > kMaxVertexAlgos) {
    throw art::Exception(art::errors::Configuration)
      << "AnaRootParser currently supports only up to " << kMaxVertexAlgos
      << " tracking algorithms, but " << GetNVertexAlgos() << " are specified."
      << "\nYou can increase kMaxVertexAlgos and recompile.";
  } // if too many trackers

  // Build my Cryostat boundaries array...Taken from Tyler Alion in Geometry Core. Should still return the same values for uBoone.
  ActiveBounds[0] = ActiveBounds[2] = ActiveBounds[4] = DBL_MAX;
  ActiveBounds[1] = ActiveBounds[3] = ActiveBounds[5] = -DBL_MAX;
  // assume single cryostats
  auto const* geom = lar::providerFrom<geo::Geometry>();
  for (geo::TPCGeo const& TPC: geom->IterateTPCs()) {
    // get center in world coordinates
    double origin[3] = {0.};
    double center[3] = {0.};
    TPC.LocalToWorld(origin, center);
    double tpcDim[3] = {TPC.HalfWidth(), TPC.HalfHeight(), 0.5*TPC.Length() };

    if( center[0] - tpcDim[0] < ActiveBounds[0] ) ActiveBounds[0] = center[0] - tpcDim[0];
    if( center[0] + tpcDim[0] > ActiveBounds[1] ) ActiveBounds[1] = center[0] + tpcDim[0];
    if( center[1] - tpcDim[1] < ActiveBounds[2] ) ActiveBounds[2] = center[1] - tpcDim[1];
    if( center[1] + tpcDim[1] > ActiveBounds[3] ) ActiveBounds[3] = center[1] + tpcDim[1];
    if( center[2] - tpcDim[2] < ActiveBounds[4] ) ActiveBounds[4] = center[2] - tpcDim[2];
    if( center[2] + tpcDim[2] > ActiveBounds[5] ) ActiveBounds[5] = center[2] + tpcDim[2];
  } // for all TPC
  std::cout << "Active Boundaries: "
    << "\n\tx: " << ActiveBounds[0] << " to " << ActiveBounds[1]
    << "\n\ty: " << ActiveBounds[2] << " to " << ActiveBounds[3]
    << "\n\tz: " << ActiveBounds[4] << " to " << ActiveBounds[5]
    << std::endl;
} // dune::AnaRootParser::AnaRootParser()

//-------------------------------------------------
dune::AnaRootParser::~AnaRootParser()
{
  DestroyData();
}

void dune::AnaRootParser::CreateTree(bool bClearData /* = false */) {
  if (!fTree) {
    art::ServiceHandle<art::TFileService> tfs;
    fTree = tfs->make<TTree>("anatree","analysis tree");
  }
  /*  if (!fPOT) {
      art::ServiceHandle<art::TFileService> tfs;
      fPOT = tfs->make<TTree>("pottree","pot tree");
      fPOT->Branch("pot",&SubRunData.pot,"pot/D");
      fPOT->Branch("potbnbETOR860",&SubRunData.potbnbETOR860,"potbnbETOR860/D");
      fPOT->Branch("potbnbETOR875",&SubRunData.potbnbETOR875,"potbnbETOR875/D");
      fPOT->Branch("potnumiETORTGT",&SubRunData.potnumiETORTGT,"potnumiETORTGT/D");
      }*/
  CreateData(bClearData);
  SetAddresses();
} // dune::AnaRootParser::CreateTree()


void dune::AnaRootParser::beginSubRun(const art::SubRun& sr)
{

  //  art::Handle< sumdata::POTSummary > potListHandle;
  //  //sr.getByLabel(fPOTModuleLabel,potListHandle);
  //
  //  if(sr.getByLabel(fPOTModuleLabel,potListHandle))
  //    SubRunData.pot=potListHandle->totpot;
  //  else
  //    SubRunData.pot=0.;

}

void dune::AnaRootParser::endSubRun(const art::SubRun& sr)
{

  /*  art::Handle< sumdata::POTSummary > potListHandle;
      if(sr.getByLabel(fPOTModuleLabel,potListHandle))
      SubRunData.pot=potListHandle->totpot;
      else
      SubRunData.pot=0.;

      art::Handle<sumdata::POTSummary> potSummaryHandlebnbETOR860;
      if (sr.getByLabel("beamdata","bnbETOR860",potSummaryHandlebnbETOR860)){
      SubRunData.potbnbETOR860 = potSummaryHandlebnbETOR860->totpot;
      }
      else
      SubRunData.potbnbETOR860 = 0;

      art::Handle<sumdata::POTSummary> potSummaryHandlebnbETOR875;
      if (sr.getByLabel("beamdata","bnbETOR875",potSummaryHandlebnbETOR875)){
      SubRunData.potbnbETOR875 = potSummaryHandlebnbETOR875->totpot;
      }
      else
      SubRunData.potbnbETOR875 = 0;

      art::Handle<sumdata::POTSummary> potSummaryHandlenumiETORTGT;
      if (sr.getByLabel("beamdata","numiETORTGT",potSummaryHandlenumiETORTGT)){
      SubRunData.potnumiETORTGT = potSummaryHandlenumiETORTGT->totpot;
      }
      else
      SubRunData.potnumiETORTGT = 0;

      if (fPOT) fPOT->Fill();
      */
}

void dune::AnaRootParser::analyze(const art::Event& evt)
{

  //services
  auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
//  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;

  // collect the sizes which might be needed to resize the tree data structure:

  // RawDigit
  art::Handle< std::vector<raw::RawDigit> > rawdigitVecHandle;
  std::vector<art::Ptr<raw::RawDigit> > rawdigitlist;
  if (evt.getByLabel(fRawDigitModuleLabel, rawdigitVecHandle))
    art::fill_ptr_vector(rawdigitlist, rawdigitVecHandle);

  // RecobWire
  art::Handle< std::vector<recob::Wire> > recobwireVecHandle;
  std::vector<art::Ptr<recob::Wire> > recobwirelist;
  if (evt.getByLabel(fCalDataModuleLabel, recobwireVecHandle))
    art::fill_ptr_vector(recobwirelist, recobwireVecHandle);


  // recob::Wire
//  art::Handle< std::vector<recob::Wire> > wireVecHandle;
//  evt.getByLabel(fCalDataModuleLabel,wireVecHandle);

  // * photons
  art::Handle< std::vector<sim::SimPhotonsLite> > photonHandle;
//  std::vector<art::Ptr<sim::SimPhotonsLite> > photonlist;
  evt.getByLabel(fLArG4ModuleLabel, photonHandle);
//    art::fill_ptr_vector(photonlist, photonHandle);

  // * hits
  art::Handle< std::vector<recob::Hit> > hitListHandle;
  std::vector<art::Ptr<recob::Hit> > hitlist;
  if (evt.getByLabel(fHitsModuleLabel,hitListHandle))
    art::fill_ptr_vector(hitlist, hitListHandle);

  // * clusters
  art::Handle< std::vector<recob::Cluster> > clusterListHandle;
  std::vector<art::Ptr<recob::Cluster> > clusterlist;
  if (fSaveClusterInfo){
    if (evt.getByLabel(fClusterModuleLabel,clusterListHandle))
      art::fill_ptr_vector(clusterlist, clusterListHandle);
  }

  // * flashes
  art::Handle< std::vector<recob::OpFlash> > flashListHandle;
  std::vector<art::Ptr<recob::OpFlash> > flashlist;
  if (evt.getByLabel(fOpFlashModuleLabel,flashListHandle))
    art::fill_ptr_vector(flashlist, flashListHandle);

  // * External Counters
  art::Handle< std::vector<raw::ExternalTrigger> > countListHandle;
  std::vector<art::Ptr<raw::ExternalTrigger> > countlist;
  if (evt.getByLabel(fExternalCounterModuleLabel,countListHandle))
    art::fill_ptr_vector(countlist, countListHandle);

  // * MC truth information
  art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
  std::vector<art::Ptr<simb::MCTruth> > mclist;
  if (fIsMC){
    if (evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle))
      art::fill_ptr_vector(mclist, mctruthListHandle);
  }

  // *MC truth cosmic generator information
  art::Handle< std::vector<simb::MCTruth> > mctruthcryListHandle;
  std::vector<art::Ptr<simb::MCTruth> > mclistcry;
  if (fIsMC && fSaveCryInfo){
    if (evt.getByLabel(fCryGenModuleLabel,mctruthcryListHandle)){
      art::fill_ptr_vector(mclistcry, mctruthcryListHandle);
    }
    else{
      // If we requested this info but it doesn't exist, don't use it.
      fSaveCryInfo = false;
      mf::LogError("AnaRootParser") << "Requested CRY information but none exists, hence not saving CRY information.";
    }
  }

  // ProtoDUNE beam generator information
  art::Handle< std::vector<simb::MCTruth> > mctruthprotoListHandle;
  std::vector<art::Ptr<simb::MCTruth> > mclistproto;
  if(fIsMC && fSaveProtoInfo){
    if(evt.getByLabel(fProtoGenModuleLabel,mctruthprotoListHandle)){
      art::fill_ptr_vector(mclistproto,mctruthprotoListHandle);
    }
    else{
      // If we requested this info but it doesn't exist, don't use it.
      fSaveProtoInfo = false;
      mf::LogError("AnaRootParser") << "Requested protoDUNE beam generator information but none exists, hence not saving protoDUNE beam generator information.";
    }
  }

  // *MC Shower information
  art::Handle< std::vector<sim::MCShower> > mcshowerh;
  if (fIsMC)
    evt.getByLabel(fMCShowerModuleLabel, mcshowerh);

  int nMCShowers = 0;
  if (fSaveMCShowerInfo && mcshowerh.isValid())
    nMCShowers = mcshowerh->size();

  // *MC Track information
  art::Handle< std::vector<sim::MCTrack> > mctrackh;
  if (fIsMC)
    evt.getByLabel(fMCTrackModuleLabel, mctrackh);

  int nMCTracks = 0;
  if (fSaveMCTrackInfo && mctrackh.isValid())
    nMCTracks = mctrackh->size();

  art::Ptr<simb::MCTruth> mctruthcry;
  int nCryPrimaries = 0;

  if (fSaveCryInfo){
    mctruthcry = mclistcry[0];
    nCryPrimaries = mctruthcry->NParticles();
  }

  art::Ptr<simb::MCTruth> mctruthproto;
  int nProtoPrimaries = 0;
  if( fSaveProtoInfo){
    mctruthproto = mclistproto[0];
    nProtoPrimaries = mctruthproto->NParticles();
  }

  int nPhotons=0;
  int nGeniePrimaries = 0, nGeneratorParticles = 0, nGEANTparticles = 0, nGEANTparticlesInAV = 0, nGEANTtrajectorysteps=0;

  art::Ptr<simb::MCTruth> mctruth;


  if (fIsMC) { //is MC

    //find origin
    std::vector< art::Handle< std::vector<simb::MCTruth> > > allmclists;
    evt.getManyByType(allmclists);
    for(size_t mcl = 0; mcl < allmclists.size(); ++mcl){
      art::Handle< std::vector<simb::MCTruth> > mclistHandle = allmclists[mcl];
      for(size_t m = 0; m < mclistHandle->size(); ++m){
        art::Ptr<simb::MCTruth> mct(mclistHandle, m);
        if (mct->Origin() == simb::kCosmicRay) isCosmics = true;
      }
    }
    if (fSaveCaloCosmics) isCosmics = false; //override to save calo info

    // GENIE
    if (!mclist.empty()){//at least one mc record

      //        double maxenergy = -1;
      //        int imc0 = 0;
      //        for (std::map<art::Ptr<simb::MCTruth>,double>::iterator ii=mctruthemap.begin(); ii!=mctruthemap.end(); ++ii){
      //          if ((ii->second)>maxenergy){
      //            maxenergy = ii->second;
      //            mctruth = ii->first;
      //            imc = imc0;
      //          }
      //          imc0++;
      //        }

      //imc = 0; //set imc to 0 to solve a confusion for BNB+cosmic files where there are two MCTruth
      mctruth = mclist[0];

      if (mctruth->NeutrinoSet()) nGeniePrimaries = mctruth->NParticles();
      //} //end (fSaveGenieInfo)

      const sim::ParticleList& plist = pi_serv->ParticleList();
      nGEANTparticles = plist.size();

      sim::ParticleList::const_iterator itPart = plist.begin(),
      pend = plist.end(); // iterator to pairs (track id, particle)

      std::string pri("primary");
      for(size_t iPart = 0; (iPart < plist.size()) && (itPart != pend); ++iPart)
      {
        const simb::MCParticle* pPart = (itPart++)->second;
        if (!pPart)
	{
          throw art::Exception(art::errors::LogicError)
            << "GEANT particle #" << iPart << " returned a null pointer";
        }
        nGEANTtrajectorysteps += pPart->NumberTrajectoryPoints();

	if(pPart->Mother() == 0 && pPart->Process() == pri) nGeneratorParticles++;

        TLorentzVector mcstart, mcend;
        unsigned int pstarti, pendi;
        double plen = length(*pPart, mcstart, mcend, pstarti, pendi);
        if( plen != 0) nGEANTparticlesInAV += 1;

      } // for particles

      // counting photons
      for ( auto const& pmt : (*photonHandle) )
      {
        std::map<int, int> PhotonsMap = pmt.DetectedPhotons;

        for(auto itphoton = PhotonsMap.begin(); itphoton!= PhotonsMap.end(); itphoton++)
        {
          for(int i = 0; i < itphoton->second ; i++)
          {
	    nPhotons++;
	  }
	}
      }


//      for()
//      {

//      }

    } // if have MC truth
    LOG_DEBUG("AnaRootParser") << "Expected "
    << nGEANTparticles << " GEANT particles, "
    << nGeniePrimaries << " GENIE particles";
  } // if MC

CreateData(); // tracker data is created with default constructor
if (fSaveGenieInfo){  fData->ResizeGenie(nGeniePrimaries);}
if (fSaveCryInfo){  fData->ResizeCry(nCryPrimaries);}
if (fSaveProtoInfo){  fData->ResizeProto(nProtoPrimaries);}
if (fSaveGeneratorInfo){  fData->ResizeGenerator(nGeneratorParticles);}
if (fSaveGeantInfo){  fData->ResizeGEANT(nGEANTparticles);}
if (fSaveGeantInfo && fSaveGeantInAVInfo){  fData->ResizeGEANTInAV(nGEANTparticlesInAV);}
if (fSaveGeantTrajectoryInfo){  fData->ResizeGEANTTrajectory(nGEANTtrajectorysteps);}
if (fSaveMCShowerInfo){  fData->ResizeMCShower(nMCShowers);}
if (fSaveMCTrackInfo){  fData->ResizeMCTrack(nMCTracks);}
if (fSavePhotonInfo){  fData->ResizePhotons(nPhotons);}

if(fLogLevel == 1)
{
  std::cout << "nGeneratorParticles: " << nGeneratorParticles << std::endl;
  std::cout << "nGEANTparticles: " << nGEANTparticles << std::endl;
  std::cout << "nGEANTtrajectorysteps: " << nGEANTtrajectorysteps << std::endl;
  std::cout << "nGEANTparticlesInAV: " << nGEANTparticlesInAV << std::endl;
  std::cout << "nPhotons: " << nPhotons << std::endl;
}

  fData->ClearLocalData(); // don't bother clearing tracker data yet

  const size_t NTrackers = GetNTrackers(); // number of trackers passed into fTrackModuleLabel
  const size_t NShowerAlgos = GetNShowerAlgos(); // number of shower algorithms into fShowerModuleLabel
  const size_t NChannels = rawdigitlist.size(); //number of channels holding raw waveforms
  const size_t NRecoChannels = recobwirelist.size(); //number of channels holding ROI's
  const size_t NHits     = hitlist.size(); // number of hits
  const size_t NVertexAlgos = GetNVertexAlgos(); // number of vertex algos
  const size_t NClusters = clusterlist.size(); //number of clusters
  const size_t NFlashes  = flashlist.size(); // number of flashes
  const size_t NExternCounts = countlist.size(); // number of External Counters
  // make sure there is the data, the tree and everything;
  CreateTree();

  /// transfer the run and subrun data to the tree data object
  //  fData->RunData = RunData;
  fData->SubRunData = SubRunData;

  fData->isdata = int(!fIsMC);

  //  raw trigger
  art::Handle< std::vector<raw::Trigger>> triggerListHandle;
  std::vector<art::Ptr<raw::Trigger>> triggerlist;
if (evt.getByLabel(fRawDigitModuleLabel, triggerListHandle)){  art::fill_ptr_vector(triggerlist, triggerListHandle);}

  if (triggerlist.size()){
    fData->triggernumber = triggerlist[0]->TriggerNumber();
    fData->triggertime   = triggerlist[0]->TriggerTime();
    fData->beamgatetime  = triggerlist[0]->BeamGateTime();
    fData->triggerbits   = triggerlist[0]->TriggerBits();
  }

//  vertices
std::vector< art::Handle< std::vector<recob::Vertex> > > vertexListHandle(NVertexAlgos);
std::vector< std::vector<art::Ptr<recob::Vertex> > > vertexlist(NVertexAlgos);
for (unsigned int it = 0; it < NVertexAlgos; ++it){
  if (evt.getByLabel(fVertexModuleLabel[it],vertexListHandle[it]))
    art::fill_ptr_vector(vertexlist[it], vertexListHandle[it]);
}

//  PFParticles
lar_pandora::PFParticleVector pfparticlelist;
lar_pandora::PFParticlesToClusters pfParticleToClusterMap;
lar_pandora::LArPandoraHelper::CollectPFParticles(evt, fPFParticleModuleLabel, pfparticlelist, pfParticleToClusterMap);

//  tracks
std::vector< art::Handle< std::vector<recob::Track> > > trackListHandle(NTrackers);
std::vector< std::vector<art::Ptr<recob::Track> > > tracklist(NTrackers);
for (unsigned int it = 0; it < NTrackers; ++it){
  if (evt.getByLabel(fTrackModuleLabel[it],trackListHandle[it]))
    art::fill_ptr_vector(tracklist[it], trackListHandle[it]);
}

//  showers
// It seems that sometimes the shower data product does not exist;
// in that case, we store a nullptr in place of the pointer to the vector;
// the data structure itself is owned by art::Event and we should not try
// to manage its memory
std::vector<std::vector<recob::Shower> const*> showerList;
std::vector< art::Handle< std::vector<recob::Shower> > > showerListHandle;
showerList.reserve(fShowerModuleLabel.size());
if (fSaveShowerInfo) {
  for (art::InputTag ShowerInputTag: fShowerModuleLabel) {
    art::Handle<std::vector<recob::Shower>> ShowerHandle;
    if (!evt.getByLabel(ShowerInputTag, ShowerHandle)) {
      showerList.push_back(nullptr);
      if (!bIgnoreMissingShowers) {
        throw art::Exception(art::errors::ProductNotFound)
          << "Showers with input tag '" << ShowerInputTag.encode()
          << "' were not found in the event."
          " If you really know what you are doing,"
          " set AnaRootParser's configuration parameter IgnoreMissingShowers"
          " to \"true\"; the lack of any shower set will be tolerated,"
          " and the shower list in the corresponding event will be set to"
          " a list of one shower, with an invalid ID.\n";
      } // if bIgnoreMissingShowers
      else {
        // this message is more alarming than strictly necessary; by design.
        mf::LogError("AnaRootParser")
          << "No showers found for input tag '" << ShowerInputTag.encode()
          << "' ; FILLING WITH FAKE DATA AS FOR USER'S REQUEST";
      }
    }

    else showerList.push_back(ShowerHandle.product());

    showerListHandle.push_back(ShowerHandle); // either way, put it into the handle list

  }

} // for shower input tag


// Find the simb::MCFlux objects corresponding to
// each simb::MCTruth object made by the generator with
// the label fGenieGenModuleLabel
//  art::FindOne<simb::MCFlux> find_mcflux(mctruthListHandle, evt, fGenieGenModuleLabel);

std::vector<const sim::AuxDetSimChannel*> fAuxDetSimChannels;
if (fSaveAuxDetInfo){
  evt.getView(fLArG4ModuleLabel, fAuxDetSimChannels);
}

std::vector<const sim::SimChannel*> fSimChannels;
if (fIsMC && fSaveGeantInfo){  evt.getView(fLArG4ModuleLabel, fSimChannels);}

  fData->run = evt.run();
  fData->subrun = evt.subRun();
  fData->event = evt.id().event() + fEventsPerSubrun*evt.subRun();

  art::Timestamp ts = evt.time();
  fData->evttime_seconds = ts.timeHigh();
  fData->evttime_nanoseconds = ts.timeLow();

  //copied from MergeDataPaddles.cxx
  art::Handle< raw::BeamInfo > beam;
  if (evt.getByLabel("beamdata",beam)){
    fData->beamtime = (double)beam->get_t_ms();
    fData->beamtime/=1000.; //in second
    std::map<std::string, std::vector<double>> datamap = beam->GetDataMap();
    if (datamap["E:TOR860"].size()){
      fData->potbnb = datamap["E:TOR860"][0];
    }
    if (datamap["E:TORTGT"].size()){
      fData->potnumitgt = datamap["E:TORTGT"][0];
    }
    if (datamap["E:TOR101"].size()){
      fData->potnumi101 = datamap["E:TOR101"][0];
    }
  }


//  std::cout<<detprop->NumberTimeSamples()<<" "<<detprop->ReadOutWindowSize()<<std::endl;
//  std::cout<<geom->DetHalfHeight()2<<" "<<geom->DetHalfWidth()2<<" "<<geom->DetLength()<<std::endl;
//  std::cout<<geom->Nwires(0)<<" "<<geom->Nwires(1)<<" "<<geom->Nwires(2)<<std::endl;


//RawDigit info
if (fSaveRawDigitInfo){
  fData->no_channels = (int) NChannels;
  if (NChannels>0)
  {
    fData->no_ticks = (int) rawdigitlist[0]->Samples();
    fData->no_ticksinallchannels = fData->no_channels*fData->no_ticks;
  }
  else
  {
  fData->no_ticks = 0;
  fData->no_ticksinallchannels = 0;
  }

  for (int i = 0; i < fData->no_channels ; i++) //loop over channels holding raw waveforms
  {
    fData->rawD_Channel[i] = (int) rawdigitlist[i]->Channel();
    int k=0;

    for (int j = fData->no_ticks*i; j < (fData->no_ticks*(i+1)) && j < kMaxReadoutTicksInAllChannels ; j++) //loop over ticks
    {
      fData->rawD_ADC[j] = rawdigitlist[i]->ADC(k);
      k++;
    }
  }
}

//RecobWire info
if (fSaveRecobWireInfo){
  fData->no_recochannels = (int) NRecoChannels;
  int RecoWTick=0;

  for(int i = 0; i < fData->no_recochannels; i++)  //loop over channels holding reco waveforms
  {
    fData->recoW_NTicks[i]=0;
    fData->recoW_Channel[i] = recobwirelist[i]->Channel();
    const recob::Wire::RegionsOfInterest_t& signalROI = recobwirelist[i]->SignalROI();

    for(const auto& range : signalROI.get_ranges())  //loop over ROIs
    {
      const std::vector<float>& signal = range.data();
      int NTicksInThisROI = range.end_index() - range.begin_index();

      for(int j = 0; j < NTicksInThisROI; j++) //loop over ticks
      {
        fData->recoW_Tick[RecoWTick] = j+range.begin_index();
        fData->recoW_ADC[RecoWTick] = signal.at(j);
	RecoWTick++;
      }
      fData->recoW_NTicks[i] += NTicksInThisROI;
    }
  }
  fData->no_recoticksinallchannels = RecoWTick;
}

//hit information
if (fSaveHitInfo){

  int trueID = -9999;
  float maxe = -9999;
  float purity = -9999;

  fData->no_hits = (int) NHits;
  fData->NHitsInAllTracks = (int) NHits;
  fData->no_hits_stored = TMath::Min( (int) NHits, (int) kMaxHits);
  if (NHits > kMaxHits) {
    // got this error? consider increasing kMaxHits
    // (or ask for a redesign using vectors)
    mf::LogError("AnaRootParser:limits") << "event has " << NHits
      << " hits, only kMaxHits=" << kMaxHits << " stored in tree";
  }

// trying to assign a space point to those recob::Hits t hat were assigned to a track, checking for every hit. Not really working yet.
//    art::FindManyP<recob::SpacePoint> fmsp(hitlist,evt,"pmtrack");

  auto hitResults = anab::FVectorReader<recob::Hit, 4>::create(evt, "dprawhit");
  const auto & fitParams = hitResults->vectors();

  for (size_t i = 0; i < NHits && i < kMaxHits ; ++i){//loop over hits
    fData->hit_channel[i] = hitlist[i]->Channel();
    fData->hit_tpc[i]   = hitlist[i]->WireID().TPC;
    fData->hit_view[i]   = hitlist[i]->WireID().Plane;
    fData->hit_wire[i]    = hitlist[i]->WireID().Wire;
    fData->hit_peakT[i]   = hitlist[i]->PeakTime();
    fData->hit_chargesum[i]  = hitlist[i]->SummedADC();
    fData->hit_chargeintegral[i]  = hitlist[i]->Integral();
    fData->hit_ph[i]  = hitlist[i]->PeakAmplitude();
    fData->hit_startT[i] = hitlist[i]->StartTick();
    fData->hit_endT[i] = hitlist[i]->EndTick();
    fData->hit_rms[i] = hitlist[i]->RMS();
    fData->hit_goodnessOfFit[i] = hitlist[i]->GoodnessOfFit();

    fData->hit_fitparamt0[i] = fitParams[i][0];
    fData->hit_fitparamtau1[i] = fitParams[i][1];
    fData->hit_fitparamtau2[i] = fitParams[i][2];
    fData->hit_fitparamampl[i] = fitParams[i][3];

    fData->hit_multiplicity[i] = hitlist[i]->Multiplicity();

    //quantities from backtracker are different from the real value in MC

    if( fIsMC )
      HitsBackTrack(  hitlist[i], trueID, maxe, purity );

    fData->hit_trueID[i] = trueID;
    fData->hit_trueEnergyMax[i] = maxe;
    fData->hit_trueEnergyFraction[i] = purity;

/*
    std::vector< art::Ptr<recob::SpacePoint> > sptv = fmsp.at(i);
        if (fmsp.at(i).size()!=0){
    std::cout << "sptv[i]->XYZ()[0]: " << sptv[i]->XYZ()[0] << std::endl;
    std::cout << "sptv[i]->XYZ()[1]: " << sptv[i]->XYZ()[1] << std::endl;
    std::cout << "sptv[i]->XYZ()[2]: " << sptv[i]->XYZ()[2] << std::endl;
      }
    std::cout << "test3" << std::endl;
//    std::cout << "test" << std::endl;
//    art::FindManyP<recob::SpacePoint> fmsp(hitListHandle,evt,"pmtrack");
//    art::FindManyP<recob::SpacePoint> fmsp(allHits, evt, fSpacePointModuleLabel);
//    art::FindManyP<recob::SpacePoint> fmsp(hitlist,evt,"pmtrack");
    for (size_t i = 0; i < NHits && i < kMaxHits ; ++i){//loop over hits
      if (fmsp.isValid()){
        if (fmsp.at(i).size()!=0){

        std::cout << "Spacepoint in x for hit" << i << ": " << fmsp.at(i)[0]->XYZ()[0] << std::endl;
        //  fData->hit_clusterid[i] = fmcl.at(i)[0]->ID();
          //	    fData->hit_clusterKey[i] = fmcl.at(i)[0].key();
        }
        else std::cout << "No spacepoint for this hit" << std::endl;
      }
    }
    //Get space points associated with the hit
    std::vector< art::Ptr<recob::SpacePoint> > sptv = fmsp.at(hits[ipl][i]);
*/



    //std::vector<double> xyz = bt_serv->HitToXYZ(hitlist[i]);
    //when the size of simIDEs is zero, the above function throws an exception
    //and crashes, so check that the simIDEs have non-zero size before
    //extracting hit true XYZ from simIDEs
    //      if (fIsMC){
    //        std::vector<sim::IDE> ides;
    //        try{bt_serv->HitToSimIDEs(hitlist[i], ides); }
    //        catch(...){}
    //        if (ides.size()>0){
    //          std::vector<double> xyz = bt_serv->SimIDEsToXYZ(ides);
    //          fData->hit_trueX[i] = xyz[0];
    //        }
    //      }
    /*
       for (unsigned int it=0; it<fTrackModuleLabel.size();++it){
       art::FindManyP<recob::Track> fmtk(hitListHandle,evt,fTrackModuleLabel[it]);
       if (fmtk.at(i).size()!=0){
       hit_trkid[it][i] = fmtk.at(i)[0]->ID();
       }
       else
       hit_trkid[it][i] = 0;
       }
       */

/*
    if (fSaveRawDigitInfo){
      //Hit to RawDigit information
      art::FindManyP<raw::RawDigit> fmrd(hitListHandle,evt,fHitsModuleLabel);
      if (hitlist[i]->WireID().Plane==2)
      {
        int dataSize = fmrd.at(i)[0]->Samples();
        short ped = fmrd.at(i)[0]->GetPedestal();

        std::vector<short> rawadc(dataSize);
        raw::Uncompress(fmrd.at(i)[0]->ADCs(), rawadc, fmrd.at(i)[0]->Compression());
        int t0 = hitlist[i]->PeakTime() - 3*(hitlist[i]->RMS());
        if (t0<0) t0 = 0;
        int t1 = hitlist[i]->PeakTime() + 3*(hitlist[i]->RMS());
        if (t1>=dataSize) t1 = dataSize-1;
        fData->rawD_ph[i] = -1;
        fData->rawD_peakT[i] = -1;
        for (int j = t0; j<=t1; ++j){
          if (rawadc[j]-ped>fData->rawD_ph[i]){
            fData->rawD_ph[i] = rawadc[j]-ped;
            fData->rawD_peakT[i] = j;
          }
        }
        fData->rawD_charge[i] = 0;
        fData->rawD_fwhh[i] = 0;
        double mean_t = 0.0;
        double mean_t2 = 0.0;
        for (int j = t0; j<=t1; ++j){
          if (rawadc[j]-ped>=0.5*fData->rawD_ph[i]){
            ++fData->rawD_fwhh[i];
          }
          if (rawadc[j]-ped>=0.1*fData->rawD_ph[i]){
            fData->rawD_charge[i] += rawadc[j]-ped;
            mean_t += (double)j*(rawadc[j]-ped);
            mean_t2 += (double)j*(double)j*(rawadc[j]-ped);
          }
        }
        mean_t/=fData->rawD_charge[i];
        mean_t2/=fData->rawD_charge[i];
        fData->rawD_rms[i] = sqrt(mean_t2-mean_t*mean_t);
      }
    } // Save RawDigitInfo
*/
    /*      if (!evt.isRealData()&&!isCosmics){
            fData -> hit_nelec[i] = 0;
            fData -> hit_energy[i] = 0;
            const sim::SimChannel* chan = 0;
            for(size_t sc = 0; sc < fSimChannels.size(); ++sc){
            if(fSimChannels[sc]->Channel() == hitlist[i]->Channel()) chan = fSimChannels[sc];
            }
            if (chan){
            for(auto const& mapitr : chan->TDCIDEMap()){
    // loop over the vector of IDE objects.
    for(auto const& ide : mapitr.second){
    fData -> hit_nelec[i] += ide.numElectrons;
    fData -> hit_energy[i] += ide.energy;
    }
    }
    }
    }*/
  } // Loop over hits



  if (evt.getByLabel(fHitsModuleLabel,hitListHandle)){
    //Find tracks associated with hits
    art::FindManyP<recob::Track> fmtk(hitListHandle,evt,fTrackModuleLabel[0]);
    for (size_t i = 0; i < NHits && i < kMaxHits ; ++i){//loop over hits
      if (fmtk.isValid()){
        if (fmtk.at(i).size()!=0){
          fData->hit_trkid[i] = fmtk.at(i)[0]->ID();
          //	    fData->hit_trkKey[i] = fmtk.at(i)[0].key();

        }
        else
          fData->hit_trkid[i] = -1;
      }
    }
  }

  //In the case of ClusterCrawler or linecluster, use "linecluster or clustercrawler" as HitModuleLabel.
  //using cchit will not make this association. In the case of gaushit, just use gaushit
  //Not initializing clusterID to -1 since some clustering algorithms assign negative IDs!
  if (evt.getByLabel(fHitsModuleLabel,hitListHandle)){
    //Find clusters associated with hits
    art::FindManyP<recob::Cluster> fmcl(hitListHandle,evt,fClusterModuleLabel);
    for (size_t i = 0; i < NHits && i < kMaxHits ; ++i){//loop over hits
      if (fmcl.isValid()){
        if (fmcl.at(i).size()!=0){
          fData->hit_clusterid[i] = fmcl.at(i)[0]->ID();
          //	    fData->hit_clusterKey[i] = fmcl.at(i)[0].key();
        }
        else
          fData->hit_clusterid[i] = -1;
      }
    }
  }
}// end (fSaveHitInfo)

if(fSavePandoraNuVertexInfo) {
  lar_pandora::PFParticleVector particleVector;
  lar_pandora::LArPandoraHelper::CollectPFParticles(evt, fPandoraNuVertexModuleLabel, particleVector);
  lar_pandora::VertexVector vertexVector;
  lar_pandora::PFParticlesToVertices particlesToVertices;
  lar_pandora::LArPandoraHelper::CollectVertices(evt, fPandoraNuVertexModuleLabel, vertexVector, particlesToVertices);

  short nprim = 0;
  for (unsigned int n = 0; n < particleVector.size(); ++n) {
    const art::Ptr<recob::PFParticle> particle = particleVector.at(n);
    if(particle->IsPrimary()) nprim++;
  }

  if (nprim > kMaxVertices){
    // got this error? consider increasing kMaxClusters
    // (or ask for a redesign using vectors)
    mf::LogError("AnaRootParser:limits") << "event has " << nprim
      << " nu neutrino vertices, only kMaxVertices=" << kMaxVertices << " stored in tree";
  }

  fData->nnuvtx = nprim;

  short iv = 0;
  for (unsigned int n = 0; n < particleVector.size(); ++n) {
    const art::Ptr<recob::PFParticle> particle = particleVector.at(n);
    if(particle->IsPrimary()) {
      lar_pandora::PFParticlesToVertices::const_iterator vIter = particlesToVertices.find(particle);
      if (particlesToVertices.end() != vIter) {
        const lar_pandora::VertexVector &vertexVector = vIter->second;
        if (!vertexVector.empty()) {
          if (vertexVector.size() == 1) {
            const art::Ptr<recob::Vertex> vertex = *(vertexVector.begin());
            double xyz[3] = {0.0, 0.0, 0.0} ;
            vertex->XYZ(xyz);
            fData->nuvtxx[iv] = xyz[0];
            fData->nuvtxy[iv] = xyz[1];
            fData->nuvtxz[iv] = xyz[2];
            fData->nuvtxpdg[iv] = particle->PdgCode();
            iv++;
          }
        }
      }
    }
  }
} // save PandoraNuVertexInfo

if (fSaveClusterInfo){
  fData->nclusters = (int) NClusters;
  if (NClusters > kMaxClusters){
    // got this error? consider increasing kMaxClusters
    // (or ask for a redesign using vectors)
    mf::LogError("AnaRootParser:limits") << "event has " << NClusters
      << " clusters, only kMaxClusters=" << kMaxClusters << " stored in tree";
  }
  for(unsigned int ic=0; ic<NClusters;++ic){//loop over clusters
    art::Ptr<recob::Cluster> clusterholder(clusterListHandle, ic);
    const recob::Cluster& cluster = *clusterholder;
    fData->clusterId[ic] = cluster.ID();
    fData->clusterView[ic] = cluster.View();
    // View 3 in LArSoft corresponds to View 0 in real life (3m long channels). View 2 in LArSoft corresponds to View 1 in real life (1m long channels).
    if(fData->clusterView[ic] == 3) {fData->clusterView[ic] = 0;}
    if(fData->clusterView[ic] == 2) {fData->clusterView[ic] = 1;}
    fData->cluster_isValid[ic] = cluster.isValid();
    fData->cluster_StartCharge[ic] = cluster.StartCharge();
    fData->cluster_StartAngle[ic] = cluster.StartAngle();
    fData->cluster_EndCharge[ic] = cluster.EndCharge();
    fData->cluster_EndAngle[ic] = cluster.EndAngle();
    fData->cluster_Integral[ic] = cluster.Integral();
    fData->cluster_IntegralAverage[ic] = cluster.IntegralAverage();
    fData->cluster_SummedADC[ic] = cluster.SummedADC();
    fData->cluster_SummedADCaverage[ic] = cluster.SummedADCaverage();
    fData->cluster_MultipleHitDensity[ic] = cluster.MultipleHitDensity();
    fData->cluster_Width[ic] = cluster.Width();
    fData->cluster_NHits[ic] = cluster.NHits();
    fData->cluster_StartWire[ic] = cluster.StartWire();
    fData->cluster_StartTick[ic] = cluster.StartTick();
    fData->cluster_EndWire[ic] = cluster.EndWire();
    fData->cluster_EndTick[ic] = cluster.EndTick();

    // Transform StartWire and EndWire to channels
    if(fData->clusterView[ic] == 1){fData->cluster_StartWire[ic]+=320; fData->cluster_EndWire[ic] += 320;}

    //Cosmic Tagger information for cluster
    /*      art::FindManyP<anab::CosmicTag> fmcct(clusterListHandle,evt,fCosmicClusterTaggerAssocLabel);
            if (fmcct.isValid()){
            fData->cluncosmictags_tagger[ic]     = fmcct.at(ic).size();
            if (fmcct.at(ic).size()>0){
            if(fmcct.at(ic).size()>1)
            std::cerr << "\n Warning : more than one cosmic tag per cluster in module! assigning the first tag to the cluster" << fCosmicClusterTaggerAssocLabel;
            fData->clucosmicscore_tagger[ic] = fmcct.at(ic).at(0)->CosmicScore();
            fData->clucosmictype_tagger[ic] = fmcct.at(ic).at(0)->CosmicType();
            }
            } */
  }//end loop over clusters
}//end fSaveClusterInfo

if (fSaveFlashInfo){
  fData->no_flashes = (int) NFlashes;
  if (NFlashes > kMaxFlashes) {
    // got this error? consider increasing kMaxHits
    // (or ask for a redesign using vectors)
    mf::LogError("AnaRootParser:limits") << "event has " << NFlashes
      << " flashes, only kMaxFlashes=" << kMaxFlashes << " stored in tree";
  }
  for (size_t i = 0; i < NFlashes && i < kMaxFlashes ; ++i){//loop over hits
    fData->flash_time[i]       = flashlist[i]->Time();
    fData->flash_pe[i]         = flashlist[i]->TotalPE();
    fData->flash_ycenter[i]    = flashlist[i]->YCenter();
    fData->flash_zcenter[i]    = flashlist[i]->ZCenter();
    fData->flash_ywidth[i]     = flashlist[i]->YWidth();
    fData->flash_zwidth[i]     = flashlist[i]->ZWidth();
    fData->flash_timewidth[i]  = flashlist[i]->TimeWidth();
  }
}

if (fSaveExternCounterInfo){
  fData->no_ExternCounts = (int) NExternCounts;
  if (NExternCounts > kMaxExternCounts) {
    // got this error? consider increasing kMaxHits
    // (or ask for a redesign using vectors)
    mf::LogError("AnaRootParser:limits") << "event has " << NExternCounts
      << " External Counters, only kMaxExternCounts=" << kMaxExternCounts << " stored in tree";
  }
  for (size_t i = 0; i < NExternCounts && i < kMaxExternCounts ; ++i){//loop over hits
    fData->externcounts_time[i] = countlist[i]->GetTrigTime();
    fData->externcounts_id[i]   = countlist[i]->GetTrigID();
  }
}


// Declare object-ID-to-PFParticleID maps so we can assign hasPFParticle and PFParticleID to the tracks, showers, vertices.
std::map<Short_t, Short_t> trackIDtoPFParticleIDMap, vertexIDtoPFParticleIDMap, showerIDtoPFParticleIDMap;

//Save PFParticle information
if (fSavePFParticleInfo){
  AnaRootParserDataStruct::PFParticleDataStruct& PFParticleData = fData->GetPFParticleData();
  size_t NPFParticles = pfparticlelist.size();

  PFParticleData.SetMaxPFParticles(std::max(NPFParticles, (size_t) 1));
  PFParticleData.Clear(); // clear all the data

  PFParticleData.nPFParticles = (short) NPFParticles;

  // now set the tree addresses to the newly allocated memory;
  // this creates the tree branches in case they are not there yet
  SetPFParticleAddress();

  if (NPFParticles > PFParticleData.GetMaxPFParticles()) {
    mf::LogError("AnaRootParser:limits") << "event has " << NPFParticles
      << " PFParticles, only "
      << PFParticleData.GetMaxPFParticles() << " stored in tree";
  }

  lar_pandora::PFParticleVector neutrinoPFParticles;
  lar_pandora::LArPandoraHelper::SelectNeutrinoPFParticles(pfparticlelist, neutrinoPFParticles);
  PFParticleData.pfp_numNeutrinos = neutrinoPFParticles.size();

  for (size_t i = 0; i < std::min(neutrinoPFParticles.size(), (size_t)kMaxNPFPNeutrinos); ++i) {
    PFParticleData.pfp_neutrinoIDs[i] = neutrinoPFParticles[i]->Self();
  }

  if (neutrinoPFParticles.size() > kMaxNPFPNeutrinos)
    std::cerr << "Warning: there were " << neutrinoPFParticles.size() << " reconstructed PFParticle neutrinos; only the first " << kMaxNPFPNeutrinos << " being stored in tree" << std::endl;

  // Get a PFParticle-to-vertex map.
  lar_pandora::VertexVector allPfParticleVertices;
  lar_pandora::PFParticlesToVertices pfParticleToVertexMap;
  lar_pandora::LArPandoraHelper::CollectVertices(evt, fPFParticleModuleLabel, allPfParticleVertices, pfParticleToVertexMap);

  // Get a PFParticle-to-track map.
  lar_pandora::TrackVector allPfParticleTracks;
  lar_pandora::PFParticlesToTracks pfParticleToTrackMap;
  lar_pandora::LArPandoraHelper::CollectTracks(evt, fPFParticleModuleLabel, allPfParticleTracks, pfParticleToTrackMap);

  // Get a PFParticle-to-shower map.
  lar_pandora::ShowerVector allPfParticleShowers;
  lar_pandora::PFParticlesToShowers pfParticleToShowerMap;
  lar_pandora::LArPandoraHelper::CollectShowers(evt, fPFParticleModuleLabel, allPfParticleShowers, pfParticleToShowerMap);

  for (size_t i = 0; i < NPFParticles && i < PFParticleData.GetMaxPFParticles() ; ++i){
    PFParticleData.pfp_selfID[i] = pfparticlelist[i]->Self();
    PFParticleData.pfp_isPrimary[i] = (Short_t)pfparticlelist[i]->IsPrimary();
    PFParticleData.pfp_numDaughters[i] = pfparticlelist[i]->NumDaughters();
    PFParticleData.pfp_parentID[i] = pfparticlelist[i]->Parent();
    PFParticleData.pfp_pdgCode[i] = pfparticlelist[i]->PdgCode();
    PFParticleData.pfp_isNeutrino[i] = lar_pandora::LArPandoraHelper::IsNeutrino(pfparticlelist[i]);

    // Set the daughter IDs.
    std::vector<size_t> daughterIDs = pfparticlelist[i]->Daughters();

    for (size_t j = 0; j < daughterIDs.size(); ++j)
      PFParticleData.pfp_daughterIDs[i][j] = daughterIDs[j];

    // Set the vertex ID.
    auto vertexMapIter = pfParticleToVertexMap.find(pfparticlelist[i]);
    if (vertexMapIter != pfParticleToVertexMap.end()) {
      lar_pandora::VertexVector pfParticleVertices = vertexMapIter->second;

      if (pfParticleVertices.size() > 1)
        std::cerr << "Warning: there was more than one vertex found for PFParticle with ID " << pfparticlelist[i]->Self() << ", storing only one" << std::endl;

      if (pfParticleVertices.size() > 0) {
        PFParticleData.pfp_vertexID[i] = pfParticleVertices.at(0)->ID();
        vertexIDtoPFParticleIDMap.insert(std::make_pair(pfParticleVertices.at(0)->ID(), pfparticlelist[i]->Self()));
      }
    }
    else
      std::cerr << "Warning: there was no vertex found for PFParticle with ID " << pfparticlelist[i]->Self() << std::endl;

    if (lar_pandora::LArPandoraHelper::IsTrack(pfparticlelist[i])){
      PFParticleData.pfp_isTrack[i] = 1;

      // Set the track ID.
      auto trackMapIter = pfParticleToTrackMap.find(pfparticlelist[i]);
      if (trackMapIter != pfParticleToTrackMap.end()) {
        lar_pandora::TrackVector pfParticleTracks = trackMapIter->second;

        if (pfParticleTracks.size() > 1)
          std::cerr << "Warning: there was more than one track found for PFParticle with ID " << pfparticlelist[i]->Self() << std::endl;

        if (pfParticleTracks.size() > 0) {
          PFParticleData.pfp_trackID[i] = pfParticleTracks.at(0)->ID();
          trackIDtoPFParticleIDMap.insert(std::make_pair(pfParticleTracks.at(0)->ID(), pfparticlelist[i]->Self()));
        }
      }
      else
        std::cerr << "Warning: there was no track found for track-like PFParticle with ID " << pfparticlelist[i]->Self() << std::endl;
    }
    else
      PFParticleData.pfp_isTrack[i] = 0;

    if (lar_pandora::LArPandoraHelper::IsShower(pfparticlelist[i])) {
      PFParticleData.pfp_isShower[i] = 1;
      // Set the shower ID.
      auto showerMapIter = pfParticleToShowerMap.find(pfparticlelist[i]);
      if (showerMapIter != pfParticleToShowerMap.end()) {
        lar_pandora::ShowerVector pfParticleShowers = showerMapIter->second;

        if (pfParticleShowers.size() > 1)
          std::cerr << "Warning: there was more than one shower found for PFParticle with ID " << pfparticlelist[i]->Self() << std::endl;

        if (pfParticleShowers.size() > 0) {
          PFParticleData.pfp_showerID[i] = pfParticleShowers.at(0)->ID();
          showerIDtoPFParticleIDMap.insert(std::make_pair(pfParticleShowers.at(0)->ID(), pfparticlelist[i]->Self()));
        }
      }
      else
        std::cerr << "Warning: there was no shower found for shower-like PFParticle with ID " << pfparticlelist[i]->Self() << std::endl;
    }
    else
      PFParticleData.pfp_isShower[i] = 0;

    // Set the cluster IDs.
    auto clusterMapIter = pfParticleToClusterMap.find(pfparticlelist[i]);
    if (clusterMapIter != pfParticleToClusterMap.end()) {
      lar_pandora::ClusterVector pfParticleClusters = clusterMapIter->second;
      PFParticleData.pfp_numClusters[i] = pfParticleClusters.size();

      for (size_t j = 0; j < pfParticleClusters.size(); ++j)
        PFParticleData.pfp_clusterIDs[i][j] = pfParticleClusters[j]->ID();
    }
    //else
    //  std::cerr << "Warning: there were no clusters found for PFParticle with ID " << pfparticlelist[i]->Self() << std::endl;
  }
} // if fSavePFParticleInfo

if (fSaveShowerInfo){

  // fill data from all the shower algorithms
  for (size_t iShowerAlgo = 0; iShowerAlgo < NShowerAlgos; ++iShowerAlgo) {
    AnaRootParserDataStruct::ShowerDataStruct& ShowerData
      = fData->GetShowerData(iShowerAlgo);
    std::vector<recob::Shower> const* pShowers = showerList[iShowerAlgo];
    art::Handle< std::vector<recob::Shower> > showerHandle = showerListHandle[iShowerAlgo];

    if (pShowers){
      FillShowers(ShowerData, *pShowers, fSavePFParticleInfo, showerIDtoPFParticleIDMap);

      if(fMVAPIDShowerModuleLabel[iShowerAlgo].size()){
        art::FindOneP<anab::MVAPIDResult> fmvapid(showerHandle, evt, fMVAPIDShowerModuleLabel[iShowerAlgo]);
        if(fmvapid.isValid()){
          for(unsigned int iShower=0;iShower<showerHandle->size();++iShower){
            const art::Ptr<anab::MVAPIDResult> pid = fmvapid.at(iShower);
            ShowerData.shwr_pidmvamu[iShower] = pid->mvaOutput.at("muon");
            ShowerData.shwr_pidmvae[iShower] = pid->mvaOutput.at("electron");
            ShowerData.shwr_pidmvapich[iShower] = pid->mvaOutput.at("pich");
            ShowerData.shwr_pidmvaphoton[iShower] = pid->mvaOutput.at("photon");
            ShowerData.shwr_pidmvapr[iShower] = pid->mvaOutput.at("proton");
          }
        } // fmvapid.isValid()
      }
    }
    else ShowerData.MarkMissing(fTree); // tree should reflect lack of data
  } // for iShowerAlgo

} // if fSaveShowerInfo

//track information for multiple trackers
if (fSaveTrackInfo) {

  // Get Geometry
  art::ServiceHandle<geo::Geometry> geomhandle;

  for (unsigned int iTracker=0; iTracker < NTrackers; ++iTracker){
    AnaRootParserDataStruct::TrackDataStruct& TrackerData = fData->GetTrackerData(iTracker);

    size_t NTracks = tracklist[iTracker].size();

    // allocate enough space for this number of tracks (but at least for one of them!)
    TrackerData.SetMaxTracks(std::max(NTracks, (size_t) 1));
    TrackerData.Clear(); // clear all the data

    TrackerData.ntracks = (int) NTracks;
    /*
    //First loop over tracks to get array sizes (number of hits in each track)
    for(size_t iTrk=0; iTrk < NTracks; ++iTrk){//loop over tracks
    art::FindMany<anab::Calorimetry> fmcaltemp(trackListHandle[iTracker], evt, fCalorimetryModuleLabel[iTracker]);
    if (fmcaltemp.isValid()){
    std::vector<const anab::Calorimetry*> calostemp = fmcaltemp.at(iTrk);
    for (size_t ical = 0; ical<calostemp.size(); ++ical){
    if (!calostemp[ical]) continue;
    if (!calostemp[ical]->PlaneID().isValid) continue;
    int planenumtemp = calostemp[ical]->PlaneID().Plane;
    if (planenumtemp<0||planenumtemp>2) continue;
    //For now make the second argument as 13 for muons.
    //	    const size_t NHitsTemp = calostemp[ical] -> dEdx().size();
    //	    TrackerData.NHitsInAllTracks+=(int) NHitsTemp;
    //	    fData->NHitsInAllTracks+=(int) NHitsTemp;
    } // for calorimetry info
    } // if has calorimetry info
    }//loop over tracks
    */

    // now set the tree addresses to the newly allocated memory;
    // this creates the tree branches in case they are not there yet
    SetTrackerAddresses(iTracker);
    if (NTracks > TrackerData.GetMaxTracks()) {
      // got this error? it might be a bug,
      // since we are supposed to have allocated enough space to fit all tracks
      mf::LogError("AnaRootParser:limits") << "event has " << NTracks
        << " " << fTrackModuleLabel[iTracker] << " tracks, only "
        << TrackerData.GetMaxTracks() << " stored in tree";
    }

    //call the track momentum algorithm that gives you momentum based on track range
    trkf::TrackMomentumCalculator trkm;
    trkm.SetMinLength(50); //change the minimal track length requirement to 50 cm

    recob::MCSFitResult res;

    int HitIterator=0;
    int HitIterator2=0;

    int trueID = -9999;
    float maxe = -9999.;
    float purity = -9999.;

    for(size_t iTrk=0; iTrk < NTracks; ++iTrk){//loop over tracks

      //save t0 from reconstructed flash track matching for every track
      art::FindManyP<anab::T0> fmt0(trackListHandle[iTracker],evt,fFlashT0FinderLabel[iTracker]);
      if (fmt0.isValid()){
        if(fmt0.at(iTrk).size()>0){
          if(fmt0.at(iTrk).size()>1)
            std::cerr << "\n Warning : more than one cosmic tag per track in module! assigning the first tag to the track" << fFlashT0FinderLabel[iTracker];
          TrackerData.trkflashT0[iTrk] = fmt0.at(iTrk).at(0)->Time();
        }
      }

      //save t0 from reconstructed flash track matching for every track
      art::FindManyP<anab::T0> fmmct0(trackListHandle[iTracker],evt,fMCT0FinderLabel[iTracker]);
      if (fmmct0.isValid()){
        if(fmmct0.at(iTrk).size()>0){
          if(fmmct0.at(iTrk).size()>1)
            std::cerr << "\n Warning : more than one cosmic tag per track in module! assigning the first tag to the cluster" << fMCT0FinderLabel[iTracker];
          TrackerData.trktrueT0[iTrk] = fmmct0.at(iTrk).at(0)->Time();
        }
      }

      //Cosmic Tagger information
      art::FindManyP<anab::CosmicTag> fmct(trackListHandle[iTracker],evt,fCosmicTaggerAssocLabel[iTracker]);
      if (fmct.isValid()){
        TrackerData.trkncosmictags_tagger[iTrk]     = fmct.at(iTrk).size();
        if (fmct.at(iTrk).size()>0){
          if(fmct.at(iTrk).size()>1)
            std::cerr << "\n Warning : more than one cosmic tag per track in module! assigning the first tag to the track" << fCosmicTaggerAssocLabel[iTracker];
          TrackerData.trkcosmicscore_tagger[iTrk] = fmct.at(iTrk).at(0)->CosmicScore();
          TrackerData.trkcosmictype_tagger[iTrk] = fmct.at(iTrk).at(0)->CosmicType();
        }
      }

      //Containment Tagger information
      art::FindManyP<anab::CosmicTag> fmcnt(trackListHandle[iTracker],evt,fContainmentTaggerAssocLabel[iTracker]);
      if (fmcnt.isValid()){
        TrackerData.trkncosmictags_containmenttagger[iTrk]     = fmcnt.at(iTrk).size();
        if (fmcnt.at(iTrk).size()>0){
          if(fmcnt.at(iTrk).size()>1)
            std::cerr << "\n Warning : more than one containment tag per track in module! assigning the first tag to the track" << fContainmentTaggerAssocLabel[iTracker];
          TrackerData.trkcosmicscore_containmenttagger[iTrk] = fmcnt.at(iTrk).at(0)->CosmicScore();
          TrackerData.trkcosmictype_containmenttagger[iTrk] = fmcnt.at(iTrk).at(0)->CosmicType();
        }
      }

      //Flash match compatibility information
      //Unlike CosmicTagger, Flash match doesn't assign a cosmic tag for every track. For those tracks, AnaRootParser initializes them with -999 or -9999
      art::FindManyP<anab::CosmicTag> fmbfm(trackListHandle[iTracker],evt,fFlashMatchAssocLabel[iTracker]);
      if (fmbfm.isValid()){
        TrackerData.trkncosmictags_flashmatch[iTrk] = fmbfm.at(iTrk).size();
        if (fmbfm.at(iTrk).size()>0){
          if(fmbfm.at(iTrk).size()>1)
            std::cerr << "\n Warning : more than one cosmic tag per track in module! assigning the first tag to the track" << fFlashMatchAssocLabel[iTracker];
          TrackerData.trkcosmicscore_flashmatch[iTrk] = fmbfm.at(iTrk).at(0)->CosmicScore();
          TrackerData.trkcosmictype_flashmatch[iTrk] = fmbfm.at(iTrk).at(0)->CosmicType();
          //std::cout<<"\n"<<evt.event()<<"\t"<<iTrk<<"\t"<<fmbfm.at(iTrk).at(0)->CosmicScore()<<"\t"<<fmbfm.at(iTrk).at(0)->CosmicType();
        }
      }

      art::Ptr<recob::Track> ptrack(trackListHandle[iTracker], iTrk);
      const recob::Track& track = *ptrack;

      TVector3 pos, dir_start, dir_end, end;
      TVector3 dir_start_flipped, dir_end_flipped;
      double tlen = 0., mom = 0.;
      int TrackID = -1;

      int ntraj = 0;
      //we need to use Bezier methods for Bezier tracks
      if (fTrackModuleLabel[iTracker].find("beziertracker")!=std::string::npos) {
        trkf::BezierTrack btrack(*ptrack);
        ntraj = btrack.NSegments();
        if(ntraj > 0) {
          double xyz[3];
          btrack.GetTrackPoint(0,xyz);
          pos.SetXYZ(xyz[0],xyz[1],xyz[2]);
          btrack.GetTrackDirection(0,xyz);
          dir_start.SetXYZ(xyz[0],xyz[1],xyz[2]);
          btrack.GetTrackDirection(1,xyz);
          dir_end.SetXYZ(xyz[0],xyz[1],xyz[2]);
          btrack.GetTrackPoint(1,xyz);
          end.SetXYZ(xyz[0],xyz[1],xyz[2]);

          tlen = btrack.GetLength();
          if (btrack.GetTrajectory().NPoints() > 0)
            mom = btrack.GetTrajectory().StartMomentum();
          // fill bezier track reco branches
          TrackID = iTrk;  //bezier has some screwed up track IDs
        }
      }
      else {   //use the normal methods for other kinds of tracks
        ntraj = track.NumberTrajectoryPoints();
        if (ntraj > 0) {
          pos       = track.Vertex();
          dir_start = track.VertexDirection();
          dir_end   = track.EndDirection();
          end       = track.End();

          dir_start_flipped.SetXYZ(dir_start.Z(), dir_start.Y(), dir_start.X());
          dir_end_flipped.SetXYZ(dir_end.Z(), dir_end.Y(), dir_end.X());

          tlen        = length(track);
          if(track.NumberTrajectoryPoints() > 0)
            mom = track.VertexMomentum();
          // fill non-bezier-track reco branches
          TrackID = track.ID();
        }
      }

      if (ntraj > 0) {
        double theta_xz = std::atan2(dir_start.X(), dir_start.Z());
        double theta_yz = std::atan2(dir_start.Y(), dir_start.Z());
        double dpos = bdist(pos);  // FIXME - Passing an uncorrected position....
        double dend = bdist(end);  // FIXME - Passing an uncorrected position....

        TrackerData.trkId[iTrk]                 = TrackID;
        TrackerData.trkstartx[iTrk]             = pos.X();
        TrackerData.trkstarty[iTrk]             = pos.Y();
        TrackerData.trkstartz[iTrk]             = pos.Z();
        TrackerData.trkstartd[iTrk]		= dpos;
        TrackerData.trkendx[iTrk]		= end.X();
        TrackerData.trkendy[iTrk]		= end.Y();
        TrackerData.trkendz[iTrk]		= end.Z();
        TrackerData.trkendd[iTrk]		= dend;
        TrackerData.trkstarttheta[iTrk]		= (180.0/3.14159)*dir_start_flipped.Theta();
        TrackerData.trkstartphi[iTrk]		= (180.0/3.14159)*dir_start_flipped.Phi();
        TrackerData.trkstartdirectionx[iTrk]	= dir_start.X();
        TrackerData.trkstartdirectiony[iTrk]	= dir_start.Y();
        TrackerData.trkstartdirectionz[iTrk]	= dir_start.Z();

	if(fLogLevel == 2)
	{
	  std::cout << std::endl;
	  std::cout << "start.X(): " << pos.X() << "\t" << "start.Y(): " << pos.Y() << "\t" << "start.Z(): " << pos.Z() << std::endl;
	  std::cout << "end.X(): " << end.X() << "\t" << "end.Y(): " << end.Y() << "\t" << "end.Z(): " << end.Z() << std::endl;
	  std::cout << "dir_start.X(): " << dir_start.X() << "\t" << "dir_start.Y(): " << dir_start.Y() << "\t" << "dir_start.Z(): " << dir_start.Z() << std::endl;
	  std::cout << "dir_end.X(): " << dir_end.X() << "\t" << "dir_end.Y(): " << dir_end.Y() << "\t" << "dir_end.Z(): " << dir_end.Z() << std::endl;
	  std::cout << "dir_start_flipped.Theta(): " << (180.0/3.14159)*dir_start_flipped.Theta() << "\t" << "dir_start_flipped.Phi(): " << (180.0/3.14159)*dir_start_flipped.Phi() << std::endl;
	  std::cout << "dir_end_flipped.Theta(): " << (180.0/3.14159)*dir_end_flipped.Theta() << "\t" << "dir_end_flipped.Phi(): " << (180.0/3.14159)*dir_end_flipped.Phi() << std::endl;
	  std::cout << std::endl;
	}

        TrackerData.trkendtheta[iTrk]	  = (180.0/3.14159)*dir_end_flipped.Theta();
        TrackerData.trkendphi[iTrk]	  = (180.0/3.14159)*dir_end_flipped.Phi();
        TrackerData.trkenddirectionx[iTrk] 	  = dir_end.X();
        TrackerData.trkenddirectiony[iTrk] 	  = dir_end.Y();
        TrackerData.trkenddirectionz[iTrk] 	  = dir_end.Z();

        TrackerData.trkthetaxz[iTrk]  	     = theta_xz;
        TrackerData.trkthetayz[iTrk]  	     = theta_yz;
        TrackerData.trkmom[iTrk]	     = mom;
	TrackerData.trkchi2PerNDF[iTrk]	     = track.Chi2PerNdof();
	TrackerData.trkNDF[iTrk]	     = track.Ndof();
        TrackerData.trklen[iTrk]	     = tlen;
        TrackerData.trklenstraightline[iTrk] = sqrt(pow(pos.X()-end.X(),2) + pow(pos.Y()-end.Y(),2) + pow(pos.Z()-end.Z(),2));
        TrackerData.trkmomrange[iTrk] 	     = trkm.GetTrackMomentum(tlen,13);
        TrackerData.trkmommschi2[iTrk]	     = trkm.GetMomentumMultiScatterChi2(ptrack);
        TrackerData.trkmommsllhd[iTrk]	     = trkm.GetMomentumMultiScatterLLHD(ptrack);

        //uBoone MCS
        res = fMCSFitter.fitMcs(*ptrack);
        TrackerData.trkmommscmic[iTrk]	  = res.bestMomentum();
        TrackerData.trkmommscfwd[iTrk]	  = res.fwdMomentum();
        TrackerData.trkmommscbwd[iTrk]	  = res.bwdMomentum();
        TrackerData.trkmommscllfwd[iTrk]  = res.fwdLogLikelihood();
        TrackerData.trkmommscllbwd[iTrk]  = res.bwdLogLikelihood();


        if (fSavePFParticleInfo) {
          auto mapIter = trackIDtoPFParticleIDMap.find(TrackID);
          if (mapIter != trackIDtoPFParticleIDMap.end()) {
            // This track has a corresponding PFParticle.
            TrackerData.trkhasPFParticle[iTrk] = 1;
            TrackerData.trkPFParticleID[iTrk] = mapIter->second;
          }
          else
            TrackerData.trkhasPFParticle[iTrk] = 0;
        }

      } // if we have trajectory

      // find vertices associated with this track
      /*
         art::FindMany<recob::Vertex> fmvtx(trackListHandle[iTracker], evt, fVertexModuleLabel[iTracker]);
         if(fmvtx.isValid()) {
         std::vector<const recob::Vertex*> verts = fmvtx.at(iTrk);
      // should have two at most
      for(size_t ivx = 0; ivx < verts.size(); ++ivx) {
      verts[ivx]->XYZ(xyz);
      // find the vertex in TrackerData to get the index
      short theVtx = -1;
      for(short jvx = 0; jvx < TrackerData.nvtx; ++jvx) {
      if(TrackerData.vtx[jvx][2] == xyz[2]) {
      theVtx = jvx;
      break;
      }
      } // jvx
      // decide if it should be assigned to the track Start or End.
      // A simple dz test should suffice
      if(fabs(xyz[2] - TrackerData.trkstartz[iTrk]) <
      fabs(xyz[2] - TrackerData.trkendz[iTrk])) {
      TrackerData.trksvtxid[iTrk] = theVtx;
      } else {
      TrackerData.trkevtxid[iTrk] = theVtx;
      }
      } // vertices
      } // fmvtx.isValid()
      */


      /* //commented out because now have several Vertices
         Float_t minsdist = 10000;
         Float_t minedist = 10000;
         for (int ivx = 0; ivx < NVertices && ivx < kMaxVertices; ++ivx){
         Float_t sdist = sqrt(pow(TrackerData.trkstartx[iTrk]-VertexData.vtxx[ivx],2)+
         pow(TrackerData.trkstarty[iTrk]-VertexData.vtxy[ivx],2)+
         pow(TrackerData.trkstartz[iTrk]-VertexData.vtxz[ivx],2));
         Float_t edist = sqrt(pow(TrackerData.trkendx[iTrk]-VertexData.vtxx[ivx],2)+
         pow(TrackerData.trkendy[iTrk]-VertexData.vtxy[ivx],2)+
         pow(TrackerData.trkendz[iTrk]-VertexData.vtxz[ivx],2));
         if (sdist<minsdist){
         minsdist = sdist;
         if (minsdist<10) TrackerData.trksvtxid[iTrk] = ivx;
         }
         if (edist<minedist){
         minedist = edist;
         if (minedist<10) TrackerData.trkevtxid[iTrk] = ivx;
         }
         }*/

      // find particle ID info
      art::FindMany<anab::ParticleID> fmpid(trackListHandle[iTracker], evt, fParticleIDModuleLabel[iTracker]);
      if(fmpid.isValid()) {
        std::vector<const anab::ParticleID*> pids = fmpid.at(iTrk);
        //if(pids.size() > 1) {
        //mf::LogError("AnaRootParser:limits")
        //<< "the " << fTrackModuleLabel[iTracker] << " track #" << iTrk
        //<< " has " << pids.size()
        //<< " set of ParticleID variables. Only one stored in the tree";
        //}
        for (size_t ipid = 0; ipid < pids.size(); ++ipid){
          if (!pids[ipid]->PlaneID().isValid) continue;
          int planenum = pids[ipid]->PlaneID().Plane;
          if (planenum<0||planenum>2) continue;
          TrackerData.trkpidpdg[iTrk][planenum] = pids[ipid]->Pdg();
          TrackerData.trkpidchi[iTrk][planenum] = pids[ipid]->MinChi2();
          TrackerData.trkpidchipr[iTrk][planenum] = pids[ipid]->Chi2Proton();
          TrackerData.trkpidchika[iTrk][planenum] = pids[ipid]->Chi2Kaon();
          TrackerData.trkpidchipi[iTrk][planenum] = pids[ipid]->Chi2Pion();
          TrackerData.trkpidchimu[iTrk][planenum] = pids[ipid]->Chi2Muon();
          TrackerData.trkpidpida[iTrk][planenum] = pids[ipid]->PIDA();
        }
      } // fmpid.isValid()
      if(fMVAPIDTrackModuleLabel[iTracker].size()){
        art::FindOneP<anab::MVAPIDResult> fmvapid(trackListHandle[iTracker], evt, fMVAPIDTrackModuleLabel[iTracker]);
        if(fmvapid.isValid()) {
          const art::Ptr<anab::MVAPIDResult> pid = fmvapid.at(iTrk);
          TrackerData.trkpidmvamu[iTrk] = pid->mvaOutput.at("muon");
          TrackerData.trkpidmvae[iTrk] = pid->mvaOutput.at("electron");
          TrackerData.trkpidmvapich[iTrk] = pid->mvaOutput.at("pich");
          TrackerData.trkpidmvaphoton[iTrk] = pid->mvaOutput.at("photon");
          TrackerData.trkpidmvapr[iTrk] = pid->mvaOutput.at("proton");
        } // fmvapid.isValid()
      }


        art::FindManyP<recob::Hit, recob::TrackHitMeta> fmthm(trackListHandle[iTracker], evt, "pmtrack");

//      if (fmthm.isValid()){
	auto vhit = fmthm.at(iTrk);
	auto vmeta = fmthm.data(iTrk);

        TrackerData.NHitsPerTrack[iTrk] = vhit.size();
      	art::FindManyP<recob::SpacePoint> fmspts(vhit, evt, "pmtrack");

	int NHitsView0 = 0;
	int NHitsView1 = 0;

	if(fLogLevel == 2)
	{
	  std::cout << "track.NumberTrajectoryPoints(): " << track.NumberTrajectoryPoints() << std::endl;
	  std::cout << "track.NPoints(): " << track.NPoints() << std::endl;
	  std::cout << "vhit.size(): " << vhit.size() << std::endl;
	  std::cout << "vmeta.size(): " << vmeta.size() << std::endl;
	  std::cout << "fmspts.size(): " << fmspts.size() << std::endl;
	}

  	  for (unsigned int h = 0; h < vhit.size(); h++)
  	  {
	    //corrected pitch
  	    double angleToVert = geomhandle->WireAngleToVertical(vhit[h]->View(), vhit[h]->WireID().TPC, vhit[h]->WireID().Cryostat) - 0.5*::util::pi<>();
  	    const TVector3& dir = tracklist[iTracker][iTrk]->DirectionAtPoint(h);
  	    const TVector3& loc = tracklist[iTracker][iTrk]->LocationAtPoint(h);
  	    double cosgamma = std::abs(std::sin(angleToVert)*dir.Y() + std::cos(angleToVert)*dir.Z());


	    TrackerData.hittrklocaltrackdirectionx[HitIterator2] = dir.X();
	    TrackerData.hittrklocaltrackdirectiony[HitIterator2] = dir.Y();
	    TrackerData.hittrklocaltrackdirectionz[HitIterator2] = dir.Z();


	    //XYZ
	    std::vector< art::Ptr<recob::SpacePoint> > sptv = fmspts.at(h);
	    TrackerData.hittrkx[HitIterator2] = sptv[0]->XYZ()[0];
	    TrackerData.hittrky[HitIterator2] = sptv[0]->XYZ()[1];
	    TrackerData.hittrkz[HitIterator2] = sptv[0]->XYZ()[2];

      	    TVector3 dir_hit_flipped;
            dir_hit_flipped.SetXYZ(dir.Z(), dir.Y(), dir.X());

	    TrackerData.hittrklocaltrackdirectiontheta[HitIterator2] = (180.0/3.14159)*dir_hit_flipped.Theta();
	    TrackerData.hittrklocaltrackdirectionphi[HitIterator2] = (180.0/3.14159)*dir_hit_flipped.Phi();

	    //dx
    	    if(vhit[h]->WireID().Plane == 0) TrackerData.hittrkpitchC[HitIterator2] = std::abs(geomhandle->WirePitch()/( sin(dir_hit_flipped.Theta())*sin(dir_hit_flipped.Phi()) ));
    	    if(vhit[h]->WireID().Plane == 1) TrackerData.hittrkpitchC[HitIterator2] = std::abs(geomhandle->WirePitch()/( sin(dir_hit_flipped.Theta())*cos(dir_hit_flipped.Phi()) ));

	    TrackerData.hittrkds[HitIterator2] = vmeta[h]->Dx();

	    if(fLogLevel == 2)
	    {
	      std::cout << "pos.X(): " << sptv[0]->XYZ()[0] << "\t" << "pos.Y(): " << sptv[0]->XYZ()[1] << "\t" << "pos.Z(): " << sptv[0]->XYZ()[2] << std::endl;
	      std::cout << "pos2.X(): " << loc.X() << "\t" << "pos2.Y(): " << loc.Y() << "\t" << "pos2.Z(): " << loc.Z() << std::endl;
	      std::cout << "dir.X(): " << dir.X() << "\t" << "dir.Y(): " << dir.Y() << "\t" << "dir.Z(): " << dir.Z() << std::endl;
	      std::cout << "dir_hit_flipped.Theta(): " << (180.0/3.14159)*dir_hit_flipped.Theta() << "\t" << "dir_hit_flipped.Phi(): " << (180.0/3.14159)*dir_hit_flipped.Phi() << std::endl;
	      std::cout << "vmeta[h]->Dx(): " << vmeta[h]->Dx() << std::endl;
	      std::cout << "Dx corrected pitch old: " << geomhandle->WirePitch()/cosgamma << std::endl;
	      std::cout << "Dx corrected pitch new: " << TrackerData.hittrkpitchC[HitIterator2] << std::endl;
	      std::cout << "view: " << vhit[h]->WireID().Plane << std::endl;
	    }

	    //hit variables
	    TrackerData.hittrkchannel[HitIterator2] = vhit[h]->Channel();
	    TrackerData.hittrktpc[HitIterator2] = vhit[h]->WireID().TPC;
	    TrackerData.hittrkview[HitIterator2] = vhit[h]->WireID().Plane;
	    TrackerData.hittrkwire[HitIterator2] = vhit[h]->WireID().Wire;
	    TrackerData.hittrkpeakT[HitIterator2] = vhit[h]->PeakTime();
	    TrackerData.hittrkchargeintegral[HitIterator2] = vhit[h]->Integral();
	    TrackerData.hittrkph[HitIterator2] = vhit[h]->PeakAmplitude();
	    TrackerData.hittrkchargesum[HitIterator2] = vhit[h]->SummedADC();
	    TrackerData.hittrkstarT[HitIterator2] = vhit[h]->StartTick();
	    TrackerData.hittrkendT[HitIterator2] = vhit[h]->EndTick();
	    TrackerData.hittrkrms[HitIterator2] = vhit[h]->RMS();
	    TrackerData.hittrkgoddnessofFit[HitIterator2] = vhit[h]->GoodnessOfFit();
	    TrackerData.hittrkmultiplicity[HitIterator2] = vhit[h]->Multiplicity();

      //quantities from backtracker are different from the real value in MC
      if( fIsMC )
        HitsBackTrack(  vhit[h], trueID, maxe, purity );

      TrackerData.hittrktrueID[HitIterator2] = trueID;
      TrackerData.hittrktrueEnergyMax[HitIterator2] = maxe;
      TrackerData.hittrktrueEnergyFraction[HitIterator2] = purity;

            HitIterator2++;

	    if(vhit[h]->WireID().Plane == 0) NHitsView0++;
	    if(vhit[h]->WireID().Plane == 1) NHitsView1++;
	  }
        TrackerData.ntrkhitsperview[iTrk][0] = NHitsView0;
        TrackerData.ntrkhitsperview[iTrk][1] = NHitsView1;

//      }

/*
	std::cout << "tracklist[iTracker][iTrk]->NumberTrajectoryPoints(): " << tracklist[iTracker][iTrk]->NumberTrajectoryPoints() << std::endl;
      for(size_t itp = 0; itp < tracklist[iTracker][iTrk]->NumberTrajectoryPoints(); ++itp)
      {
        const TVector3& pos = tracklist[iTracker][iTrk]->LocationAtPoint(itp);
      }
*/
//
//      art::FindManyP<recob::Hit> fmht(trackListHandle[iTracker], evt, fTrackModuleLabel);
//      std::vector< art::Ptr<recob::Hit> > allHits = fmht.at(trkIter);
//      art::FindManyP<recob::SpacePoint> fmspts(allHits, evt, fSpacePointModuleLabel);
//
/*
      art::FindManyP<recob::SpacePoint> fmspts(vhit, evt, "pmtrack");
	for (size_t h = 0; h < vhit.size(); ++h)
	{
	  std::vector< art::Ptr<recob::SpacePoint> > sptv = fmspts.at(h);

	  for (size_t j = 0; j < sptv.size(); ++j)
	  {
	    std::cout << "sptv[j]->XYZ()[0]: " << sptv[j]->XYZ()[0] << std::endl;
	    std::cout << "sptv[j]->XYZ()[1]: " << sptv[j]->XYZ()[1] << std::endl;
	    std::cout << "sptv[j]->XYZ()[2]: " << sptv[j]->XYZ()[2] << std::endl;
	    std::cout << "sptv[j]->ErrXYZ()[0]: " << sptv[j]->ErrXYZ()[0] << std::endl;
	    std::cout << "sptv[j]->ErrXYZ()[1]: " << sptv[j]->ErrXYZ()[1] << std::endl;
	    std::cout << "sptv[j]->ErrXYZ()[2]: " << sptv[j]->ErrXYZ()[2] << std::endl;
	  }
	}
*/
//


/////////////////////
      art::FindMany<anab::Calorimetry> fmcal(trackListHandle[iTracker], evt, fCalorimetryModuleLabel[iTracker]);
      if (fmcal.isValid()){
        std::vector<const anab::Calorimetry*> calos = fmcal.at(iTrk);
        if (calos.size() > TrackerData.GetMaxPlanesPerTrack(iTrk)) {
          // if you get this message, there is probably a bug somewhere since
          // the calorimetry planes should be 3.
          mf::LogError("AnaRootParser:limits")
            << "the " << fTrackModuleLabel[iTracker] << " track #" << iTrk
            << " has " << calos.size() << " planes for calorimetry , only "
            << TrackerData.GetMaxPlanesPerTrack(iTrk) << " stored in tree";
        }

        //	  for (size_t ical = 0; ical<calos.size(); ++ical){
        for (size_t ical = calos.size() - 1; ical <= 1 ; --ical){  //reverse order so that we access info on plane 0 (which is view 0 in real life) first
          if (!calos[ical]) continue;
          if (!calos[ical]->PlaneID().isValid) continue;
          int planenum = calos[ical]->PlaneID().Plane;
          if (planenum<0||planenum>1) continue; //only two planes (0 and 1) in dual phase
          TrackerData.trkke[iTrk][planenum]    = calos[ical]->KineticEnergy();
          TrackerData.trkrange[iTrk][planenum] = calos[ical]->Range();
          //For now make the second argument as 13 for muons.
          TrackerData.trkpitchc[iTrk][planenum]= calos[ical] -> TrkPitchC();
          const size_t NHits = calos[ical] -> dEdx().size();
//          TrackerData.ntrkhitsperview[iTrk][planenum] = (int) NHits;
          if (NHits > TrackerData.GetMaxHitsPerTrack(iTrk, planenum)) {
            // if you get this error, you'll have to increase kMaxTrackHits
            mf::LogError("AnaRootParser:limits")
              << "the " << fTrackModuleLabel[iTracker] << " track #" << iTrk
              << " has " << NHits << " hits on calorimetry plane #" << planenum
              <<", only "
              << TrackerData.GetMaxHitsPerTrack(iTrk, planenum) << " stored in tree";
          }
          if (!isCosmics){
            for(size_t iTrkHit = 0; iTrkHit < NHits && iTrkHit < TrackerData.GetMaxHitsPerTrack(iTrk, planenum); ++iTrkHit) {
              TrackerData.trkdedx[iTrk][planenum][iTrkHit]  = (calos[ical] -> dEdx())[iTrkHit];
              TrackerData.trkdqdx[iTrk][planenum][iTrkHit]  = (calos[ical] -> dQdx())[iTrkHit];
              TrackerData.trkresrg[iTrk][planenum][iTrkHit] = (calos[ical] -> ResidualRange())[iTrkHit];
              TrackerData.trktpc[iTrk][planenum][iTrkHit]   = (calos[ical] -> PlaneID()).TPC;
              const auto& TrkPos = (calos[ical] -> XYZ())[iTrkHit];
              auto& TrkXYZ = TrackerData.trkxyz[iTrk][planenum][iTrkHit];
              TrkXYZ[0] = TrkPos.X();
              TrkXYZ[1] = TrkPos.Y();
              TrkXYZ[2] = TrkPos.Z();
              //TrackerData.trkdedx2.push_back((calos[ical] -> dEdx())[iTrkHit]);
              TrackerData.trkdedx2[HitIterator] = (calos[ical] -> dEdx())[iTrkHit];
              TrackerData.trkdqdx2[HitIterator] = (calos[ical] -> dQdx())[iTrkHit];
              TrackerData.trktpc2[HitIterator] = (calos[ical] -> PlaneID()).TPC;
              TrackerData.trkx2[HitIterator] = TrkPos.X();
              TrackerData.trky2[HitIterator] = TrkPos.Y();
              TrackerData.trkz2[HitIterator] = TrkPos.Z();
              HitIterator++;
            } // for track hits
          }
        } // for calorimetry info


	//best plane
        if(TrackerData.ntrkhitsperview[iTrk][0] > TrackerData.ntrkhitsperview[iTrk][1] && TrackerData.ntrkhitsperview[iTrk][0] > TrackerData.ntrkhitsperview[iTrk][2]) TrackerData.trkpidbestplane[iTrk] = 0;
        else if(TrackerData.ntrkhitsperview[iTrk][1] > TrackerData.ntrkhitsperview[iTrk][0] && TrackerData.ntrkhitsperview[iTrk][1] > TrackerData.ntrkhitsperview[iTrk][2]) TrackerData.trkpidbestplane[iTrk] = 1;
        else if(TrackerData.ntrkhitsperview[iTrk][2] > TrackerData.ntrkhitsperview[iTrk][0] && TrackerData.ntrkhitsperview[iTrk][2] > TrackerData.ntrkhitsperview[iTrk][1]) TrackerData.trkpidbestplane[iTrk] = 2;
        else if(TrackerData.ntrkhitsperview[iTrk][2] == TrackerData.ntrkhitsperview[iTrk][0] && TrackerData.ntrkhitsperview[iTrk][2] > TrackerData.ntrkhitsperview[iTrk][1]) TrackerData.trkpidbestplane[iTrk] = 2;
        else if(TrackerData.ntrkhitsperview[iTrk][2] == TrackerData.ntrkhitsperview[iTrk][1] && TrackerData.ntrkhitsperview[iTrk][2] > TrackerData.ntrkhitsperview[iTrk][0]) TrackerData.trkpidbestplane[iTrk] = 2;
        else if(TrackerData.ntrkhitsperview[iTrk][1] == TrackerData.ntrkhitsperview[iTrk][0] && TrackerData.ntrkhitsperview[iTrk][1] > TrackerData.ntrkhitsperview[iTrk][2]) TrackerData.trkpidbestplane[iTrk] = 0;
        else if(TrackerData.ntrkhitsperview[iTrk][1] == TrackerData.ntrkhitsperview[iTrk][0] && TrackerData.ntrkhitsperview[iTrk][1] == TrackerData.ntrkhitsperview[iTrk][2]) TrackerData.trkpidbestplane[iTrk] = 2;

        // FIXME - Do i want to add someway to work out the best TPC???....
      } // if has calorimetry info
/*
      //track truth information
      if (fIsMC){
        //get the hits on each plane
        art::FindManyP<recob::Hit>      fmht(trackListHandle[iTracker], evt, fTrackModuleLabel[iTracker]);
        std::vector< art::Ptr<recob::Hit> > allHits = fmht.at(iTrk);
        std::vector< art::Ptr<recob::Hit> > hits[kNplanes];

        for(size_t ah = 0; ah < allHits.size(); ++ah){
//          if ( allHits[ah]->WireID().Plane >= 0 &&  // always true
             if( allHits[ah]->WireID().Plane <  3){
            hits[allHits[ah]->WireID().Plane].push_back(allHits[ah]);
          }
        }
        for (size_t ipl = 0; ipl < 3; ++ipl){
          double maxe = 0;
          HitsBackTrack(hits[ipl],TrackerData.trkidtruth[iTrk][ipl],TrackerData.trkpurtruth[iTrk][ipl],maxe);
          //std::cout<<"\n"<<iTracker<<"\t"<<iTrk<<"\t"<<ipl<<"\t"<<trkidtruth[iTracker][iTrk][ipl]<<"\t"<<trkpurtruth[iTracker][iTrk][ipl]<<"\t"<<maxe;
          if (TrackerData.trkidtruth[iTrk][ipl]>0){
            const art::Ptr<simb::MCTruth> mc = pi_serv->TrackIdToMCTruth_P(TrackerData.trkidtruth[iTrk][ipl]);
            TrackerData.trkorigin[iTrk][ipl] = mc->Origin();
            const simb::MCParticle *particle = pi_serv->TrackIdToParticle_P(TrackerData.trkidtruth[iTrk][ipl]);
            double tote = 0;
            const std::vector<const sim::IDE*> vide=bt_serv->TrackIdToSimIDEs_Ps(TrackerData.trkidtruth[iTrk][ipl]);
            for (auto ide: vide) {
              tote += ide->energy;
            }
            TrackerData.trkpdgtruth[iTrk][ipl] = particle->PdgCode();
            TrackerData.trkefftruth[iTrk][ipl] = maxe/(tote/kNplanes); //tote include both induction and collection energies
            //std::cout<<"\n"<<trkpdgtruth[iTracker][iTrk][ipl]<<"\t"<<trkefftruth[iTracker][iTrk][ipl];
          }
        }

        double maxe = 0;
        HitsBackTrack(allHits,TrackerData.trkg4id[iTrk],TrackerData.trkpurity[iTrk],maxe);
        if (TrackerData.trkg4id[iTrk]>0){
          const art::Ptr<simb::MCTruth> mc = pi_serv->TrackIdToMCTruth_P(TrackerData.trkg4id[iTrk]);
          TrackerData.trkorig[iTrk] = mc->Origin();
        }
        if (allHits.size()){
          std::vector<art::Ptr<recob::Hit> > all_hits;
          art::Handle<std::vector<recob::Hit> > hithandle;
          float totenergy = 0;
          if (evt.get(allHits[0].id(), hithandle)){
            art::fill_ptr_vector(all_hits, hithandle);
            for(size_t h = 0; h < all_hits.size(); ++h){

              art::Ptr<recob::Hit> hit = all_hits[h];
              std::vector<sim::IDE> ides;
              //bt_serv->HitToSimIDEs(hit,ides);
              std::vector<sim::TrackIDE> eveIDs = bt_serv->HitToEveTrackIDEs(hit);

              for(size_t e = 0; e < eveIDs.size(); ++e){
                //std::cout<<h<<" "<<e<<" "<<eveIDs[e].trackID<<" "<<eveIDs[e].energy<<" "<<eveIDs[e].energyFrac<<std::endl;
                if (eveIDs[e].trackID==TrackerData.trkg4id[iTrk]) totenergy += eveIDs[e].energy;
              }
            }
          }
          if (totenergy) TrackerData.trkcompleteness[iTrk] = maxe/totenergy;
        }
      }//end if (fIsMC)
*/
    }//end loop over track
      //std::cout << "HitIterator: " << HitIterator << std::endl;
  }//end loop over track module labels
}// end (fSaveTrackInfo)

  /*trkf::TrackMomentumCalculator trkm;
    std::cout<<"\t"<<trkm.GetTrackMomentum(200,2212)<<"\t"<<trkm.GetTrackMomentum(-10, 13)<<"\t"<<trkm.GetTrackMomentum(300,-19)<<"\n";
    */

  //Save Vertex information for multiple algorithms
  if (fSaveVertexInfo){
    for (unsigned int iVertexAlg=0; iVertexAlg < NVertexAlgos; ++iVertexAlg){
      AnaRootParserDataStruct::VertexDataStruct& VertexData = fData->GetVertexData(iVertexAlg);

      size_t NVertices = vertexlist[iVertexAlg].size();

      VertexData.SetMaxVertices(std::max(NVertices, (size_t) 1));
      VertexData.Clear(); // clear all the data

      VertexData.nvtx = (short) NVertices;

      // now set the tree addresses to the newly allocated memory;
      // this creates the tree branches in case they are not there yet
      SetVertexAddresses(iVertexAlg);
      if (NVertices > VertexData.GetMaxVertices()) {
        // got this error? it might be a bug,
        // since we are supposed to have allocated enough space to fit all tracks
        mf::LogError("AnaRootParser:limits") << "event has " << NVertices
          << " " << fVertexModuleLabel[iVertexAlg] << " tracks, only "
          << VertexData.GetMaxVertices() << " stored in tree";
      }

      for (size_t i = 0; i < NVertices && i < kMaxVertices ; ++i){//loop over hits
        VertexData.vtxId[i] = vertexlist[iVertexAlg][i]->ID();
        Double_t xyz[3] = {};
        vertexlist[iVertexAlg][i] -> XYZ(xyz);
        VertexData.vtxx[i] = xyz[0];
        VertexData.vtxy[i] = xyz[1];
        VertexData.vtxz[i] = xyz[2];

        if (fSavePFParticleInfo) {
          auto mapIter = vertexIDtoPFParticleIDMap.find(vertexlist[iVertexAlg][i]->ID());
          if (mapIter != vertexIDtoPFParticleIDMap.end()) {
            // This vertex has a corresponding PFParticle.
            VertexData.vtxhasPFParticle[i] = 1;
            VertexData.vtxPFParticleID[i] = mapIter->second;
          }
          else
            VertexData.vtxhasPFParticle[i] = 0;
        }

        // find PFParticle ID info
        art::FindMany<recob::PFParticle> fmPFParticle(vertexListHandle[iVertexAlg], evt, fPFParticleModuleLabel);
        if(fmPFParticle.isValid()) {
          std::vector<const recob::PFParticle*> pfparticles = fmPFParticle.at(i);
          if(pfparticles.size() > 1)
            std::cerr << "Warning: more than one associated PFParticle found for a vertex. Only one stored in tree." << std::endl;
          if (pfparticles.size() == 0)
            VertexData.vtxhasPFParticle[i] = 0;
          else {
            VertexData.vtxhasPFParticle[i] = 1;
            VertexData.vtxPFParticleID[i] = pfparticles.at(0)->Self();
          }
        } // fmPFParticle.isValid()
      }
    }
  }


  //mc truth information
  if (fIsMC){

    if (fSaveCryInfo){

      //store cry (cosmic generator information)
      fData->mcevts_truthcry = mclistcry.size();
      fData->cry_no_primaries = nCryPrimaries;
      //fData->cry_no_primaries;
      for(Int_t iPartc = 0; iPartc < mctruthcry->NParticles(); ++iPartc){
        const simb::MCParticle& partc(mctruthcry->GetParticle(iPartc));
        fData->cry_primaries_pdg[iPartc]=partc.PdgCode();
        fData->cry_Eng[iPartc]=partc.E();
        fData->cry_Px[iPartc]=partc.Px();
        fData->cry_Py[iPartc]=partc.Py();
        fData->cry_Pz[iPartc]=partc.Pz();
        fData->cry_P[iPartc]=partc.P();
        fData->cry_StartPointx[iPartc] = partc.Vx();
        fData->cry_StartPointy[iPartc] = partc.Vy();
        fData->cry_StartPointz[iPartc] = partc.Vz();
        fData->cry_StartPointt[iPartc] = partc.T();
        fData->cry_status_code[iPartc]=partc.StatusCode();
        fData->cry_mass[iPartc]=partc.Mass();
        fData->cry_trackID[iPartc]=partc.TrackId();
        fData->cry_ND[iPartc]=partc.NumberDaughters();
        fData->cry_mother[iPartc]=partc.Mother();
      } // for cry particles
    }// end fSaveCryInfo

    // Save the protoDUNE beam generator information
    if(fSaveProtoInfo){
      fData->proto_no_primaries = nProtoPrimaries;
      for(Int_t iPartp = 0; iPartp < nProtoPrimaries; ++iPartp){
        const simb::MCParticle& partp(mctruthproto->GetParticle(iPartp));

        fData->proto_isGoodParticle[iPartp] = (partp.Process() == "primary");
        fData->proto_vx[iPartp] = partp.Vx();
        fData->proto_vy[iPartp] = partp.Vy();
        fData->proto_vz[iPartp] = partp.Vz();
        fData->proto_t[iPartp] = partp.T();
        fData->proto_px[iPartp] = partp.Px();
        fData->proto_py[iPartp] = partp.Py();
        fData->proto_pz[iPartp] = partp.Pz();
        fData->proto_momentum[iPartp] = partp.P();
        fData->proto_energy[iPartp] = partp.E();
        fData->proto_pdg[iPartp] = partp.PdgCode();
      }
    }

    //save neutrino interaction information
    fData->mcevts_truth = mclist.size();
    if (fData->mcevts_truth > 0){//at least one mc record
      if (fSaveGenieInfo){
        int neutrino_i = 0;
        for(unsigned int iList = 0; (iList < mclist.size()) && (neutrino_i < kMaxTruth) ; ++iList){
          if (mclist[iList]->NeutrinoSet()){
            fData->nuPDG_truth[neutrino_i]  = mclist[iList]->GetNeutrino().Nu().PdgCode();
            fData->ccnc_truth[neutrino_i]   = mclist[iList]->GetNeutrino().CCNC();
            fData->mode_truth[neutrino_i]   = mclist[iList]->GetNeutrino().Mode();
            fData->Q2_truth[neutrino_i]     = mclist[iList]->GetNeutrino().QSqr();
            fData->W_truth[neutrino_i]      = mclist[iList]->GetNeutrino().W();
            fData->X_truth[neutrino_i]      = mclist[iList]->GetNeutrino().X();
            fData->Y_truth[neutrino_i]      = mclist[iList]->GetNeutrino().Y();
            fData->hitnuc_truth[neutrino_i] = mclist[iList]->GetNeutrino().HitNuc();
            fData->enu_truth[neutrino_i]    = mclist[iList]->GetNeutrino().Nu().E();
            fData->nuvtxx_truth[neutrino_i] = mclist[iList]->GetNeutrino().Nu().Vx();
            fData->nuvtxy_truth[neutrino_i] = mclist[iList]->GetNeutrino().Nu().Vy();
            fData->nuvtxz_truth[neutrino_i] = mclist[iList]->GetNeutrino().Nu().Vz();
            if (mclist[iList]->GetNeutrino().Nu().P()){
              fData->nu_dcosx_truth[neutrino_i] = mclist[iList]->GetNeutrino().Nu().Px()/mclist[iList]->GetNeutrino().Nu().P();
              fData->nu_dcosy_truth[neutrino_i] = mclist[iList]->GetNeutrino().Nu().Py()/mclist[iList]->GetNeutrino().Nu().P();
              fData->nu_dcosz_truth[neutrino_i] = mclist[iList]->GetNeutrino().Nu().Pz()/mclist[iList]->GetNeutrino().Nu().P();
            }
            fData->lep_mom_truth[neutrino_i] = mclist[iList]->GetNeutrino().Lepton().P();
            if (mclist[iList]->GetNeutrino().Lepton().P()){
              fData->lep_dcosx_truth[neutrino_i] = mclist[iList]->GetNeutrino().Lepton().Px()/mclist[iList]->GetNeutrino().Lepton().P();
              fData->lep_dcosy_truth[neutrino_i] = mclist[iList]->GetNeutrino().Lepton().Py()/mclist[iList]->GetNeutrino().Lepton().P();
              fData->lep_dcosz_truth[neutrino_i] = mclist[iList]->GetNeutrino().Lepton().Pz()/mclist[iList]->GetNeutrino().Lepton().P();
            }

            //flux information
            //
            // Double-check that a simb::MCFlux object is associated with the
            // current simb::MCTruth object. For GENIE events, these should
            // always accompany each other. Other generators (e.g., MARLEY) may
            // create simb::MCTruth objects without corresponding simb::MCFlux
            // objects. -- S. Gardiner
            /*	    if (find_mcflux.isValid()) {
                    auto flux_maybe_ref = find_mcflux.at(iList);
                    if (flux_maybe_ref.isValid()) {
                    auto flux_ref = flux_maybe_ref.ref();
                    fData->vx_flux[neutrino_i]        = flux_ref.fvx;
                    fData->vy_flux[neutrino_i]        = flux_ref.fvy;
                    fData->vz_flux[neutrino_i]        = flux_ref.fvz;
                    fData->pdpx_flux[neutrino_i]      = flux_ref.fpdpx;
                    fData->pdpy_flux[neutrino_i]      = flux_ref.fpdpy;
                    fData->pdpz_flux[neutrino_i]      = flux_ref.fpdpz;
                    fData->ppdxdz_flux[neutrino_i]    = flux_ref.fppdxdz;
                    fData->ppdydz_flux[neutrino_i]    = flux_ref.fppdydz;
                    fData->pppz_flux[neutrino_i]      = flux_ref.fpppz;

                    fData->ptype_flux[neutrino_i]      = flux_ref.fptype;
                    fData->ppvx_flux[neutrino_i]       = flux_ref.fppvx;
                    fData->ppvy_flux[neutrino_i]       = flux_ref.fppvy;
                    fData->ppvz_flux[neutrino_i]       = flux_ref.fppvz;
                    fData->muparpx_flux[neutrino_i]    = flux_ref.fmuparpx;
                    fData->muparpy_flux[neutrino_i]    = flux_ref.fmuparpy;
                    fData->muparpz_flux[neutrino_i]    = flux_ref.fmuparpz;
                    fData->mupare_flux[neutrino_i]     = flux_ref.fmupare;

                    fData->tgen_flux[neutrino_i]     = flux_ref.ftgen;
                    fData->tgptype_flux[neutrino_i]  = flux_ref.ftgptype;
                    fData->tgppx_flux[neutrino_i]    = flux_ref.ftgppx;
                    fData->tgppy_flux[neutrino_i]    = flux_ref.ftgppy;
                    fData->tgppz_flux[neutrino_i]    = flux_ref.ftgppz;
                    fData->tprivx_flux[neutrino_i]   = flux_ref.ftprivx;
                    fData->tprivy_flux[neutrino_i]   = flux_ref.ftprivy;
                    fData->tprivz_flux[neutrino_i]   = flux_ref.ftprivz;

                    fData->dk2gen_flux[neutrino_i]   = flux_ref.fdk2gen;
                    fData->gen2vtx_flux[neutrino_i]   = flux_ref.fgen2vtx;

                    fData->tpx_flux[neutrino_i]    = flux_ref.ftpx;
                    fData->tpy_flux[neutrino_i]    = flux_ref.ftpy;
                    fData->tpz_flux[neutrino_i]    = flux_ref.ftpz;
                    fData->tptype_flux[neutrino_i] = flux_ref.ftptype;
                    } // flux_maybe_ref.isValid()
                    } // find_mcflux.isValid() */
            neutrino_i++;
          }//mclist is NeutrinoSet()
        }//loop over mclist

        if (mctruth->NeutrinoSet()){
          //genie particles information
          fData->genie_no_primaries = mctruth->NParticles();

          size_t StoreParticles = std::min((size_t) fData->genie_no_primaries, fData->GetMaxGeniePrimaries());
          if (fData->genie_no_primaries > (int) StoreParticles) {
            // got this error? it might be a bug,
            // since the structure should have enough room for everything
            mf::LogError("AnaRootParser:limits") << "event has "
              << fData->genie_no_primaries << " MC particles, only "
              << StoreParticles << " stored in tree";
          }
          for(size_t iPart = 0; iPart < StoreParticles; ++iPart){
            const simb::MCParticle& part(mctruth->GetParticle(iPart));
            fData->genie_primaries_pdg[iPart]=part.PdgCode();
            fData->genie_Eng[iPart]=part.E();
            fData->genie_Px[iPart]=part.Px();
            fData->genie_Py[iPart]=part.Py();
            fData->genie_Pz[iPart]=part.Pz();
            fData->genie_P[iPart]=part.P();
            fData->genie_status_code[iPart]=part.StatusCode();
            fData->genie_mass[iPart]=part.Mass();
            fData->genie_trackID[iPart]=part.TrackId();
            fData->genie_ND[iPart]=part.NumberDaughters();
            fData->genie_mother[iPart]=part.Mother();
          } // for particle
          //const simb::MCNeutrino& nu(mctruth->GetNeutrino());
        } //if neutrino set
      }// end (fSaveGenieInfo)

      //Extract MC Shower information and fill the Shower branches
      if (fSaveMCShowerInfo){
        fData->no_mcshowers = nMCShowers;
        size_t shwr = 0;
        for(std::vector<sim::MCShower>::const_iterator imcshwr = mcshowerh->begin();
            imcshwr != mcshowerh->end(); ++imcshwr) {
          const sim::MCShower& mcshwr = *imcshwr;
          fData->mcshwr_origin[shwr]          = mcshwr.Origin();
          fData->mcshwr_pdg[shwr]	      = mcshwr.PdgCode();
          fData->mcshwr_TrackId[shwr]	      = mcshwr.TrackID();
          fData->mcshwr_Process[shwr]	      = mcshwr.Process();
          fData->mcshwr_startX[shwr]          = mcshwr.Start().X();
          fData->mcshwr_startY[shwr]          = mcshwr.Start().Y();
          fData->mcshwr_startZ[shwr]          = mcshwr.Start().Z();
          fData->mcshwr_endX[shwr]            = mcshwr.End().X();
          fData->mcshwr_endY[shwr]            = mcshwr.End().Y();
          fData->mcshwr_endZ[shwr]            = mcshwr.End().Z();
          if (mcshwr.DetProfile().E()!= 0){
            fData->mcshwr_isEngDeposited[shwr] = 1;
            fData->mcshwr_CombEngX[shwr]        = mcshwr.DetProfile().X();
            fData->mcshwr_CombEngY[shwr]        = mcshwr.DetProfile().Y();
            fData->mcshwr_CombEngZ[shwr]        = mcshwr.DetProfile().Z();
            fData->mcshwr_CombEngPx[shwr]       = mcshwr.DetProfile().Px();
            fData->mcshwr_CombEngPy[shwr]       = mcshwr.DetProfile().Py();
            fData->mcshwr_CombEngPz[shwr]       = mcshwr.DetProfile().Pz();
            fData->mcshwr_CombEngE[shwr]        = mcshwr.DetProfile().E();
            fData->mcshwr_dEdx[shwr]            = mcshwr.dEdx();
            fData->mcshwr_StartDirX[shwr]       = mcshwr.StartDir().X();
            fData->mcshwr_StartDirY[shwr]       = mcshwr.StartDir().Y();
            fData->mcshwr_StartDirZ[shwr]       = mcshwr.StartDir().Z();
          }
          else
            fData->mcshwr_isEngDeposited[shwr] = 0;
          fData->mcshwr_Motherpdg[shwr]       = mcshwr.MotherPdgCode();
          fData->mcshwr_MotherTrkId[shwr]     = mcshwr.MotherTrackID();
          fData->mcshwr_MotherProcess[shwr]   = mcshwr.MotherProcess();
          fData->mcshwr_MotherstartX[shwr]    = mcshwr.MotherStart().X();
          fData->mcshwr_MotherstartY[shwr]    = mcshwr.MotherStart().Y();
          fData->mcshwr_MotherstartZ[shwr]    = mcshwr.MotherStart().Z();
          fData->mcshwr_MotherendX[shwr]      = mcshwr.MotherEnd().X();
          fData->mcshwr_MotherendY[shwr]      = mcshwr.MotherEnd().Y();
          fData->mcshwr_MotherendZ[shwr]      = mcshwr.MotherEnd().Z();
          fData->mcshwr_Ancestorpdg[shwr]     = mcshwr.AncestorPdgCode();
          fData->mcshwr_AncestorTrkId[shwr]   = mcshwr.AncestorTrackID();
          fData->mcshwr_AncestorProcess[shwr] = mcshwr.AncestorProcess();
          fData->mcshwr_AncestorstartX[shwr]  = mcshwr.AncestorStart().X();
          fData->mcshwr_AncestorstartY[shwr]  = mcshwr.AncestorStart().Y();
          fData->mcshwr_AncestorstartZ[shwr]  = mcshwr.AncestorStart().Z();
          fData->mcshwr_AncestorendX[shwr]    = mcshwr.AncestorEnd().X();
          fData->mcshwr_AncestorendY[shwr]    = mcshwr.AncestorEnd().Y();
          fData->mcshwr_AncestorendZ[shwr]    = mcshwr.AncestorEnd().Z();
          ++shwr;
        }
        fData->mcshwr_Process.resize(shwr);
        fData->mcshwr_MotherProcess.resize(shwr);
        fData->mcshwr_AncestorProcess.resize(shwr);
      }//End if (fSaveMCShowerInfo){

      //Extract MC Track information and fill the Shower branches
      if (fSaveMCTrackInfo){
        fData->no_mctracks = nMCTracks;
        size_t trk = 0;
        for(std::vector<sim::MCTrack>::const_iterator imctrk = mctrackh->begin();imctrk != mctrackh->end(); ++imctrk) {
          const sim::MCTrack& mctrk = *imctrk;
          TLorentzVector tpcstart, tpcend, tpcmom;
          double plen = driftedLength(mctrk, tpcstart, tpcend, tpcmom);
          fData->mctrk_origin[trk]          = mctrk.Origin();
          fData->mctrk_pdg[trk]	            = mctrk.PdgCode();
          fData->mctrk_TrackId[trk]	    = mctrk.TrackID();
          fData->mctrk_Process[trk]	    = mctrk.Process();
          fData->mctrk_startX[trk]          = mctrk.Start().X();
          fData->mctrk_startY[trk]          = mctrk.Start().Y();
          fData->mctrk_startZ[trk]          = mctrk.Start().Z();
          fData->mctrk_endX[trk]            = mctrk.End().X();
          fData->mctrk_endY[trk]            = mctrk.End().Y();
          fData->mctrk_endZ[trk]            = mctrk.End().Z();
          fData->mctrk_Motherpdg[trk]       = mctrk.MotherPdgCode();
          fData->mctrk_MotherTrkId[trk]     = mctrk.MotherTrackID();
          fData->mctrk_MotherProcess[trk]   = mctrk.MotherProcess();
          fData->mctrk_MotherstartX[trk]    = mctrk.MotherStart().X();
          fData->mctrk_MotherstartY[trk]    = mctrk.MotherStart().Y();
          fData->mctrk_MotherstartZ[trk]    = mctrk.MotherStart().Z();
          fData->mctrk_MotherendX[trk]      = mctrk.MotherEnd().X();
          fData->mctrk_MotherendY[trk]      = mctrk.MotherEnd().Y();
          fData->mctrk_MotherendZ[trk]      = mctrk.MotherEnd().Z();
          fData->mctrk_Ancestorpdg[trk]     = mctrk.AncestorPdgCode();
          fData->mctrk_AncestorTrkId[trk]   = mctrk.AncestorTrackID();
          fData->mctrk_AncestorProcess[trk] = mctrk.AncestorProcess();
          fData->mctrk_AncestorstartX[trk]  = mctrk.AncestorStart().X();
          fData->mctrk_AncestorstartY[trk]  = mctrk.AncestorStart().Y();
          fData->mctrk_AncestorstartZ[trk]  = mctrk.AncestorStart().Z();
          fData->mctrk_AncestorendX[trk]    = mctrk.AncestorEnd().X();
          fData->mctrk_AncestorendY[trk]    = mctrk.AncestorEnd().Y();
          fData->mctrk_AncestorendZ[trk]    = mctrk.AncestorEnd().Z();

          fData->mctrk_len_drifted[trk]       = plen;

          if (plen != 0){
            fData->mctrk_startX_drifted[trk] = tpcstart.X();
            fData->mctrk_startY_drifted[trk] = tpcstart.Y();
            fData->mctrk_startZ_drifted[trk] = tpcstart.Z();
            fData->mctrk_endX_drifted[trk]   = tpcend.X();
            fData->mctrk_endY_drifted[trk]   = tpcend.Y();
            fData->mctrk_endZ_drifted[trk]   = tpcend.Z();
            fData->mctrk_p_drifted[trk]      = tpcmom.Vect().Mag();
            fData->mctrk_px_drifted[trk]     = tpcmom.X();
            fData->mctrk_py_drifted[trk]     = tpcmom.Y();
            fData->mctrk_pz_drifted[trk]     = tpcmom.Z();
          }
          ++trk;
        }

        fData->mctrk_Process.resize(trk);
        fData->mctrk_MotherProcess.resize(trk);
        fData->mctrk_AncestorProcess.resize(trk);
      }//End if (fSaveMCTrackInfo){


      //Photon particles information
      if (fSavePhotonInfo){
	int photoncounter=0;
        for ( auto const& pmt : (*photonHandle) )
        {
          std::map<int, int> PhotonsMap = pmt.DetectedPhotons;

          for(auto itphoton = PhotonsMap.begin(); itphoton!= PhotonsMap.end(); itphoton++)
          {
            for(int iphotonatthistime = 0; iphotonatthistime < itphoton->second ; iphotonatthistime++)
            {
	      fData->photons_time[photoncounter]=itphoton->first;
	      fData->photons_channel[photoncounter]=pmt.OpChannel;
	      photoncounter++;
	    }
	  }
        }
	fData->numberofphotons=photoncounter;
      }

      //Generator particles information
      if (fSaveGeneratorInfo){
        const sim::ParticleList& plist = pi_serv->ParticleList();

        size_t generator_particle=0;
        sim::ParticleList::const_iterator itPart = plist.begin(),
        pend = plist.end(); // iterator to pairs (track id, particle)

        for(size_t iPart = 0; (iPart < plist.size()) && (itPart != pend); ++iPart){

          const simb::MCParticle* pPart = (itPart++)->second;
          if (!pPart) {
            throw art::Exception(art::errors::LogicError)
              << "GEANT particle #" << iPart << " returned a null pointer";
          }

          if (iPart < fData->GetMaxGeneratorparticles()) {

            std::string pri("primary");
	    if( pPart->Mother() == 0 && pPart->Process() == pri )
	    {
              fData->TrackId[generator_particle]=pPart->TrackId();
              fData->pdg[generator_particle]=pPart->PdgCode();
              fData->status[generator_particle] = pPart->StatusCode();
              fData->Eng[generator_particle]=pPart->E();
              fData->Mass[generator_particle]=pPart->Mass();
              fData->Px[generator_particle]=pPart->Px();
              fData->Py[generator_particle]=pPart->Py();
              fData->Pz[generator_particle]=pPart->Pz();
              fData->P[generator_particle]=pPart->Momentum().Vect().Mag();
              fData->StartPointx[generator_particle]=pPart->Vx();
              fData->StartPointy[generator_particle]=pPart->Vy();
              fData->StartPointz[generator_particle]=pPart->Vz();
              fData->StartT[generator_particle] = pPart->T();

      	      TVector3 momentum_start_flipped;
              momentum_start_flipped.SetXYZ(pPart->Pz(), pPart->Py(), pPart->Px());

	      fData->theta[generator_particle] = (180.0/3.14159)*momentum_start_flipped.Theta();
	      fData->phi[generator_particle] = (180.0/3.14159)*momentum_start_flipped.Phi();
	    }
            ++generator_particle;
          }
          else if (iPart == fData->GetMaxGeneratorparticles()) {
            // got this error? it might be a bug,
            // since the structure should have enough room for everything
            mf::LogError("AnaRootParser:limits") << "event has "
              << plist.size() << " MC particles, only "
              << fData->GetMaxGeneratorparticles() << " will be stored in tree";
          }
        } // for particles
        fData->generator_list_size = generator_particle;
        fData->processname.resize(generator_particle);
        LOG_DEBUG("AnaRootParser")
          << "Counted "
          << fData->generator_list_size << " Generator particles";
      }//fSaveGeneratorInfo




      //GEANT particles information
      if (fSaveGeantInfo){

        const sim::ParticleList& plist = pi_serv->ParticleList();

        std::string pri("primary");
        int primary=0;
        int active = 0;
        size_t geant_particle=0;
        sim::ParticleList::const_iterator itPart = plist.begin(),
          pend = plist.end(); // iterator to pairs (track id, particle)

        // helper map track ID => index
        std::map<int, size_t> TrackIDtoIndex;
        std::vector<int> gpdg;
        std::vector<int> gmother;
        for(size_t iPart = 0; (iPart < plist.size()) && (itPart != pend); ++iPart){
          const simb::MCParticle* pPart = (itPart++)->second;
          if (!pPart) {
            throw art::Exception(art::errors::LogicError)
              << "GEANT particle #" << iPart << " returned a null pointer";
          }

          bool isPrimary = pPart->Process() == pri;
          int TrackID = pPart->TrackId();
          TrackIDtoIndex.emplace(TrackID, iPart);
          gpdg.push_back(pPart->PdgCode());
          gmother.push_back(pPart->Mother());
          if (iPart < fData->GetMaxGEANTparticles()) {
//            if (pPart->E()<fG4minE&&(!isPrimary)) continue;
            if (isPrimary) ++primary;

            TLorentzVector mcstart, mcend;
            unsigned int pstarti, pendi;
            double plen = length(*pPart, mcstart, mcend, pstarti, pendi);
            bool isActive = plen != 0;

//            TLorentzVector mcstartdrifted, mcenddrifted;
//	      unsigned int pstartdriftedi, penddriftedi;
//            double plendrifted = driftedLength(*pPart, mcstartdrifted, mcenddrifted, pstartdriftedi, penddriftedi);
//            bool isDrifted = plendrifted!= 0;

	    if(fLogLevel == 3)
	    {
	      std::cout << "pPart->TrackId():" << pPart->TrackId() << std::endl;
	      std::cout << "pPart->Mother():" << pPart->Mother() << std::endl;
	      std::cout << "pPart->PdgCode():" << pPart->PdgCode() << std::endl;
	      std::cout << std::endl;
	    }


            fData->process_primary[geant_particle] = int(isPrimary);
            fData->processname[geant_particle]= pPart->Process();
            fData->Mother[geant_particle]=pPart->Mother();
            fData->TrackId[geant_particle]=pPart->TrackId();
            fData->pdg[geant_particle]=pPart->PdgCode();
            fData->status[geant_particle] = pPart->StatusCode();
            fData->Eng[geant_particle]=pPart->E();
            fData->EndE[geant_particle]=pPart->EndE();
            fData->Mass[geant_particle]=pPart->Mass();
            fData->Px[geant_particle]=pPart->Px();
            fData->Py[geant_particle]=pPart->Py();
            fData->Pz[geant_particle]=pPart->Pz();
            fData->P[geant_particle]=pPart->Momentum().Vect().Mag();
            fData->StartPointx[geant_particle]=pPart->Vx();
            fData->StartPointy[geant_particle]=pPart->Vy();
            fData->StartPointz[geant_particle]=pPart->Vz();
            fData->StartT[geant_particle] = pPart->T();
            fData->EndPointx[geant_particle]=pPart->EndPosition()[0];
            fData->EndPointy[geant_particle]=pPart->EndPosition()[1];
            fData->EndPointz[geant_particle]=pPart->EndPosition()[2];
            fData->EndT[geant_particle] = pPart->EndT();

	    fData->NTrajectoryPointsPerParticle[geant_particle] = pPart->NumberTrajectoryPoints();

	    //fData->theta[geant_particle] = pPart->Momentum().Theta();
	    //fData->phi[geant_particle] = pPart->Momentum().Phi();
	    //Change definition of theta and phi (swap x and z coordinate since x is "up" in dual phase)
      	    TVector3 momentum_start_flipped;
            momentum_start_flipped.SetXYZ(pPart->Pz(), pPart->Py(), pPart->Px());

	    fData->theta[geant_particle] = (180.0/3.14159)*momentum_start_flipped.Theta();
	    fData->phi[geant_particle] = (180.0/3.14159)*momentum_start_flipped.Phi();

            fData->theta_xz[geant_particle] = std::atan2(pPart->Px(), pPart->Pz());
            fData->theta_yz[geant_particle] = std::atan2(pPart->Py(), pPart->Pz());
//            fData->pathlen_drifted[geant_particle]  = plendrifted;
            fData->NumberDaughters[geant_particle]=pPart->NumberDaughters();
            fData->inTPCActive[geant_particle] = int(isActive);
//            fData->inTPCDrifted[geant_particle] = int(isDrifted);
            art::Ptr<simb::MCTruth> const& mc_truth = pi_serv->ParticleToMCTruth_P(pPart);
            if (mc_truth){
              fData->origin[geant_particle] = mc_truth->Origin();
              fData->MCTruthIndex[geant_particle] = mc_truth.key();
            }

            if (isActive && fSaveGeantInAVInfo){
              fData->pathlen_tpcAV[active]  = plen;
	      fData->TrackId_tpcAV[active] =  pPart->TrackId();
	      fData->PDGCode_tpcAV[active] = pPart->PdgCode();

              fData->StartPointx_tpcAV[active] = mcstart.X();
              fData->StartPointy_tpcAV[active] = mcstart.Y();
              fData->StartPointz_tpcAV[active] = mcstart.Z();
              fData->StartT_tpcAV[active] = mcstart.T();
              fData->StartE_tpcAV[active] = pPart->E(pstarti);
              fData->StartP_tpcAV[active] = pPart->P(pstarti);
              fData->StartPx_tpcAV[active] = pPart->Px(pstarti);
              fData->StartPy_tpcAV[active] = pPart->Py(pstarti);
              fData->StartPz_tpcAV[active] = pPart->Pz(pstarti);
              fData->EndPointx_tpcAV[active] = mcend.X();
              fData->EndPointy_tpcAV[active] = mcend.Y();
              fData->EndPointz_tpcAV[active] = mcend.Z();
              fData->EndT_tpcAV[active] = mcend.T();
              fData->EndE_tpcAV[active] = pPart->E(pendi);
              fData->EndP_tpcAV[active] = pPart->P(pendi);
              fData->EndPx_tpcAV[active] = pPart->Px(pendi);
              fData->EndPy_tpcAV[active] = pPart->Py(pendi);
              fData->EndPz_tpcAV[active] = pPart->Pz(pendi);

	      //Change definition of theta and phi (swap x and z coordinate since x is "up" in dual phase)
      	      TVector3 momentum_start_tpcAv_flipped;
              momentum_start_tpcAv_flipped.SetXYZ(pPart->Pz(pstarti), pPart->Py(pstarti), pPart->Px(pstarti));
      	      TVector3 momentum_end_tpcAv_flipped;
              momentum_end_tpcAv_flipped.SetXYZ(pPart->Pz(pendi), pPart->Py(pendi), pPart->Px(pendi));

              fData->thetastart_tpcAV[active] = (180.0/3.14159)*momentum_start_tpcAv_flipped.Theta();
              fData->phistart_tpcAV[active] = (180.0/3.14159)*momentum_start_tpcAv_flipped.Phi();

              fData->thetaend_tpcAV[active] = (180.0/3.14159)*momentum_end_tpcAv_flipped.Theta();
              fData->phiend_tpcAV[active] = (180.0/3.14159)*momentum_end_tpcAv_flipped.Phi();

	      active++;
            }

/*
            if (isDrifted){
              fData->StartPointx_drifted[geant_particle] = mcstartdrifted.X();
              fData->StartPointy_drifted[geant_particle] = mcstartdrifted.Y();
              fData->StartPointz_drifted[geant_particle] = mcstartdrifted.Z();
              fData->StartT_drifted[geant_particle] = mcstartdrifted.T();
              fData->StartE_drifted[geant_particle] = pPart->E(pstartdriftedi);
              fData->StartP_drifted[geant_particle] = pPart->P(pstartdriftedi);
              fData->StartPx_drifted[geant_particle] = pPart->Px(pstartdriftedi);
              fData->StartPy_drifted[geant_particle] = pPart->Py(pstartdriftedi);
              fData->StartPz_drifted[geant_particle] = pPart->Pz(pstartdriftedi);
              fData->EndPointx_drifted[geant_particle] = mcenddrifted.X();
              fData->EndPointy_drifted[geant_particle] = mcenddrifted.Y();
              fData->EndPointz_drifted[geant_particle] = mcenddrifted.Z();
              fData->EndT_drifted[geant_particle] = mcenddrifted.T();
              fData->EndE_drifted[geant_particle] = pPart->E(penddriftedi);
              fData->EndP_drifted[geant_particle] = pPart->P(penddriftedi);
              fData->EndPx_drifted[geant_particle] = pPart->Px(penddriftedi);
              fData->EndPy_drifted[geant_particle] = pPart->Py(penddriftedi);
              fData->EndPz_drifted[geant_particle] = pPart->Pz(penddriftedi);
            }
*/
            //access auxiliary detector parameters
            if (fSaveAuxDetInfo) {
              unsigned short nAD = 0; // number of cells that particle hit

              // find deposit of this particle in each of the detector cells
              for (const sim::AuxDetSimChannel* c: fAuxDetSimChannels) {

                // find if this cell has a contribution (IDE) from this particle,
                // and which one
                const std::vector<sim::AuxDetIDE>& setOfIDEs = c->AuxDetIDEs();
                // using a C++ "lambda" function here; this one:
                // - sees only TrackID from the current scope
                // - takes one parameter: the AuxDetIDE to be tested
                // - returns if that IDE belongs to the track we are looking for
                std::vector<sim::AuxDetIDE>::const_iterator iIDE
                  = std::find_if(
                      setOfIDEs.begin(), setOfIDEs.end(),
                      [TrackID](const sim::AuxDetIDE& IDE){ return IDE.trackID == TrackID; }
                      );
                if (iIDE == setOfIDEs.end()) continue;

                // now iIDE points to the energy released by the track #i (TrackID)

                // look for IDE with matching trackID
                // find trackIDs stored in setOfIDEs with the same trackID, but negative,
                // this is an untracked particle who's energy should be added as deposited by this original trackID
                float totalE = 0.; // total energy deposited around by the GEANT particle in this cell
                for(const auto& adtracks: setOfIDEs) {
                  if( fabs(adtracks.trackID) == TrackID )
                    totalE += adtracks.energyDeposited;
                } // for

                // fill the structure
                if (nAD < kMaxAuxDets) {
                  fData->AuxDetID[geant_particle][nAD] = c->AuxDetID();
                  fData->entryX[geant_particle][nAD]   = iIDE->entryX;
                  fData->entryY[geant_particle][nAD]   = iIDE->entryY;
                  fData->entryZ[geant_particle][nAD]   = iIDE->entryZ;
                  fData->entryT[geant_particle][nAD]   = iIDE->entryT;
                  fData->exitX[geant_particle][nAD]    = iIDE->exitX;
                  fData->exitY[geant_particle][nAD]    = iIDE->exitY;
                  fData->exitZ[geant_particle][nAD]    = iIDE->exitZ;
                  fData->exitT[geant_particle][nAD]    = iIDE->exitT;
                  fData->exitPx[geant_particle][nAD]   = iIDE->exitMomentumX;
                  fData->exitPy[geant_particle][nAD]   = iIDE->exitMomentumY;
                  fData->exitPz[geant_particle][nAD]   = iIDE->exitMomentumZ;
                  fData->CombinedEnergyDep[geant_particle][nAD] = totalE;
                }
                ++nAD;
              } // for aux det sim channels
              fData->NAuxDets[geant_particle] = nAD;

              if (nAD > kMaxAuxDets) {
                // got this error? consider increasing kMaxAuxDets
                mf::LogError("AnaRootParser:limits")
                  << "particle #" << iPart
                  << " touches " << nAD << " auxiliary detector cells, only "
                  << kMaxAuxDets << " of them are saved in the tree";
              } // if too many detector cells
            } // if (fSaveAuxDetInfo)

            ++geant_particle;
          }
          else if (iPart == fData->GetMaxGEANTparticles()) {
            // got this error? it might be a bug,
            // since the structure should have enough room for everything
            mf::LogError("AnaRootParser:limits") << "event has "
              << plist.size() << " MC particles, only "
              << fData->GetMaxGEANTparticles() << " will be stored in tree";
          }
        } // for particles

        fData->geant_list_size_in_tpcAV = active;
        fData->no_primaries = primary;
        fData->geant_list_size = geant_particle;
        fData->processname.resize(geant_particle);
        LOG_DEBUG("AnaRootParser")
          << "Counted "
          << fData->geant_list_size << " GEANT particles ("
          << fData->geant_list_size_in_tpcAV << " in AV), "
          << fData->no_primaries << " primaries, "
          << fData->genie_no_primaries << " GENIE particles";

        FillWith(fData->MergedId, 0);

        // for each particle, consider all the direct ancestors with the same
        // PDG ID, and mark them as belonging to the same "group"
        // (having the same MergedId)


        int currentMergedId = 1;
        for(size_t iPart = 0; iPart < geant_particle; ++iPart){
        // if the particle already belongs to a group, don't bother
        if (fData->MergedId[iPart]) continue;
        // the particle starts its own group
        fData->MergedId[iPart] = currentMergedId;
        int currentMotherTrackId = fData->Mother[iPart];
        while (currentMotherTrackId > 0) {
        if (TrackIDtoIndex.find(currentMotherTrackId)==TrackIDtoIndex.end()) break;
        unsigned int gindex = TrackIDtoIndex[currentMotherTrackId];
        if (gindex>=plist.size()) break;
        // if the mother particle is of a different type,
        // don't bother with iPart ancestry any further
        if (gpdg[gindex]!=fData->pdg[iPart]) break;
        if (TrackIDtoIndex.find(currentMotherTrackId)!=TrackIDtoIndex.end()){
        unsigned int igeantMother = TrackIDtoIndex[currentMotherTrackId];
        if (igeantMother<geant_particle){
        fData->MergedId[igeantMother] = currentMergedId;
        }
        }
        currentMotherTrackId = gmother[gindex];
        }
        ++currentMergedId;
        }// for merging check

      } // if (fSaveGeantInfo)

      //GEANT trajectory information
      if (fSaveGeantTrajectoryInfo){

        const sim::ParticleList& plist = pi_serv->ParticleList();
        sim::ParticleList::const_iterator itPart = plist.begin(),
          pend = plist.end(); // iterator to pairs (track id, particle)
	int trajpointcounter = 0;

        for(size_t iPart = 0; (iPart < plist.size()) && (itPart != pend); ++iPart)
	{
          const simb::MCParticle* pPart = (itPart++)->second;
          if (!pPart)
	  {
            throw art::Exception(art::errors::LogicError)
              << "GEANT particle #" << iPart << " returned a null pointer";
          }

          unsigned int numberTrajectoryPoints = pPart->NumberTrajectoryPoints();
          for(unsigned int trajpoint=0; trajpoint < numberTrajectoryPoints; ++trajpoint)
          {
            const TLorentzVector& trajpointPosition= pPart->Position(trajpoint);
            //const TLorentzVector& trajpointMomentum= pPart->Momentum(trajpoint);

	    fData->TrajTrackId[trajpointcounter] = pPart->TrackId();
	    fData->TrajPDGCode[trajpointcounter] = pPart->PdgCode();

	    fData->TrajX[trajpointcounter] = trajpointPosition.X();
	    fData->TrajY[trajpointcounter] = trajpointPosition.Y();
	    fData->TrajZ[trajpointcounter] = trajpointPosition.Z();
	    fData->TrajT[trajpointcounter] = pPart->T(trajpoint);
	    fData->TrajE[trajpointcounter] = pPart->E(trajpoint);
	    fData->TrajP[trajpointcounter] = pPart->P(trajpoint);
	    fData->TrajPx[trajpointcounter] = pPart->Px(trajpoint);
	    fData->TrajPy[trajpointcounter] = pPart->Py(trajpoint);
	    fData->TrajPz[trajpointcounter] = pPart->Pz(trajpoint);

      	    TVector3 trajpointMomentum_flipped;
            trajpointMomentum_flipped.SetXYZ(pPart->Pz(trajpoint), pPart->Py(trajpoint), pPart->Px(trajpoint));

            fData->TrajTheta[trajpointcounter] = (180.0/3.14159)*trajpointMomentum_flipped.Theta();
            fData->TrajPhi[trajpointcounter] = (180.0/3.14159)*trajpointMomentum_flipped.Phi();

	    if(fLogLevel == 4)
	    {
	    std::cout << std::endl;
	    std::cout << "trajpointcounter: " << trajpointcounter << std::endl;
	    std::cout << "fData->TrajX[trajpointcounter]: " << fData->TrajX[trajpointcounter] << std::endl;
	    std::cout << "fData->TrajY[trajpointcounter]: " << fData->TrajY[trajpointcounter] << std::endl;
	    std::cout << "fData->TrajZ[trajpointcounter]: " << fData->TrajZ[trajpointcounter] << std::endl;
	    std::cout << "fData->TrajT[trajpointcounter]: " << fData->TrajT[trajpointcounter] << std::endl;
	    std::cout << "fData->TrajE[trajpointcounter]: " << fData->TrajE[trajpointcounter] << std::endl;
	    std::cout << "fData->TrajP[trajpointcounter]: " << fData->TrajP[trajpointcounter] << std::endl;
	    std::cout << "fData->TrajPx[trajpointcounter]: " << fData->TrajPx[trajpointcounter] << std::endl;
	    std::cout << "fData->TrajPy[trajpointcounter]: " << fData->TrajPy[trajpointcounter] << std::endl;
	    std::cout << "fData->TrajPz[trajpointcounter]: " << fData->TrajPz[trajpointcounter] << std::endl;
	    std::cout << "fData->TrajTheta[trajpointcounter]: " << fData->TrajTheta[trajpointcounter] << std::endl;
	    std::cout << "fData->TrajPhi[trajpointcounter]: " << fData->TrajPhi[trajpointcounter] << std::endl;
	    }

	    trajpointcounter++;
	  }
        } // for particles

	fData->geant_trajectory_size = trajpointcounter;

      }// if (fSaveGeantTrajectoryInfo)

    }//if (mcevts_truth)
    }//if (fIsMC){

    fData->taulife = detprop->ElectronLifetime();

    fTree->Fill();

    if (mf::isDebugEnabled()) {
      // use mf::LogDebug instead of LOG_DEBUG because we reuse it in many lines;
      // thus, we protect this part of the code with the line above
      mf::LogDebug logStream("AnaRootParserStructure");
      logStream
        << "Tree data structure contains:"
        << "\n - " << fData->no_hits << " hits (" << fData->GetMaxHits() << ")"
        << "\n - " << fData->NHitsInAllTracks << " hits (" << fData->GetMaxHits() << ")"
        << "\n - " << fData->genie_no_primaries << " genie primaries (" << fData->GetMaxGeniePrimaries() << ")"
        << "\n - " << fData->geant_list_size << " GEANT particles (" << fData->GetMaxGEANTparticles() << "), "
        << fData->no_primaries << " primaries"
        << "\n - " << fData->geant_list_size_in_tpcAV << " GEANT particles in AV "
        << "\n - " << ((int) fData->kNTracker) << " trackers:"
        ;

      size_t iTracker = 0;
      for (auto tracker = fData->TrackData.cbegin();
          tracker != fData->TrackData.cend(); ++tracker, ++iTracker
          ) {
        logStream
          << "\n -> " << tracker->ntracks << " " << fTrackModuleLabel[iTracker]
          << " tracks (" << tracker->GetMaxTracks() << ")"
          ;
        for (int iTrk = 0; iTrk < tracker->ntracks; ++iTrk) {
          logStream << "\n    [" << iTrk << "] "<< tracker->ntrkhitsperview[iTrk][0];
          for (size_t ipl = 1; ipl < tracker->GetMaxPlanesPerTrack(iTrk); ++ipl)
            logStream << " + " << tracker->ntrkhitsperview[iTrk][ipl];
          logStream << " hits (" << tracker->GetMaxHitsPerTrack(iTrk, 0);
          for (size_t ipl = 1; ipl < tracker->GetMaxPlanesPerTrack(iTrk); ++ipl)
            logStream << " + " << tracker->GetMaxHitsPerTrack(iTrk, ipl);
          logStream << ")";
        } // for tracks
      } // for trackers
    } // if logging enabled

    // if we don't use a permanent buffer (which can be huge),
    // delete the current buffer, and we'll create a new one on the next event

    if (!fUseBuffer) {
      LOG_DEBUG("AnaRootParserStructure") << "Freeing the tree data structure";
      DestroyData();
    }
    } // dune::AnaRootParser::analyze()


    void dune::AnaRootParser::FillShower( AnaRootParserDataStruct::ShowerDataStruct& showerData, size_t iShower,
        recob::Shower const& shower, const bool fSavePFParticleInfo,
        const std::map<Short_t, Short_t> &showerIDtoPFParticleIDMap
        ) const {

      showerData.showerID[iShower]        = shower.ID();
      showerData.shwr_bestplane[iShower]  = shower.best_plane();
      showerData.shwr_length[iShower]     = shower.Length();

      TVector3 const& dir_start = shower.Direction();
      showerData.shwr_startdcosx[iShower] = dir_start.X();
      showerData.shwr_startdcosy[iShower] = dir_start.Y();
      showerData.shwr_startdcosz[iShower] = dir_start.Z();

      TVector3 const& pos_start = shower.ShowerStart();
      showerData.shwr_startx[iShower]     = pos_start.X();
      showerData.shwr_starty[iShower]     = pos_start.Y();
      showerData.shwr_startz[iShower]     = pos_start.Z();

      if (fSavePFParticleInfo) {
        auto mapIter = showerIDtoPFParticleIDMap.find(shower.ID());
        if (mapIter != showerIDtoPFParticleIDMap.end()) {
          // This vertex has a corresponding PFParticle.
          showerData.shwr_hasPFParticle[iShower] = 1;
          showerData.shwr_PFParticleID[iShower] = mapIter->second;
        }
        else
          showerData.shwr_hasPFParticle[iShower] = 0;
      }

      if (shower.Energy().size() == kNplanes)
        std::copy_n
          (shower.Energy().begin(),    kNplanes, &showerData.shwr_totEng[iShower][0]);
      if (shower.dEdx().size() == kNplanes)
        std::copy_n
          (shower.dEdx().begin(),      kNplanes, &showerData.shwr_dedx[iShower][0]);
      if (shower.MIPEnergy().size() == kNplanes)
        std::copy_n
          (shower.MIPEnergy().begin(), kNplanes, &showerData.shwr_mipEng[iShower][0]);

    } // dune::AnaRootParser::FillShower()


    void dune::AnaRootParser::FillShowers( AnaRootParserDataStruct::ShowerDataStruct& showerData,
        std::vector<recob::Shower> const& showers, const bool fSavePFParticleInfo,
        const std::map<Short_t, Short_t> &showerIDtoPFParticleIDMap
        ) const {

      const size_t NShowers = showers.size();

      //
      // prepare the data structures, the tree and the connection between them
      //

      // allocate enough space for this number of showers
      // (but at least for one of them!)
      showerData.SetMaxShowers(std::max(NShowers, (size_t) 1));
      showerData.Clear(); // clear all the data

      // now set the tree addresses to the newly allocated memory;
      // this creates the tree branches in case they are not there yet
      showerData.SetAddresses(fTree);
      if (NShowers > showerData.GetMaxShowers()) {
        // got this error? it might be a bug,
        // since we are supposed to have allocated enough space to fit all showers
        mf::LogError("AnaRootParser:limits") << "event has " << NShowers
          << " " << showerData.Name() << " showers, only "
          << showerData.GetMaxShowers() << " stored in tree";
      }

      //
      // now set the data
      //
      // set the record of the number of showers
      // (data structures are already properly resized)
      showerData.nshowers = (Short_t) NShowers;

      // set all the showers one by one
      for (size_t i = 0; i < NShowers; ++i) FillShower(showerData, i, showers[i], fSavePFParticleInfo, showerIDtoPFParticleIDMap);

    } // dune::AnaRootParser::FillShowers()



    void dune::AnaRootParser::HitsBackTrack( art::Ptr<recob::Hit> const& hit, int & trackid, float & maxe, float & purity){

      trackid = -1;
      purity = -1;

      art::ServiceHandle<cheat::BackTrackerService> bt_serv;

      std::map<int,double> trkide;
      std::vector<sim::IDE> ides;

      //bt_serv->HitToSimIDEs(hit,ides);
      std::vector<sim::TrackIDE> eveIDs = bt_serv->HitToEveTrackIDEs(hit);

      for(size_t e = 0; e < eveIDs.size(); ++e){
          //std::cout<<" "<<e<<" "<<eveIDs[e].trackID<<" "<<eveIDs[e].energy<<" "<<eveIDs[e].energyFrac<<std::endl;
          trkide[eveIDs[e].trackID] += eveIDs[e].energy;
      }

      maxe = 0;
      double tote = 0;
      for (std::map<int,double>::iterator ii = trkide.begin(); ii!=trkide.end(); ++ii){
        tote += ii->second;
        if ((ii->second)>maxe){
          maxe = ii->second;
          trackid = ii->first;
        }
      }

      //std::cout << "the total energy of this reco track is: " << tote << std::endl;
      if (tote>0){
        purity = maxe/tote;
      }
    }

    // Calculate distance to boundary.
    double dune::AnaRootParser::bdist(const TVector3& pos)
    {
      double d1 = -ActiveBounds[0] + pos.X();
      double d2 =  ActiveBounds[1] - pos.X();
      double d3 = -ActiveBounds[2] + pos.Y();
      double d4 =  ActiveBounds[3] - pos.Y();
      double d5 = -ActiveBounds[4] + pos.Z();
      double d6 =  ActiveBounds[5] - pos.Z();

      double Xmin = std::min(d1, d2);
      double result = 0;
      if (Xmin<0) result = std::min(std::min(std::min( d3, d4), d5), d6);         // - FIXME Passing uncorrected hits means X positions are very wrong ( outside of ActiveVolume )
      else result = std::min(std::min(std::min(std::min(Xmin, d3), d4), d5), d6);
      if (result<-1) result = -1;
      /*
         std::cout << "\n" << std::endl;
         std::cout << "-"<<ActiveBounds[0]<<" + "<<pos.X()<< " = " << d1 << "\n"
         <<      ActiveBounds[1]<<" - "<<pos.X()<< " = " << d2 << "\n"
         << "-"<<ActiveBounds[2]<<" + "<<pos.Y()<< " = " << d3 << "\n"
         <<      ActiveBounds[3]<<" - "<<pos.Y()<< " = " << d4 << "\n"
         << "-"<<ActiveBounds[4]<<" + "<<pos.Z()<< " = " << d5 << "\n"
         <<      ActiveBounds[5]<<" - "<<pos.Z()<< " = " << d6 << "\n"
         << "And the Minimum is " << result << std::endl;
         */
      return result;
    }

    // Length of reconstructed track, trajectory by trajectory.
    double dune::AnaRootParser::length(const recob::Track& track)
    {
      double result = 0.;
      TVector3 disp = track.LocationAtPoint(0);
      int n = track.NumberTrajectoryPoints();

      for(int i = 1; i < n; ++i) {
        const TVector3& pos = track.LocationAtPoint(i);
        disp -= pos;
        result += disp.Mag();
        disp = pos;
      }
      return result;
    }


    double dune::AnaRootParser::driftedLength(const sim::MCTrack& mctrack, TLorentzVector& tpcstart, TLorentzVector& tpcend, TLorentzVector& tpcmom){
      auto const* geom = lar::providerFrom<geo::Geometry>();
      auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();

      //compute the drift x range
      double vDrift = detprop->DriftVelocity()*1e-3; //cm/ns
      double xrange[2] = {DBL_MAX, -DBL_MAX };
      for (unsigned int c=0; c<geom->Ncryostats(); ++c) {
        for (unsigned int t=0; t<geom->NTPC(c); ++t) {
          double Xat0 = detprop->ConvertTicksToX(0,0,t,c);
          double XatT = detprop->ConvertTicksToX(detprop->NumberTimeSamples(),0,t,c);
          xrange[0] = std::min(std::min(Xat0, XatT), xrange[0]);
          xrange[1] = std::max(std::max(Xat0, XatT), xrange[1]);
        }
      }

      double result = 0.;
      TVector3 disp;
      bool first = true;

      for(auto step: mctrack) {
        // check if the particle is inside a TPC
        if (step.X() >= ActiveBounds[0] && step.X() <= ActiveBounds[1] &&
            step.Y() >= ActiveBounds[2] && step.Y() <= ActiveBounds[3] &&
            step.Z() >= ActiveBounds[4] && step.Z() <= ActiveBounds[5] ){
          // Doing some manual shifting to account for
          // an interaction not occuring with the beam dump
          // we will reconstruct an x distance different from
          // where the particle actually passed to to the time
          // being different from in-spill interactions
          double newX = step.X()+(step.T()*vDrift);
          if (newX < xrange[0] || newX > xrange[1]) continue;

          TLorentzVector pos(newX,step.Y(),step.Z(),step.T());
          if(first){
            tpcstart = pos;
            tpcmom = step.Momentum();
            first = false;
          }
          else {
            disp -= pos.Vect();
            result += disp.Mag();
          }
          disp = pos.Vect();
          tpcend = pos;
        }
      }
      return result;
    }

    // Length of MC particle, trajectory by trajectory (with the manual shifting for x correction)
    double dune::AnaRootParser::driftedLength(const simb::MCParticle& p, TLorentzVector& start, TLorentzVector& end, unsigned int &starti, unsigned int &endi)
    {
      auto const* geom = lar::providerFrom<geo::Geometry>();
      auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();

      //compute the drift x range
      double vDrift = detprop->DriftVelocity()*1e-3; //cm/ns
      double xrange[2] = {DBL_MAX, -DBL_MAX };
      for (unsigned int c=0; c<geom->Ncryostats(); ++c) {
        for (unsigned int t=0; t<geom->NTPC(c); ++t) {
          double Xat0 = detprop->ConvertTicksToX(0,0,t,c);
          double XatT = detprop->ConvertTicksToX(detprop->NumberTimeSamples(),0,t,c);
          xrange[0] = std::min(std::min(Xat0, XatT), xrange[0]);
          xrange[1] = std::max(std::max(Xat0, XatT), xrange[1]);
        }
      }

      double result = 0.;
      TVector3 disp;
      bool first = true;

      for(unsigned int i = 0; i < p.NumberTrajectoryPoints(); ++i) {
        // check if the particle is inside a TPC
        if (p.Vx(i) >= ActiveBounds[0] && p.Vx(i) <= ActiveBounds[1] &&
            p.Vy(i) >= ActiveBounds[2] && p.Vy(i) <= ActiveBounds[3] &&
            p.Vz(i) >= ActiveBounds[4] && p.Vz(i) <= ActiveBounds[5]){
          // Doing some manual shifting to account for
          // an interaction not occuring with the beam dump
          // we will reconstruct an x distance different from
          // where the particle actually passed to to the time
          // being different from in-spill interactions
          double newX = p.Vx(i)+(p.T(i)*vDrift);
          if (newX < xrange[0] || newX > xrange[1]) continue;
          TLorentzVector pos(newX,p.Vy(i),p.Vz(i),p.T());
          if(first){
            start = pos;
            starti=i;
            first = false;
          }
          else {
            disp -= pos.Vect();
            result += disp.Mag();
          }
          disp = pos.Vect();
          end = pos;
          endi = i;
        }
      }
      return result;
    }

    // Length of MC particle, trajectory by trajectory (with out the manual shifting for x correction)
    double dune::AnaRootParser::length(const simb::MCParticle& p, TLorentzVector& start, TLorentzVector& end, unsigned int &starti, unsigned int &endi)
    {
      double result = 0.;
      TVector3 disp;
      bool first = true;

      for(unsigned int i = 0; i < p.NumberTrajectoryPoints(); ++i) {
        // check if the particle is inside a TPC
        if (p.Vx(i) >= ActiveBounds[0] && p.Vx(i) <= ActiveBounds[1] && p.Vy(i) >= ActiveBounds[2] && p.Vy(i) <= ActiveBounds[3] && p.Vz(i) >= ActiveBounds[4] && p.Vz(i) <= ActiveBounds[5]){
          if(first){
            start = p.Position(i);
            first = false;
            starti = i;
          }else{
            disp -= p.Position(i).Vect();
            result += disp.Mag();
          }
          disp = p.Position(i).Vect();
          end = p.Position(i);
          endi = i;
        }
      }
      return result;
    }

    //......................................................................
    int dune::AnaRootParser::CountHits(const art::Event&    evt,
                              	       const art::InputTag& which,
                                       unsigned int         cryostat,
                              	       unsigned int         tpc,
                              	       unsigned int         plane)
    {
      std::vector<const recob::Hit*> temp;
      int NumberOfHitsBeforeThisPlane=0;
      evt.getView(which, temp);   //temp.size() = total number of hits for this event (number of all hits in all Cryostats, TPC's, planes and wires)
      for(size_t t = 0; t < temp.size(); ++t){
	if( temp[t]->WireID().Cryostat == cryostat&& temp[t]->WireID().TPC == tpc && temp[t]->WireID().Plane == plane ) break;
	NumberOfHitsBeforeThisPlane++;
      }
      return NumberOfHitsBeforeThisPlane;
    }

    namespace dune{

      DEFINE_ART_MODULE(AnaRootParser)

    }

    namespace lar_pandora{class LArPandoraHelper;}
