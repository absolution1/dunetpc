////////////////////////////////////////////////////////////////////////
// \file DetectorPropertiesStandard.h
//
// \brief service to contain information about detector electronics, etc
//
// \author brebel@fnal.gov
// 
// Separation of service from Detector info class:
// jpaley@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef DETINFO_DETECTORPROPERTIES_PROTODUNESP_H
#define DETINFO_DETECTORPROPERTIES_PROTODUNESP_H
// LArSoft libraries

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Principal/Run.h"


// #include "canvas/Persistency/Utilities/TypeID.h"


#include "art/Framework/Principal/Event.h"
#include "art/Persistency/Provenance/ScheduleContext.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"

#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/CoreUtils/ProviderPack.h"
#include "lardataalg/DetectorInfo/LArProperties.h"
#include "lardataalg/DetectorInfo/DetectorClocks.h"
#include "lardataalg/DetectorInfo/DetectorProperties.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "lardata/DetectorInfoServices/ServicePack.h"
// framework libraries
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
// C/C++ standard libraries
#include <set>







#include "dune/Protodune/singlephase/DetectorServices/Providers/DetectorPropertiesProtoDUNEsp.h"





///General LArSoft Utilities
namespace spdp{
  
  class DetectorPropertiesProtoDUNEsp : public detinfo::DetectorProperties {
    public:
      /// List of service providers we depend on
      using providers_type = lar::ProviderPack<
        geo::GeometryCore,
        detinfo::LArProperties,
        detinfo::DetectorClocks
        >;
       
      /// Structure for configuration parameters
      struct Configuration_t {
        using Name = fhicl::Name;
        using Comment = fhicl::Comment;
        
        fhicl::Sequence<double> Efield { Name("Efield"), Comment(
          "electric field in front of each wire plane (the last one is the big one!) [kV/cm]")
          };
    

       fhicl::Atom<bool        > fGetHVDriftfromMetaData {
          Name("GetHVDriftfromSamweb"),
          Comment("option to get HV drift field from MetaData")
        };
          fhicl::Atom<bool        > fGetReadOutWindowSizefromMetaData{
          Name("GetReadOutWindowSizefromSamweb"),
          Comment("option to get ReadoutWindowSize and NumberTimeSamples from MetaData")
        };

        fhicl::Atom<double      > Electronlifetime         {
          Name("Electronlifetime"        ),
          Comment("electron lifetime in liquid argon [us]")
        };
        fhicl::Atom<double      > Temperature              {
          Name("Temperature"             ),
          Comment("argon temperature [K]")
        };
        fhicl::Atom<double      > ElectronsToADC           {
          Name("ElectronsToADC"          ),
          Comment("conversion factor: (ADC counts)/(ionization electrons)")
        };
        fhicl::Atom<unsigned int> NumberTimeSamples        {
          Name("NumberTimeSamples"       ),
          Comment("number of TPC readout TDC clock ticks per event")
        };
        fhicl::Atom<unsigned int> ReadOutWindowSize        {
          Name("ReadOutWindowSize"       ),
          Comment("number of TPC readout TDC clock ticks per readout window")
        };
        
        // The following are not really "optional": the ones for the views which
        // are present are mandatory.
        fhicl::OptionalAtom<double      > TimeOffsetU              {
          Name("TimeOffsetU"             ),
          Comment("tick offset subtracted to to convert spacepoint coordinates to hit times on view U")
        };
        fhicl::OptionalAtom<double      > TimeOffsetV              {
          Name("TimeOffsetV"             ),
          Comment("tick offset subtracted to to convert spacepoint coordinates to hit times on view V")
        };
        fhicl::OptionalAtom<double      > TimeOffsetZ              {
          Name("TimeOffsetZ"             ),
          Comment("tick offset subtracted to to convert spacepoint coordinates to hit times on view Z")
        };
        fhicl::OptionalAtom<double      > TimeOffsetY              {
          Name("TimeOffsetY"             ),
          Comment("tick offset subtracted to to convert spacepoint coordinates to hit times on view Y")
        };
        fhicl::OptionalAtom<double      > TimeOffsetX              {
          Name("TimeOffsetX"             ),
          Comment("tick offset subtracted to to convert spacepoint coordinates to hit times on view X")
        };
        
        fhicl::Atom<double      > SternheimerA             {
          Name("SternheimerA"),
          Comment("parameter a of Sternheimer correction delta = 2log(10) x - cbar + { a (x1-x)^k } theta(x1-x), x = log10(p/m)")
        };
        fhicl::Atom<double      > SternheimerK             {
          Name("SternheimerK"),
          Comment("parameter k of Sternheimer correction delta = 2log(10) x - cbar + { a (x_1-x)^k } theta(x1-x), x = log10(p/m)")
        };
        fhicl::Atom<double      > SternheimerX0            {
          Name("SternheimerX0"),
          Comment("minimum x = log10(p/m) for the application of Sternheimer correction")
        };
        fhicl::Atom<double      > SternheimerX1            {
          Name("SternheimerX1"),
          Comment("parameter x_1 of Sternheimer correction delta = 2log(10) x - cbar + { a (x_1-x)^k } theta(x1-x), x = log10(p/m)")
        };
        fhicl::Atom<double      > SternheimerCbar             {
          Name("SternheimerCbar"),
          Comment("parameter cbar of Sternheimer correction delta = 2log(10) x - cbar + { a (x_1-x)^k } theta(x1-x), x = log10(p/m)")
        };
        fhicl::Atom<bool> SimpleBoundary { Name("SimpleBoundaryProcess" ), Comment("") };
      
      }; // Configuration_t
 
      DetectorPropertiesProtoDUNEsp();
      DetectorPropertiesProtoDUNEsp(fhicl::ParameterSet const& pset, 
                         const geo::GeometryCore* geo,
                         const detinfo::LArProperties* lp,
                         const detinfo::DetectorClocks* c,
                         std::set<std::string> const& ignore_params = {}
                         );
      /**
       * @brief Constructs the provider and sets up the dependencies
       * @param pset FHiCL parameter set for provider configuration
       * @param providers pack of providers `DetectorPropertiesStandard` depends
       *        on
       * @param ignore_params unknown configuration keys in `pset` to be
       *        tolerated
       * @see Setup()
       */
      DetectorPropertiesProtoDUNEsp(fhicl::ParameterSet const& pset,
                         providers_type providers,
                         std::set<std::string> const& ignore_params = {});
      DetectorPropertiesProtoDUNEsp(DetectorPropertiesProtoDUNEsp const&) = delete;
      virtual ~DetectorPropertiesProtoDUNEsp() = default;
      
      /**
       * @brief Configures the provider, first validating the configuration
       * @param p configuration parameter set
       * @param ignore_params parameters to be ignored (optional)
       * 
       * This method will validate the parameter set (except for the parameters
       * it's explicitly told to ignore) and extract the useful information out
       * of it.
       */
      void ValidateAndConfigure(
        fhicl::ParameterSet const& p,
        std::set<std::string> const& ignore_params = {}
        );
      
      /// Extracts the relevant configuration from the specified object
      void Configure(Configuration_t const& config);
      
      /**
       * @brief Validates the specified configuration
       * @param p configuration parameter set
       * @param ignore_params parameters to be ignored (optional)
       * @return a parsed configuration object
       * @see ValidateAndConfigure(), Configure()
       * 
       * This method will validate the parameter set (except for the parameters
       * it's explicitly told to ignore) and it returns an object ready to
       * be used with Configure().
       */
      Configuration_t ValidateConfiguration(
        fhicl::ParameterSet const& p,
        std::set<std::string> const& ignore_params = {}
        );
      bool Update(uint64_t ts); 
      bool UpdateHV(std::string filename);
      bool UpdateReadoutWindowSize(std::string filename);
      bool UpdateTemp(int run);
      bool UpdateClocks(const detinfo::DetectorClocks* clks);
      
      /**
       * @brief Sets all the providers at once
       * @param providers the pack of service providers we depend on
       * 
       * Example:
       *     
       *     lar::DetectorPropertiesStandard::providers_type providers;
       *     providers.set(lar::providerFrom<geo::Geometry>());
       *     providers.set(lar::providerFrom<detinfo::LArPropertiesService>());
       *     providers.set(lar::providerFrom<detinfo::DetectorClocksService>());
       *     detprop->Setup(providers);
       *
       */
      void Setup(providers_type providers);
        
      void SetGeometry(const geo::GeometryCore* g) { fGeo = g; }
      void SetLArProperties(const detinfo::LArProperties* lp) { fLP = lp; }
      void SetDetectorClocks(const detinfo::DetectorClocks* clks) { fClocks = clks; }
      void SetNumberTimeSamples(unsigned int nsamp) { fNumberTimeSamples=nsamp;}
      // Accessors.
      virtual double Efield(unsigned int planegap=0) const override; ///< kV/cm
      virtual double DriftVelocity(double efield=0., double temperature=0.) const override;  ///< cm/us
      
      /// dQ/dX in electrons/cm, returns dE/dX in MeV/cm.
      virtual double BirksCorrection(double dQdX) const override;
      virtual double ModBoxCorrection(double dQdX) const override;
      virtual double ElectronLifetime() const override { return fElectronlifetime; }   //< microseconds
      
      
      /**
       * @brief Returns argon density at a given temperature
       * @param temperature the temperature in kelvin
       * @return argon density in g/cm^3
       * 
       * Density is nearly a linear function of temperature.
       * See the NIST tables for details
       * Slope is between -6.2 and -6.1, intercept is 1928 kg/m^3.
       * This parameterization will be good to better than 0.5%.
       */
      virtual double Density(double temperature) const override;                          ///< g/cm^3
      
      // need to provide a definition, since the override above hides the inherited one
      virtual double Density() const override { return Density(Temperature()); }
      
      /// In kelvin.
      virtual double Temperature()                   const override { return fTemperature; }
      
      /**
       * @brief Restricted mean energy loss (dE/dx)
       * @param mom  momentum of incident particle [GeV/c]
       * @param mass mass of incident particle [GeV/c^2]
       * @param tcut maximum kinetic energy of delta rays [MeV]; 0 for unlimited
       * @return the restricted mean energy loss (dE/dx) in units of MeV/cm
       *
       * Returned value is always positive.
       * For unrestricted mean energy loss, set tcut = 0 (special case),
       * or tcut large.
       * 
       * Based on Bethe-Bloch formula as contained in particle data book.
       * Material parameters are from the configuration.
       */
      virtual double Eloss(double mom, double mass, double tcut) const override;
      
      /**
       * @brief Energy loss fluctuation (@f$ \sigma_{E}^2 / x @f$)
       * @param mom  momentum of incident particle in [GeV/c]
       * @param mass mass of incident particle [GeV/c^2]
       * @return energy loss fluctuation in MeV^2/cm
       *
       * Based on Bichsel formula referred to but not given in PDG.
       */
      virtual double ElossVar(double mom, double mass) const override;
      virtual double       SamplingRate()      const override { return fTPCClock.TickPeriod() * 1.e3; }
      virtual double       ElectronsToADC()    const override { return fElectronsToADC; }
      virtual unsigned int NumberTimeSamples() const override { return fNumberTimeSamples; }
      virtual unsigned int ReadOutWindowSize() const override { return fReadOutWindowSize; }
      virtual int          TriggerOffset()     const override;
      virtual double       TimeOffsetU()       const override{ return fTimeOffsetU; };
      virtual double       TimeOffsetV()       const override { return fTimeOffsetV; };
      virtual double       TimeOffsetZ()       const override{ return fTimeOffsetZ; };
      virtual double       TimeOffsetY()       const override{ return fTimeOffsetY; };
      virtual double       ConvertXToTicks(double X, int p, int t, int c) const override;
      virtual double       ConvertXToTicks(double X, geo::PlaneID const& planeid) const override
        { return ConvertXToTicks(X, planeid.Plane, planeid.TPC, planeid.Cryostat); }
      virtual double       ConvertTicksToX(double ticks, int p, int t, int c) const override;
      virtual double       ConvertTicksToX(double ticks, geo::PlaneID const& planeid) const override
        { return ConvertTicksToX(ticks, planeid.Plane, planeid.TPC, planeid.Cryostat); }
      virtual double       GetXTicksOffset(int p, int t, int c) const override;
      virtual double       GetXTicksOffset(geo::PlaneID const& planeid) const override
        { return GetXTicksOffset(planeid.Plane, planeid.TPC, planeid.Cryostat); }
      virtual double       GetXTicksCoefficient(int t, int c) const override;
      virtual double       GetXTicksCoefficient(geo::TPCID const& tpcid) const override
        { return GetXTicksCoefficient(tpcid.TPC, tpcid.Cryostat); }
      virtual double       GetXTicksCoefficient() const override;
      // The following methods convert between TDC counts (SimChannel time) and
      // ticks (RawDigit/Wire time).
      virtual double       ConvertTDCToTicks(double tdc) const override;
      virtual double       ConvertTicksToTDC(double ticks) const override;
      
      virtual bool SimpleBoundary()     const override  { return fSimpleBoundary; }
      
      /// Verifies that the provider is in a fully configured status
      /// @throw cet::exception (category DetectorPropertiesStandard) if not ok
      void CheckIfConfigured() const;
      
    protected:
      
     
      /// Parameters for Sternheimer density effect corrections
      struct SternheimerParameters_t {
        double a;               ///< parameter a
        double k;               ///< parameter k
        double x0;              ///< parameter x0
        double x1;              ///< parameter x1
        double cbar;            ///< parameter Cbar
      }; //  SternheimerParameters_t
      
      void         CalculateXTicksParams();
      
      // service providers we depend on;
      // in principle could be replaced by a single providerpacl_type.
      const detinfo::LArProperties* fLP;
      const detinfo::DetectorClocks* fClocks;
      const geo::GeometryCore* fGeo;
      


      bool                        fGetHVDriftfromMetaData;
      bool                        fGetReadOutWindowSizefromMetaData;
      double                         fHV_cath;   //  <KV
      std::vector<double>          fEfield;           ///< kV/cm (per inter-plane volume)
      double                         fElectronlifetime; ///< microseconds
      double                         fTemperature;      ///< kelvin
      double       fSamplingRate;      ///< in ns
      double            fElectronsToADC;    ///< conversion factor for # of ionization electrons to 1 ADC count
      unsigned int fNumberTimeSamples; ///< number of clock ticks per event
      unsigned int fReadOutWindowSize; ///< number of clock ticks per readout window
      double       fTimeOffsetU;       ///< time offset to convert spacepoint coordinates to hit times on view U
      double       fTimeOffsetV;       ///< time offset to convert spacepoint coordinates to hit times on view V
      double       fTimeOffsetZ;       ///< time offset to convert spacepoint coordinates to hit times on view Z
      double       fTimeOffsetY;       ///< time offset to convert spacepoint coordinates to hit times on view Y
      double       fTimeOffsetX;       ///< time offset to convert spacepoint coordinates to hit times on view X
      double       fHasTimeOffsetU = false; ///< whether time offset was configured for view U
      double       fHasTimeOffsetV = false; ///< whether time offset was configured for view V
      double       fHasTimeOffsetZ = false; ///< whether time offset was configured for view Z
      double       fHasTimeOffsetY = false; ///< whether time offset was configured for view Y
      double       fHasTimeOffsetX = false; ///< whether time offset was configured for view X
      
      SternheimerParameters_t fSternheimerParameters; ///< Sternheimer parameters
      
      double       fXTicksCoefficient; ///< Parameters for x<-->ticks
      std::vector<std::vector<std::vector<double> > > fXTicksOffsets;
      std::vector<std::vector<double> >               fDriftDirection;
      ::detinfo::ElecClock fTPCClock;     ///< TPC electronics clock
      
      bool fSimpleBoundary;
      /// Checks the configuration of time offsets.
      std::string CheckTimeOffsetConfigurationAfterSetup() const;
      
      /// Checks that provider configuration is complete, using setup
      /// information.
      void CheckConfigurationAfterSetup() const;
      
      /// Time-independent implementation of clock updates.
      void DoUpdateClocks();
      
    }; // class DetectorPropertiesStandard
} //namespace detinfo
#endif // DETINFO_DETECTOR_PROPERTIES_H