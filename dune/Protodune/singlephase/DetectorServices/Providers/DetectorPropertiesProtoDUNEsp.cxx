/**
 * @file lardataalg/DetectorInfo/DetectorPropertiesStandard.cxx
 * @brief Separation of service from Detector info class.
 * @author Jonathan Paley (jpaley@fnal.gov)
 */
// Framework includes
#include <cassert>
// LArSoft includes
#include "DetectorPropertiesProtoDUNEsp.h"
#include "larcorealg/CoreUtils/ProviderUtil.h" // lar::IgnorableProviderConfigKeys()
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
// Art includes
#include "fhiclcpp/make_ParameterSet.h"
// C/C++ libraries
#include <sstream> // std::ostringstream
namespace {
  
  template <typename T>
  inline T sqr(T v) { return v*v; }
  
} // local namespace
namespace detinfo{
  //--------------------------------------------------------------------
  DetectorPropertiesStandard::DetectorPropertiesStandard() :
    fLP(0), fClocks(0), fGeo(0)
  {
  }
  
  //--------------------------------------------------------------------
  DetectorPropertiesStandard::DetectorPropertiesStandard(fhicl::ParameterSet const& pset,
                                         const geo::GeometryCore* geo,
                                         const detinfo::LArProperties* lp,
                                         const detinfo::DetectorClocks* c,
                                         std::set<std::string> const& ignore_params /* = {} */
                                         ):
    fLP(lp), fClocks(c), fGeo(geo)
  {
    {
       mf::LogInfo debug("setupProvider<DetectorPropertiesStandard>");
       
       debug << "Asked to ignore " << ignore_params.size() << " keys:";
       for (auto const& key: ignore_params) debug << " '" << key << "'";
    }
    
    ValidateAndConfigure(pset, ignore_params);
    
    fTPCClock = fClocks->TPCClock();
    DoUpdateClocks();
  }
    
  //--------------------------------------------------------------------
  DetectorPropertiesStandard::DetectorPropertiesStandard(fhicl::ParameterSet const& pset,
                                         providers_type providers,
                                         std::set<std::string> const& ignore_params /* = {} */
                                         ):
    DetectorPropertiesStandard(pset,
      providers.get<geo::GeometryCore>(),
      providers.get<detinfo::LArProperties>(),
      providers.get<detinfo::DetectorClocks>(),
      ignore_params
      )
    {}
  
  //--------------------------------------------------------------------
  bool DetectorPropertiesStandard::Update(uint64_t) 
  {
    DoUpdateClocks();
    return true;
  }
  //--------------------------------------------------------------------
  bool DetectorPropertiesStandard::UpdateClocks(const detinfo::DetectorClocks* clks) 
  {
    fClocks = clks;
    
    fTPCClock = fClocks->TPCClock();
    DoUpdateClocks();
    return true;
  }
  
  //------------------------------------------------------------
  double DetectorPropertiesStandard::ConvertTDCToTicks(double tdc) const
  {
    return fClocks->TPCTDC2Tick(tdc);
  }
  
  //--------------------------------------------------------------
  double DetectorPropertiesStandard::ConvertTicksToTDC(double ticks) const
  {
    return fClocks->TPCTick2TDC(ticks);
  }
  
 
  //--------------------------------------------------------------------
  void DetectorPropertiesStandard::Configure(Configuration_t const& config) {
    
    fEfield                     = config.Efield();
    fElectronlifetime           = config.Electronlifetime();
    fTemperature                = config.Temperature();
    fElectronsToADC             = config.ElectronsToADC();
    fNumberTimeSamples          = config.NumberTimeSamples();
    fReadOutWindowSize          = config.ReadOutWindowSize();
    fHasTimeOffsetU = config.TimeOffsetU(fTimeOffsetU);
    fHasTimeOffsetV = config.TimeOffsetV(fTimeOffsetV);
    fHasTimeOffsetZ = config.TimeOffsetZ(fTimeOffsetZ);
    fHasTimeOffsetY = config.TimeOffsetY(fTimeOffsetY);
    fHasTimeOffsetX = config.TimeOffsetX(fTimeOffsetX);
    
    fSternheimerParameters.a    = config.SternheimerA();
    fSternheimerParameters.k    = config.SternheimerK();
    fSternheimerParameters.x0   = config.SternheimerX0();
    fSternheimerParameters.x1   = config.SternheimerX1();
    fSternheimerParameters.cbar = config.SternheimerCbar();
    fSimpleBoundary = config.SimpleBoundary();
    DoUpdateClocks();
    
  } // DetectorPropertiesStandard::Configure()
  
  //--------------------------------------------------------------------
  DetectorPropertiesStandard::Configuration_t
  DetectorPropertiesStandard::ValidateConfiguration(
    fhicl::ParameterSet const& p,
    std::set<std::string> const& ignore_params /* = {} */
  ) {
    std::set<std::string> ignorable_keys = lar::IgnorableProviderConfigKeys();
    ignorable_keys.insert(ignore_params.begin(), ignore_params.end());
    
    // parses and validates the parameter set:
    fhicl::Table<Configuration_t> config_table { p, ignorable_keys };
    
    return std::move(config_table());
    
  } // DetectorPropertiesStandard::ValidateConfiguration()
  
  //--------------------------------------------------------------------
  void DetectorPropertiesStandard::ValidateAndConfigure(
    fhicl::ParameterSet const& p,
    std::set<std::string> const& ignore_params /* = {} */
  ) {
    Configure(ValidateConfiguration(p, ignore_params));
  } // ValidateAndConfigure()
  
  
  //------------------------------------------------------------------------------------//
  void DetectorPropertiesStandard::Setup(providers_type providers) {
    
    SetGeometry(providers.get<geo::GeometryCore>());
    SetLArProperties(providers.get<detinfo::LArProperties>());
    SetDetectorClocks(providers.get<detinfo::DetectorClocks>());
    
    CheckConfigurationAfterSetup();
    
  } // DetectorPropertiesStandard::Setup()
  
  
  //------------------------------------------------------------------------------------//
  double DetectorPropertiesStandard::Efield(unsigned int planegap) const
  {
    if(planegap >= fEfield.size())
      throw cet::exception("DetectorPropertiesStandard") << "requesting Electric field in a plane gap that is not defined\n";
    
    return fEfield[planegap];
  }
  
  
  //------------------------------------------------
  double DetectorPropertiesStandard::Density(double temperature) const
  {
    // Default temperature use internal value.
    if(temperature == 0.)
      temperature = Temperature();
  
    double density = -0.00615*temperature + 1.928;
  
    return density;
  } // DetectorPropertiesStandard::Density()
  
  
  //----------------------------------------------------------------------------------
  // Restricted mean energy loss (dE/dx) in units of MeV/cm.
  //
  // For unrestricted mean energy loss, set tcut = 0, or tcut large.
  //
  // Arguments:
  //
  // mom  - Momentum of incident particle in GeV/c.
  // mass - Mass of incident particle in GeV/c^2.
  // tcut - Maximum kinetic energy of delta rays (MeV).
  //
  // Returned value is positive.
  //
  // Based on Bethe-Bloch formula as contained in particle data book.
  // Material parameters (stored in larproperties.fcl) are taken from
  // pdg web site http://pdg.lbl.gov/AtomicNuclearProperties/
  //
  double DetectorPropertiesStandard::Eloss(double mom, double mass, double tcut) const
  {
    // Some constants.
  
    double K = 0.307075;     // 4 pi N_A r_e^2 m_e c^2 (MeV cm^2/mol).
    double me = 0.510998918; // Electron mass (MeV/c^2).
  
    // Calculate kinematic quantities.
  
    double bg = mom / mass;           // beta*gamma.
    double gamma = sqrt(1. + bg*bg);  // gamma.
    double beta = bg / gamma;         // beta (velocity).
    double mer = 0.001 * me / mass;   // electron mass / mass of incident particle.
    double tmax = 2.*me* bg*bg / (1. + 2.*gamma*mer + mer*mer);  // Maximum delta ray energy (MeV).
  
    // Make sure tcut does not exceed tmax.
  
    if(tcut == 0. || tcut > tmax)
      tcut = tmax;
  
    // Calculate density effect correction (delta).
  
    double x = std::log10(bg);
    double delta = 0.;
    if(x >= fSternheimerParameters.x0) {
      delta = 2. * std::log(10.) * x - fSternheimerParameters.cbar;
      if(x < fSternheimerParameters.x1)
        delta += fSternheimerParameters.a * std::pow(fSternheimerParameters.x1 - x, fSternheimerParameters.k);
    }
  
    // Calculate stopping number.
  
    double B = 0.5 * std::log(2.*me*bg*bg*tcut / (1.e-12 * sqr(fLP->ExcitationEnergy())))
      - 0.5 * beta*beta * (1. + tcut / tmax) - 0.5 * delta;
  
    // Don't let the stopping number become negative.
  
    if(B < 1.)
      B = 1.;
  
    // Calculate dE/dx.
  
    double dedx = Density() * K*fLP->AtomicNumber()*B / (fLP->AtomicMass() * beta*beta);
  
    // Done.
  
    return dedx;
  } // DetectorPropertiesStandard::Eloss()
  
  //----------------------------------------------------------------------------------
  double DetectorPropertiesStandard::ElossVar(double mom, double mass) const
  {
    // Some constants.
  
    double K = 0.307075;     // 4 pi N_A r_e^2 m_e c^2 (MeV cm^2/mol).
    double me = 0.510998918; // Electron mass (MeV/c^2).
  
    // Calculate kinematic quantities.
  
    double bg = mom / mass;          // beta*gamma.
    double gamma2 = 1. + bg*bg;      // gamma^2.
    double beta2 = bg*bg / gamma2;   // beta^2.
  
    // Calculate final result.
  
    double result = gamma2 * (1. - 0.5 * beta2) * me * (fLP->AtomicNumber() / fLP->AtomicMass()) * K * Density();
    return result;
  } // DetectorPropertiesStandard::ElossVar()
  //------------------------------------------------------------------------------------//
  double DetectorPropertiesStandard::DriftVelocity(double efield, double temperature) const {
  // Drift Velocity as a function of Electric Field and LAr Temperature
  // from : W. Walkowiak, NIM A 449 (2000) 288-294
  //
  // Efield should have units of kV/cm
  // Temperature should have units of Kelvin
  // Default Efield, use internal value.
  if(efield == 0.)
    efield = Efield();
  //
  if(efield > 4.0)
    mf::LogWarning("DetectorPropertiesStandard") << "DriftVelocity Warning! : E-field value of "
                                    << efield
                                    << " kV/cm is outside of range covered by drift"
                                    << " velocity parameterization. Returned value"
                                    << " may not be correct";
  // Default temperature use internal value.
  if(temperature == 0.)
    temperature = Temperature();
  if(temperature < 87.0 || temperature > 94.0)
    mf::LogWarning("DetectorPropertiesStandard") << "DriftVelocity Warning! : Temperature value of "
                                    << temperature
                                    << " K is outside of range covered by drift velocity"
                                    << " parameterization. Returned value may not be"
                                    << " correct";
  double tshift = -87.203+temperature;
  double xFit = 0.0938163-0.0052563*tshift-0.0001470*tshift*tshift;
  double uFit = 5.18406+0.01448*tshift-0.003497*tshift*tshift-0.000516*tshift*tshift*tshift;
  double vd;
// Icarus Parameter Set, use as default
  double  P1 = -0.04640; // K^-1
  double  P2 = 0.01712;  // K^-1
  double  P3 = 1.88125;   // (kV/cm)^-1
  double  P4 =  0.99408;    // kV/cm
  double  P5 =  0.01172;   // (kV/cm)^-P6
  double  P6 =  4.20214;
  double  T0 =  105.749;  // K
      // Walkowiak Parameter Set
  double    P1W = -0.01481; // K^-1
  double  P2W = -0.0075;  // K^-1
  double   P3W =  0.141;   // (kV/cm)^-1
  double   P4W =  12.4;    // kV/cm
  double   P5W =  1.627;   // (kV/cm)^-P6
  double   P6W =  0.317;
  double   T0W =  90.371;  // K
// From Craig Thorne . . . currently not documented
// smooth transition from linear at small fields to 
//     icarus fit at most fields to Walkowiak at very high fields
   if (efield < xFit) vd=efield*uFit;
   else if (efield<0.619) { 
     vd = ((P1*(temperature-T0)+1)
               *(P3*efield*std::log(1+P4/efield) + P5*std::pow(efield,P6))
               +P2*(temperature-T0));
   }
   else if (efield<0.699) {
     vd = 12.5*(efield-0.619)*((P1W*(temperature-T0W)+1)
               *(P3W*efield*std::log(1+P4W/efield) + P5W*std::pow(efield,P6W))
               +P2W*(temperature-T0W))+
       12.5*(0.699-efield)*((P1*(temperature-T0)+1)
               *(P3*efield*std::log(1+P4/efield) + P5*std::pow(efield,P6))
               +P2*(temperature-T0));
   }
   else {
     vd = ((P1W*(temperature-T0W)+1)
               *(P3W*efield*std::log(1+P4W/efield) + P5W*std::pow(efield,P6W))
               +P2W*(temperature-T0W));     
   }
  vd /= 10.;
  return vd; // in cm/us
}
  //----------------------------------------------------------------------------------
  // The below function assumes that the user has applied the lifetime correction and
  // effective pitch between the wires (usually after 3D reconstruction). Using with
  // mean wire pitch will not give correct results.
  // parameters:
  //  dQdX in electrons/cm, charge (amplitude or integral obtained) divided by
  //         effective pitch for a given 3D track.
  // returns dEdX in MeV/cm
  double DetectorPropertiesStandard::BirksCorrection(double dQdx) const
  {
    // Correction for charge quenching using parameterization from
    // S.Amoruso et al., NIM A 523 (2004) 275
    
    double  A3t    = util::kRecombA;
    double  K3t    = util::kRecombk;                     // in KV/cm*(g/cm^2)/MeV
    double  rho    = Density();                    // LAr density in g/cm^3
    double Wion    = 1000./util::kGeVToElectrons;        // 23.6 eV = 1e, Wion in MeV/e
    double E_field  = Efield();                           // Electric Field in the drift region in KV/cm
    K3t           /= rho;                                // KV/MeV
    double dEdx    = dQdx/(A3t/Wion-K3t/E_field*dQdx);    //MeV/cm
    
    return dEdx;
  }  
  
  //----------------------------------------------------------------------------------
  // Modified Box model correction 
  double DetectorPropertiesStandard::ModBoxCorrection(double dQdx) const
  {
    // Modified Box model correction has better behavior than the Birks
    // correction at high values of dQ/dx.
    double  rho    = Density();                    // LAr density in g/cm^3
    double Wion    = 1000./util::kGeVToElectrons;        // 23.6 eV = 1e, Wion in MeV/e
    double E_field  = Efield();                           // Electric Field in the drift region in KV/cm
    double Beta    = util::kModBoxB / (rho * E_field);
    double Alpha   = util::kModBoxA;
    double dEdx = (exp(Beta * Wion * dQdx ) - Alpha) / Beta;
    
    return dEdx;
    
  }
  
  //------------------------------------------------------------------------------------//
  int  DetectorPropertiesStandard::TriggerOffset()     const 
  {
    return fTPCClock.Ticks(fClocks->TriggerOffsetTPC() * -1.);
  }
  
  
  //--------------------------------------------------------------------
  //  x<--> ticks conversion methods 
  //
  //  Ben Jones April 2012, 
  //  based on code by Herb Greenlee in SpacePointService
  //  
  
  
  //--------------------------------------------------------------------
  // Take an X coordinate, and convert to a number of ticks, the
  // charge deposit occured at t=0
 
  double DetectorPropertiesStandard::ConvertXToTicks(double X, int p, int t, int c) const
  {
    return (X / (fXTicksCoefficient * fDriftDirection.at(c).at(t)) +  fXTicksOffsets.at(c).at(t).at(p) );
  }
  
  //-------------------------------------------------------------------
  // Take a cooridnate in ticks, and convert to an x position
  // assuming event deposit occured at t=0
 
  double  DetectorPropertiesStandard::ConvertTicksToX(double ticks, int p, int t, int c) const
  {
    return (ticks - fXTicksOffsets.at(c).at(t).at(p)) * fXTicksCoefficient * fDriftDirection.at(c).at(t);  
  }
  
  //--------------------------------------------------------------------
  void DetectorPropertiesStandard::CheckIfConfigured() const
  {
    if (!fGeo) throw cet::exception(__FUNCTION__) << "Geometry is uninitialized!";
    if (!fLP) throw cet::exception(__FUNCTION__) << "LArPropertiesStandard is uninitialized!";
    if (!fClocks) throw cet::exception(__FUNCTION__) << "DetectorClocks is uninitialized!";
  }
  
  
  //--------------------------------------------------------------------
  // Recalculte x<-->ticks conversion parameters from detector constants
  
  void DetectorPropertiesStandard::CalculateXTicksParams()
  {
    CheckIfConfigured();
    
    double samplingRate   = SamplingRate();
    double efield         = Efield();
    double temperature    = Temperature();
    double driftVelocity  = DriftVelocity(efield, temperature);
    
    fXTicksCoefficient    = 0.001 * driftVelocity * samplingRate;
    double triggerOffset  = TriggerOffset();
    fXTicksOffsets.clear();
    fXTicksOffsets.resize(fGeo->Ncryostats());
    fDriftDirection.clear();
    fDriftDirection.resize(fGeo->Ncryostats());
    for(size_t cstat = 0; cstat < fGeo->Ncryostats(); ++cstat){
      fXTicksOffsets[cstat].resize(fGeo->Cryostat(cstat).NTPC());
      fDriftDirection[cstat].resize(fGeo->Cryostat(cstat).NTPC());
      for(size_t tpc = 0; tpc < fGeo->Cryostat(cstat).NTPC(); ++tpc) {
        const geo::TPCGeo& tpcgeom = fGeo->Cryostat(cstat).TPC(tpc);
        const double dir((tpcgeom.DriftDirection() == geo::kNegX) ? +1.0 :-1.0);
        fDriftDirection[cstat][tpc] = dir;
        int nplane = tpcgeom.Nplanes();
        fXTicksOffsets[cstat][tpc].resize(nplane, 0.);
        for(int plane = 0; plane < nplane; ++plane) {
          const geo::PlaneGeo& pgeom = tpcgeom.Plane(plane);
          
          
          // Get field in gap between planes
          double efieldgap[3];
          double driftVelocitygap[3];
          double fXTicksCoefficientgap[3];
          for (int igap = 0; igap<3; ++igap){
            efieldgap[igap] = Efield(igap);
            driftVelocitygap[igap] = DriftVelocity(efieldgap[igap], temperature);
            fXTicksCoefficientgap[igap] = 0.001 * driftVelocitygap[igap] * samplingRate;
          }
          
          // Calculate geometric time offset.
          // only works if xyz[0]<=0
          const double* xyz = tpcgeom.PlaneLocation(0);
          
          fXTicksOffsets[cstat][tpc][plane] = -xyz[0]/(dir * fXTicksCoefficient) + triggerOffset;
          if (nplane==3){
            /*
         |    ---------- plane = 2 (collection)
         |                      Coeff[2]
         |    ---------- plane = 1 (2nd induction)
         |                      Coeff[1]
         |    ---------- plane = 0 (1st induction) x = xyz[0]
         |                      Coeff[0]
         |    ---------- x = 0
         V     For plane = 0, t offset is -xyz[0]/Coeff[0]
         x   */
            for (int ip = 0; ip < plane; ++ip){
              fXTicksOffsets[cstat][tpc][plane] += tpcgeom.PlanePitch(ip,ip+1)/fXTicksCoefficientgap[ip+1];
            }
          }          
          else if (nplane==2){ ///< special case for ArgoNeuT
            /*
         |    ---------- plane = 1 (collection)
         |                      Coeff[2]
         |    ---------- plane = 0 (2nd induction) x = xyz[0]
         |    ---------- x = 0, Coeff[1]
         V    ---------- first induction plane
         x                      Coeff[0]
For plane = 0, t offset is pitch/Coeff[1] - (pitch+xyz[0])/Coeff[0]
                         = -xyz[0]/Coeff[0] - pitch*(1/Coeff[0]-1/Coeff[1])
            */
            for (int ip = 0; ip < plane; ++ip){
              fXTicksOffsets[cstat][tpc][plane] += tpcgeom.PlanePitch(ip,ip+1)/fXTicksCoefficientgap[ip+2];
            }
            fXTicksOffsets[cstat][tpc][plane] -= tpcgeom.PlanePitch()*(1/fXTicksCoefficient-1/fXTicksCoefficientgap[1]);
          }
          
          // Add view dependent offset
          // FIXME the offset should be plane-dependent
          geo::View_t view = pgeom.View();
          switch (view) {
            case geo::kU:
              fXTicksOffsets[cstat][tpc][plane] += fTimeOffsetU;
              break;
            case geo::kV:
              fXTicksOffsets[cstat][tpc][plane] += fTimeOffsetV;
              break;
            case geo::kZ:
              fXTicksOffsets[cstat][tpc][plane] += fTimeOffsetZ;
              break;
            case geo::kY:
              fXTicksOffsets[cstat][tpc][plane] += fTimeOffsetY;
              break;
            case geo::kX:
              fXTicksOffsets[cstat][tpc][plane] += fTimeOffsetX;
              break;
            default:
              throw cet::exception(__FUNCTION__) << "Bad view = " << view << "\n" ;
          } // switch
        }
      }
    }
  }
  //--------------------------------------------------------------------
  // Get scale factor for x<-->ticks
  double DetectorPropertiesStandard::GetXTicksCoefficient(int t, int c) const
  {
    return fXTicksCoefficient * fDriftDirection.at(c).at(t);
  }
  //--------------------------------------------------------------------
  // Get scale factor for x<-->ticks
  double DetectorPropertiesStandard::GetXTicksCoefficient() const
  {
    return fXTicksCoefficient;
  }
  //--------------------------------------------------------------------
  //  Get offset for x<-->ticks
  double DetectorPropertiesStandard::GetXTicksOffset(int p, int t, int c) const
  {
    return fXTicksOffsets.at(c).at(t).at(p);        
  }
  
  //--------------------------------------------------------------------
  std::string DetectorPropertiesStandard::CheckTimeOffsetConfigurationAfterSetup
    () const
  {
    
    std::ostringstream errors;
    auto const views = fGeo->Views();
    
    if ((views.count(geo::kU) != 0) != fHasTimeOffsetU) {
      if (fHasTimeOffsetU)
        errors << "TimeOffsetU has been specified, but no U view is present.\n";
      else
        errors << "TimeOffsetU missing for view U.\n";
    }
    if ((views.count(geo::kV) != 0) != fHasTimeOffsetV) {
      if (fHasTimeOffsetV)
        errors << "TimeOffsetV has been specified, but no V view is present.\n";
      else
        errors << "TimeOffsetV missing for view Z.\n";
    }
    if ((views.count(geo::kZ) != 0) != fHasTimeOffsetZ) {
      if (fHasTimeOffsetZ)
        errors << "TimeOffsetZ has been specified, but no Z view is present.\n";
      else
        errors << "TimeOffsetZ missing for view Z.\n";
    }
    if ((views.count(geo::kY) != 0) != fHasTimeOffsetY) {
      if (fHasTimeOffsetY)
        errors << "TimeOffsetY has been specified, but no Y view is present.\n";
      else
        errors << "TimeOffsetY missing for view Y.\n";
    }
    if ((views.count(geo::kX) != 0) != fHasTimeOffsetX) {
      if (fHasTimeOffsetX)
        errors << "TimeOffsetX has been specified, but no X view is present.\n";
      else
        errors << "TimeOffsetX missing for view X.\n";
    }
    
    return errors.str();
    
  } // DetectorPropertiesStandard::CheckTimeOffsetConfigurationAfterSetup()
  
  //--------------------------------------------------------------------
  void DetectorPropertiesStandard::CheckConfigurationAfterSetup() const {
    
    std::string errors;
    
    errors += CheckTimeOffsetConfigurationAfterSetup();
    
    if (!errors.empty()) {
      throw cet::exception("DetectorPropertiesStandard")
        << "Detected configuration errors: \n" << errors;
    }
    
  } // DetectorPropertiesStandard::CheckConfigurationAfterSetup()
  
  //--------------------------------------------------------------------
  void DetectorPropertiesStandard::DoUpdateClocks() 
  {
    CalculateXTicksParams();
  }
  //--------------------------------------------------------------------
  
  
} // namespace