/**Adapted for ProtoDUNE by Owen Goodwin (ogoodwin@fnal.gov) from
 * lardataalg/DetectorInfo/DetectorPropertiesStandard.cxx
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
#include <map>
#include <string>

#include "art_root_io/RootDB/SQLite3Wrapper.h"
#include "IFDH_service.h"


#include "dune/Protodune/singlephase/DetectorServices/Services/DetectorPropertiesServiceProtoDUNEsp.h"
#include "lardataalg/DetectorInfo/LArProperties.h"
#include "larcore/Geometry/Geometry.h"



#include "lardata/DetectorInfoServices/ServicePack.h" // lar::extractProviders()
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"







namespace {

  template <typename T>
  inline T sqr(T v) { return v*v; }

} // local namespace
namespace spdp{
  //--------------------------------------------------------------------
  DetectorPropertiesProtoDUNEsp::DetectorPropertiesProtoDUNEsp() :
    fLP(0), fGeo(0)
  {

  }

  //--------------------------------------------------------------------
  DetectorPropertiesProtoDUNEsp::DetectorPropertiesProtoDUNEsp(fhicl::ParameterSet const& pset,
                                                               const geo::GeometryCore* geo,
                                                               const detinfo::LArProperties* lp,
                                                               std::set<std::string> const& ignore_params /* = {} */
                                                               ):
    fLP(lp), fGeo(geo)
  {
    {
      mf::LogInfo debug("setupProvider<DetectorPropertiesStandard>");

      debug << "Asked to ignore " << ignore_params.size() << " keys:";
      for (auto const& key: ignore_params) debug << " '" << key << "'";
    }

    ValidateAndConfigure(pset, ignore_params);

  }

  //--------------------------------------------------------------------
  DetectorPropertiesProtoDUNEsp::DetectorPropertiesProtoDUNEsp(fhicl::ParameterSet const& pset,
                                                               providers_type providers,
                                                               std::set<std::string> const& ignore_params /* = {} */
                                                               ):
    DetectorPropertiesProtoDUNEsp(pset,
                                  providers.get<geo::GeometryCore>(),
                                  providers.get<detinfo::LArProperties>(),
                                  ignore_params
                                  )
  {}

  bool DetectorPropertiesProtoDUNEsp::UpdateTemp(int run)
  {


    if(fUseRunDependentTemperature){ // updata the temperature based on run number for data

      if ((run > 5903) & (run < 6930)){ //first runs after Nov 17 2019 and March 1 2019, where tempreture average was lower
        fTemperature = 87.36;
      }
      else if (run >= 6930)
      {
        fTemperature = 87.65;
      }
    }


    return true;
  }

  bool DetectorPropertiesProtoDUNEsp::UpdateReadoutWindowSize(detinfo::DetectorClocksData const& clockData,
                                                              std::string filename){

    bool retVal=false;

    if(fGetReadOutWindowSizefromMetaData){
      art::ServiceHandle<ifdh_ns::IFDH> ifdh;
      auto metadata=ifdh->getMetadata(filename);
      retVal = true;
      std::string window_str="DUNE_data.readout_window: ";
      if(metadata.find(window_str)!=std::string::npos){
        int n1 = metadata.find(window_str);
        n1 += window_str.length();
        int n2 = metadata.find("\n", n1);
        double window=std::stod(metadata.substr(n1, n2-n1)); //milliseconds
        double ticks=window*1000/(clockData.TPCClock().TickPeriod()); //sampling rate 2Mhz
        fNumberTimeSamples=ticks;
        fReadOutWindowSize=ticks;
        std::cout<<"Setting ReadOutWindowSize and NumberTimeSamples as : "<<fReadOutWindowSize<<" (Ticks), Value retrieved from Metadata"<<std::endl;
      }
      else{
        std::cout<<"Run has no specific readout window size entry in Metadata, default to fcl configured value : "<<fReadOutWindowSize<<" (Ticks)"<<std::endl;
      }
    }

    return retVal;

  }






  bool DetectorPropertiesProtoDUNEsp::UpdateHV(std::string filename)
  {



    bool retVal=false;

    if(fGetHVDriftfromMetaData){
      retVal = true;

      art::ServiceHandle<ifdh_ns::IFDH> ifdh;
      auto metadata=ifdh->getMetadata(filename);

      int run = 0;
      std::string run_str="Runs: ";
      if(metadata.find(run_str)!=std::string::npos){
        int n1 = metadata.find(run_str);
        n1 += run_str.length();
        int n2 = metadata.find(".", n1);
        run = std::stoi(metadata.substr(n1, n2-n1));
        std::cout<<"Run number from metadata: "<<run<<std::endl;
      }

      std::string hv_str="detector.hv_value: ";
      if(metadata.find(hv_str)!=std::string::npos){
        int n1 = metadata.find(hv_str);
        n1 += hv_str.length();
        int n2 = metadata.find("\n", n1);

        fHV_cath=std::stod(metadata.substr(n1, n2-n1));
        std::cout<<"Using HV on cathode as: "<<fHV_cath<<"KV,  Value retreived from samweb MetaData"<<std::endl;
      }
      else{
        std::cout<<"HV Metadata entry not found"<<std::endl;
        std::cout<<"Run has no nomninal HV value entered in Metadata, default to 180KV"<<std::endl;
        fHV_cath=180;
      }



      double Gplane_bias=(-0.665); //
      double Uplane_bias=(-0.370);
      double Vplane_bias=0;
      double Xplane_bias=(0.820);


      //double Gplane_bias=(-0.665*(fHV_cath/180)); // this is scaling, wasnt actually done.
      //double Uplane_bias=(-0.370*(fHV_cath/180));
      //double Vplane_bias=0;
      //double Xplane_bias=0.820*(fHV_cath/180);
      fEfield={fHV_cath/360,std::fabs(Gplane_bias-Uplane_bias)/0.47625,std::fabs(Uplane_bias-Vplane_bias)/0.47625,std::fabs(Vplane_bias-Xplane_bias)/0.47625};
      //fEfield[0]=fHV_cath/360;
      //fEfield[1]=std::fabs(Gplane_bias-Uplane_bias)/0.475;
      //fEfield[2]=std::fabs(Uplane_bias-Vplane_bias)/0.475;
      //fEfield[3]=std::fabs(Vplane_bias-Xplane_bias)/0.475;
      if (std::abs(fHV_cath-180) < 1e-6){
        // Use the corrected E field from Flavio, Francesco and Stefania
        // Ref: https://indico.fnal.gov/event/20939/contribution/1/material/slides/0.pdf
        if (run < 6725){ //first run in Feb 8, 2019
          fEfield[0] = 0.4867;
        }
        else{
          fEfield[0] = 0.4995;
        }
      }
      std::cout<<"Calculated E field in 4 plane gaps as: "<<fEfield[0]<<","<<fEfield[1]<<","<<fEfield[2]<<","<<fEfield[3]<<std::endl;
    }//End GetHVDriftfrom MetaData if


    return retVal;
  }




  //--------------------------------------------------------------------
  detinfo::DetectorPropertiesData
  DetectorPropertiesProtoDUNEsp::DataFor(const detinfo::DetectorClocksData& clockData) const
  {
    return CalculateXTicksParams(clockData);
  }

  //--------------------------------------------------------------------
  void DetectorPropertiesProtoDUNEsp::Configure(Configuration_t const& config) {

    fEfield                     = config.Efield();
    fGetHVDriftfromMetaData    = config.fGetHVDriftfromMetaData();
    fGetReadOutWindowSizefromMetaData = config.fGetReadOutWindowSizefromMetaData();
    fUseRunDependentTemperature = config.fUseRunDependentTemperature();
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

  } // DetectorPropertiesStandard::Configure()

  //--------------------------------------------------------------------
  DetectorPropertiesProtoDUNEsp::Configuration_t
  DetectorPropertiesProtoDUNEsp::ValidateConfiguration(
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
  void DetectorPropertiesProtoDUNEsp::ValidateAndConfigure(
                                                           fhicl::ParameterSet const& p,
                                                           std::set<std::string> const& ignore_params /* = {} */
                                                           ) {
    Configure(ValidateConfiguration(p, ignore_params));
  } // ValidateAndConfigure()


  //------------------------------------------------------------------------------------//
  void DetectorPropertiesProtoDUNEsp::Setup(providers_type providers) {

    SetGeometry(providers.get<geo::GeometryCore>());
    SetLArProperties(providers.get<detinfo::LArProperties>());

    CheckConfigurationAfterSetup();

  } // DetectorPropertiesStandard::Setup()


  //------------------------------------------------------------------------------------//
  double DetectorPropertiesProtoDUNEsp::Efield(unsigned int planegap) const
  {
    if(planegap >= fEfield.size())
      throw cet::exception("DetectorPropertiesStandard") << "requesting Electric field in a plane gap that is not defined\n";


    return fEfield[planegap];
  }


  //------------------------------------------------
  double DetectorPropertiesProtoDUNEsp::Density(double temperature) const
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
  double DetectorPropertiesProtoDUNEsp::Eloss(double mom, double mass, double tcut) const
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
  double DetectorPropertiesProtoDUNEsp::ElossVar(double mom, double mass) const
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
  double DetectorPropertiesProtoDUNEsp::DriftVelocity(double efield, double temperature) const {
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
  double DetectorPropertiesProtoDUNEsp::BirksCorrection(double dQdx) const
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
  double DetectorPropertiesProtoDUNEsp::ModBoxCorrection(double dQdx) const
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

  //--------------------------------------------------------------------
  void DetectorPropertiesProtoDUNEsp::CheckIfConfigured() const
  {
    if (!fGeo) throw cet::exception(__FUNCTION__) << "Geometry is uninitialized!";
    if (!fLP) throw cet::exception(__FUNCTION__) << "LArPropertiesStandard is uninitialized!";
  }


  //--------------------------------------------------------------------
  // Recalculte x<-->ticks conversion parameters from detector constants

  detinfo::DetectorPropertiesData
  DetectorPropertiesProtoDUNEsp::CalculateXTicksParams(detinfo::DetectorClocksData const& clockData) const
  {
    CheckIfConfigured();

    double samplingRate   = sampling_rate(clockData);
    double efield         = Efield();
    double temperature    = Temperature();
    double driftVelocity  = DriftVelocity(efield, temperature);

    double const x_ticks_coefficient = 0.001 * driftVelocity * samplingRate;
    double triggerOffset  = trigger_offset(clockData);
    std::vector<std::vector<std::vector<double>>> x_ticks_offsets(fGeo->Ncryostats());
    std::vector<std::vector<double>> drift_direction(fGeo->Ncryostats());
    for(size_t cstat = 0; cstat < fGeo->Ncryostats(); ++cstat){
      x_ticks_offsets[cstat].resize(fGeo->Cryostat(cstat).NTPC());
      drift_direction[cstat].resize(fGeo->Cryostat(cstat).NTPC());
      for(size_t tpc = 0; tpc < fGeo->Cryostat(cstat).NTPC(); ++tpc) {
        const geo::TPCGeo& tpcgeom = fGeo->Cryostat(cstat).TPC(tpc);
        const double dir((tpcgeom.DriftDirection() == geo::kNegX) ? +1.0 :-1.0);
        drift_direction[cstat][tpc] = dir;
        int nplane = tpcgeom.Nplanes();
        x_ticks_offsets[cstat][tpc].resize(nplane, 0.);

        for(int plane = 0; plane < nplane; ++plane) {
          const geo::PlaneGeo& pgeom = tpcgeom.Plane(plane);
          const double* xyz = tpcgeom.PlaneLocation(0);
          x_ticks_offsets[cstat][tpc][plane] = -xyz[0]/(dir * x_ticks_coefficient) + triggerOffset;

          // Add view dependent offset
          // FIXME the offset should be plane-dependent
          geo::View_t view = pgeom.View();
          switch (view) {
          case geo::kU:
            x_ticks_offsets[cstat][tpc][plane] += fTimeOffsetU;
            break;
          case geo::kV:
            x_ticks_offsets[cstat][tpc][plane] += fTimeOffsetV;
            break;
          case geo::kZ:
            x_ticks_offsets[cstat][tpc][plane] += fTimeOffsetZ;
            break;
          case geo::kY:
            x_ticks_offsets[cstat][tpc][plane] += fTimeOffsetY;
            break;
          case geo::kX:
            x_ticks_offsets[cstat][tpc][plane] += fTimeOffsetX;
            break;
          default:
            throw cet::exception(__FUNCTION__) << "Bad view = " << view << "\n" ;
          } // switch
        }
      }
    }
    return detinfo::DetectorPropertiesData{*this, x_ticks_coefficient, move(x_ticks_offsets), move(drift_direction)};
  }

  //--------------------------------------------------------------------
  std::string DetectorPropertiesProtoDUNEsp::CheckTimeOffsetConfigurationAfterSetup
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
  void DetectorPropertiesProtoDUNEsp::CheckConfigurationAfterSetup() const {

    std::string errors;

    errors += CheckTimeOffsetConfigurationAfterSetup();

    if (!errors.empty()) {
      throw cet::exception("DetectorPropertiesStandard")
        << "Detected configuration errors: \n" << errors;
    }

  } // DetectorPropertiesStandard::CheckConfigurationAfterSetup()

} // namespace
