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
#include <map>

namespace {
  
  template <typename T>
  inline T sqr(T v) { return v*v; }
  
} // local namespace
namespace spdp{
  //--------------------------------------------------------------------
  DetectorPropertiesProtoDUNEsp::DetectorPropertiesProtoDUNEsp() :
    fLP(0), fClocks(0), fGeo(0)
  {
  }
  
  //--------------------------------------------------------------------
  DetectorPropertiesProtoDUNEsp::DetectorPropertiesProtoDUNEsp(fhicl::ParameterSet const& pset,
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
  DetectorPropertiesProtoDUNEsp::DetectorPropertiesProtoDUNEsp(fhicl::ParameterSet const& pset,
                                         providers_type providers,
                                         std::set<std::string> const& ignore_params /* = {} */
                                         ):
    DetectorPropertiesProtoDUNEsp(pset,
      providers.get<geo::GeometryCore>(),
      providers.get<detinfo::LArProperties>(),
      providers.get<detinfo::DetectorClocks>(),
      ignore_params
      )
    {}
  
  //--------------------------------------------------------------------
  bool DetectorPropertiesProtoDUNEsp::Update(uint64_t)
  {
    DoUpdateClocks();


    // auto *tpchv = lar::providerFrom<slowctrls::TPCHVServiceProtoDUNE>();

    // float t = 1542662093;
    // std::string name="NP04_DCS_01:Heinz_V";
    // double rawHV = tpchv->GetValue(name, t);
    // std::cout<<rawHV<<std::endl;

  
    // auto *tpchv = lar::providerFrom<slowctrls::TPCHVServiceProtoDUNE>();
    //slowctrls::SlowControls* tpchv_copy = tpchv;
    // float t = 1542662093;
    // double rawHV = tpchv->GetValue("NP04_DCS_01:Heinz_V", t);
    
    // //double rawCurr = tpchv->GetValue("NP04_DCS_01:Heinz_V",t);
    // std::cout<<rawHV<<std::endl;
    //std::cout<<rawCurr<<std::endl;
    return true;
  }

  bool DetectorPropertiesProtoDUNEsp::UpdateHV(uint64_t run) 
  {
    
    bool retVal = true;
   
    
 
std::map<uint64_t,  double > HVmap = {
  //Good Beam run list
 {5141, 180}, {5143, 180}, {5145, 180}, {5146, 160}, 
 {5152, 180}, {5158, 180}, {5174, 180}, {5181, 180}, 
 {5185, 180}, {5190, 180}, {5194, 180}, {5199, 180}, 
 {5203, 180}, {5204, 180}, {5205, 150}, {5209, 150}, 
 {5211, 150}, {5212, 130}, {5213, 130}, {5815, 180}, 
 {5216, 130}, {5219, 180}, {5225, 180}, {5235, 180}, 
 {5240, 180}, {5244, 180}, {5249, 140}, {5250, 160}, 
 {5254, 160}, {5257, 140}, {5258, 140}, {5259, 140}, 
 {5260, 140}, {5261, 140}, {5267, 140}, {5276, 140}, 
 {5282, 140}, {5283, 140}, {5284, 140}, {5287, 140}, 
 {5290, 140}, {5293, 140}, {5298, 140}, {5301, 140}, 
 {5303, 140}, {5304, 140}, {5308, 180}, {5311, 180}, 
 {5313, 175}, {5315, 180}, {5338, 180}, {5341, 160}, 
 {5387, 180}, {5423, 180}, {5424, 180}, {5426, 180}, 
 {5455, 180}, {5456, 180}, {5457, 180}, {5458, 180}, 
 {5460, 180}, {5809, 180}, {5810, 180}, {5814, 180}, 
 {5816, 180}, {5817, 180}, {5842, 180}, {5843, 180}, 
 {5844, 180}, {5429, 180}, {5430, 180}, {5431, 180}, 
 {5432, 180}, {5433, 180}, {5434, 180}, {5437, 180}, 
 {5438, 180}, {5439, 180}, {5441, 180}, {5442, 180}, 
 {5449, 180}, {5450, 180}, {5451, 180}, {5452, 180}, 
 {5818, 180}, {5819, 180}, {5824, 180}, {5758, 180}, 
 {5759, 180}, {5760, 180}, {5762, 180}, {5765, 180}, 
 {5766, 180}, {5768, 180}, {5769, 180}, {5770, 180}, 
 {5771, 180}, {5772, 180}, {5773, 180}, {5774, 180}, 
 {5775, 180}, {5776, 180}, {5777, 180}, {5778, 180}, 
 {5779, 180}, {5780, 180}, {5783, 180}, {5784, 180}, 
 {5785, 180}, {5786, 180}, {5788, 180}, {5791, 180}, 
 {5792, 180}, {5794, 180}, {5796, 180}, {5797, 180}, 
 {5825, 180}, {5826, 180}, {5827, 180}, {5831, 180}, 
 {5833, 180}, {5836, 180}, {5837, 180}, {5834, 180}, 
 {5835, 180}, {5838, 180}, {5839, 180}, {5840, 180}, 
 {5841, 180},

  //Post Beam runs

  {5847, 180}, {5848, 180}, {5849, 180}, {5850, 180}, 
  {5851, 180}, {5852, 162}, {5853, 162}, {5854, 162}, 
  {5855, 144}, {5856, 126}, {5857, 108}, {5858, 90}, 
  {5859, 72}, {5860, 54}, {5861, 36}, {5863, 18}, 
  {5877, 0}, {5880, 0}, {5886, 0}, {5887, 0}, 
  {5888, 0}, {5889, 0}, {5890, 0}, {5891, 0}, 
  {5892, 0}, {5893, 0}, {5894, 0}, {5895, 0}, 
  {5896, 0}, {5897, 0}, {5898, 0}, {5899, 0}, 
  {5900, 0}, {5901, 0}, {5903, 0}, {5906, 180}, 
  {5907, 180}, {5910, 180}, {5912, 180}, {5924, 180}, 
  {5925, 180}, {5926, 180}, {5927, 180}, {5929, 180}, 
  {5930, 180}, {5931, 180}, {5932, 180}, {5933, 180}, 
  {5934, 180}, {5935, 180}, {5936, 180}, {5937, 180}, 
  {5940, 180}, {5941, 180}, {5944, 180}, {5945, 180}, 
  {5946, 180}, {5947, 180}, {5948, 180}, {5949, 180}, 
  {5950, 180}, {5951, 180}, {5952, 180}, {5953, 180}, 
  {5954, 180}, {5955, 180}, {5956, 180}, {5957, 180}, 
  {5958, 180}, {5960, 180}, {5961, 180}, {5962, 180}, 
  {5963, 180}, {5964, 180}, {5965, 180}, {5968, 180}, 
  {5969, 180}, {5970, 180}, {5971, 180}, {5972, 180}, 
  {5973, 180}, {5974, 180}, {5978, 180}, {5979, 180}, 
  {5980, 180}, {5981, 180}, {5982, 180}, {5983, 180}, 
  {5984, 180}, {5995, 0}, {5996, 0}, {5997, 0}, 
  {5998, 0}, {5999, 0}, {6000, 0}, {6001, 0}, 
  {6002, 0}, {6003, 0}, {6004, 0}, {6005, 0}, 
  {6006, 0}, {6012, 0}, {6013, 0}, {6014, 0}, 
  {6015, 0}, {6020, 0}, {6022, 0}, {6023, 0}, 
  {6024, 0}, {6025, 0}, {6027, 0}, {6028, 0}, 
  {6029, 0}, {6030, 0}, {6032, 0}, {6033, 0}, 
  {6034, 0}, {6035, 0}, {6036, 0}, {6037, 0}, 
  {6038, 0}, {6039, 0}, {6040, 0}, {6041, 0}, 
  {6042, 0}, {6043, 0}, {6045, 0}, {6046, 0}, 
  {6068, 0}, {6071, 180}, {6073, 180}, {6078, 180}, 
  {6079, 180}, {6081, 180}, {6082, 180}, {6083, 180}, 
  {6084, 180}, {6085, 180}, {6086, 180}, {6087, 180}, 
  {6088, 180}, {6089, 180}, {6090, 180}, {6091, 180}, 
  {6092, 180}, {6093, 180}, {6094, 180}, {6095, 180}, 
  {6096, 180}, {6097, 180}, {6098, 180}

};


   if(HVmap.find(run) == HVmap.end()) {
    std::cout<<"Run has no nomninal HV value entered, default to 180KV"<<std::endl;
    fHV_cath=180;
  }  
  else
    {fHV_cath=HVmap.at(run);
     std::cout<<"Using HV on cathode as: "<<fHV_cath<<"KV, \n Value retrived from table of runs"<<std::endl; 
  }



    
    double Gplane_bias=(-0.665*(fHV_cath/180));
    double Uplane_bias=(-0.370*(fHV_cath/180));
    double Vplane_bias=0;
    double Xplane_bias=0.820*(fHV_cath/180);
    fEfield={fHV_cath/360,std::fabs(Gplane_bias-Uplane_bias)/0.475,std::fabs(Uplane_bias-Vplane_bias)/0.475,std::fabs(Vplane_bias-Xplane_bias)/0.475};
    //fEfield[0]=fHV_cath/360;
    //fEfield[1]=std::fabs(Gplane_bias-Uplane_bias)/0.475;
    //fEfield[2]=std::fabs(Uplane_bias-Vplane_bias)/0.475;
    //fEfield[3]=std::fabs(Vplane_bias-Xplane_bias)/0.475;
    std::cout<<"Calculated E field in 4 plane gaps:"<<fEfield[0]<<","<<fEfield[1]<<","<<fEfield[2]<<","<<fEfield[3]<<std::endl;
    
    
    return retVal;
  }




  //--------------------------------------------------------------------
  bool DetectorPropertiesProtoDUNEsp::UpdateClocks(const detinfo::DetectorClocks* clks) 
  {
    fClocks = clks;
    
    fTPCClock = fClocks->TPCClock();
    DoUpdateClocks();
    return true;
  }
  
  //------------------------------------------------------------
  double DetectorPropertiesProtoDUNEsp::ConvertTDCToTicks(double tdc) const
  {
    return fClocks->TPCTDC2Tick(tdc);
  }
  
  //--------------------------------------------------------------
  double DetectorPropertiesProtoDUNEsp::ConvertTicksToTDC(double ticks) const
  {
    return fClocks->TPCTick2TDC(ticks);
  }
  
 
  //--------------------------------------------------------------------
  void DetectorPropertiesProtoDUNEsp::Configure(Configuration_t const& config) {
    
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
    SetDetectorClocks(providers.get<detinfo::DetectorClocks>());
    
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
  
  //------------------------------------------------------------------------------------//
  int  DetectorPropertiesProtoDUNEsp::TriggerOffset()     const 
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
 
  double DetectorPropertiesProtoDUNEsp::ConvertXToTicks(double X, int p, int t, int c) const
  {
    return (X / (fXTicksCoefficient * fDriftDirection.at(c).at(t)) +  fXTicksOffsets.at(c).at(t).at(p) );
  }
  
  //-------------------------------------------------------------------
  // Take a cooridnate in ticks, and convert to an x position
  // assuming event deposit occured at t=0
 
  double  DetectorPropertiesProtoDUNEsp::ConvertTicksToX(double ticks, int p, int t, int c) const
  {
    return (ticks - fXTicksOffsets.at(c).at(t).at(p)) * fXTicksCoefficient * fDriftDirection.at(c).at(t);  
  }
  
  //--------------------------------------------------------------------
  void DetectorPropertiesProtoDUNEsp::CheckIfConfigured() const
  {
    if (!fGeo) throw cet::exception(__FUNCTION__) << "Geometry is uninitialized!";
    if (!fLP) throw cet::exception(__FUNCTION__) << "LArPropertiesStandard is uninitialized!";
    if (!fClocks) throw cet::exception(__FUNCTION__) << "DetectorClocks is uninitialized!";
  }
  
  
  //--------------------------------------------------------------------
  // Recalculte x<-->ticks conversion parameters from detector constants
  
  void DetectorPropertiesProtoDUNEsp::CalculateXTicksParams()
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
  double DetectorPropertiesProtoDUNEsp::GetXTicksCoefficient(int t, int c) const
  {
    return fXTicksCoefficient * fDriftDirection.at(c).at(t);
  }
  //--------------------------------------------------------------------
  // Get scale factor for x<-->ticks
  double DetectorPropertiesProtoDUNEsp::GetXTicksCoefficient() const
  {
    return fXTicksCoefficient;
  }
  //--------------------------------------------------------------------
  //  Get offset for x<-->ticks
  double DetectorPropertiesProtoDUNEsp::GetXTicksOffset(int p, int t, int c) const
  {
    return fXTicksOffsets.at(c).at(t).at(p);        
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
  
  //--------------------------------------------------------------------
  void DetectorPropertiesProtoDUNEsp::DoUpdateClocks() 
  {
    CalculateXTicksParams();
  }
  //--------------------------------------------------------------------
  
  
} // namespace