// test_GeoRopChannelGroupService.cxx

// David Adams
// October 2016
//
// Test GeoRopChannelGroupService.

#include "../GeoRopChannelGroupService.h"
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include "dune/ArtSupport/ArtServiceHelper.h"

using std::string;
using std::cout;
using std::endl;
using std::istringstream;
using std::ostringstream;
using std::ofstream;
using std::setw;
using std::vector;
using art::ServiceHandle;

#undef NDEBUG
#include <cassert>

int test_GeoRopChannelGroupService(string sgeo) {
  const string myname = "test_GeoRopChannelGroupService: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  const string line = "-----------------------------";

  cout << myname << line << endl;
  cout << myname << "Create fcl file." << endl;
  string fname = "test_GeoRopChannelGroupService.fcl";
  ofstream fout(fname.c_str());
  fout << "#include \"services_dune.fcl\"" << endl;
  fout << "services: @local::" << sgeo << endl;
  fout << "services.ChannelGroupService: {" << endl;
  fout << "  service_provider: \"GeoRopChannelGroupService\"" << endl;
  fout << "}" << endl;

  cout << myname << line << endl;
  cout << myname << "Fetch art service helper." << endl;
  ArtServiceHelper& ash = ArtServiceHelper::instance();
  ash.print();

  cout << line << endl;
  cout << myname << "Read service configurations." << endl;
  assert( ash.addServices(fname, true) == 0 );

  cout << myname << line << endl;
  cout << myname << "Load services." << endl;
  assert( ash.loadServices() == 1 );

  cout << myname << line << endl;
  cout << myname << "Fetch ChannelGroupService." << endl;
  ServiceHandle<ChannelGroupService> hcgs;
  hcgs->print(cout, myname);

  cout << myname << line << endl;
  cout << myname << "Geo fcl: " << sgeo << endl;

  cout << myname << line << endl;
  unsigned int nrop = hcgs->size();
  cout << myname << "Check ROP count: " << nrop << endl;
  if ( sgeo == "dunefd_services" ) {
    assert( nrop == 600 );
  } else {
    assert( nrop == 16 );
  }

  cout << myname << line << endl;
  cout << myname << "ROP names and channels:" << endl;
  for ( unsigned int irop=0; irop<nrop; ++irop ) {
    string name = hcgs->name(irop);
    const vector<ChannelGroupService::Index> chans = hcgs->channels(irop);
    cout << myname << "  " << setw(10) << name << ": ["
         << setw(6) << chans.front() << ", "
         << setw(6) << chans.back() << "]" << endl;
    assert( name.size() );
    assert( name != "NoSuchRop" );
    assert( chans.size() );
  }

  cout << myname << line << endl;
  cout << myname << "Done." << endl;
  return 0;
}

int main(int argc, char* argv[]) {
  string sgeo = "dune35t_basic_services";
  if ( argc > 1 ) {
    sgeo = argv[1];
    if ( sgeo == "-h" ) {
      cout << "Usage: " << argv[0] << endl;
      cout << "       " << argv[0] << " dune35t_basic_services" << endl;
      cout << "       " << argv[0] << " dunefd_services" << endl;
      return 0;
    }
  }
  return test_GeoRopChannelGroupService(sgeo);
  ArtServiceHelper::close();
  return 0;
}


//******************************************************************************
//*** Big configuration chunks
//***

// this is a copy of the "standard" LArProperties configuration
extern const std::string LArPropertiesServiceConfigurationString { R"cfg(
 service_provider: "LArPropertiesServiceStandard"

 # For following parameters, see http://pdg.lbl.gov/AtomicNuclearProperties/
 RadiationLength:  19.55   # g/cm^2
 AtomicNumber:     18      # Ar atomic number.
 AtomicMass:       39.948  # Ar atomic mass (g/mol).
 ExcitationEnergy: 188.0   # Ar mean excitation energy (eV).

# realistic Argon 39 decays
# Argon39DecayRate: 0.00141 # decays per cm^3 per second.  Assumes 1.01 Bq/kg and a density of 1.396 g/cc
# switch them off for faster simulation
 Argon39DecayRate: 0.0

 # Optical properties	
 # Fast and slow scintillation emission spectra, from [J Chem Phys vol 91 (1989) 1469]
 FastScintEnergies:    [ 6.0,  6.7,  7.1,  7.4,  7.7, 7.9,  8.1,  8.4,  8.5,  8.6,  8.8,  9.0,  9.1,  9.4,  9.8,  10.4,  10.7]
 SlowScintEnergies:    [ 6.0,  6.7,  7.1,  7.4,  7.7, 7.9,  8.1,  8.4,  8.5,  8.6,  8.8,  9.0,  9.1,  9.4,  9.8,  10.4,  10.7]
 FastScintSpectrum:    [ 0.0,  0.04, 0.12, 0.27, 0.44, 0.62, 0.80, 0.91, 0.92, 0.85, 0.70, 0.50, 0.31, 0.13, 0.04,  0.01, 0.0]
 SlowScintSpectrum:    [ 0.0,  0.04, 0.12, 0.27, 0.44, 0.62, 0.80, 0.91, 0.92, 0.85, 0.70, 0.50, 0.31, 0.13, 0.04,  0.01, 0.0]
 ScintResolutionScale: 1.     # resolution factor used by G4 scintillation
 ScintFastTimeConst:   6.     # fast scintillation time constant (ns)
 ScintSlowTimeConst:   1590.  # slow scintillation time constant (ns)
 ScintBirksConstant:   0.069  # birks constant for LAr (1/MeV cm)
 ScintYield:           24000. # total scintillation yield (ph/Mev)         
 ScintPreScale:        0.03   # Scale factor to reduce number of photons simulated
                              # Later QE should be corrected for this scale
 ScintYieldRatio:      0.3    # fast / slow scint ratio (needs revisitting)
 ScintByParticleType:  false  # whether to use different yields and
                              # quenching per particle in fast op sim
 EnableCerenkovLight: true    # whether to switch on cerenkov light (slow)




 # Scintillation yields and fast/slow ratios per particle type 
 MuonScintYield:          24000
 MuonScintYieldRatio:     0.23
 PionScintYield:          24000
 PionScintYieldRatio:     0.23 
 ElectronScintYield:      20000
 ElectronScintYieldRatio: 0.27
 KaonScintYield:          24000
 KaonScintYieldRatio:     0.23
 ProtonScintYield:        19200
 ProtonScintYieldRatio:   0.29
 AlphaScintYield:         16800
 AlphaScintYieldRatio:    0.56


 # Refractive index as a function of energy (eV) from arXiv:1502.04213v1
 RIndexEnergies: [ 1.900,  2.934,  3.592,  5.566,  6.694,  7.540,  8.574,  9.044,  9.232,  9.420,  9.514,  9.608,  9.702,  9.796,  9.890,  9.984,  10.08,  10.27,  10.45,  10.74,  10.92 ]
 RIndexSpectrum: [ 1.232,  1.236,  1.240,  1.261,  1.282,  1.306,  1.353,  1.387,  1.404,  1.423,  1.434,  1.446,  1.459,  1.473,  1.488,  1.505,  1.524,  1.569,  1.627,  1.751,  1.879 ]

 # absorption length as function of energy
 AbsLengthEnergies: [ 4,     5,     6,     7,     8,     9,     10,    11   ]       
 AbsLengthSpectrum: [ 2000., 2000., 2000., 2000., 2000., 2000., 2000., 2000.] 

 # Rayleigh scattering length (cm) @ 90K as a function of energy (eV) from arXiv:1502.04213
 RayleighEnergies: [   2.80,   3.00,   3.50,   4.00,  5.00,  6.00,  7.00,  8.00,  8.50,  9.00,  9.20,  9.40,  9.50,  9.60,  9.70,  9.80,  9.90,  10.0,  10.2,  10.4,  10.6, 10.8 ]
 RayleighSpectrum: [ 47923., 35981., 18825., 10653., 3972., 1681., 750.9, 334.7, 216.8, 135.0, 109.7, 88.06, 78.32, 69.34, 61.06, 53.46, 46.50, 40.13, 28.91, 19.81, 12.61, 7.20 ]

 # Surface reflectivity data - vector of energy spectrum per
 #   surface type
 ReflectiveSurfaceEnergies:           [ 7, 9, 10 ]             
 ReflectiveSurfaceNames:            [ "STEEL_STAINLESS_Fe7Cr2Ni" ]  
 ReflectiveSurfaceReflectances:     [ [ 0.25, 0.25, 0.25 ] ]        
 ReflectiveSurfaceDiffuseFractions: [ [ 0.5,  0.5,  0.5  ] ]        
)cfg"};


// this is a copy of the "standard" DetectorProperties configuration
extern const std::string DetectorPropertiesServiceConfigurationString { R"cfg(
service_provider: "DetectorPropertiesServiceStandard"

# Drift properties 
SternheimerA:     0.1956  # Ar Sternheimer parameter a.
SternheimerK:     3.0000  # Ar Sternheimer parameter k.
SternheimerX0:    0.2000  # Ar Sternheimer parameter x0.
SternheimerX1:    3.0000  # Ar Sternheimer parameter x0.
SternheimerCbar:  5.2146  # Ar Sternheimer parameter Cbar.

Temperature:       87
Electronlifetime:  3.0e3
Efield:           [0.5,0.666,0.8]  #(predicted for microBooNE)
ElectronsToADC:    6.8906513e-3 # 1fC = 43.008 ADC counts for DUNE fd
NumberTimeSamples: 4492         # drift length/drift velocity*sampling rate = (359.4 cm)/(0.16 cm/us)*(2 MHz)
ReadOutWindowSize: 4492         # drift length/drift velocity*sampling rate = (359.4 cm)/(0.16 cm/us)*(2 MHz)
TimeOffsetU:       0.
TimeOffsetV:       0.
TimeOffsetZ:       0.

SimpleBoundaryProcess: true  #enable opticalBoundaryProcessSimple instead of G4 default

)cfg"};
