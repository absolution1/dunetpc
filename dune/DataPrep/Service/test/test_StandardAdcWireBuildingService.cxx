// test_StandardAdcWireBuildingService.cxx
//
// David Adams
// May 2016
//
// Test StandardAdcWireBuildingService.
//

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "lardataobj/RawData/RawDigit.h"
#include "dune/ArtSupport/ArtServiceHelper.h"
#include "dune/DuneInterface/AdcWireBuildingService.h"

#undef NDEBUG
#include <cassert>

using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::istringstream;
using std::ofstream;
using std::setw;
using std::setprecision;
using std::fixed;
using art::ServiceHandle;
using recob::Wire;

typedef vector<unsigned int> IndexVector;

//**********************************************************************

int test_StandardAdcWireBuildingService(int a_LogLevel =1) {
  const string myname = "test_StandardAdcWireBuildingService: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  cout << myname << "Create top-level FCL." << endl;
  string fclfile = "test_StandardAdcWireBuildingService.fcl";
  ofstream fout(fclfile.c_str());
  fout << "#include \"services_dune.fcl\"" << endl;
  fout << "services:      @local::dune35t_services" << endl;
  fout << "services.AdcWireBuildingService: {" << endl;
  fout << "  service_provider: StandardAdcWireBuildingService" << endl;
  fout << "  LogLevel:       " << a_LogLevel << endl;
  fout << "}" << endl;
  fout.close();

  cout << myname << "Fetch art service helper." << endl;
  ArtServiceHelper& ash = ArtServiceHelper::instance();
  ash.print();

  if ( ash.serviceStatus() == 0 ) {
    cout << myname << line << endl;
    cout << myname << "Add supporting services." << endl;
    assert( ash.addService("ExptGeoHelperInterface", fclfile, true) == 0 );
    assert( ash.addService("Geometry", fclfile, true) == 0 );

    cout << myname << line << endl;
    cout << myname << "Add wire building service." << endl;
    assert( ash.addService("AdcWireBuildingService", fclfile, true) == 0 );
    ash.print();

    cout << myname << line << endl;
    cout << myname << "Load services." << endl;
    assert( ash.loadServices() == 1 );
    ash.print();
  }

  cout << myname << line << endl;
  float fac = 250.0;
  cout << myname << "Create a response function." << endl;
  const int nsigshape = 10;
  AdcSignalVector sigshape;
  for ( unsigned int isig=0; isig<5; ++isig ) sigshape.push_back(fac*isig);
  for ( unsigned int isig=0; isig<5; ++isig ) sigshape.push_back(-sigshape.at(5-isig-1));
  assert( sigshape.size() == nsigshape );

  cout << myname << line << endl;
  cout << myname << "Create digits." << endl;
  const unsigned int ndig = 4;
  const unsigned int ntrk = 3;
  cout << myname << "  Digit count: " << ndig << endl;
  cout << myname << "  Track count: " << ntrk << endl;
  std::vector<raw::RawDigit> digits;
  digits.reserve(ndig);
  unsigned int sigoff[ntrk][ndig]; 
  sigoff[0][0] = 10;
  sigoff[0][1] = 11;
  sigoff[0][2] = 12;
  sigoff[0][3] = 13;
  sigoff[1][0] = 30;
  sigoff[1][1] = 32;
  sigoff[1][2] = 34;
  sigoff[1][3] = 36;
  sigoff[2][0] = 60;
  sigoff[2][1] = 59;
  sigoff[2][2] = 58;
  sigoff[2][3] = 57;
  unsigned int nsig = 80;
  float ped[ndig] = {2000.2, 2000.4, 2001.0, 1999.5};
  unsigned int chan[ndig] = {120, 121, 122, 123};
  AdcSignalVector sigsin(nsig, 0.0);
  AdcChannelDataMap acds;
  for ( unsigned int idig=0; idig<ndig; ++idig ) {
    cout << myname << "  Digit " << idig << endl;
    for ( int isig=0; isig<10; ++isig ) {
      sigsin[sigoff[0][idig]+isig] += sigshape[isig];
      sigsin[sigoff[1][idig]+isig] += sigshape[isig];
      sigsin[sigoff[2][idig]+isig] += sigshape[isig];
    }
    assert(sigsin.size() == nsig);
    AdcCountVector adcsin;
    for ( unsigned int isig=0; isig<nsig; ++isig) {
      AdcSignal sig = sigsin[isig] + ped[idig];
      AdcCount adc = 0.0;
      if ( sig > 0.0 ) adc = int(sig+0.5);
      if ( adc > 4095 ) adc = 4095;
      adcsin.push_back(adc);
    }
    assert(adcsin.size() == nsig);
    raw::RawDigit dig(chan[idig], nsig, adcsin, raw::kNone);
    dig.SetPedestal(ped[idig]);
    cout << myname << "      Compressed size: " << dig.NADC() << endl;
    cout << myname << "    Uncompressed size: " << dig.Samples() << endl;
    cout << myname << "             Pedestal: " << dig.GetPedestal() << endl;
    cout << myname << "              Channel: " << dig.Channel() << endl;
    assert(dig.Samples() == nsig);
    assert(dig.Channel() == chan[idig]);
    assert(dig.GetPedestal() == ped[idig]);
    AdcRoiVector rois;
    rois.push_back(AdcRoi(sigoff[0][idig], sigoff[0][idig]+nsigshape));
    rois.push_back(AdcRoi(sigoff[1][idig], sigoff[1][idig]+nsigshape));
    rois.push_back(AdcRoi(sigoff[2][idig], sigoff[2][idig]+nsigshape));
    cout << myname << "    ROI count: " << rois.size() << endl;
    for ( const AdcRoi& roi : rois ) {
      cout << myname << "      (" << roi.first << ", " << roi.second << ")" << endl;
    }
    cout << myname << "    ROI count: " << rois.size() << endl;
    AdcChannelData acd;
    acd.samples = sigsin;
    acd.channel = chan[idig];
    acd.rois = rois;
    digits.push_back(std::move(dig));
    acd.digit = &digits.back();
    const AdcChannelData& acdMoved = acds[chan[idig]] = std::move(acd);
    cout << myname << "    Moved data channel: " << chan[idig] << " ?= " << acdMoved.channel
         << " ?= " << acdMoved.digit->Channel() << endl;
  }

  cout << myname << line << endl;
  cout << myname << "Fetch wire building service." << endl;
  ServiceHandle<AdcWireBuildingService> hwib;
  hwib->print(cout, myname);

  cout << myname << line << endl;
  cout << myname << "Build wires." << endl;
  vector<Wire> wires;
  wires.reserve(ndig);
  for ( AdcChannelDataMap::value_type& iacd : acds ) {
    AdcChannelData& acd = iacd.second;
    assert( acd.channel == iacd.first );
    assert( acd.digit->Channel() == iacd.first );
    int rstat = hwib->build(acd, &wires);
    cout << myname << "Channel: " << acd.channel
         << ", # samples: " << acd.samples.size()
         << ", # ROI: " << acd.rois.size()
         << ", wire: " << acd.wire
         << ", " << &(wires.back()) << endl;
    assert( rstat == 0 );
    cout << myname << "  Size of Wire: " << sizeof(wires[0]) << endl;
    for ( unsigned int iwir=0; iwir<wires.size(); ++iwir ) {
      cout << myname << "  Wire " << iwir << ": " << &(wires[iwir])
           << " is for channel " << wires[iwir].Channel() << endl;
    }
  }

  cout << myname << line << endl;
  cout << myname << "Check wires." << endl;
  cout << myname << "Wire count: " << wires.size() << endl;
  assert( wires.size() == ndig );
  for ( unsigned int idig=0; idig<ndig; ++idig ) {
    const Wire& wire = wires[idig];
    const AdcChannelData& acd = acds[wire.Channel()];
    cout << myname << "  Wire " << idig << endl;
    cout << myname << "    Channel " << wire.Channel() << endl;
    cout << myname << "    Signal count: " << wire.NSignal() << endl;
    cout << myname << "    ROI count: " << wire.SignalROI().n_ranges() << endl;
    cout << myname << "    Wire: " << acd.wire << " ?= " << &(wires[idig]) << endl;
    assert( wire.Channel() == acd.channel );
    assert( wire.NSignal() == acd.samples.size() );
    assert( wire.SignalROI().size() == acd.samples.size() );
    assert( wire.SignalROI().n_ranges() == acd.rois.size() );
  }

  cout << myname << line << endl;
  cout << myname << "Done." << endl;
  return 0;
}

//**********************************************************************

int main(int argc, char* argv[]) {
  int a_LogLevel = 1;
  if ( argc > 1 ) {
    istringstream ssarg(argv[1]);
    ssarg >> a_LogLevel;
  }
  return test_StandardAdcWireBuildingService(a_LogLevel);
}

//**********************************************************************
