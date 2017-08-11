// TestDigit.h

// David Adams
// August 2017
//
// Class to construct a digit for use in testing.

#ifndef TestDigit_H
#define TestDigit_H

class TestDigit {
public:

  // Ctor.
  TestDigit(int i=0);

  // Data.
  std::string line = "-----------------------------";
  AdcCount lowbits = 63;
  AdcCount highbits = 4095 - 63;
  AdcIndex channel = -1;
  AdcSignal pedestal = 0.0;
  const raw::RawDigit* pdig = nullptr;
  AdcIndex nsig = 0;
  AdcSignalVector sigsin;
  AdcFlagVector expflags;

private:

  std::unique_ptr<raw::RawDigit> digitPtr;

};

//**********************************************************************

TestDigit::TestDigit(int idig) {
  const string myname = "TestDigit::ctor: ";
  cout << myname << line << endl;
  cout << myname << "Creating raw digit " << idig << endl;
  float fac = 250.0;
  pedestal = 2000.2;
  for ( int i=0; i<10; ++i ) sigsin.push_back(fac*i);
  for ( int i=10; i>=0; --i ) sigsin.push_back(fac*i);
  for ( int i=19; i>=0; --i ) sigsin.push_back(-sigsin[i]);
  if ( idig == 2 ) {   // Add a gap and then repeat the waveform.
    for ( unsigned int ksig = 0; ksig<30; ++ksig ) {
      sigsin.push_back(pedestal);
    }
    for ( unsigned int isig=0; isig<nsig; ++isig) {
      sigsin.push_back(sigsin[isig]);
    }
  }
  nsig = sigsin.size();
  assert(sigsin.size() == nsig);
  AdcCountVector adcsin;
  unsigned int isig_stucklo = 5;
  unsigned int isig_stuckhi = 15;
  for ( unsigned int isig=0; isig<nsig; ++isig) {
    AdcSignal sig = sigsin[isig] + pedestal;
    AdcCount adc = 0.0;
    if ( sig > 0.0 ) adc = int(sig+0.5);
    if ( adc > 4095 ) adc = 4095;
    AdcCount adchigh = adc & highbits;
    if ( isig == isig_stucklo ) adc = adchigh;           // Stuck low bits to zero.
    if ( isig == isig_stuckhi ) adc = adchigh + lowbits; // Stuck low bits to one.
    adcsin.push_back(adc);
  }
  assert(adcsin.size() == nsig);
  channel = 123;
  digitPtr.reset(new raw::RawDigit(channel, nsig, adcsin, raw::kNone));
  raw::RawDigit& dig = *digitPtr;
  dig.SetPedestal(pedestal);
  cout << myname << "    Compressed size: " << dig.NADC() << endl;
  cout << myname << "  Uncompressed size: " << dig.Samples() << endl;
  cout << myname << "           Pedestal: " << dig.GetPedestal() << endl;
  cout << myname << "            Channel: " << dig.Channel() << endl;
  assert(dig.Samples() == nsig);
  assert(dig.Channel() == channel);
  assert(dig.GetPedestal() == pedestal);

  cout << myname << line << endl;
  cout << myname << "Create the expected flag vector." << endl;
  expflags.resize(nsig, AdcGood);
  for ( unsigned int isig=0; isig<nsig; ++isig) {
    AdcCount adc = adcsin[isig];
    AdcCount adclow = adc & lowbits;
    if ( adc <= 0 ) expflags[isig] = AdcUnderflow;
    else if ( adc >= 4095 ) expflags[isig] = AdcOverflow;
    else if ( adclow == 0 ) expflags[isig] = AdcStuckOff;
    else if ( adclow == lowbits ) expflags[isig] = AdcStuckOn;
  }
  pdig = digitPtr.get();
}

//**********************************************************************

#endif
