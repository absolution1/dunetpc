// InterpolatingAdcMitigationService_service.cc

#include "InterpolatingAdcMitigationService.h"
#include <iostream>
#include "art/Framework/Services/Registry/ServiceHandle.h"

using std::string;
using std::cout;
using std::endl;

//**********************************************************************

InterpolatingAdcMitigationService::
InterpolatingAdcMitigationService(fhicl::ParameterSet const& pset, art::ActivityRegistry&)
: m_LogLevel(1) {
  const string myname = "InterpolatingAdcMitigationService::ctor: ";
  pset.get_if_present<int>("LogLevel", m_LogLevel);
  m_SkipUnderflows = pset.get<bool>("SkipUnderflows");
  m_SkipOverflows = pset.get<bool>("SkipOverflows");
  m_MaxConsecutiveSamples = pset.get<int>("MaxConsecutiveSamples");
  m_MaxConsecutiveFlag = pset.get<int>("MaxConsecutiveFlag");
  print(cout, myname);
}

//**********************************************************************

int InterpolatingAdcMitigationService::
update(AdcChannelData& data) const {
  const string myname = "InterpolatingAdcMitigationService:update: ";
  if ( m_LogLevel >= 3 ) {
    cout << myname << "Entering..." << endl;
    cout << myname << "          Channel: " << data.channel << endl;
    cout << myname << "Input vector size: " << data.samples.size() << endl;
  }
  AdcSignalVector& sigs = data.samples;
  AdcFlagVector& flags = data.flags;
  unsigned int isigFirst = sigs.size();   // First sample in the bad sample sequence.
  bool useMax = m_MaxConsecutiveSamples >= 0;
  unsigned int maxsig = m_MaxConsecutiveSamples;
  unsigned int nupdated = 0;
  // Loop over samples.
  for ( unsigned int isig=0; isig<sigs.size(); ++isig ) {
    AdcFlag flag = flags[isig];
    unsigned int isLast = isig == sigs.size()-1;
    // Set isBad which indicates if this is value that is flagged to be updated.
    bool isBad = false;
    if ( flag != AdcGood ) {
      isBad = true;
      if      ( m_SkipUnderflows && flag == AdcUnderflow ) isBad = false;
      else if ( m_SkipOverflows  && flag ==  AdcOverflow ) isBad = false;
    }
    // Record if this is the first sample in a new bad sequence.
    if ( isBad && isigFirst > isig ) isigFirst = isig;
    // Check if we have found the end of a sequence of bad samples.
    // This is indicated with a nonzero value for isigLast.
    // We will update samples in the range [isigFirst, isigLast].
    bool endBadSequence = false;
    unsigned int isigLast = 0;
    if ( !isBad && isigFirst < isig ) {
      endBadSequence = true;
      isigLast = isig - 1;
    }
    if ( isLast && isBad ) {
      endBadSequence = true;
      isigLast = isig;
    }
    // Update the current sequence of bad samples.
    if ( endBadSequence ) {
      unsigned int nsig = isigLast - isigFirst + 1;
      bool tooMany = useMax && nsig > maxsig;
      if ( m_LogLevel > 2 ) {
        cout << myname << " Updating at sample " << isig << ":"
             << " range=[" << isigFirst << "," << isigLast << "],"
             << " tooMany=" << tooMany << ", isLast=" << isLast << endl;
      }
      // For too many samples or beginning or end of data, use MaxConsecutiveFlag.
      if ( tooMany || isigFirst == 0 || isLast ) {
        if ( m_MaxConsecutiveFlag == 1 ) {
          for ( unsigned isig=isigFirst; isig<=isigLast; ++isig ) {
            sigs[isig] = 0.0;
            flags[isig] = AdcSetFixed;
            ++nupdated;
          }
        }
      // Otherwise, interpolate: sig_i = a*i + b
      } else {
        unsigned int isig1 = isigFirst - 1;
        unsigned int isig2 = isigLast + 1;
        double sig1 = sigs[isig1];
        double sig2 = sigs[isigLast+1];
        double fac = 1.0/(isig2 - isig1);
        double a = fac*(sig2 - sig1);
        double b = fac*(isig2*sig1 - isig1*sig2);
        for ( unsigned isig=isigFirst; isig<=isigLast; ++isig ) {
          sigs[isig] = a*isig + b;
          flags[isig] = AdcInterpolated;
          ++nupdated;
        }
      }
      isigFirst = sigs.size();
    }
  }
  if ( m_LogLevel >= 2 ) cout << myname << "Channel " << data.channel
                              << ": # updated/total = " << nupdated << "/" << sigs.size() << endl;
  return 0;
}

//**********************************************************************

std::ostream& InterpolatingAdcMitigationService::
print(std::ostream& out, std::string prefix) const {
  out << prefix << "InterpolatingAdcMitigationService:"                   << endl;
  out << prefix << "               LogLevel: " << m_LogLevel              << endl;
  out << prefix << "         SkipUnderflows: " << m_SkipUnderflows        << endl;
  out << prefix << "          SkipOverflows: " << m_SkipOverflows         << endl;
  out << prefix << "  MaxConsecutiveSamples: " << m_MaxConsecutiveSamples << endl;
  out << prefix << "     MaxConsecutiveFlag: " << m_MaxConsecutiveFlag    << endl;
  return out;
}

//**********************************************************************

DEFINE_ART_SERVICE_INTERFACE_IMPL(InterpolatingAdcMitigationService, AdcMitigationService)

//**********************************************************************
