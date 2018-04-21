// AdcDataDumper_tool.cc

#include "AdcDataDumper.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

using std::string;
using std::cout;
using std::endl;
using std::ofstream;
using std::ostream;
using std::ostringstream;
using std::setw;

using Index = unsigned int;

//**********************************************************************
// Local definitions.
//**********************************************************************

namespace {

template<class T>
T divide(double sum, Index nval) {
  cout << "AdcDataDumper::divide: Invalid type." << endl;
  return 0.0;
}

template<>
AdcSignal divide<AdcSignal>(double sum, Index nval) {
  return sum/nval;
}

template<>
AdcCount divide<AdcCount>(double sum, Index nval) {
  return AdcCount(sum/nval + 0.5);
}

/*  Same type as AdcCount
template<>
AdcFlag divide<AdcFlag>(double sum, Index nval) {
  if ( sum == 0.0 ) return 0;
  return 100;
}
*/

//**********************************************************************

// Apply limits and rebinning to an input vector.

template<class V>
class DisplayVector {
public:
  V vout;
  DisplayVector(Index first, Index a_rebin, Index maxn, const V& vin) {
    Index rebin = a_rebin == 0 ? 1 : a_rebin;
    for ( Index isam1=first; isam1<vin.size(); isam1+=rebin ) {
      typename V::value_type binval = vin[isam1];
      if ( rebin > 1 ) {
        double binsum = 0.0;
        Index isam2 = isam1 + rebin;
        if ( isam2 > vin.size() ) break;
        for ( Index isam=isam1; isam<isam2; ++isam ) {
          binsum += vin[isam];
        }
        binval = divide<typename V::value_type>(binsum, rebin);
      }
      vout.push_back(binval);
      if ( vout.size() >= maxn ) break;
    }
  }
};

//**********************************************************************

char charThresh(double val, double thresh) {
  if ( val > thresh ) return '+';
  if ( val < -thresh ) return '-';
  return ' ';
}

}  // end unnamed namespace

//**********************************************************************
// Class methods.
//**********************************************************************

AdcDataDumper::AdcDataDumper(fhicl::ParameterSet const& ps)
: m_FileName(ps.get<string>("FileName")),
  m_Prefix(ps.get<string>("Prefix")),
  m_NewFile(ps.get<bool>("NewFile")),
  m_ShowChannelCount(ps.get<bool>("ShowChannelCount")),
  m_ShowTickCounts(ps.get<bool>("ShowTickCounts")),
  m_ShowRaw(ps.get<bool>("ShowRaw")),
  m_ShowPrepared(ps.get<bool>("ShowPrepared")),
  m_ShowFirst(ps.get<unsigned int>("ShowFirst")),
  m_ShowRebin(ps.get<unsigned int>("ShowRebin")),
  m_ShowMax(ps.get<unsigned int>("ShowMax")),
  m_ShowThreshold(ps.get<float>("ShowThreshold")),
  m_ShowOpt(ps.get<unsigned int>("ShowOpt")),
  m_pout(nullptr) {
  if ( m_FileName.size() == 0 ) m_pout = &cout;
  else if ( ! m_NewFile ) m_pout = new ofstream(m_FileName.c_str());
}

//**********************************************************************

AdcDataDumper::~AdcDataDumper() {
  if ( m_pout != nullptr && m_pout != &cout ) {
    delete m_pout;
    m_pout = nullptr;
  }
}

//**********************************************************************

DataMap AdcDataDumper::viewMap(const AdcChannelDataMap& acds) const {
  DataMap ret;
  ostream* pout = m_pout;
  bool newfile = pout == nullptr;
  // If file is not already set, build the file name and open it.
  if ( newfile ) {
    if ( ! m_NewFile ) return ret.setStatus(1);
    string fname = m_FileName;
    string::size_type npos = string::npos;
    string::size_type ipos = fname.find("%PAT%");
    ipos = fname.find("%CHAN1%");
    if ( ipos != npos ) {
      string srep = "NOCHAN";
      if ( acds.size() ) {
        ostringstream ssrep;
        ssrep << acds.begin()->first;
        srep = ssrep.str();
      }
      while ( ipos != npos ) {
        fname.replace(ipos, 7, srep);
        ipos = fname.find("%CHAN1%", ipos+7);
      }
    }
    ipos = fname.find("%CHAN2%");
    if ( ipos != npos ) {
      string srep = "NOCHAN";
      if ( acds.size() ) {
        ostringstream ssrep;
        ssrep << acds.rbegin()->first;
        string srep = ssrep.str();
        while ( ipos != npos ) {
          fname.replace(ipos, 7, srep);
          ipos = fname.find("%CHAN2%", ipos+7);
        }
      }
    }
    pout = new ofstream(fname.c_str());
  }
  if ( pout == nullptr ) return ret.setStatus(2);
  ostream& out = *pout;
  string pre = m_Prefix + ":";
  if ( m_ShowChannelCount ) out << pre << "  Channel count: " << acds.size() << endl;
  Index wcha = 6;
  Index wcou = 0;
  if ( m_ShowRaw || m_ShowPrepared ) {
    out << pre << " Values are displayed starting at tick " << m_ShowFirst;
    if ( m_ShowRebin > 1 ) out << " with rebinning of " << m_ShowRebin;
    else out << " without rebinning";
    out << endl;
  }
  for ( const AdcChannelDataMap::value_type& iacd : acds ) {
    const AdcChannelData& acd = iacd.second;
    ostringstream sschanpre;
    sschanpre << pre << setw(wcha) << acd.channel << ":";
    string chanpre = sschanpre.str();
    sschanpre.str("");
    sschanpre << pre << setw(wcha+1) << " ";
    string nochanpre = sschanpre.str();
    if ( m_ShowTickCounts ) {
      out << chanpre;
      out << " nraw=" << setw(wcou) << acd.raw.size();
      out << " nsam=" << setw(wcou) << acd.samples.size();
      out << " nflg=" << setw(wcou) << acd.flags.size();
      out << " nsig=" << setw(wcou) << acd.signal.size();
      out << " nroi=" << setw(wcou) << acd.rois.size();
      out << endl;
      chanpre = nochanpre;
    }
    if ( m_ShowRaw ) {
      DisplayVector<AdcCountVector> dvec(m_ShowFirst, m_ShowRebin, m_ShowMax, acd.raw);
      out << chanpre << " Raw:";
      if ( m_ShowOpt == 2 ) {
        out << " |";
        for ( AdcCount val : dvec.vout ) out << setw(1) << charThresh(val, m_ShowThreshold);
        out << "|";
      } else {
        for ( AdcCount val : dvec.vout ) out << setw(6) << val;
      }
      out << endl;
      chanpre = nochanpre;
    }
    if ( m_ShowPrepared ) {
      DisplayVector<AdcSignalVector> dvec(m_ShowFirst, m_ShowRebin, m_ShowMax, acd.samples);
      out << chanpre << " Prp:";
      if ( m_ShowOpt == 2 ) {
        out << " |";
        for ( AdcCount val : dvec.vout ) out << setw(1) << charThresh(val, m_ShowThreshold);
        out << "|";
      } else {
        for ( AdcCount val : dvec.vout ) out << setw(6) << int(round(val));
      }
      out << endl;
      chanpre = nochanpre;
    }
  }
  if ( newfile ) delete pout;
  return ret;
}

//**********************************************************************
