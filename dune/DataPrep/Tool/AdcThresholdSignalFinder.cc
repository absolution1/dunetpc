#include "AdcThresholdSignalFinder.h"

#include "dune/ArtSupport/DuneToolManager.h"
#include "dune/DuneInterface/Data/RunData.h"
#include "dune/DuneInterface/Tool/RunDataTool.h"
#include "dune/DuneCommon/Utility/RootParFormula.h"
#include <iostream>

using std::string;
using std::cout;
using std::endl;

//**********************************************************************
// Local defintitions.
//**********************************************************************

namespace {

string boolToString(bool val) {
  return val ? "true" : "false";
}

}  // end unnamed namespace

//**********************************************************************
// Class methods.
//**********************************************************************

AdcThresholdSignalFinder::AdcThresholdSignalFinder(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")),
  m_Threshold(new RootParFormula("Threshold", ps.get<Name>("Threshold"))),
  m_BinsBefore(ps.get<unsigned int>("BinsBefore")),
  m_BinsAfter(ps.get<unsigned int>("BinsAfter")),
  m_FlagPositive(ps.get<bool>("FlagPositive")),
  m_FlagNegative(ps.get<bool>("FlagNegative")),
  m_prdtool(nullptr) {
  const string myname = "AdcThresholdSignalFinder::ctor: ";
  DuneToolManager* pdtm = DuneToolManager::instance();
  if ( pdtm == nullptr ) {
    cout << myname << "ERROR: Unable to retrieve tool manager." << endl;
  } else {
    string stnam = "runDataTool";
    if ( m_Threshold->npar() ) {
      m_prdtool = pdtm->getShared<RunDataTool>(stnam);
      if ( m_prdtool == nullptr ) {
        cout << myname << "ERROR: RunDataTool " << stnam
             << " not found. Scale factor formula will not be evaluated." << endl;
      } else {
        cout << myname << "RunDataTool retrieved." << endl;
      }
    }
  }
  if ( m_LogLevel >= 1 ) {
    cout << myname << "Configuration parameters:" << endl;
    cout << myname << "      LogLevel: " << m_LogLevel << endl;
    cout << myname << "     Threshold: " << m_Threshold->formulaString() << endl;
    cout << myname << "    BinsBefore: " << m_BinsBefore << endl;
    cout << myname << "     BinsAfter: " << m_BinsAfter << endl;
    cout << myname << "  FlagPositive: " << boolToString(m_FlagPositive) << endl;
    cout << myname << "  FlagNegative: " << boolToString(m_FlagNegative) << endl;
  }
}

//**********************************************************************

DataMap AdcThresholdSignalFinder::update(AdcChannelData& acd) const {
  const string myname = "AdcThresholdSignalFinder::update: ";
  DataMap ret;
  AdcIndex nsam = acd.samples.size();
  if ( nsam == 0 ) {
    cout << myname << "ERROR: No samples found in channel " << acd.channel() << endl;
    acd.signal.clear();
    acd.rois.clear();
    return ret.setStatus(1);
  }
  if ( m_prdtool != nullptr ) {
    RunData rdat = m_prdtool->runData(acd.run());
    if ( ! rdat.isValid() ) cout << myname << "WARNING: RunData not found." << endl;
    else rdat.setFormulaPars(*m_Threshold);
  }
  if ( ! m_Threshold->ready() ) {
    cout << myname << "WARNING: Using default scale factor " << m_Threshold->defaultEval() << endl;
  } 
  float thr = m_Threshold->eval();
  if ( m_LogLevel >= 2 ) cout << myname << "Finding ROIs for channel " << acd.channel() << endl;
  AdcIndex nsamlo = m_BinsBefore;
  AdcIndex nsamhi = m_BinsAfter;
  acd.signal.clear();
  acd.signal.resize(nsam, false);
  AdcIndex nbinAbove = 0;
  AdcIndex isamUnknown = 0;  // First sample which is not known to be inside or ouside a ROI.
  for ( AdcIndex isam=0; isam<nsam; ++isam ) {
    float val = acd.samples[isam];
    bool keep = ( m_FlagPositive && val >  thr ) ||
                ( m_FlagNegative && val < -thr );
    if ( keep ) {
      ++nbinAbove;
      AdcIndex jsam1 = isam > nsamlo ? isam - nsamlo : 0;
      if ( jsam1 < isamUnknown ) jsam1 = isamUnknown;
      AdcIndex jsam2 = isam + nsamhi + 1;
      if ( jsam2 > nsam ) jsam2 = nsam;
      if ( m_LogLevel >= 4 ) cout << myname << "Trigger: sample[" << isam
                                  << "] = " << val << " ==> range: ["
                                  << jsam1 << ", " << jsam2 << ")" << endl;
      for ( AdcIndex jsam=jsam1; jsam<jsam2; ++jsam ) acd.signal[jsam] = true;
      isamUnknown = jsam2;
    }
  }
  acd.roisFromSignal();
  if ( m_LogLevel >= 3 ) {
    cout << myname << "  # ticks above threshold: " << nbinAbove << endl;
    cout << myname << "  # ROI: " << acd.rois.size() << endl;
  }
  ret.setFloat("threshold", thr);
  ret.setInt("nThresholdBins", nbinAbove);
  ret.setInt("nroi", acd.rois.size());
  return ret;
}

//**********************************************************************

DataMap AdcThresholdSignalFinder::view(const AdcChannelData& acd) const {
  AdcChannelData acdtmp;
  acdtmp.samples = acd.samples;
  return update(acdtmp);
}

//**********************************************************************
