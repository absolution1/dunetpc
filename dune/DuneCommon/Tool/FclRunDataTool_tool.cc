// FclRunDataTool_tool.cc

#include "FclRunDataTool.h"
#include "dune/DuneCommon/StringManipulator.h"
#include "fhiclcpp/intermediate_table.h"
#include "fhiclcpp/make_ParameterSet.h"
#include "TString.h"
#include "TSystem.h"
#include <iostream>

using std::cout;
using std::endl;
using std::string;

//**********************************************************************

namespace {

using Name = FclRunDataTool::Name;
using Index = RunData::Index;
using IndexVector = RunData::IndexVector;

int parseFcl(string path, string fclname, RunData& rdat) {
  TString ts(fclname.c_str());
  gSystem->FindFile(path.c_str(), ts);
  string pfname = ts.Data();
  if ( pfname.size() == 0 ) return 1;
  fhicl::intermediate_table tbl;
  cet::filepath_maker policy;
  fhicl::parse_document(pfname, policy, tbl);
  fhicl::ParameterSet ps;
  fhicl::make_ParameterSet(tbl, ps);
  ps.get_if_present<Index>("run", rdat.accessRun());
  ps.get_if_present<string>("cryostat", rdat.accessCryostat());
  ps.get_if_present<IndexVector>("apas", rdat.accessApas());
  ps.get_if_present<float>("gain", rdat.accessGain());
  ps.get_if_present<float>("shaping", rdat.accessShaping());
  ps.get_if_present<float>("leakage", rdat.accessLeakage());
  ps.get_if_present<float>("hvfrac", rdat.accessHvfrac());
  ps.get_if_present<Index>("pulserAmplitude", rdat.accessPulserAmplitude());
  ps.get_if_present<Index>("pulserSource", rdat.accessPulserSource());
  ps.get_if_present<Index>("pulserPeriod", rdat.accessPulserPeriod());
  ps.get_if_present<Name>("phaseGroup", rdat.accessPhaseGroup());
  ps.get_if_present<IndexVector>("phases", rdat.accessPhases());
  return 0;
}

}  // end unnamed namespace

//**********************************************************************

FclRunDataTool::FclRunDataTool(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<Index>("LogLevel")),
  m_FileNames(ps.get<NameVector>("FileNames")) {
  const Name myname = "FclRunDataTool::ctor: ";
  m_fclPath = gSystem->Getenv("FHICL_FILE_PATH");
  if ( m_LogLevel ) {
    cout << "Configuration:" << endl;
    cout << "   LogLevel: " << m_LogLevel << endl;
    cout << "  FileNames: ";
    if ( m_FileNames.size() == 0 ) {
      cout << "<empty>";
    }
    for ( Name fname : m_FileNames ) {
      cout << "\n    " << fname;
    }
    cout << endl;
  }
}

//**********************************************************************

RunData FclRunDataTool::runData(Index run, Index subRun) const {
  const Name myname = "FclRunDataTool::runData: ";
  RunData rdat;
  bool useSubRun = subRun;
  if ( m_LogLevel >= 2 ) {
    cout << myname << "Fetching tool for run " << run;
    if ( useSubRun ) cout << ", subrun " << subRun;
    cout << endl;
  }
  for ( Name fname : m_FileNames ) {
    StringManipulator sman(fname);
    sman.replaceFixedWidth("%RUN%", run, 6);
    sman.replaceFixedWidth("%SUBRUN%", subRun, 6);
    if ( parseFcl(m_fclPath, fname, rdat) ) {
      if ( m_LogLevel >= 3 ) cout << "Unable to find/read " << fname << endl;
    } else {
      if ( m_LogLevel >= 3 ) cout << "  Read " << fname << endl;
    }
  }
  return rdat;
}

//**********************************************************************
