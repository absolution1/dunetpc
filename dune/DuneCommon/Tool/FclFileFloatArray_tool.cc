// FclFileFloatArray_tool.cc

#include "FclFileFloatArray.h"
#include "fhiclcpp/intermediate_table.h"
#include "fhiclcpp/make_ParameterSet.h"
#include "art/Utilities/make_tool.h"
#include <iostream>
#include "TSystem.h"

using std::cout;
using std::endl;

//**********************************************************************

FclFileFloatArray::FclFileFloatArray(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<Index>("LogLevel")),
  m_FileNames(ps.get<NameVector>("FileNames")),
  m_DefaultValue(ps.get<float>("DefaultValue")) {
  const Name myname = "FclFileFloatArray::ctor: ";
  if ( m_LogLevel ) {
    cout << "Configuration:" << endl;
    cout << "      LogLevel: " << m_LogLevel << endl;
    cout << "     FileNames: ";
    if ( m_FileNames.size() == 0 ) {
      cout << "<empty>";
    }
    for ( Name fname : m_FileNames ) {
      cout << "\n       " << fname;
    }
    cout << endl;
    cout << "  DefaultValue: " << m_DefaultValue << endl;
  }
  init();
}

//**********************************************************************

void FclFileFloatArray::init() {
  const Name myname = "FclFileFloatArray::init: ";
  bool first = true;
  Name fclPath = gSystem->Getenv("FHICL_FILE_PATH");
  for ( Name fname : m_FileNames ) {
    TString ts(fname.c_str());
    gSystem->FindFile(fclPath.c_str(), ts);
    Name pfname = ts.Data();
    if ( pfname.size() == 0 ) {
      cout << myname << "WARNING: Unable to find file " << fname << endl;
      continue;
    }
    cet::filepath_maker policy;
    auto ps = fhicl::ParameterSet::make(fhicl::parse_document(pfname, policy));
    auto pfat = art::make_tool<FloatArrayTool>(ps);
    if ( m_LogLevel >= 2 ) cout << myname << "Reading tool " << pfat->label() << endl;
    if ( first ) {
      m_unit = pfat->unit();
      first = false;
    } else {
      if ( pfat->unit() != unit() ) {
        cout << myname << "WARNING: Inconsistent units: " << pfat->unit() << " != " << m_unit << endl;
        cout << myname << "Skipping tool at " << fname << endl;
        continue;
      }
    }
    pfat->fill(m_vals, defaultValue());
  }
}

//**********************************************************************

DEFINE_ART_CLASS_TOOL(FclFileFloatArray)
