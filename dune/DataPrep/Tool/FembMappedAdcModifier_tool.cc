// FembMappedAdcModifier_tool.cc

#include "FembMappedAdcModifier.h"
#include "dune/ArtSupport/DuneToolManager.h"
#include "TSystem.h"
#include <iostream>
#include <sstream>

using std::string;
using std::cout;
using std::endl;
using std::ostringstream;

using Name = std::string;
using ToolPtr = std::unique_ptr<AdcChannelTool>;

//**********************************************************************

FembMappedAdcModifier::FembMappedAdcModifier(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")),
  m_ToolBase(ps.get<Name>("ToolBase")),
  m_DirName(ps.get<Name>("DirName")) {
  const string myname = "FembMappedAdcModifier::ctor: ";
  if ( m_LogLevel >= 1 ) {
    cout << myname << "      Log level: " << m_LogLevel << endl;
    cout << myname << "      Tool base: " << m_ToolBase << endl;
    cout << myname << "        Fcl dir: " << m_DirName << endl;
  }
}

//**********************************************************************

DataMap FembMappedAdcModifier::view(const AdcChannelData& acd) const {
  DataMap result;
  AdcChannelData acdtmp(acd);
  return update(acdtmp);
}

//**********************************************************************

DataMap FembMappedAdcModifier::update(AdcChannelData& acd) const {
  const string myname = "FembMappedAdcModifier::update: ";
  DataMap res;
  AdcChannel ifmb = acd.fembID;
  ostringstream sstool;
  sstool << m_ToolBase;
  if  ( ifmb == AdcChannelData::badIndex ) sstool << "Default";
  else sstool << ifmb;
  Name toolName = sstool.str();
  DuneToolManager* pdtm = nullptr;
  std::unique_ptr<DuneToolManager> pdtmManaged;
  if ( m_DirName.size() ) {
    Name fclfile = m_DirName + "/" + toolName + ".fcl";
    if ( gSystem->AccessPathName(fclfile.c_str()) ) {
      cout << myname << "ERROR: No such file: " << fclfile << endl;
      return res.setStatus(2000);
    }
    pdtmManaged.reset(new DuneToolManager(fclfile));
    pdtm = pdtmManaged.get();
  } else {
    pdtm = DuneToolManager::instance();
  }
  ToolPtr ptool = pdtm->getPrivate<AdcChannelTool>(toolName);
  if ( ptool == nullptr ) {
    if ( m_LogLevel >= 2 ) {
      cout << myname << "ERROR: Unable find tool " << toolName << endl;
    }
    return res.setStatus(1000);
  }
  return ptool->update(acd);
}

//**********************************************************************
