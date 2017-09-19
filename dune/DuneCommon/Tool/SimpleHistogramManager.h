// SimpleHistogramManager.h
//
// David Adams
// August 2017
//
// Tool that provides access to and management of histograms.

#ifndef SimpleHistogramManager_H
#define SimpleHistogramManager_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/HistogramManager.h"
#include <map>

class SimpleHistogramManager : public HistogramManager {

public:

  // Ctor from FCL.
  SimpleHistogramManager(fhicl::ParameterSet const& ps);

  // Dtor.
  ~SimpleHistogramManager() override;

  int manage(TH1* ph) override;

  NameVector names() const override;

  TH1* get(Name hname) const override;

  TH1* last() const override;

  int release(Name hname) override;

private:

  // Parameters.
  int m_LogLevel;

  // Local data.
  NameVector m_names;
  std::map<Name, TH1*> m_hists;

};

DEFINE_ART_CLASS_TOOL(SimpleHistogramManager)

#endif
