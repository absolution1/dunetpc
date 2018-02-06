// TestTool.h

// David Adams
// August 2017
//
// Example tool for testing.
//
// Configuration:
//   Label:  String associated with tool.

#ifndef TestTool_H
#define TestTool_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/AdcChannelViewer.h"
#include <iostream>

class TestTool {

public:

  TestTool(fhicl::ParameterSet const& ps);

  ~TestTool();

  std::string label() { ++m_count; return m_Label; }

private:

  // Configuration data.
  std::string m_Label;

  // Call count.
  int m_count;

};

DEFINE_ART_CLASS_TOOL(TestTool)

#endif
