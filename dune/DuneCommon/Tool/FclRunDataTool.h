// FclRunDataTool.h
//
// David Adams
// May 2018
//
// Fcl-based tool to find conditions data for a run.
// 
// A search is done along the usual fcl path for each of the files in FileNames
// afer making the following substituions in those names:
//   %RUN% --> run number
//   %SUBRUN% --> run number
// with the run and subrun sustitutions padded with leading zeros to attain
// a length of (at least) six.
//
// Parameters:
//   LogLevel - Message logging level (0=none, 1=ctor, 2=each call, ...)
//   FileNames - Vector of file name patterns

#ifndef FclRunDataTool_H
#define FclRunDataTool_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/RunDataTool.h"
#include <vector>

class FclRunDataTool : public RunDataTool {

public:

  using Name = std::string;
  using NameVector = std::vector<Name>;

  // Ctor.
  FclRunDataTool(fhicl::ParameterSet const& ps);

  // Dtor.
  ~FclRunDataTool() override =default;

  // Return run data.
  RunData runData(Index run, Index subRun) const override;

private:

  // Parameters.
  Index m_LogLevel;
  NameVector m_FileNames;

  Name m_fclPath;

};

DEFINE_ART_CLASS_TOOL(FclRunDataTool)

#endif
