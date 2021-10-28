// FclFileFloatArray.h
//
// David Adams
// May 2018
//
// Fcl-based tool to give access to build a float array from an array of
// fcl files, each  specifying a FloatArrayTool.
//
// Parameters:
//   LogLevel - Message logging level (0=none, 1=ctor, 2=each call, ...)
//   FileNames - Vector of file name patterns
//   DefaultValue - Value assigned to values that are out of range

#ifndef FclFileFloatArray_H
#define FclFileFloatArray_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/FloatArrayTool.h"
#include <vector>

class FclFileFloatArray : public FloatArrayTool {

public:

  using Name = std::string;
  using NameVector = std::vector<Name>;

  // Ctor.
  FclFileFloatArray(fhicl::ParameterSet const& ps);

  // Dtor.
  ~FclFileFloatArray() override =default;

  // Return the size.
  Index size() const override { return m_vals.size(); }

  // Return the values.
  const FloatVector& values() const override { return m_vals; }

  // Return the default value.
  float defaultValue() const override { return m_DefaultValue; }

  // Return the units.
  Name unit() const override { return m_unit; }

private:

  // Parameters.
  Index m_LogLevel;
  NameVector m_FileNames;
  float m_DefaultValue;

  Name m_unit;
  FloatVector m_vals;

  // Read the values.
  void init();

};


#endif
