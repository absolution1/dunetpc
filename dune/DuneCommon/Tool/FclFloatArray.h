// FclFloatArray.h
//
// David Adams
// January 2019
//
// Tool that returns a float array stored as fcl.
// 
// Parameters:
//   LogLevel - Message logging level:
//              0: No messages
//              1: Show config in ctor with a few values
//              2: Show config in ctor with all values
//   DefaultValue - Value to return when an index is out of range.
//   Offset - Index of the first value
//   Label - string label
//   Values - Array values.

#ifndef FclFloatArray_H
#define FclFloatArray_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/FloatArrayTool.h"
#include <vector>

class FclFloatArray : public FloatArrayTool {

public:

  using Name = std::string;

  // Ctor.
  FclFloatArray(fhicl::ParameterSet const& ps);

  // Dtor.
  ~FclFloatArray() override =default;

  // Return the offset.
  Index offset() const override { return m_Offset; }
  
  // Return the label.
  std::string label() const override { return m_Label; }
 
  // Return the units.
  Name unit() const override { return m_Unit; }

  // Return the values without offset.
  const FloatVector& values() const override { return m_Values; }

  // Retrurn the default value.
  float defaultValue() const override { return m_DefaultValue; }

private:

  // Parameters.
  Index m_LogLevel;
  float m_DefaultValue;
  Index m_Offset;
  Name m_Label;
  Name m_Unit;
  FloatVector m_Values;


};


#endif
