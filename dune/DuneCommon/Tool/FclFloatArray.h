// FclFloatArray.h
//
// David Adams
// January 2019
//
// Tool that returns a float array stored in a fcl.
// 
// Parameters:
//   LogLevel - Message logging level:
//              0: No messages
//              1: Show config in ctor with a few values
//              2: Show config in ctor with all values
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
 
  // Return the values without offset.
  const FloatVector& values() const override { return m_Values; }

private:

  // Parameters.
  Index m_LogLevel;
  Index m_Offset;
  std::string m_Label;
  FloatVector m_Values;

};

DEFINE_ART_CLASS_TOOL(FclFloatArray)

#endif
