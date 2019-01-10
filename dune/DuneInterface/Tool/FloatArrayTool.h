// FloatArrayTool.h
//
// David Adams
// January 2019
//
// Interface for a tool providing access to an array
// of floating values. It also holds an offset, i.e.
// the index of the first value, and a label.

#ifndef FloatArrayTool_H
#define FloatArrayTool_H

#include "dune/DuneInterface/Data/RunData.h"
#include <vector>

class FloatArrayTool {

public:

  using Index = unsigned int;
  using FloatVector = std::vector<float>;

  virtual ~FloatArrayTool() =default;

  // Return the offset.
  virtual Index offset() const =0;

  // Return the label.
  virtual std::string label() const =0;

  // Return the values without offset.
  virtual const FloatVector& values() const =0;

  // Return the number of values.
  virtual Index size() const { return values().size(); }

  // Insert the values into a supplied vector.
  // Values are place at the offset position.
  // If needed, the size of the input vector is increased.
  // Existing values may be overwritten.
  // New values below the offset are set to valdef.
  // Returns 0 for success.
  virtual int fill(FloatVector& vals, float valdef =0.0) const {
    Index len = offset() + size();
    if ( len > vals.size() ) vals.resize(len, valdef);
    std::copy(values().begin(), values().end(), vals.begin() + offset());
    return 0;
  }

};

#endif
