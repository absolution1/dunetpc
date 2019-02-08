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

#include <vector>

class FloatArrayTool {

public:

  using Index = unsigned int;
  using FloatVector = std::vector<float>;

  virtual ~FloatArrayTool() =default;

  // Return the offset.
  virtual Index offset() const { return 0; }

  // Return the label.
  virtual std::string label() const { return ""; }

  // Return the units for the values.
  virtual std::string unit() const { return ""; }

  // Return the values without offset.
  virtual const FloatVector& values() const =0;

  // Return the number of values.
  virtual Index size() const { return values().size(); }

  // Return the default value, i.e. the value returned when
  // a caller requests a value that is out of range.
  virtual float defaultValue() const { return 0.0; }

  // Return if an index is in range.
  virtual bool inRange(Index ival) {
    if ( ival < offset() ) return false;
    if ( ival + offset() >= size() ) return false;
    return true;
  }

  // Return the value for an (offset) index.
  virtual float value(Index ival) {
    if ( ! inRange(ival) ) return defaultValue();
    return values()[ival + offset()];
  }

  // Insert the values into a supplied vector.
  // Values are place at the offset position.
  // If needed, the size of the input vector is increased.
  // Existing values may be overwritten.
  // New values below the offset are set to valdef.
  // Returns 0 for success.
  virtual int fill(FloatVector& vals, float valdef) const {
    Index len = offset() + size();
    if ( len > vals.size() ) vals.resize(len, valdef);
    std::copy(values().begin(), values().end(), vals.begin() + offset());
    return 0;
  }

};

#endif
