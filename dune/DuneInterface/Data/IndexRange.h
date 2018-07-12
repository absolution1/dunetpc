// IndexRange.h
//
// David Adams
// July 2018
//
// Data structure that holds a named and labeled contiguous
// range of indices. The intention is that the name is unique
// (in some context) and can be used as an identifier and to
// construct other (e.g. histogram) names. The label is expected
// to be used to construct titles and labels for histograms,
// plots, etc.

#ifndef IndexRange_H
#define IndexRange_H

#include <string>

class IndexRange {

public:

  using Index = unsigned int;
  using Name = std::string;

  // Data.
  Name name;         // Name
  Name label;        // Label
  Index begin =0;    // First index
  Index end =0;      // Last + 1 index

  // Ctors.
  IndexRange() { };
  IndexRange(Index len) : end(len) { }
  IndexRange(Index a_begin, Index a_end) : begin(a_begin), end(a_end) { }

  // Length of the range.
  Index size() const { return end>begin ? end - begin : 0; }

  // Other helpers.
  bool isValid() const { return end > begin; }
  Index first() const { return begin; }
  Index last() const { return end > 0 ? end - 1 : 0; }
  
};

#endif
