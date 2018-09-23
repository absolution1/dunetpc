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
//
// An arbitrarry number of labesl (op to four  in the ctor) may
// be assigned to each range.

#ifndef IndexRange_H
#define IndexRange_H

#include <string>

class IndexRange {

public:

  using Index = unsigned int;
  using Name = std::string;
  using NameVector = std::vector<Name>;

  // Data.
  Name name;         // Name
  NameVector labels; // Label
  Index begin =0;    // First index
  Index end =0;      // Last + 1 index

  // Ctors.
  IndexRange() { };
  IndexRange(Index len) : end(len) { }
  IndexRange(Index a_begin, Index a_end) : begin(a_begin), end(a_end) { }
  IndexRange(Name a_name, Index a_begin, Index a_end,
             Name lab0="", Name lab1="", Name lab2="", Name lab3="")
  : name(a_name), begin(a_begin), end(a_end) {
    if ( lab0.size() ) setLabel(0, lab0);
    if ( lab1.size() ) setLabel(1, lab1);
    if ( lab2.size() ) setLabel(2, lab2);
    if ( lab3.size() ) setLabel(3, lab3);
  }

  // Length of the range.
  Index size() const { return end>begin ? end - begin : 0; }

  // Other helpers.
  bool isValid() const { return end > begin; }
  Index first() const { return begin; }
  Index last() const { return end > 0 ? end - 1 : 0; }
  std::string rangeString() const {
    if ( ! isValid() ) return "[]";
    return "[" + std::to_string(begin) + ", " + std::to_string(end) + ")";
  }
  void setLabel(Name lab) { setLabel(0, lab); }
  void setLabel(Index ilab, Name lab) {
    labels.resize(ilab+1, "");
    labels[ilab] = lab;
  }
  Name label(Index ilab =0) const {
    if ( ilab >= labels.size() ) return "";
    return labels[ilab];
  }
  bool contains(Index ival) const { return ival >= begin && ival < end; }
  
};

#endif
