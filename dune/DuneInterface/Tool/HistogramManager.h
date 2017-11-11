// HistogramManager.h
//
// David Adams
// August 2017
//
// Interface for a tool that manages histograms.
//
// Ih histogram has an associated directory, then the tool
// should unset that association and take mangement, i.e. ensure
// the histogram is deleted at the appropriate time.

#ifndef HistogramManager_H
#define HistogramManager_H

#include <string>
#include <vector>

class TH1;

class HistogramManager {

public:

  using Name = std::string;
  using NameVector = std::vector<Name>;

  virtual ~HistogramManager() =default;

  // Take management of a histogram.
  virtual int manage(TH1* ph) =0;

  // Return the names of histograms managed here.
  virtual NameVector names() const =0;

  // Return a histogram.
  virtual TH1* get(Name hname) const =0;

  // Return the last histogram managed here.
  virtual TH1* last() const =0;

  // Delete the named histogram.
  // Defult argument "* releases all histograms.
  virtual int release(Name hname ="*") =0;

};

#endif
