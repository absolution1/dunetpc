// IndexVectorMapTool.h

// David Adams
// June 2018
//
// Interface for tools that map one index to another, e.g.
// online to offline channel or vice versa.

#ifndef IndexVectorMapTool_H
#define IndexVectorMapTool_H

#include <vector>
#include <algorithm>

class IndexVectorMapTool {

public:

  using Index = unsigned int;
  using IndexVector = std::vector<Index>;

  static Index badIndex() { return -1; }

  virtual ~IndexVectorMapTool() =default;

  // Return the vector for an index.
  virtual const IndexVector& get(Index idx) const =0;

  // Return if the vector for an index contains a value.
  virtual bool contains(Index idx, Index val) {
    const IndexVector vec = get(idx);
    return find(vec.begin(), vec.end(), val) != vec.end();
  }

};

#endif
