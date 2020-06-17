// IndexVectorMapTool.h

// David Adams
// June 2018
// June 2020: Change get to return by balue instead of ref.
//
// Interface for tools that map one index to a vector of indices.

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
  virtual IndexVector get(Index idx) const =0;

  // Return if the vector for an index contains a value.
  virtual bool contains(Index idx, Index val) {
    IndexVector vec = get(idx);
    return find(vec.begin(), vec.end(), val) != vec.end();
  }

};

#endif
