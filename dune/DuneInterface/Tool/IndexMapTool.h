// IndexMapTool.h

// David Adams
// June 2018
//
// Interface for tools that map one index to another, e.g.
// online to offline channel or vice versa.

#ifndef IndexMapTool_H
#define IndexMapTool_H

class IndexMapTool {

public:

  using Index = unsigned int;

  static Index badIndex() { return -1; }

  virtual ~IndexMapTool() =default;

  virtual Index get(Index idx) const =0;

};

#endif
