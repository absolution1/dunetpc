// IndexLabelTool.h

// David Adams
// June 2018
//
// Interface for tools maps and index to a string.

#ifndef IndexLabelTool_H
#define IndexLabelTool_H

class IndexLabelTool {

public:

  using Index = unsigned int;
  using Name = std::string;

  virtual ~IndexLabelTool() =default;

  virtual Name get(Index idx) const =0;

};

#endif
