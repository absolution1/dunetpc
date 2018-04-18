// ITestTool.h

// David Adams
// April 2018
//
// Example tool interface for testing.

#ifndef ITestTool_H
#define ITestTool_H

#include <string>

class ITestTool {

public:

  virtual ~ITestTool() =default;

  virtual std::string mytype() const { return "ITestTool"; }

  virtual std::string label() =0;

private:

};

#endif
