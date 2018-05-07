// test_make_tool.cxx

#include "dune/ArtSupport/Tool/TestTool.h"
#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/intermediate_table.h"
#include "fhiclcpp/make_ParameterSet.h"
#include "art/Utilities/make_tool.h"
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>

using std::string;
using std::cout;
using std::endl;
using std::ofstream;

#undef NDEBUG
#include <cassert>

//**********************************************************************

using std::string;
using std::cout;
using std::endl;
using std::ofstream;
using std::ostream;
using std::ostringstream;
using std::setw;
using std::fixed;
using std::setprecision;

using Index = unsigned int;


//**********************************************************************

class INotTestTool {

public:

  virtual ~INotTestTool() =default;

  virtual std::string mytype() { return "INotTestTool"; }

  virtual std::string label() { return "BadLabel"; }

  virtual std::string extraLabel() { return "BadExtraLabel"; }

private:

};


//**********************************************************************

int test_make_tool(bool doCrash, bool useExistingFcl) {
  const string myname = "test_make_tool: ";
  cout << myname << "Starting test" << endl;
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";
  string scfg;

  cout << myname << line << endl;
  string fclfile = "test_make_tool.fcl";
  if ( ! useExistingFcl ) {
    cout << myname << "Creating top-level FCL." << endl;
    ofstream fout(fclfile.c_str());
    fout << "tools: {" << endl;
    fout << "  mytool1: {" << endl;
    fout << "    tool_type: TestTool" << endl;
    fout << "    Label: \"Tool 1\"" << endl;
    fout << "  }" << endl;
    fout << "  mytool2: {" << endl;
    fout << "    tool_type: TestTool" << endl;
    fout << "    Label: \"Tool 2\"" << endl;
    fout << "  }" << endl;
    fout << "}" << endl;
    fout.close();
  } else {
    cout << myname << "Using existing top-level FCL." << endl;
  }

  cet::filepath_lookup policy("FHICL_FILE_PATH");
  fhicl::intermediate_table tbl;
  fhicl::parse_document(fclfile, policy, tbl);
  fhicl::ParameterSet psTop;
  fhicl::make_ParameterSet(tbl, psTop);
  fhicl::ParameterSet pstools = psTop.get<fhicl::ParameterSet>("tools");
  fhicl::ParameterSet pstool1 = pstools.get<fhicl::ParameterSet>("mytool1");
  fhicl::ParameterSet pstool2 = pstools.get<fhicl::ParameterSet>("mytool2");

  cout << myname << line << endl;
  cout << "Retrieve with the correct interface type." << endl;
  std::unique_ptr<ITestTool> pt1 = art::make_tool<ITestTool>(pstool1);
  cout << myname << "Tool pointer: " << pt1.get() << endl;
  cout << myname << "Tool type: " << pt1->mytype() << endl;
  assert( pt1->mytype() == "TestTool" );
  cout << myname << "Tool label: " << pt1->label() << endl;
  assert( pt1->label() == "Tool 1" );

  cout << myname << line << endl;
  cout << "Retrieve with the incorrect interface type." << endl;
  bool exceptionRaised = false;
  bool nullReturn = false;
  std::unique_ptr<INotTestTool> pt2;
  try {
    pt2 = art::make_tool<INotTestTool>(pstool2);
  } catch (...) {
    exceptionRaised = true;
    nullReturn = pt2 == nullptr;
    pt2.reset(nullptr);
  }
  if ( pt2 ) {
    cout << myname << "Tool retrieval succeeded. This is not good." << endl;
    cout << myname << "Tool pointer: " << pt2.get() << endl;
    //assert( pt2.get() == nullptr );
    cout << myname << "Tool type: " << pt2->mytype() << endl;
    //assert( pt2->mytype() == "BadTestTool" );
    cout << myname << "Tool label: " << pt2->label() << endl;
    assert( pt2->label() == "Tool 2" );
    if ( doCrash )  {
      cout << myname << "Trying to fetch extra label..." << endl;
      cout << myname << "Tool extra label: " << pt2->extraLabel() << endl;
    }
  } else {
    cout << myname << "Tool retrieval failed. As it should." << endl;
    if ( exceptionRaised ) cout << myname << "  Exception was raised." << endl;
    if ( nullReturn ) cout << myname << "  Tool pointer was null." << endl;
  }

  cout << myname << line << endl;
  cout << myname << "Done." << endl;
  return 0;
}

//**********************************************************************

int main(int narg, char** carg) {
  bool doCrash = narg > 1;
  bool keepFcl = narg > 2;
  int rstat = test_make_tool(doCrash, keepFcl);
  return rstat;
}
