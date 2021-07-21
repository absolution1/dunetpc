// test_RootParFormula.cxx

// David Adams
// July 2021.
//
// This is a test and demonstration for RootParFormula.

#undef NDEBUG

#include "../RootParFormula.h"
#include <string>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <vector>
#include "TH1F.h"
#include "TCanvas.h"

using std::string;
using std::cout;
using std::endl;
using std::setw;
using std::vector;

using Index = unsigned int;
using DVec = std::vector<double>;

//**********************************************************************

bool floatcheck(double x1, double x2) {
  return fabs(x1-x2) < 1.e-4;
}

//**********************************************************************

int test_RootParFormula() {
  const string myname = "test_RootParFormula: ";
  cout << myname << "Starting test" << endl;
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";
  string scfg;

  cout << myname << line << endl;
  cout << myname << "Create formula and check." << endl;
  string snam = "myform";
  string sform = "[offset] + [scale]*x";
  RootParFormula rpf(snam, sform);
  cout << myname << "Name: " << rpf.name() << endl;
  assert( rpf.name() == snam );
  cout << myname << "Formula: " << rpf.formulaString() << endl;
  assert( rpf.formulaString() == sform );
  cout << myname << "# params: " << rpf.npar() << endl;
  assert( rpf.npar() == 2 );
  cout << myname << "Params: [";
  bool first = true;
  for ( string nam : rpf.pars() ) {
    if ( first ) first = false;
    else cout << ", ";
    cout << nam;
  }
  cout << "]" << endl;
  assert( rpf.pars().size() == 2 );
  assert( rpf.pars()[0] == "offset" );
  assert( rpf.pars()[1] == "scale" );
  cout << myname << "# vars: " << rpf.nvar() << endl;
  assert( rpf.nvar() == 1 );
  assert( rpf.setPars().size() == 0 );
  assert( rpf.unsetPars().size() == 2 );
  assert( rpf.resetPars().size() == 0 );
  assert( rpf.pars().size() == 2 );
  assert( ! rpf.ready() );

  cout << myname << line << endl;
  cout << myname << "Set parameters." << endl;
  assert( rpf.unsetPars().size() == 2 );
  assert( rpf.resetPars().size() == 0 );
  assert( ! rpf.ready() );
  rpf.setParValue("offset", 10.0);
  assert( rpf.unsetPars().size() == 1 );
  assert( rpf.resetPars().size() == 0 );
  assert( ! rpf.ready() );
  rpf.setParValue("scale", 3.0);
  assert( rpf.unsetPars().size() == 0 );
  assert( rpf.resetPars().size() == 0 );
  assert( rpf.ready() );
  rpf.setParValue("scale", 2.0);
  assert( rpf.unsetPars().size() == 0 );
  assert( rpf.resetPars().size() == 1 );
  assert( rpf.ready() );

  cout << myname << line << endl;
  cout << myname << "Check eval." << endl;
  DVec vals = { 0.0, 1.0, 2.0 };
  DVec evcs = { 10.0, 12.0, 14.0 };
  for ( Index ival=0; ival<vals.size(); ++ival ) {
    double eval = rpf.eval(vals[ival]);
    cout << myname << "  " << vals[ival] << ": " << eval << endl;
    assert( floatcheck(eval, evcs[ival]) );
    DVec vars = { vals[ival] };
    double eval2 = rpf.eval(vars);
    assert( floatcheck(eval2, evcs[ival]) );
  }

  cout << myname << line << endl;
  cout << myname << "Set default eval." << endl;
  double def = rpf.defaultEval();
  assert( floatcheck(def, 0.0 ) );
  assert( rpf.setDefaultEval(99.0) == 0 );
  def = rpf.defaultEval();
  assert( floatcheck(def, 99.0) );

  cout << myname << line << endl;
  cout << myname << "Reset pars." << endl;
  assert( rpf.ready() );
  assert( floatcheck(rpf.eval(0.0), 10.0) );
  assert( rpf.unsetParValues() == 0 );
  assert( ! rpf.ready() );
  assert( floatcheck(rpf.eval(0.0), 99.0) );
  rpf.setParValue("offset", 10.0);
  assert( ! rpf.ready() );
  rpf.setParValue("scale", 3.0);
  assert( rpf.ready() );
  assert( floatcheck(rpf.eval(0.0), 10.0) );

  cout << myname << line << endl;
  cout << myname << "Done." << endl;

  return 0;
}

//**********************************************************************

int main() {
  return test_RootParFormula();
}

//**********************************************************************
