// test_MyClass.cxx

// David Adams
// September 2019
//
// This is a test and demonstration for MyClass.

#undef NDEBUG

#include "../MyClass.h"
#include <string>
#include <iostream>
#include <fstream>
#include <cassert>
#include "TFile.h"

using std::string;
using std::cout;
using std::endl;
using std::ofstream;

//**********************************************************************

int test_MyClass() {
  const string myname = "test_MyClass: ";
  cout << myname << "Starting test" << endl;
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";
  string scfg;

  cout << myname << line << endl;
  cout << myname << "Create object." << endl;
  MyClass obj1;
  obj1.myint = 123;
  obj1.myfloat = 246.78;

  string rfnam = "test_MyClass.root";
  cout << myname << line << endl;
  cout << myname << "Write object to " << rfnam << "." << endl;
  TFile* prout = TFile::Open(rfnam.c_str(), "RECREATE");
  prout->WriteObject(&obj1, "myobj1");
  prout->Write();
  delete prout;

  cout << myname << line << endl;
  cout << myname << "Read object from " << rfnam << "." << endl;
  TFile* prin = TFile::Open(rfnam.c_str(), "READ");
  MyClass* pobj = nullptr;
  prin->GetObject("myobj1", pobj);
  delete prin;

  cout << myname << line << endl;
  cout << myname << "Check read object." << endl;
  assert( pobj->myint == obj1.myint );
  assert( pobj->myfloat == obj1.myfloat );
  delete pobj;

  cout << myname << line << endl;
  string snam = "test_MyClass.C";
  cout << myname << "Creating test script " << snam << endl;
  ofstream fout(snam.c_str());
  fout << "void test_MyClass() {" << endl;
  fout << "  TFile* pfin = TFile::Open(\"" << rfnam << "\");" << endl;
  fout << "  MyClass* pobj;" << endl;
  fout << "  pfin->GetObject(\"myobj1\", pobj);" << endl;
  fout << "  cout << \"  myint = \" << pobj->myint << endl;" << endl;
  fout << "  cout << \"  myfloat = \" << pobj->myfloat << endl;" << endl;
  fout << "  delete pobj;" << endl;
  fout << "}" << endl;
  cout << "Check root file with \"root.exe " << snam << "\"" << endl;

  cout << myname << line << endl;
  cout << myname << "Done." << endl;
  return 0;
}

//**********************************************************************

int main() {
  return test_MyClass();
}

//**********************************************************************
