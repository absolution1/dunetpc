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
#include <vector>
#include "TFile.h"
#include "TH1F.h"

using std::string;
using std::cout;
using std::endl;
using std::ofstream;
using std::vector;

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
  MyClass* pobj1 = new MyClass;
  MyClass& obj1 = *pobj1;
  obj1.myint = 123;
  obj1.myfloat = 246.78;
  TH1F h0("histo", "histo", 10, 0, 10);
  h0.SetDirectory(nullptr);
  for ( int i=0; i<5; ++i ) h0.Fill(5+0.3*i);
  obj1.myobj = h0.Clone("myhisto");
  TH1* ph1 = dynamic_cast<TH1*>(obj1.myobj);
  assert( ph1 != nullptr );

  string rfnam = "test_MyClass.root";
  cout << myname << line << endl;
  cout << myname << "Write object to " << rfnam << "." << endl;
  TFile* prout = TFile::Open(rfnam.c_str(), "RECREATE");
  assert( obj1.readCount == 0 );
  prout->WriteObject(&obj1, "myobj1");
  assert( obj1.readCount == 0 );
  prout->Write();
  delete prout;

  cout << myname << line << endl;
  cout << myname << "Read object from " << rfnam << "." << endl;
  TFile* prin = TFile::Open(rfnam.c_str(), "READ");
  MyClass* pobj2 = nullptr;
  prin->GetObject("myobj1", pobj2);

  vector<string> msgs = {"file open", "file closed" };
  for ( string msg : msgs ) {
    cout << myname << line << endl;
    cout << myname << "Check read object with input " << msg << endl;
    cout << "Read count: " << pobj2->readCount << endl;
    assert( pobj2->myint == obj1.myint );
    assert( pobj2->myfloat == obj1.myfloat );
    assert( pobj2->readCount == 1 );
    cout << "TObject pointer: " << pobj2->myobj << endl;
    assert( pobj2->myobj != nullptr );
    pobj2->myobj->Print();
    TH1* ph2 = dynamic_cast<TH1*>(pobj2->myobj);
    assert( ph2 != nullptr );
    assert( ph2 != ph1 );
    assert( ph2->GetEntries() == ph1->GetEntries() );
    if ( prin != nullptr ) {
      cout << myname << "Closing input file." << endl;
      delete prin;
      prin = 0;
    }
  }

  cout << myname << line << endl;
  cout << myname << "Deleting read object." << endl;
  delete pobj2;

  cout << myname << line << endl;
  cout << myname << "Deleting original object." << endl;
  delete pobj1;

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
