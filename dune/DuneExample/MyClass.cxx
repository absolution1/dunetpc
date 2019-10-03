// MyClass.cxx

#include "MyClass.h"
#include <iostream>
#include <string>
#include "TBuffer.h"
#include "TClass.h"

using std::cout;
using std::endl;
using std::string;

MyClass::MyClass() {
  cout << "Creating object." << endl;
}

MyClass::~MyClass() {
  cout << "Deleting object." << endl;
  delete myobj;
}

void MyClass::Streamer(TBuffer& buf) {
  const string myname = "MyClass: ";
  //TClass* pclass = MyClass::Class();  // With TObject?
  TClass* pclass = TClass::GetClass("MyClass");
  if ( buf.IsReading() ) {
    if ( pclass == nullptr ) {
      cout << myname << "Dictionary not found for read." << endl;
      return;
    }
    cout << myname << "Reading object." << endl;
    pclass->ReadBuffer(buf, this);
    ++readCount;
  } else {
    if ( pclass == nullptr ) {
      cout << myname << "Dictionary not found for write." << endl;
      return;
    }
    cout << myname << "Writing object." << endl;
    pclass->WriteBuffer(buf, this);
  }
}
