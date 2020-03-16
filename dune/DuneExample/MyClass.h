// MyClass.h
//
// David Adams
// September 2019
//
// Example Root class.

#ifndef MyClass_H
#define MyClass_H

#include "TObject.h"

class MyClass {

public:

  int myint =0;
  float myfloat = 0.0;
  TObject* myobj = nullptr;
  unsigned int readCount =0;

  MyClass();
  ~MyClass();
 
  // Custom streamer.
  void Streamer(TBuffer& rbuf);

};

#endif
