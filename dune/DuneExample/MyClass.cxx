// MyClass.cxx

#include "MyClass.h"
#include <iostream>

using std::cout;
using std::endl;

MyClass::MyClass() {
  cout << "Creating object." << endl;
}

MyClass::~MyClass() {
  cout << "Deleting object." << endl;
}
