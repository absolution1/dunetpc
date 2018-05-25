// StringManipulator.
//
// David Adams
// August 2017
//
// Class to manipulate a string.

#include <string>
#include <sstream>

#ifndef StringManipulator_H
#define StringManipulator_H

class StringManipulator {

public:

  using Index = unsigned int;

  // Fetch a fill value for an output stream as follows:
  //   '0' - zero or positive integers
  //   '-' - negative integers
  template<typename T>
  static char
  getFill(typename std::enable_if<std::is_integral<T>::value, const T&>::type val) {
    return val < 0 ? '-' : '0';
  }
  //   '_' - otherwise
  template<typename T>
  static char
  getFill(typename std::enable_if<!std::is_integral<T>::value, const T&>::type val) {
    return '_';
  }

  // Ctor.
  StringManipulator(std::string& strin) : m_str(strin) { }

  // Return the manipulated string.
  const std::string& string() { return m_str; }

  // Replace all occurences of a substring with a value.
  // Returns the number of replacements (<0 for error).
  template<typename T>
  int replace(std::string substr, const T& xsub);

  // Replace all occurences of a substring with a value
  // with fixed width.
  // Returns the number of replacements (<0 for error).
  template<typename T>
  int replaceFixedWidth(std::string substr, const T& xsub, Index width);

private:

  std::string& m_str;

};

//**********************************************************************

template<typename T>
int StringManipulator::replace(std::string ssubout, const T& xsubin) {
  std::string ssubin;
  bool havesub = false;
  std::string::size_type ipos = m_str.find(ssubout);
  int count = 0;
  while ( ipos != std::string::npos ) {
    if ( ! havesub ) {
      std::ostringstream sssubin;
      sssubin << xsubin;
      ssubin = sssubin.str();
    }
    std::string::size_type lout = ssubout.size();
    m_str.replace(ipos, lout, ssubin);
    ipos = m_str.find(ssubout, ipos+lout);
    ++count;
  }
  return count;
}

//**********************************************************************

template<typename T>
int StringManipulator::
replaceFixedWidth(std::string ssubout, const T& xsubin, Index width) {
  std::string ssubin;
  char cfill = getFill<T>(xsubin);
  bool havesub = false;
  std::string::size_type ipos = m_str.find(ssubout);
  int count = 0;
  while ( ipos != std::string::npos ) {
    if ( ! havesub ) {
      std::ostringstream sssubin;
      sssubin.fill(cfill);
      sssubin.width(width);
      sssubin << xsubin;
      ssubin = sssubin.str();
    }
    std::string::size_type lout = ssubout.size();
    m_str.replace(ipos, lout, ssubin);
    ipos = m_str.find(ssubout, ipos+lout);
    ++count;
  }
  return count;
}

//**********************************************************************

#endif
