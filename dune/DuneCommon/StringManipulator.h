// StringManipulator.
//
// David Adams
// August 2017
//
// Class to manipulate a string.

#include <string>
#include <sstream>

class StringManipulator {

public:

  // Ctor.
  StringManipulator(std::string& strin) : m_str(strin) { }

  // Return the manipulated string.
  const std::string& string() { return m_str; }

  // Replace all occurences of a substring with a value.
  // Returns the number of replacements (<0 for error).
  template<typename T>
  int replace(std::string substr, const T& xsub);

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
