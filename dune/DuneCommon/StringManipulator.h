// StringManipulator.h
//
// David Adams
// August 2017
//
// Class to manipulate a string.

#include <string>
#include <sstream>
#include <vector>

#ifndef StringManipulator_H
#define StringManipulator_H

class StringManipulator {

public:

  using Index = unsigned int;
  using StringVector = std::vector<std::string>;

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

  // Ctor from a string.
  // If copy is false, this object will modify the passed string.
  // If true, the modified copy may be retrieved with str().
  StringManipulator(std::string& strin, bool copy)
  : m_str(copy ? m_strCopy : strin), m_strCopy(strin) { }

  // Ctor from a const string.
  explicit StringManipulator(const std::string& sin) : m_str(m_strCopy), m_strCopy(sin) { }

  // Ctor from a C-string.
  explicit StringManipulator(const char* chin) : m_str(m_strCopy), m_strCopy(chin) { }

  // Set/get log level.
  void setLogLevel(Index logLevel) { m_LogLevel = logLevel; }
  Index logLevel() const { return m_LogLevel; }

  // Return the manipulated string.
  const std::string& str() const { return m_str; }

  // Return the manipulated C string.
  const char* c_str() const { return m_str.c_str(); }

  // Return the size of the string.
  std::string::size_type size() const { return str().size(); }

  // Replace all occurences of a substring with a value.
  // Returns the number of replacements (<0 for error).
  template<typename T>
  int replace(std::string substr, const T& xsub);

  // Replace all occurences of a substring with a value
  // with fixed width.
  // Returns the number of replacements (<0 for error).
  template<typename T>
  int replaceFixedWidth(std::string substr, const T& xsub, Index width);

  // Split a string into substrings bounded by the ends and by any of the
  // Iff fullSplit is true, empty strings are included.
  // E.g. for seps = "/,", "who,am//I?/" -->
  //   {"who", "am", "I?"}         for fullSplit = false and
  //   {"who", "am", "", "I?", ""} for fullSplit = true
  const std::vector<std::string>& split(std::string seps, bool fullSplit =false);

  // Return the strings obtained by iterating over parts in each split region.
  // e.g. "R. {A,M}. Nixon" --> "R. A. Nixon", "R. M. Nixon".
  // The >2 split pattern spat specifies the splitting:
  //   first char flags the start of a split region
  //   last char flags the end of a split region
  //   remaining chars are separators (seps in split(...)
  const StringVector& patternSplit(std::string spat);

  // Return the strings from the last split.
  const StringVector& splits() const { return m_splits; }

  // Conversions to aritmentic types.
  // These fail if the string is empty or has extra characters following.

  // Return if string only contains digits (ddd..., where d = 0, 1, ..., 9).
  bool isDigits() const;

  // Return if the string is an unsigned int (ddd... or +ddd...).
  bool isUnsignedInt() const;

  // Return if the string is a signed int (ddd..., +ddd... or -ddd...).
  bool isInt() const;

  // Return if the string is a float.
  bool isFloat() const;

  // Convert string to int. Return valbad if not int.
  int toInt(int valbad =0) const;

  // Convert string to unsigned int. Return valbad if not unsigned int.
  unsigned int toUnsignedInt(unsigned int valbad =0) const;

  // Convert string to float. Return valbad if not int.
  float toFloat(float valbad =0) const;

private:

  int m_LogLevel = 0;
  std::string& m_str;
  std::string m_strCopy;
  StringVector m_splits;

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
