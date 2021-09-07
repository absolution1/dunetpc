// StringManipulator.h

#include "dune/DuneCommon/Utility/StringManipulator.h"
#include <iostream>
#include <sstream>
#include <cctype>
#include <cstdlib>

using std::string;
using std::cout;
using std::endl;
using std::ostringstream;

using StringVector = StringManipulator::StringVector;
using StringVV = std::vector<StringVector>;

//**********************************************************************

const StringVector& StringManipulator::split(string seps, bool fullSplit) {
  m_splits.clear();
  if ( seps.size() == 0 ) {
    m_splits.push_back(m_str);
    return m_splits;
  }
  // Replace all separators with the first,
  string str = m_str;
  char csep = seps[0];
  for ( char ch : seps.substr(1) ) {
    for ( char& rch : str ) if ( rch == ch ) rch = csep;
  }
  // Split.
  string word;
  bool isSep = false;
  for ( char ch : m_str ) {
    isSep = ch == csep;
    if ( isSep ) {
      if ( fullSplit || word.size() ) {
        m_splits.push_back(word);
        word = "";
      }
    } else {
      word += ch;
    }
  }
  if ( word.size() || (fullSplit && isSep) ) m_splits.push_back(word);
  return m_splits;
}

//**********************************************************************

const StringVector& StringManipulator::patternSplit(string spat) {
  const string myname = "StringManipulator::patternSplit: ";
  string word;
  bool inPat = false;
  m_splits.clear();
  if ( spat.size() < 2 ) return m_splits;
  char chStart = spat[0];
  char chEnd = spat[spat.size() - 1 ];
  string seps = spat.substr(1, spat.size() - 2 );
  // Build the sequence of words and word vectors used to create the
  // full vector of strings.
  StringVV wordvv;
  for ( char ch : m_str ) {
    // In gap between patterns.
    if ( ! inPat ) {
      if ( ch == chStart ) {
        wordvv.push_back(StringVector(1, word));
        word = "";
        inPat = true;
      } else {
        word += ch;
      }
    // In a pattern.
    } else {
      if ( ch == chEnd ) {
        StringManipulator sman(word,false);
        wordvv.push_back(sman.split(seps));
        word = "";
        inPat = false;
      } else {
        word += ch;
      }
    }
  }
  if ( word.size() ) wordvv.push_back(StringVector(1, word));
  if ( logLevel() >= 1 ) {
    cout << myname << "Word sequence count: " << wordvv.size() << endl;
  }
  // Construct the vector of strings from the wors sequence.
  StringVector& strs1(m_splits);
  strs1.push_back("");
  for ( const StringVector& wordvec : wordvv ) {
    StringVector strs2;
    for ( string str1 : strs1 ) {
      for ( string word : wordvec ) {
        strs2.push_back(str1 + word);
      }
    }
    strs1 = strs2;
  }
  if ( logLevel() >= 1 ) {
    cout << myname << "Split count: " << m_splits.size() << endl;
  }
  return m_splits;
}
  
//**********************************************************************

string StringManipulator::
floatToString(float val, int prec, bool trunc, string sdot, string smin) {
  ostringstream ssout;
  ssout.setf(std::ios::fixed);
  ssout.precision(prec);
  ssout << val;
  string sout = ssout.str();
  if ( trunc && sout.size() ) {
    string::size_type ipos = sout.size() - 1;
    while ( ipos > 0 ) {
      if ( sout[ipos] == '.' ) {
        --ipos;
        break;
      } else if ( sout[ipos] == '0' ) --ipos;
      else break;
    }
    sout = sout.substr(0, ipos + 1);
    if ( sdot.size() && sdot != "." ) {
      for ( ipos=0; ipos<sout.size(); ++ipos ) {
        if ( sout[ipos] == '.' ) {
          sout.replace(ipos, 1, sdot);
          ipos += sdot.size() - 1;
        }
      }
    }
    if ( smin.size() && smin != "." ) {
      for ( ipos=0; ipos<sout.size(); ++ipos ) {
        if ( sout[ipos] == '-' ) {
          sout.replace(ipos, 1, smin);
          ipos += smin.size() - 1;
        }
      }
    }
  }
  return sout;
}

//**********************************************************************

bool StringManipulator::isDigits() const {
  if ( str().empty() ) return false;
  for ( char c : str() ) if ( ! std::isdigit(c) ) return false;
  return true;
}

//**********************************************************************

bool StringManipulator::isInt() const {
  if ( str().empty() ) return false;
  char c = str()[0];
  string schk = (c == '+' || c == '-') ? str().substr(1) : str();
  return StringManipulator(schk).isDigits();
}

//**********************************************************************

bool StringManipulator::isUnsignedInt() const {
  if ( str().empty() ) return false;
  char c = str()[0];
  string schk = c == '+' ? str().substr(1) : str();
  return StringManipulator(schk).isDigits();
}

//**********************************************************************

bool StringManipulator::isFloat() const {
  if ( str().size() == 0 ) return false;
  char* pch = nullptr;
  strtof(c_str(), &pch);
  return pch[0] == '\0';
}

//**********************************************************************

int StringManipulator::toInt(int badval) const {
  if ( ! isInt() ) return badval;
  return std::strtol(c_str(), nullptr, 10);
}

//**********************************************************************

unsigned int StringManipulator::toUnsignedInt(unsigned int badval) const {
  if ( ! isUnsignedInt() ) return badval;
  return std::strtoul(c_str(), nullptr, 10);
}

//**********************************************************************

float StringManipulator::toFloat(float badval) const {
  if ( str().size() == 0 ) return badval;
  char* pch = nullptr;
  float val = strtof(c_str(), &pch);
  return pch[0] == '\0' ? val : badval;
}

//**********************************************************************
