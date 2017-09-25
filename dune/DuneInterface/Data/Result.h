// Result.h

#ifndef Result_H
#define Result_H

// Class to store results of a calculation in maps indexed
// by name. Different types of results are supported:
//     Int: int
//   Float: double
//    Hist: TH1*
//
// Supports automatic conversion from int:
//   Result myfun() { return 3; }
// and automatic conversion to int or bool.
//   Result res = myfun();
//   if ( res ) cout << "Failed with error " << res.status() << endl;
//   if ( res > 10 ) cout << "Failed badly." << endl;
//   if ( ! res ) cout << "Success!" << endl;
// Explicit comparisons require conversion:
//   if ( res.status() == 2 ) cout << "Error 2: File not found." << endl;
//   if ( int(res) == 3 ) cout << "Error three" << endl;
// Example draw of histogram:
//   if ( res.hasHist("myhist") ) res.getHist("myhist")->Draw();

#include <map>
#include <string>

class TH1;

class Result {

public:

  using Name = std::string;
  using IntMap = std::map<Name,int>;
  using FloatMap = std::map<Name,double>;
  using HistMap = std::map<Name,TH1*>;

  // Ctor.
  Result(int stat =0) : m_stat(stat) { }

  // Status.
  // Typically 0 indicates success.
  int status() const { return m_stat; }
  operator int() const { return m_stat; }
  explicit operator bool() const { return m_stat; }

  // Return if a result is stored for a given name.
  bool haveInt(Name name)   const { return maphas<int>(m_ints, name); }
  bool haveFloat(Name name) const { return maphas<double>(m_flts, name); }
  bool haveHist(Name name)  const { return maphas<TH1*>(m_hsts, name); }

  // Return a result for a given name.
  // The indicated default is returned if a value is not stored.
  int getInt(Name name, int def =0) const { return mapget<int>(m_ints, name, def); }
  double getFloat(Name name, double def =0.0) const { return mapget<double>(m_flts, name, def); }
  TH1* getHist(Name name, TH1* def =nullptr) const { return mapget<TH1*>(m_hsts, name, def); }

  // Return the maps holding the results.
  const IntMap& getIntMap() const { return m_ints; }
  const FloatMap& getFloatMap() const { return m_flts; }
  const HistMap& getHistMap() const { return m_hsts; }
  
private:

  int m_stat;
  IntMap m_ints;
  FloatMap m_flts;
  HistMap m_hsts;

  template<typename T>
  bool maphas(const std::map<Name,T>& vals, Name name) const {
    typename std::map<Name,T>::const_iterator ival = vals.find(name);
    return ival != vals.end();
  }

  template<typename T>
  T mapget(const std::map<Name,T>& vals, Name name, T def) const {
    typename std::map<Name,T>::const_iterator ival = vals.find(name);
    if ( ival == vals.end() ) return def;
    return ival->second;
  }

};

#endif
