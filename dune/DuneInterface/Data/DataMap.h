// DataMap.h

#ifndef DataMap_H
#define DataMap_H

// Class to store results of a calculation in maps indexed
// by name. Different types of results are supported:
//     Int: int
//   Float: double
//    Hist: TH1*
//
// Supports automatic conversion from int:
//   DataMap myfun() { return 3; }
// and automatic conversion to int or bool.
//   DataMap res = myfun();
//   if ( res ) cout << "Failed with error " << res.status() << endl;
//   if ( res > 10 ) cout << "Failed badly." << endl;
//   if ( ! res ) cout << "Success!" << endl;
// Explicit comparisons require conversion:
//   if ( res.status() == 2 ) cout << "Error 2: File not found." << endl;
//   if ( int(res) == 3 ) cout << "Error three" << endl;
// Example draw of histogram:
//   if ( res.hasHist("myhist") ) res.getHist("myhist")->Draw();
//
// This class does not manage (i.e. delete histograms). Its filler should ensure
// tha any stored histograms have appropriate lifetime.

#include <map>
#include <string>
#include <iostream>

class TH1;

class DataMap {

public:

  using Name = std::string;
  using IntMap = std::map<Name,int>;
  using FloatMap = std::map<Name,double>;
  using HistMap = std::map<Name,TH1*>;

  // Ctor from status flag.
  explicit DataMap(int stat =0) : m_stat(stat) { }

  // Setters.
  DataMap& setStatus(int stat) { m_stat = stat; return *this; }
  void setInt(Name name, int val) { m_ints[name] = val; }
  void setFloat(Name name, double val) { m_flts[name] = val; }
  void setHist(Name name, TH1* ph) { m_hsts[name] = ph; }

  // Extend this map with another.
  void extend(const DataMap& rhs) {
    m_stat += rhs.status();
    mapextend<int>(m_ints, rhs.m_ints);
    mapextend<double>(m_flts, rhs.m_flts);
    mapextend<TH1*>(m_hsts, rhs.m_hsts);
  }
  DataMap& operator+=(const DataMap& rhs) { extend(rhs); return *this; }

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
  
  // Display contents.
  void print() const {
    using std::cout;
    using std::endl;
    cout << "Status: " << status() << endl;
    if ( m_ints.size() ) {
      cout << "Integers:" << endl;
      for ( typename IntMap::value_type ient : m_ints ) {
        cout << "  " << ient.first << ": " << ient.second << endl;
      }
    } else {
      cout << "No integers." << endl;
    }
    if ( m_flts.size() ) {
      cout << "Floats:" << endl;
      for ( typename FloatMap::value_type ient : m_flts ) {
        cout << "  " << ient.first << ": " << ient.second << endl;
      }
    } else {
      cout << "No floats." << endl;
    }
    if ( m_hsts.size() ) {
      cout << "Histograms:" << endl;
      for ( typename HistMap::value_type ient : m_hsts ) {
        TH1* ph = ient.second;
        //Name hnam = ph == nullptr ? "NULL" : ph->GetName();
        cout << "  " << ient.first << ": " << ph << endl;
      }
    } else {
      cout << "No histograms." << endl;
    }
  }

private:

  int m_stat;
  IntMap m_ints;
  FloatMap m_flts;
  HistMap m_hsts;

  // Return if a map has key.
  template<typename T>
  bool maphas(const std::map<Name,T>& vals, Name name) const {
    typename std::map<Name,T>::const_iterator ival = vals.find(name);
    return ival != vals.end();
  }

  // Return a map's value for a key or default if they is not present.
  template<typename T>
  T mapget(const std::map<Name,T>& vals, Name name, T def) const {
    typename std::map<Name,T>::const_iterator ival = vals.find(name);
    if ( ival == vals.end() ) return def;
    return ival->second;
  }

  // Extend or overwrite the entries in a map with those from another.
  template<typename T>
  void mapextend(std::map<Name,T>& lhs, const std::map<Name,T>& rhs) {
    for ( typename std::map<Name,T>::value_type ient : rhs ) {
      Name key = ient.first;
      T val = ient.second;
      lhs[key] = val;
    }
  }

};

#endif
