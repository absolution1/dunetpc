// DataMap.h

#ifndef DataMap_H
#define DataMap_H

// Class to store results of a calculation in maps indexed
// by name. Different types of results are supported:
//           Int: int
//         Float: double
//   FloatVector: vector<float>
//          Hist: TH1*
//         Graph: TGraph*
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
// that any stored histograms have appropriate lifetime.

#include <vector>
#include <map>
#include <string>
#include <iostream>

#include "TH1.h"
#include "TGraph.h"

class DataMap {

public:

  using Name = std::string;
  using FloatVector = std::vector<float>;
  using IntMap = std::map<Name,int>;
  using FloatMap = std::map<Name,double>;
  using FloatVectorMap = std::map<Name,FloatVector>;
  using HistVector = std::vector<TH1*>;
  using HistMap = std::map<Name,TH1*>;
  using HistVectorMap = std::map<Name,HistVector>;
  using GraphPtr = std::shared_ptr<TGraph>;
  using GraphMap = std::map<Name,GraphPtr>;
  using SharedHistVector = std::vector<std::shared_ptr<TH1>>;

  // Ctor from status flag.
  explicit DataMap(int stat =0) : m_stat(stat) { }

  // Dtor.
  ~DataMap() {
    const std::string myname = "DataMap::dtor: ";
    //std::cout << myname << "@" << this << std::endl;
  }

  // Setters.
  // If own is true, the result has owns the histogram(s).
  // If results are copied to other results, then this ownership
  // is shared so that the histogram is deleted when the last of
  // these results is deleted.
  // NOTE: Code should be fixed to remove shared pointers when a
  // histogram is overwritten (i.e. another added with the same name).
  // Graphs are always owned in this same sense.
  DataMap& setStatus(int stat) { m_stat = stat; return *this; }
  void setInt(Name name, int val) { m_ints[name] = val; }
  void setFloat(Name name, double val) { m_flts[name] = val; }
  void setFloatVector(Name name, const FloatVector& val) { m_fltvecs[name] = val; }
  void setHist(Name name, TH1* ph, bool own =false) {
    m_hsts[name] = ph;
    if ( own ) {
      ph->SetDirectory(nullptr);
      m_sharedHsts.push_back(std::shared_ptr<TH1>(ph));
    }
  }
  void setHistVector(Name name, const HistVector& hsts, bool own =false) {
    const std::string myname = "DataMap::setHistVector: ";
    m_hstvecs[name] = hsts;
    if ( own ) {
      for ( TH1* ph : hsts ) {
        bool manage = true;
        for ( std::shared_ptr<TH1>& phShared : m_sharedHsts ) {
          if ( phShared.get() == ph ) {
            std::cout << myname << "Ignoring duplicate attempt to manage a histogram." << std::endl;
            manage = false;
            break;
          }
        }
        if ( manage ) {
          ph->SetDirectory(nullptr);
          m_sharedHsts.push_back(std::shared_ptr<TH1>(ph));
        }
      }
    }
  }
  void setGraph(Name name, TGraph* pg) {
    m_grfs[name] = GraphPtr(pg);
  }

  // Extend this map with another.
  void extend(const DataMap& rhs) {
    m_stat += rhs.status();
    mapextend<int>(m_ints, rhs.m_ints);
    mapextend<double>(m_flts, rhs.m_flts);
    mapextend<FloatVector>(m_fltvecs, rhs.m_fltvecs);
    mapextend<TH1*>(m_hsts, rhs.m_hsts);
    mapextend<HistVector>(m_hstvecs, rhs.m_hstvecs);
    m_sharedHsts.insert(m_sharedHsts.end(), rhs.m_sharedHsts.begin(), rhs.m_sharedHsts.end());
    mapextend<GraphPtr>(m_grfs, rhs.m_grfs);
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
  bool haveFloatVector(Name name) const { return maphas<FloatVector>(m_fltvecs, name); }
  bool haveHist(Name name)  const { return maphas<TH1*>(m_hsts, name); }
  bool haveHistVector(Name name)  const { return maphas<HistVector>(m_hstvecs, name); }
  bool haveGraph(Name name)  const { return maphas<GraphPtr>(m_grfs, name); }

  // Return a result for a given name.
  // The indicated default is returned if a value is not stored.
  int getInt(Name name, int def =0) const { return mapget<int>(m_ints, name, def); }
  double getFloat(Name name, double def =0.0) const { return mapget<double>(m_flts, name, def); }
  const FloatVector& getFloatVector(Name name) const { return mapgetobj<FloatVector>(m_fltvecs, name); }
  TH1* getHist(Name name, TH1* def =nullptr) const { return mapget<TH1*>(m_hsts, name, def); }
  const HistVector& getHistVector(Name name) const { return mapgetobj<HistVector>(m_hstvecs, name); }
  TGraph* getGraph(Name name) const { return mapgetobj<GraphPtr>(m_grfs, name).get(); }

  // Return the maps holding the results.
  const IntMap& getIntMap() const { return m_ints; }
  const FloatMap& getFloatMap() const { return m_flts; }
  const FloatVectorMap& getFloatVectorMap() const { return m_fltvecs; }
  const HistMap& getHistMap() const { return m_hsts; }
  const HistVectorMap& getHistVectorMap() const { return m_hstvecs; }
  const GraphMap& getGraphMap() const { return m_grfs; }

  // Return the named histograms in a vector.
  HistVector getHists() const {
    HistVector hsts;
    for ( HistMap::value_type ihst : m_hsts ) hsts.push_back(ihst.second);
    return hsts;
  }
  
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
    }
    if ( m_flts.size() ) {
      cout << "Floats:" << endl;
      for ( typename FloatMap::value_type ient : m_flts ) {
        cout << "  " << ient.first << ": " << ient.second << endl;
      }
    }
    if ( m_fltvecs.size() ) {
      cout << "Float vectors:" << endl;
      unsigned int maxflt = 20;
      for ( typename FloatVectorMap::value_type ient : m_fltvecs ) {
        unsigned nflt = ient.second.size();
        cout << "  " << ient.first << "[" << nflt << "]:";
        for ( unsigned int iflt=0; iflt<nflt; ++iflt ) {
          if ( iflt ) cout << ",";
          cout << " " << ient.second[iflt];
          if ( iflt > maxflt ) {
            cout << ", ..." << endl;
            break;
          }
        }
        cout << endl;
      }
    }
    if ( m_hsts.size() ) {
      cout << "Histograms:" << endl;
      for ( typename HistMap::value_type ient : m_hsts ) {
        TH1* ph = ient.second;
        //Name hnam = ph == nullptr ? "NULL" : ph->GetName();
        cout << "  " << ient.first << ": " << ph << endl;
      }
    }
    if ( m_hstvecs.size() ) {
      cout << "Histogram vectors:" << endl;
      for ( typename HistVectorMap::value_type ient : m_hstvecs ) {
        unsigned nhst = ient.second.size();
        cout << "  " << ient.first << "[" << nhst << "]" << endl;
      }
    }
    if ( m_grfs.size() ) {
      cout << "Graphs:" << endl;
      for ( typename GraphMap::value_type ient : m_grfs ) {
        GraphPtr pg = ient.second;
        //Name hnam = ph == nullptr ? "NULL" : ph->GetName();
        cout << "  " << ient.first << ": " << pg.get() << endl;
      }
    }
  }

private:

  int m_stat;
  IntMap m_ints;
  FloatMap m_flts;
  FloatVectorMap m_fltvecs;
  HistMap m_hsts;
  HistVectorMap m_hstvecs;
  SharedHistVector m_sharedHsts;
  GraphMap m_grfs;

  // Return if a map has key.
  template<typename T>
  bool maphas(const std::map<Name,T>& vals, Name name) const {
    typename std::map<Name,T>::const_iterator ival = vals.find(name);
    return ival != vals.end();
  }

  // Return a map's value for a key or default if key is not present.
  template<typename T>
  T mapget(const std::map<Name,T>& vals, Name name, T def) const {
    typename std::map<Name,T>::const_iterator ival = vals.find(name);
    if ( ival == vals.end() ) return def;
    return ival->second;
  }

  // Return a map's value for a key or default ctor object if key is not present.

  template<typename T>
  const T& mapgetobj(const std::map<Name,T>& vals, Name name) const {
    static T def;
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
