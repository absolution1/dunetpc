// DataMap.h

#ifndef DataMap_H
#define DataMap_H

// Class to store results of a calculation in maps indexed
// by name. Different types of results are supported:
//           Int: int
//         Float: float
//        String: std::string
//     IntVector: vector<Int>
//   FloatVector: vector<Float>
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
// This class manages (i.e. deletes) its graphs and optionally manages
// its histograms. Any copies of the result share in that management.
// If a histogram is not managed, the caller should ensure that it has
// appropriate lifetime.

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <sstream>

#include "TH1.h"
#include "TGraph.h"

class DataMap {

public:

  using Name = std::string;
  using String = std::string;
  using StringMap = std::map<Name, String>;
  using IntVector = std::vector<int>;
  using IntMap = std::map<Name,int>;
  using IntVectorMap = std::map<Name,IntVector>;
  using Float = float;
  using FloatVector = std::vector<float>;
  using FloatMap = std::map<Name,Float>;
  using FloatVectorMap = std::map<Name,FloatVector>;
  using HistVector = std::vector<TH1*>;
  using HistMap = std::map<Name,TH1*>;
  using HistVectorMap = std::map<Name,HistVector>;
  using GraphPtr = std::shared_ptr<TGraph>;
  using GraphMap = std::map<Name,GraphPtr>;
  using SharedHistPtr = std::shared_ptr<TH1>;
  using SharedHistVector = std::vector<SharedHistPtr>;

public:

  // Set/fetch debug level.
  //   0 - silent
  //   1 - report when a non-existent name is requested
  //   2 - report all access
  static int dbg(int setDbg =-1) {
    static int val = 0;
    if ( setDbg >= 0 ) {
      val = setDbg;
    }
    return val;
  }

  // Class to manage map entry objects.
  template<typename T>
  struct MapEntry {
    using Value    = T;
    using Map      = typename std::map<Name,Value>;
    using Iterator = typename Map::const_iterator;
    using Entry    = typename Map::value_type;
    const Entry& m_ent;
    MapEntry(Iterator ient) : m_ent(*ient) { }
    std::string toString() const {
      std::ostringstream sout;
      sout << m_ent.first << ": " << m_ent.second;
      return sout.str();
    }
  };

  // Specialization to handle vectors.
  template<typename T>
  struct MapEntry<std::vector<T>> {
    using Value    = typename std::vector<T>;
    using Map      = typename std::map<Name,Value>;
    using Iterator = typename Map::const_iterator;
    using Entry    = typename Map::value_type;
    const Entry& m_ent;
    MapEntry(Iterator ient) : m_ent(*ient) { }
    std::string toString() const {
      std::ostringstream sout;
      sout << m_ent.first << "[" << m_ent.second.size() << "]";
      return sout.str();
    }
  };

public:

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
  void setIntVector(Name name, const IntVector& val) { m_intvecs[name] = val; }
  void setFloat(Name name, float val) { m_flts[name] = val; }
  void setFloatVector(Name name, const FloatVector& val) { m_fltvecs[name] = val; }
  void setString(Name name, String val) { m_strs[name] = val; }
  void setHist(Name name, TH1* ph, bool own =false) {
    m_hsts[name] = ph;
    if ( own && ph != nullptr ) {
      ph->SetDirectory(nullptr);
      m_sharedHsts.push_back(std::shared_ptr<TH1>(ph));
    }
  }
  void setHist(TH1* ph, bool own =false) {
    if ( ph != nullptr ) setHist(ph->GetName(), ph, own);
  }
  void setHist(std::shared_ptr<TH1> ph) {
    m_hsts[ph->GetName()] = ph.get();
    m_sharedHsts.push_back(ph);
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
  void setHistVector(Name name, const SharedHistVector& hsts) {
    const std::string myname = "DataMap::setHistVector: ";
    m_hstvecs[name].clear();
    for ( SharedHistPtr ph : hsts ) {
      m_hstvecs[name].push_back(ph.get());
      m_sharedHsts.push_back(ph);
    }
  }
  void setGraph(Name name, TGraph* pg) {
    m_grfs[name] = GraphPtr(pg);
  }
  void setGraph(TGraph* pg) {
    if ( pg != nullptr ) setGraph(pg->GetName(), pg);
  }

  // Extend this map with another.
  void extend(const DataMap& rhs) {
    m_stat += rhs.status();
    mapextend<int>(m_ints, rhs.m_ints);
    mapextend<IntVector>(m_intvecs, rhs.m_intvecs);
    mapextend<Float>(m_flts, rhs.m_flts);
    mapextend<FloatVector>(m_fltvecs, rhs.m_fltvecs);
    mapextend<String>(m_strs, rhs.m_strs);
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
  bool haveIntVector(Name name) const { return maphas<IntVector>(m_intvecs, name); }
  bool haveFloat(Name name) const { return maphas<Float>(m_flts, name); }
  bool haveFloatVector(Name name) const { return maphas<FloatVector>(m_fltvecs, name); }
  bool haveString(Name name) const { return maphas<String>(m_strs, name); }
  bool haveHist(Name name)  const { return maphas<TH1*>(m_hsts, name); }
  bool haveHistVector(Name name)  const { return maphas<HistVector>(m_hstvecs, name); }
  bool haveGraph(Name name)  const { return maphas<GraphPtr>(m_grfs, name); }

  // Return a result for a given name.
  // The indicated default is returned if a value is not stored.
  int getInt(Name name, int def =0) const { return mapget<int>(m_ints, name, def); }
  const IntVector& getIntVector(Name name) const { return mapgetobj<IntVector>(m_intvecs, name); }
  Float getFloat(Name name, Float def =0.0) const { return mapget<Float>(m_flts, name, def); }
  const FloatVector& getFloatVector(Name name) const { return mapgetobj<FloatVector>(m_fltvecs, name); }
  String getString(Name name, String def ="") const { return mapget<String>(m_strs, name, def); }
  TH1* getHist(Name name, TH1* def =nullptr) const { return mapget<TH1*>(m_hsts, name, def); }
  const HistVector& getHistVector(Name name) const { return mapgetobj<HistVector>(m_hstvecs, name); }
  TGraph* getGraph(Name name) const { return mapgetobj<GraphPtr>(m_grfs, name).get(); }

  // Return the maps holding the results.
  const IntMap& getIntMap() const { return m_ints; }
  const IntVectorMap& getIntVectorMap() const { return m_intvecs; }
  const FloatMap& getFloatMap() const { return m_flts; }
  const FloatVectorMap& getFloatVectorMap() const { return m_fltvecs; }
  const StringMap& getStringMap() const { return m_strs; }
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
  void print(std::string prefix ="") const {
    using std::cout;
    using std::endl;
    cout << prefix << "Status: " << status() << endl;
    if ( m_ints.size() ) {
      cout << prefix << "Integers:" << endl;
      for ( typename IntMap::value_type ient : m_ints ) {
        cout << prefix << "  " << ient.first << ": " << ient.second << endl;
      }
    }
    if ( m_intvecs.size() ) {
      cout << prefix << "Int vectors:" << endl;
      unsigned int maxval = 20;
      for ( typename IntVectorMap::value_type ient : m_intvecs ) {
        unsigned nval = ient.second.size();
        cout << prefix << "  " << ient.first << "[" << nval << "]:";
        for ( unsigned int ival=0; ival<nval; ++ival ) {
          if ( ival ) cout << ",";
          cout << " " << ient.second[ival];
          if ( ival > maxval ) {
            cout << ", ...";
            break;
          }
        }
        cout << endl;
      }
    }
    if ( m_flts.size() ) {
      cout << prefix << "Floats:" << endl;
      for ( typename FloatMap::value_type ient : m_flts ) {
        cout << prefix << "  " << ient.first << ": " << ient.second << endl;
      }
    }
    if ( m_fltvecs.size() ) {
      cout << prefix << "Float vectors:" << endl;
      unsigned int maxval = 20;
      for ( typename FloatVectorMap::value_type ient : m_fltvecs ) {
        unsigned nval = ient.second.size();
        cout << prefix << "  " << ient.first << "[" << nval << "]:";
        for ( unsigned int ival=0; ival<nval; ++ival ) {
          if ( ival ) cout << ",";
          cout << " " << ient.second[ival];
          if ( ival > maxval ) {
            cout << ", ...";
            break;
          }
        }
        cout << endl;
      }
    }
    if ( m_strs.size() ) {
      cout << prefix << "Strings:" << endl;
      for ( typename StringMap::value_type ient : m_strs ) {
        cout << prefix << "  " << ient.first << ": " << ient.second << endl;
      }
    }
    if ( m_hsts.size() ) {
      cout << prefix << "Histograms:" << endl;
      for ( typename HistMap::value_type ient : m_hsts ) {
        TH1* ph = ient.second;
        //Name hnam = ph == nullptr ? "NULL" : ph->GetName();
        cout << prefix << "  " << ient.first << ": " << ph << endl;
      }
    }
    if ( m_hstvecs.size() ) {
      cout << prefix << "Histogram vectors:" << endl;
      for ( typename HistVectorMap::value_type ient : m_hstvecs ) {
        unsigned nhst = ient.second.size();
        cout << prefix << "  " << ient.first << "[" << nhst << "]" << endl;
      }
    }
    if ( m_grfs.size() ) {
      cout << prefix << "Graphs:" << endl;
      for ( typename GraphMap::value_type ient : m_grfs ) {
        GraphPtr pg = ient.second;
        //Name hnam = ph == nullptr ? "NULL" : ph->GetName();
        cout << prefix << "  " << ient.first << ": " << pg.get() << endl;
      }
    }
  }

private:

  int m_stat;
  IntMap m_ints;
  IntVectorMap m_intvecs;
  FloatMap m_flts;
  FloatVectorMap m_fltvecs;
  StringMap m_strs;
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
    const std::string myname = "DataMap::mapget: ";
    typename std::map<Name,T>::const_iterator ival = vals.find(name);
    if ( ival == vals.end() ) {
      if ( dbg() ) std::cout << myname << "Unknown entry: " << name << std::endl;
      return def;
    }
    //if ( dbg() == 2 ) std::cout << myname << ient->first << std::endl;
    if ( dbg() == 2 ) std::cout << myname << MapEntry<T>(ival).toString() << std::endl;
    return ival->second;
  }

  // Return a map's value for a key or default ctor object if key is not present.
  template<typename T>
  const T& mapgetobj(const std::map<Name,T>& vals, Name name) const {
    const std::string myname = "DataMap::mapgetobj: ";
    static T def;
    typename std::map<Name,T>::const_iterator ival = vals.find(name);
    if ( ival == vals.end() ) {
      if ( dbg() ) std::cout << myname << "Unknown entry: " << name << std::endl;
      return def;
    }
    //if ( dbg() == 2 ) std::cout << myname << ient->first << std::endl;
    if ( dbg() == 2 ) std::cout << myname << MapEntry<T>(ival).toString() << std::endl;
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
