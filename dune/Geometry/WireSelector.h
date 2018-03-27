// WireSelector.h
//
// David Adams
// March 2018
//
// Class to select wires and return measured coordinate and
// readout channel for each selected wire.
//
// The planes are selected automatically when the this object is constructed
// and after calls to any of the plane selection methods. Call fillData()
// to construct the vector of wires.

#ifndef WireSelector_H
#define WireSelector_H

#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include <vector>
#include <map>

namespace geo {
class GeometryCore;
}

class WireSelector {

public:

  using Index = unsigned int;
  using IndexVector = std::vector<Index>;
  using GeometryCore = geo::GeometryCore;
  using View = geo::View_t;
  using PlaneID = geo::PlaneID;
  using PlaneIDVector = std::vector<PlaneID>;
  
  // Class that describes a wire and its associated drift volume.
  // x is the position of the wire in the drift direction, i.e. perpendicular to the plane
  // y is the global coordinate along the wire
  // z is the global coordinate the wire number measures, i.e. perpendicular to the wire in the wire plane.
  struct WireInfo {
    float x;
    float y;
    float z;
    float driftMax;  // < 0 if TPC volume is x < w_wire
    float length;
    float pitch;
    Index channel;    // Channel that reads out this wire.
    WireInfo() : y(0.0), z(0.0), length(0.0), pitch(0.0), channel(0) {}
    WireInfo(float a_x, float a_y, float a_z,
             float a_driftMax, float a_length, float a_pitch,
             Index a_channel)
      : x(a_x), y(a_y), z(a_z), driftMax(a_driftMax), length(a_length), pitch(a_pitch), channel(a_channel) { }
    float x1() const { return driftMax > 0.0 ? x : x + driftMax; }
    float x2() const { return driftMax < 0.0 ? x : x + driftMax; }
    float y1() const { return y - 0.5*length; }
    float y2() const { return y + 0.5*length; }
    float z1() const { return z - 0.5*pitch; }
    float z2() const { return z + 0.5*pitch; }
  };

  using WireInfoVector = std::vector<WireInfo>;
  using WireInfoMap = std::multimap<Index, const WireInfo*>;

  // Ctor from a cryostat.
  explicit WireSelector(Index icry =0);

  // Const methods.

  // Return the geometry.
  const GeometryCore* geometry() const { return m_pgeo; }

  // Return the selected cryostats.
  const IndexVector& cryostats() const { return m_icrys; }

  // View if a view selection was made.
  View view() const { return m_view; }

  // Wire angle if a wire angle selection was made.
  double wireAngle() const { return m_wireAngle; }
  double wireAngleTolerance() const { return m_wireAngleTolerance; }

  // TPC drift range.
  double driftMin() const { return m_driftMin; }
  double driftMax() const { return m_driftMax; }

  // Return the selected wire planes. Use geo::kUnknown to select all.
  const PlaneIDVector& planeIDs() const { return m_pids; }

  // Return the wire data.
  bool haveData() const { return m_haveData; };
  const WireInfoVector& data() const { return m_data; }
  const WireInfoMap& dataMap() const { return m_datamap; }

  // Non-const methods.

  // Fetch the geometry (if needed).
  const GeometryCore* geometry();

  // Specify the list of cryostats.
  // If this is not called, only crystat 0 is used.
  void selectCryostats(const IndexVector& a_icrys);

  // Select planes with the specified geometry view (geo::kU, geo::kV, geo::kZ, ...).
  // Use geo::kUnknown to select all planes.
  // Vector of selected planes is updated.
  void selectView(View view);

  // Select planes within tol of a specified wire angle.
  // Selects all if tol > pi.
  // Vector of selected planes is updated.
  void selectWireAngle(double wireAngle, double tol =0.001);

  // Select planes in TPCs with max drift distance in the range [dmin, dmax).
  // Vector of selected planes is updated.
  void selectDrift(double dmin, double dmax =1.e20);

  // Select planes using the current criteria.
  // This is called automatically by ctors and plane selectors.
  void selectPlanes();

  // Returns the wire data after building it if it is not already present.
  const WireInfoVector& fillData();

  // Returns the channel-mapped wire data after building it if it is not already present.
  const WireInfoMap& fillDataMap();

  // Clear the wire data.
  void clearData();

private:

  const GeometryCore* m_pgeo =nullptr;
  IndexVector m_icrys;
  View m_view =geo::kUnknown;
  double m_wireAngle =999.;
  double m_wireAngleTolerance =999.;
  double m_driftMin = 0.0;
  double m_driftMax = 1.e20;
  PlaneIDVector m_pids;
  bool m_haveData =false;
  WireInfoVector m_data;
  WireInfoMap m_datamap;

};

#endif
