// ChannelGeo.h
//
// David Adams
// August 2018
//
// Return geometry info for a channel using the current
// geometry service.

#ifndef ChannelGeo_H
#define ChannelGeo_H

#include <utility>
#include <vector>
#include "TVector3.h"

namespace geo {
  class GeometryCore;
}

class ChannelGeo {

public:

  using Index = unsigned int;
  using Point = TVector3;
  using EndPoints = std::pair<Point, Point>;
  using EndPointsVector = std::vector<EndPoints>;

  // Ctor from channel number and geometry service.
  ChannelGeo(Index icha, const geo::GeometryCore* pgeo);

  // Ctor from channel number using the current geometry service.
  ChannelGeo(Index icha);

  // Return the channel number.
  Index channel() const { return m_icha; }

  // Is this a valid channel?
  bool isValid() const { return m_valid; }

  // Return the number of wires (wire segments).
  Index nWires() const { return m_ends.size(); }

  // Return the wire endpoints for this channel.
  const EndPointsVector& wires() const { return m_ends; }

  // Return the top or bottom endpoint for this channel.
  const Point& top() const { return m_top; }
  const Point& bottom() const { return m_bot; }

private:

  Index m_icha;
  bool m_valid;
  EndPointsVector m_ends;
  Point m_top;
  Point m_bot;

};

#endif
