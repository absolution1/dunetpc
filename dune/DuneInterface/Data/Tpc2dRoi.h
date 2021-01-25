// Tpc2dRoi.h
//
// David Adams
// January 2021
//
// Class that holds the prepared data for a TPC 2D ROI (region of interest),
// i.e. an array of floats indexed by channel and tick.
//
// The data may be read with data().value(..) and written wiht data().setValue(..).

#ifndef Tpc2dRoi_H
#define Tpc2dRoi_H

#include "dune/DuneInterface/Data/Real2dData.h"

class Tpc2dRoi {

public:

  using Index = unsigned int;
  using DataArray = Real2dData<float>;
  using IndexArray = DataArray::IndexArray;

  // Ctor for an ROI with no size.
  Tpc2dRoi() : m_sampleOffset(0), m_channelOffset(0) { }

  // Ctor for an ROI with no size.
  Tpc2dRoi(Index ncha, Index nsam, Index icha0, Index isam0 =0)
  : m_data({ncha, nsam}), m_sampleOffset(isam0), m_channelOffset(icha0) { }

  // Return the offset for the first tick wrt to a common reference.
  // Step size is one data tick assumed the same for all channels in the plane.
  Index sampleOffset() const { return m_sampleOffset; }

  // Return the number of ticks.
  Index sampleSize() const { return data().nSamples()[1]; }

  // Return the first channel.
  Index channelOffset() const { return m_channelOffset; }

  // Return the number of channels.
  Index channelSize() const { return data().nSamples()[0]; }

  // Return the data array.
  const DataArray& data() const { return m_data; }
  DataArray& data() { return m_data; }

private:

  Real2dData<float> m_data;
  Index m_sampleOffset;
  Index m_channelOffset;

};

#endif
