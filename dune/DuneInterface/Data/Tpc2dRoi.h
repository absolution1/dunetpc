// Tpc2dRoi.h
//
// David Adams
// January 2021
//
// Class that holds the prepared data for a rectangular TPC 2D ROI (region of interest),
// i.e. an array of floats indexed by channel and tick.
//
// The data may be read with data().value(..) and written wiht data().setValue(..).

#ifndef Tpc2dRoi_H
#define Tpc2dRoi_H

#include "dune/DuneInterface/Data/Real2dData.h"
#include "dune/DuneInterface/Data/FftwReal2dDftData.h"
#include <memory>

class Tpc2dRoi {

public:

  using Index = unsigned int;
  using LongIndex = unsigned long int;
  using DataArray = Real2dData<float>;
  using IndexArray = DataArray::IndexArray;
  using Dft = FftwDouble2dDftData;
  using DftPtr = std::unique_ptr<Dft>;

  // Ctor for an ROI with no size.
  Tpc2dRoi() : m_sampleOffset(0), m_channelOffset(0) { }

  // Ctor for an ROI with size.
  Tpc2dRoi(Index ncha, Index nsam, Index icha0, LongIndex isam0 =0)
  : m_data({ncha, nsam}), m_sampleOffset(isam0), m_channelOffset(icha0) { }

  // Return the offset for the first tick wrt to a common reference.
  // Step size is one data tick assumed the same for all channels in the plane.
  LongIndex sampleOffset() const { return m_sampleOffset; }

  // Return the number of ticks.
  Index sampleSize() const { return data().nSamples()[1]; }

  // Return the first channel.
  Index channelOffset() const { return m_channelOffset; }

  // Return the number of channels.
  Index channelSize() const { return data().nSamples()[0]; }

  // Return the data array.
  const DataArray& data() const { return m_data; }
  DataArray& data() { return m_data; }

  // Return a value.
  // Returns valdef for indices out of range.
  float value(Index icha, LongIndex itck, float valdef =0.0) const {
    if ( icha < channelOffset() ) return valdef;
    if ( itck < sampleOffset() ) return valdef;
    Index chk = 0;
    DataArray::IndexArray idxs;
    idxs[0] = icha - channelOffset();
    idxs[1] = itck - sampleOffset();
    float val = data().value(idxs, &chk);
    if ( chk ) return valdef;
    return val;
  }

  // Return a pointer to DFT. Null if undefined.
  Dft* dft() { return m_pdft.get(); }
  const Dft* dft() const { return m_pdft.get(); }

  // Reset the DFT data pointer. Existing data is deleted.
  // This class now owns that DFT data.
  void resetDft(Dft* pdft) { m_pdft.reset(pdft); }

private:

  Real2dData<float> m_data;
  LongIndex m_sampleOffset;
  Index m_channelOffset;
  DftPtr m_pdft;

};

#endif
