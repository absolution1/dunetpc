// SignalShapingService

// David Adams
// January 2017

// Interface for the class that performs convolution and deconvolution.
//
// I have simple adopted the interface used in SignalShapingServiceDUNE.
// A method to return the expected size for the signal vector has been added.

#ifndef SignalShapingService_H
#define SignalShapingService_H

#include <vector>

namespace util {
class SignalShaping;
}

class SignalShapingService {
public:

  using Channel = unsigned int;
  using FloatVector = std::vector<float>;
  using DoubleVector = std::vector<double>;
  using DoubleVectorVector = std::vector<DoubleVector>;

  virtual ~SignalShapingService() { }

  virtual DoubleVectorVector GetNoiseFactVec() const =0;
  virtual double GetASICGain(Channel channel) const =0;
  virtual double GetShapingTime(Channel channel) const =0; 
  virtual double GetRawNoise(Channel channel) const =0;
  virtual double GetDeconNoise(Channel channel) const =0;
  virtual double GetDeconNorm() const =0;
  virtual unsigned int GetSignalSize() const =0;   // The number of ticks for convolution.

  // Accessors.
  virtual int FieldResponseTOffset(Channel channel) const =0;
  virtual const util::SignalShaping& SignalShaping(Channel channel) const =0;

  // Do convolution calcution (for simulation).
  virtual void Convolute(Channel channel,  FloatVector& sigs) const =0;
  virtual void Convolute(Channel channel, DoubleVector& sigs) const =0;

  // Do deconvolution calcution (for reconstruction).
  virtual void Deconvolute(Channel channel,  FloatVector& sigs) const =0;
  virtual void Deconvolute(Channel channel, DoubleVector& sigs) const =0;
    
};

#endif
