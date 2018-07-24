// TickModTreeData.h

#ifndef TickModTreeData_H
#define TickModTreeData_H

#include "Rtypes.h"

class StickyCodeMetrics;

class TickModTreeData {

public:

  using UShort = unsigned short;
  using Index = unsigned int;
  using Float = float;

  static Index badIndex() { return 999999; }

  Index run      =badIndex();
  Index chan     =badIndex();
  Index femb     =badIndex();
  Index fembChan =badIndex();
  Index itkm     =badIndex();
  // Sticky code data.
  Index nsample;
  Index maxAdc;
  Index maxAdc2;
  Float meanAdc;
  Float meanAdc2;
  Float maxFraction;
  Float zeroFraction;
  Float oneFraction;
  Float highFraction;
  Float fitMean;
  Float fitSigma;
  Float fitExcess;

  // Ctor.
  TickModTreeData();

  // Clear the data.
  void clear();

  // Fill sticky code data.
  void fill(const StickyCodeMetrics& scm);

  ClassDefNV(TickModTreeData, 1);

};


#endif
