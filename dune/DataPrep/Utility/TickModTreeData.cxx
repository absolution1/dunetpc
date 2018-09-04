// TickModTreeData.cxx

#include "TickModTreeData.h"
#include "dune/DataPrep/Utility/StickyCodeMetrics.h"

#include "TTree.h"

//**********************************************************************

TickModTreeData::TickModTreeData() {
  clear();
}

//**********************************************************************

void TickModTreeData::clear() {
  run      = badIndex();
  chan     = badIndex();
  femb     = badIndex();
  fembChan = badIndex();
  itkm     = badIndex();
  nsample = 0;
  maxAdc = badIndex();
  maxAdc2 = badIndex();
  meanAdc = -1.0;
  meanAdc2 = -1.0;
  maxFraction = -1.0;
  zeroFraction = -1.0;
  oneFraction = -1.0;
  highFraction = -1.0;
  fitStatus = -1;
  fitMean = -1.0;
  fitSigma = -1.0;
  fitExcess = -1.0;
}

//**********************************************************************

void TickModTreeData::fill(const StickyCodeMetrics& scm) {
  nsample      = scm.nsample();
  maxAdc       = scm.maxAdc();
  maxAdc2      = scm.maxAdc2();
  meanAdc      = scm.meanAdc();
  meanAdc2     = scm.meanAdc2();
  maxFraction  = scm.maxFraction();
  zeroFraction = scm.zeroFraction();
  oneFraction  = scm.oneFraction();
  highFraction = scm.highFraction();
  fitStatus    = scm.fitStatus();
  fitMean      = scm.fitMean();
  fitSigma     = scm.fitSigma();
  fitExcess    = scm.fitExcess();
}

//**********************************************************************

void TickModTreeData::createBranches(TTree* ptree) {
  ptree->Branch("run", &run);
  ptree->Branch("chan", &chan);
}

//**********************************************************************
