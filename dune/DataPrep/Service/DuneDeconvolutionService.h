// DuneDeconvolutionService.h
//
// David Adams
// July 2016
//
// Implementation of service that uses the DUNE signal shaping service (SSS) to
// carry out deconvolution.
//
// Note that SSS has a fixed FFT size and the data is temporarily padded to that
// size starting from the first sample.
//
// Code is copied from and should be equivalent to that in
// dunetpc/dune/CalData/CalWireDUNE35t_module.cc.
//
// Configuration:
//   LogLevel - message logging level: 0=none, 1=initialization, 2+=every event

#ifndef DuneDeconvolutionService_H
#define DuneDeconvolutionService_H

#include "dune/DuneInterface/AdcDeconvolutionService.h"

class DuneDeconvolutionService : public AdcDeconvolutionService {

public:

  DuneDeconvolutionService(fhicl::ParameterSet const& pset, art::ActivityRegistry&);

  int update(AdcChannelData& data) const;

  std::ostream& print(std::ostream& out =std::cout, std::string prefix ="") const;

private:

  // Configuration parameters.
  int  m_LogLevel;

};

DECLARE_ART_SERVICE_INTERFACE_IMPL(DuneDeconvolutionService, AdcDeconvolutionService, LEGACY)

#endif
