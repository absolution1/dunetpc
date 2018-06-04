// quietHistFit.cxx

#include "dune/DuneCommon/quietHistFit.h"

#include "TF1.h"
#include "TH1F.h"

int quietHistFit(TH1* ph, std::string fname, std::string fopt) {
  // Block Root info message for new Canvas produced in fit.
  int levelSave = gErrorIgnoreLevel;
  gErrorIgnoreLevel = 1001;
  // Block non-default (e.g. art) from handling the Root "error".
  // We switch to the Root default handler while making the call to Print.
  ErrorHandlerFunc_t pehSave = nullptr;
  ErrorHandlerFunc_t pehDefault = DefaultErrorHandler;
  if ( GetErrorHandler() != pehDefault ) {
    pehSave = SetErrorHandler(pehDefault);
  }
  int fstat = ph->Fit(fname.c_str(), fopt.c_str());
  if ( pehSave != nullptr ) SetErrorHandler(pehSave);
  gErrorIgnoreLevel = levelSave;
  return fstat;
}

int quietHistFit(TH1* ph, TF1* pf, std::string fopt) {
  // Block Root info message for new Canvas produced in fit.
  int levelSave = gErrorIgnoreLevel;
  gErrorIgnoreLevel = 1001;
  // Block non-default (e.g. art) from handling the Root "error".
  // We switch to the Root default handler while making the call to Print.
  ErrorHandlerFunc_t pehSave = nullptr;
  ErrorHandlerFunc_t pehDefault = DefaultErrorHandler;
  if ( GetErrorHandler() != pehDefault ) {
    pehSave = SetErrorHandler(pehDefault);
  }
  int fstat = ph->Fit(pf, fopt.c_str());
  if ( pehSave != nullptr ) SetErrorHandler(pehSave);
  gErrorIgnoreLevel = levelSave;
  return fstat;
}
