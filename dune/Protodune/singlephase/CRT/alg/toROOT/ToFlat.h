//File: flat.h
//Brief: Converts Triggers to a flat TTree.  Useful when Sarah wants to access CRT output 
//       in python.  
//Usage: Provide one or more EBuilder output files on
// the command line.  Files should have been produced by EBuilder.
//
//       fromEBuilder/toROOTFlat <one EBuilder file name> [another EBuilder file name]...
//
//Author: Andrew Olivier aolivier@ur.rochester.edu

#ifndef CRT_TOFLAT_H
#define CRT_TOFLAT_H

//ROOT includes
#include "TTree.h"

//crt-core includes
#include "dunetpc/dune/Protodune/singlephase/CRT/data/CRTTrigger.h"

//c++ includes
#include <iostream>

namespace CRT
{
  class ToFlat
  {
    public:
      
      template <class TFS> //TFS is an object with an interface like ART's TFileService
      ToFlat(TFS& fileService) 
      {
        fTree = fileService->template make<TTree>("EBuilder", "CRT::Triggers");

        fTree->Branch("NTriggers", &fNTriggers, "NTriggers/l");
        fTree->Branch("Modules", &fModules, "Modules[NTriggers]/s");
        fTree->Branch("ADCs", &fADCs, "ADCs[NTriggers][64]/s"); //Second index is explicitly channel number
        fTree->Branch("Timestamps", &fTimestamps, "Timestamps[NTriggers]/l");
      }

      virtual ~ToFlat(); 

      void AnalyzeEvent(const std::vector<CRT::Trigger>& triggers);

    protected:
      static constexpr size_t NChannels = 64; //There are 64 strips = channels in each CRT module
      static constexpr size_t MaxTriggers = 32*2; //I don't expect many duplicate modules in an event.  Just in case, leave space for 
                                                  //up to 2 of each module in a single event.  

      //TTree Branches
      //TODO: Writing variable-sized arrays now, so these could probably be vectors
      size_t fNTriggers;
      short fADCs[MaxTriggers][NChannels];
      unsigned long long fTimestamps[MaxTriggers];
      short fModules[MaxTriggers];
    
      TTree* fTree; //TTree to which ROOT output will be written

    private:
      void ResetBranches();
  };
}

#endif //CRT_TOFLAT_H
