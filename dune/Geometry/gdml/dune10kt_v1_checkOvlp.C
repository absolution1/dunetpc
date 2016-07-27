#include <vector>
#include <string>
#include <sstream>
#include <iostream>

#include "TSystem.h"
#include "TGeoManager.h"

void dune10kt_v1_checkOvlp(bool nowires = true)
{
  gSystem->Load("libGeom");
  gSystem->Load("libGdml");

  std::vector<std::string> geoms;
  geoms.push_back("dune10kt_v1");
  geoms.push_back("dune10kt_v1_workspace");
  geoms.push_back("dune10kt_v1_45deg");
  geoms.push_back("dune10kt_v1_45deg_workspace");
  geoms.push_back("dune10kt_v1_3mmpitch");
  geoms.push_back("dune10kt_v1_3mmpitch_workspace");

  for(unsigned int i=0; i<geoms.size(); i++){
    std::cout << "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" << std::endl;

    std::stringstream file;
    file << geoms[i];
    if( nowires ) file << "_nowires";
    file << ".gdml";

    std::cout << "Checking " << file.str() << "\n" << std::endl;

    TGeoManager::Import(file.str().c_str());
    gGeoManager->GetTopNode();
    gGeoManager->CheckOverlaps(1e-5);
    gGeoManager->PrintOverlaps();
  }

}

