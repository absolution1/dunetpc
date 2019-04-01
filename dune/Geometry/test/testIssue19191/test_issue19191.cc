#include "larcore/Geometry/Geometry.h"
#include "dune/ArtSupport/ArtServiceHelper.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include <iostream>

using namespace std;

int main() {
//  ArtServiceHelper& ash = ArtServiceHelper::instance();
//  ash.addServices("test_issue19191.fcl", true);
  
  art::ServiceHandle<geo::Geometry> geom;
  double xlo, ylo, zlo;
  double xhi, yhi, zhi;
  geom->WorldBox(&xlo, &ylo, &zlo, &xhi, &yhi, &zhi);
  cout << "WorldBox: " << xlo << "," << xhi << endl;
  return 0;
}

