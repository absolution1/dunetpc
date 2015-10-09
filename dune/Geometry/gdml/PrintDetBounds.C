
#include <vector>
#include <iostream>
#include "TGeoNode.h"
#include "TGeoBBox.h"
#include "TGeoVolume.h"
#include "TGeoMatrix.h"
#include "TGeoManager.h"
#include "TSystem.h"

class Bounds
{

private:
  TGeoNode* fWorld;
  TGeoNode* fDetEnc;
  TGeoHMatrix* fDetMatrix;

  ///////////// SELECT GEOMETRY HERE ///////////////
  //
     const char* fGeomPath = "dune10kt_v1_nowires.gdml";
     //const char * fGeomPath = "lbne10kt_nowires.gdml";
  //
  //////////////////////////////////////////////////

public: 

  Bounds();
  ~Bounds();
  void FindDetEnclosure( std::vector<TGeoNode*>& path,
			 unsigned int depth);

};

//......................................................................
Bounds::Bounds(){

  // Import Geometry File
  gSystem->Load("libGeom");
  gSystem->Load("libGdml");
  TGeoManager::Import(fGeomPath);


  // Find Det Enclosure
  std::vector<TGeoNode*> path(8);
  path[0] = gGeoManager->GetTopNode();
  fWorld = path[0];
  this->FindDetEnclosure(path, 0); // this sets fDetEnc and fDetMatrix
   

  // Get the Det Enc's center in world coordinates
  double worldDet[4];
  double local[4] = {0.,0.,0.,1.};
  fDetMatrix->LocalToMaster(local,worldDet);


  // Get xyz dimensions and bounds of Det Enc (indexed xyz=012)
  TGeoVolume *detvol = fDetEnc->GetVolume();
  double DetHalfDim[3] = { ((TGeoBBox*)detvol->GetShape())->GetDX(),
			   ((TGeoBBox*)detvol->GetShape())->GetDY(),
			   ((TGeoBBox*)detvol->GetShape())->GetDZ()  };

  double DetBounds[6] = { worldDet[0]-DetHalfDim[0],    // 0: -x
			  worldDet[0]+DetHalfDim[0],    // 1: +x
			  worldDet[1]-DetHalfDim[1],    // 2: -y
			  worldDet[1]+DetHalfDim[1],    // 3: +y
			  worldDet[2]-DetHalfDim[2],    // 4: -z
			  worldDet[2]+DetHalfDim[2]  }; // 5: +z

  // Get xyz dimensions of World (indexed xyz=123)
  TGeoVolume *worvol = fWorld->GetVolume();
  double WorldHalfDim[3] = { ((TGeoBBox*)worvol->GetShape())->GetDX(),
			     ((TGeoBBox*)worvol->GetShape())->GetDY(),
			     ((TGeoBBox*)worvol->GetShape())->GetDZ()  };

  // Print info
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "Units: cm" << std::endl;
  std::cout << "Coordinates:  \t\tx makes a RHS \ty is vertical \tz is beam direction" << std::endl;
  std::cout << "World Dimensions: \tx = " << 2*WorldHalfDim[0]
	    << "\ty = " <<  2*WorldHalfDim[1]
	    << "\tz = " <<  2*WorldHalfDim[2] << std::endl;
  std::cout << "Det Dimensions: \tx = " << 2*DetHalfDim[0]
	    << "\ty = " <<  2*DetHalfDim[1]
	    << "\tz = " <<  2*DetHalfDim[2] << std::endl;
  std::cout << "Det Center: \t\tx = " << worldDet[0]
	    << "\t\ty = " <<  worldDet[1]
	    << "\tz = " <<  worldDet[2] << std::endl; 
  // DetEnc center should be this (for dune10kt_v1_nowires):
  //<position name="posDetEnclosure" unit="cm" x="0" y="120.75" z="2777.99875"/>
  std::cout << std::endl;
  std::cout << "Detector Enclosure Boundaries (everything outside is rock):" << std::endl;
  std::cout << "\tx: " << DetBounds[0] << "\tto\t"<< DetBounds[1] << std::endl;
  std::cout << "\ty: " << DetBounds[2] << "\tto\t"<< DetBounds[3] << std::endl;
  std::cout << "\tz: " << DetBounds[4] << "\tto\t"<< DetBounds[5] << std::endl;
  std::cout << std::endl;
  std::cout << "Amount of rock on each side:" << std::endl;
  std::cout << "\t-x: " << WorldHalfDim[0]+DetBounds[0] << std::endl;
  std::cout << "\t+x: " << WorldHalfDim[0]-DetBounds[1] << std::endl;
  std::cout << "\t-y: " << WorldHalfDim[1]+DetBounds[2] << std::endl;
  std::cout << "\t+y: " << WorldHalfDim[1]-DetBounds[3] << std::endl;
  std::cout << "\t-z: " << WorldHalfDim[2]+DetBounds[4] << std::endl;
  std::cout << "\t+z: " << WorldHalfDim[2]-DetBounds[5] << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;

}

//......................................................................
void Bounds::FindDetEnclosure( std::vector<TGeoNode*>& path,
			       unsigned int depth)
{

  const char* nm = path[depth]->GetName();
  if( (strncmp(nm, "volDetEnclosure", 15) == 0) ){

     fDetEnc = path[depth];

     fDetMatrix = new TGeoHMatrix(*fDetEnc->GetMatrix());

     return;
  }

  //explore the next layer down
  unsigned int deeper = depth+1;
  if(deeper >= path.size()){
    std::cout << "exceeded maximum TGeoNode depth\n";
  }

  const TGeoVolume *v = path[depth]->GetVolume();
  int nd = v->GetNdaughters();
  for(int i = 0; i < nd; ++i){
    path[deeper] = v->GetNode(i);
    this->FindDetEnclosure(path, deeper);
  }

}

//......................................................................
Bounds::~Bounds(){
  //delete fWorld;
  //delete fDetEnc;
}


void PrintDetBounds(){
  new Bounds();
}
