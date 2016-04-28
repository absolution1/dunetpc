
#include <vector>

typedef struct _drawopt 
{
  const char* volume;
  int         color;
} drawopt;

 int showAuxDets = 1;
 int showWirePlanes = 1;
 int showLongActive = 1;
 int showShortActive = 1;
 int showCathode = 1;
 int showGasAr = 1;
 int showCryoShell = 1;
 int showCryoPadding = 1;
 int showTrenchAndDirt = 0;
 bool OnlyDrawAPAs = false;

/////////////////////////////////////////////////
/////////////////////////////////////////////////
dune35t_geo(TString volName="")
{
  gSystem->Load("libGeom");
  gSystem->Load("libGdml");

  TGeoManager::Import("dune35t4apa_v6_nowires.gdml");
  gGeoManager->DefaultColors();

  //for (int i=0;; ++i) 
  //  {
  //    if (optuboone[i].volume==0) break;
  //    gGeoManager->FindVolumeFast(optuboone[i].volume)->SetLineColor(optuboone[i].color);
  //  }
  TList* mat = gGeoManager->GetListOfMaterials();
  TIter next(mat);
  TObject *obj;
  while (obj = next()) {
    obj->Print();
  }




 if(OnlyDrawAPAs){
   showWirePlanes = 1;
   showAuxDets = 0;
   showLongActive = 0;
   showShortActive = 0;
   showCathode = 0;
   showGasAr = 0;
   showCryoShell = 0;
   showCryoPadding = 0;
   showTrenchAndDirt = 0;
 }
 
gGeoManager->FindVolumeFast("volAPAFrameYSide-0")->SetTransparency(0);
gGeoManager->FindVolumeFast("volAPAFrameYSide-0")->SetLineColor(14);
gGeoManager->FindVolumeFast("volAPAFrameZSide-0")->SetTransparency(0);
gGeoManager->FindVolumeFast("volAPAFrameZSide-0")->SetLineColor(14);
gGeoManager->FindVolumeFast("volAPAFrameYSide-1")->SetTransparency(0);
gGeoManager->FindVolumeFast("volAPAFrameYSide-1")->SetLineColor(14);
gGeoManager->FindVolumeFast("volAPAFrameZSide-1")->SetTransparency(0);
gGeoManager->FindVolumeFast("volAPAFrameZSide-1")->SetLineColor(14);
gGeoManager->FindVolumeFast("volAPAFrameYSide-2")->SetTransparency(0);
gGeoManager->FindVolumeFast("volAPAFrameYSide-2")->SetLineColor(14);
gGeoManager->FindVolumeFast("volAPAFrameZSide-2")->SetTransparency(0);
gGeoManager->FindVolumeFast("volAPAFrameZSide-2")->SetLineColor(14);
gGeoManager->FindVolumeFast("volAPAFrameYSide-3")->SetTransparency(0);
gGeoManager->FindVolumeFast("volAPAFrameYSide-3")->SetLineColor(14);
gGeoManager->FindVolumeFast("volAPAFrameZSide-3")->SetTransparency(0);
gGeoManager->FindVolumeFast("volAPAFrameZSide-3")->SetLineColor(14);

 
 gGeoManager->GetVolume("volNeckSteelShell")->SetLineColor(19);
 gGeoManager->GetVolume("volNeckSteelShell")->SetVisibility(showCryoShell);
 gGeoManager->GetVolume("volNeckSteelShell")->SetTransparency(25);
 gGeoManager->GetVolume("volTopSteelShell")->SetLineColor(19);
 gGeoManager->GetVolume("volTopSteelShell")->SetVisibility(showCryoShell);
 gGeoManager->GetVolume("volTopSteelShell")->SetTransparency(25);
 gGeoManager->GetVolume("volBotSteelShell")->SetLineColor(19);
 gGeoManager->GetVolume("volBotSteelShell")->SetVisibility(showCryoShell);
 gGeoManager->GetVolume("volBotSteelShell")->SetTransparency(25);

 gGeoManager->GetVolume("volCathode")->SetLineColor(kOrange+2);
 gGeoManager->GetVolume("volCathode")->SetVisibility(showCathode);
 gGeoManager->GetVolume("volCathode")->SetTransparency(70);
 gGeoManager->GetVolume("volCPATubeYSide")->SetVisibility(showCathode);
 gGeoManager->GetVolume("volCPATubeYSide")->SetTransparency(70);
 gGeoManager->GetVolume("volCPATubeZSide")->SetVisibility(showCathode);
 gGeoManager->GetVolume("volCPATubeZSide")->SetTransparency(70);

 gGeoManager->GetVolume("volGaseousArgon")->SetLineColor(kYellow-7);
 gGeoManager->GetVolume("volGaseousArgon")->SetVisibility(showGasAr);
 gGeoManager->GetVolume("volGaseousArgon")->SetTransparency(85);


 gGeoManager->FindVolumeFast("volTPCActiveSmallestLongDrift")->SetVisibility(showLongActive);
 gGeoManager->FindVolumeFast("volTPCActiveSmallestLongDrift")->SetTransparency(31);
 gGeoManager->FindVolumeFast("volTPCActiveSmallestLongDrift")->SetLineColor(3);
 gGeoManager->FindVolumeFast("volTPCActiveMidLongDrift")->SetVisibility(showLongActive);
 gGeoManager->FindVolumeFast("volTPCActiveMidLongDrift")->SetTransparency(31);
 gGeoManager->FindVolumeFast("volTPCActiveMidLongDrift")->SetLineColor(3);
 gGeoManager->FindVolumeFast("volTPCActiveLargestLongDriftUpstream")->SetVisibility(showLongActive);
 gGeoManager->FindVolumeFast("volTPCActiveLargestLongDriftUpstream")->SetTransparency(31);
 gGeoManager->FindVolumeFast("volTPCActiveLargestLongDriftUpstream")->SetLineColor(3);
 gGeoManager->FindVolumeFast("volTPCActiveLargestLongDriftDownstream")->SetVisibility(showLongActive);
 gGeoManager->FindVolumeFast("volTPCActiveLargestLongDriftDownstream")->SetTransparency(31);
 gGeoManager->FindVolumeFast("volTPCActiveLargestLongDriftDownstream")->SetLineColor(3);
 gGeoManager->FindVolumeFast("volTPCActiveSmallestShortDrift")->SetVisibility(showShortActive);
 gGeoManager->FindVolumeFast("volTPCActiveSmallestShortDrift")->SetTransparency(31);
 gGeoManager->FindVolumeFast("volTPCActiveSmallestShortDrift")->SetLineColor(3);
 gGeoManager->FindVolumeFast("volTPCActiveMidShortDrift")->SetVisibility(showShortActive);
 gGeoManager->FindVolumeFast("volTPCActiveMidShortDrift")->SetTransparency(31);
 gGeoManager->FindVolumeFast("volTPCActiveMidShortDrift")->SetLineColor(3);
 gGeoManager->FindVolumeFast("volTPCActiveLargestShortDriftUpstream")->SetVisibility(showShortActive);
 gGeoManager->FindVolumeFast("volTPCActiveLargestShortDriftUpstream")->SetTransparency(31);
 gGeoManager->FindVolumeFast("volTPCActiveLargestShortDriftUpstream")->SetLineColor(3);
 gGeoManager->FindVolumeFast("volTPCActiveLargestShortDriftDownstream")->SetVisibility(showShortActive);
 gGeoManager->FindVolumeFast("volTPCActiveLargestShortDriftDownstream")->SetTransparency(31);
 gGeoManager->FindVolumeFast("volTPCActiveLargestShortDriftDownstream")->SetLineColor(3);


 gGeoManager->FindVolumeFast("volTPCPlaneUSmallestLongDrift")->SetVisibility(showWirePlanes);
 gGeoManager->FindVolumeFast("volTPCPlaneUSmallestLongDrift")->SetTransparency(1);
 gGeoManager->FindVolumeFast("volTPCPlaneUSmallestLongDrift")->SetLineColor(kMagenta);
 gGeoManager->FindVolumeFast("volTPCPlaneUMidLongDrift")->SetVisibility(showWirePlanes);
 gGeoManager->FindVolumeFast("volTPCPlaneUMidLongDrift")->SetTransparency(1);
 gGeoManager->FindVolumeFast("volTPCPlaneUMidLongDrift")->SetLineColor(kMagenta);
 gGeoManager->FindVolumeFast("volTPCPlaneULargestLongDriftUpstream")->SetVisibility(showWirePlanes);
 gGeoManager->FindVolumeFast("volTPCPlaneULargestLongDriftUpstream")->SetTransparency(1);
 gGeoManager->FindVolumeFast("volTPCPlaneULargestLongDriftUpstream")->SetLineColor(kMagenta);
 gGeoManager->FindVolumeFast("volTPCPlaneULargestLongDriftDownstream")->SetVisibility(showWirePlanes);
 gGeoManager->FindVolumeFast("volTPCPlaneULargestLongDriftDownstream")->SetTransparency(1);
 gGeoManager->FindVolumeFast("volTPCPlaneULargestLongDriftDownstream")->SetLineColor(kMagenta);
 gGeoManager->FindVolumeFast("volTPCPlaneUSmallestShortDrift")->SetVisibility(showWirePlanes);
 gGeoManager->FindVolumeFast("volTPCPlaneUSmallestShortDrift")->SetTransparency(1);
 gGeoManager->FindVolumeFast("volTPCPlaneUSmallestShortDrift")->SetLineColor(kMagenta);
 gGeoManager->FindVolumeFast("volTPCPlaneUMidShortDrift")->SetVisibility(showWirePlanes);
 gGeoManager->FindVolumeFast("volTPCPlaneUMidShortDrift")->SetTransparency(1);
 gGeoManager->FindVolumeFast("volTPCPlaneUMidShortDrift")->SetLineColor(kMagenta);
 gGeoManager->FindVolumeFast("volTPCPlaneULargestShortDriftUpstream")->SetVisibility(showWirePlanes);
 gGeoManager->FindVolumeFast("volTPCPlaneULargestShortDriftUpstream")->SetTransparency(1);
 gGeoManager->FindVolumeFast("volTPCPlaneULargestShortDriftUpstream")->SetLineColor(kMagenta);
 gGeoManager->FindVolumeFast("volTPCPlaneULargestShortDriftDownstream")->SetVisibility(showWirePlanes);
 gGeoManager->FindVolumeFast("volTPCPlaneULargestShortDriftDownstream")->SetTransparency(1);
 gGeoManager->FindVolumeFast("volTPCPlaneULargestShortDriftDownstream")->SetLineColor(kMagenta);


 gGeoManager->GetVolume("volFoamSouth")->SetLineColor(46);
 gGeoManager->GetVolume("volFoamSouth")->SetVisibility(showCryoPadding);
 gGeoManager->GetVolume("volFoamSouth")->SetTransparency(85);
 gGeoManager->GetVolume("volFoamNorth")->SetLineColor(46);
 gGeoManager->GetVolume("volFoamNorth")->SetVisibility(showCryoPadding);
 gGeoManager->GetVolume("volFoamNorth")->SetTransparency(85);
 gGeoManager->GetVolume("volFoamNorthNeck")->SetLineColor(46);
 gGeoManager->GetVolume("volFoamNorthNeck")->SetVisibility(showCryoPadding);
 gGeoManager->GetVolume("volFoamNorthNeck")->SetTransparency(85);
 gGeoManager->GetVolume("volFoamEastWest")->SetLineColor(46);
 gGeoManager->GetVolume("volFoamEastWest")->SetVisibility(showCryoPadding);
 gGeoManager->GetVolume("volFoamEastWest")->SetTransparency(85);
 gGeoManager->GetVolume("volFoamEastWestNeck")->SetLineColor(46);
 gGeoManager->GetVolume("volFoamEastWestNeck")->SetVisibility(showCryoPadding);
 gGeoManager->GetVolume("volFoamEastWestNeck")->SetTransparency(85);
 gGeoManager->GetVolume("volFoamBottom")->SetLineColor(46);
 gGeoManager->GetVolume("volFoamBottom")->SetVisibility(showCryoPadding);
 gGeoManager->GetVolume("volFoamBottom")->SetTransparency(85);
 gGeoManager->GetVolume("volFoamTop")->SetLineColor(46);
 gGeoManager->GetVolume("volFoamTop")->SetVisibility(showCryoPadding);
 gGeoManager->GetVolume("volFoamTop")->SetTransparency(85);

 gGeoManager->GetVolume("volBottomConcreteShell")->SetLineColor(19); 
 gGeoManager->GetVolume("volBottomConcreteShell")->SetVisibility(showCryoPadding);
 gGeoManager->GetVolume("volBottomConcreteShell")->SetTransparency(25);
 gGeoManager->GetVolume("volNeckConcreteShell")->SetLineColor(19);
 gGeoManager->GetVolume("volNeckConcreteShell")->SetVisibility(showCryoPadding);
 gGeoManager->GetVolume("volNeckConcreteShell")->SetTransparency(25);


 gGeoManager->GetVolume("volTrenchBottomConcreteShell")->SetLineColor(19); 
 gGeoManager->GetVolume("volTrenchBottomConcreteShell")->SetVisibility(showTrenchAndDirt);
 gGeoManager->GetVolume("volTrenchBottomConcreteShell")->SetTransparency(30);
 gGeoManager->GetVolume("volTrenchTopConcrete")->SetLineColor(19); 
 gGeoManager->GetVolume("volTrenchTopConcrete")->SetVisibility(showTrenchAndDirt);
 gGeoManager->GetVolume("volTrenchTopConcrete")->SetTransparency(70);
 gGeoManager->GetVolume("volDirtWithHole")->SetLineColor(kOrange+9); 
 gGeoManager->GetVolume("volDirtWithHole")->SetVisibility(showTrenchAndDirt);
 gGeoManager->GetVolume("volDirtWithHole")->SetTransparency(70);
 gGeoManager->GetVolume("volBerm")->SetLineColor(kOrange+9); 
 gGeoManager->GetVolume("volBerm")->SetVisibility(showTrenchAndDirt);
 gGeoManager->GetVolume("volBerm")->SetTransparency(70);

 TGeoIterator next1(gGeoManager->GetTopVolume());
 TGeoNode *node = 0;
 
 while ( (node=(TGeoNode*)next1()) ){
   const char* nm = node->GetName();
   if( (strncmp(nm, "volAuxDet", 9) == 0) ){
     node->GetVolume()->SetLineColor(kRed-3); 
     node->GetVolume()->SetVisibility(showAuxDets);
     node->GetVolume()->SetTransparency(30);
   }
 }
//    std::vector<const TGeoNode*> path(8);
//    path[0] = gGeoManager->GetTopNode();
//    RunThroughVolumeTree(path, 0);


 //gGeoManager->GetTopNode();
 gGeoManager->CheckOverlaps(1e-5,"d");
 gGeoManager->PrintOverlaps();
 gGeoManager->SetMaxVisNodes(70000);


  //gGeoManager->GetTopVolume()->Draw("ogl");
  //gGeoManager->FindVolumeFast("volDetEnclosure")->Draw("ogl");
  gGeoManager->FindVolumeFast("volWorld")->Draw("ogl");
  //gGeoManager->FindVolumeFast("volWorld")->Draw("");
  //gGeoManager->FindVolumeFast("volTPCLargestShortDriftUpstream")->Draw("ogl");



  TFile *tf = new TFile("dune.root", "RECREATE");
 
  gGeoManager->Write();

  tf->Close();
}



  //......................................................................
  void RunThroughVolumeTree( std::vector<const TGeoNode*>& path,
			     unsigned int depth)
  {
    const char* nm = path[depth]->GetName();
    
    std::cout << nm; // << std::endl;

    if( (strncmp(nm, "volAuxDet", 9) == 0) ){
      path[depth]->GetVolume()->SetLineColor(kRed-3); 
      path[depth]->GetVolume()->SetVisibility(showAuxDets);
      path[depth]->GetVolume()->SetTransparency(30);
      return;
    }
    
    //explore the next layer down
    unsigned int deeper = depth+1;
    if(deeper >= path.size()){
      std::cout << "exceeded maximum TGeoNode depth" << std::endl;
    }
    
    const TGeoVolume *v = path[depth]->GetVolume();
    int nd = v->GetNdaughters();
    std::cout << ": " << nd << " daughters " << std::endl;
    for(int i = 0; i < nd; ++i){
      path[deeper] = v->GetNode(i); // what is wrong with this line??
      RunThroughVolumeTree(path, deeper);
    }
    
  }
