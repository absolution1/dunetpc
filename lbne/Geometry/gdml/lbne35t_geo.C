typedef struct _drawopt 
{
  const char* volume;
  int         color;
} drawopt;


/////////////////////////////////////////////////
/////////////////////////////////////////////////
lbne35t_geo(TString volName="")
{
  gSystem->Load("libGeom");
  gSystem->Load("libGdml");

  //TGeoManager::Import("lbne35t4apa_v3.gdml");
  TGeoManager::Import("lbne35t4apa_v3_nowires.gdml");
  gGeoManager->DefaultColors();

  drawopt optuboone[] = {
// color in volumes later
    {"volCathode",	kOrange-6},
    {0, 0}
  };

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

 int showLongActive = 1;
 int showShortActive = 1;
 int showCathode = 1;
 int showGasAr = 1;
 int showCryoShell = 1;

 bool OnlyDrawAPAs = false;


 if(OnlyDrawAPAs){
   int showLongActive = 0;
   int showShortActive = 0;
   int showCathode = 0;
   int showGasAr = 0;
   int showCryoShell = 0;
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
 //gGeoManager->GetVolume("volCPATubeYSide")->SetLineColor(kOrange+2);
 gGeoManager->GetVolume("volCPATubeYSide")->SetVisibility(showCathode);
 gGeoManager->GetVolume("volCPATubeYSide")->SetTransparency(70);
 //gGeoManager->GetVolume("volCPATubeZSide")->SetLineColor(kOrange+2);
 gGeoManager->GetVolume("volCPATubeZSide")->SetVisibility(showCathode);
 gGeoManager->GetVolume("volCPATubeZSide")->SetTransparency(70);

 gGeoManager->GetVolume("volGaseousArgon")->SetLineColor(kYellow-7);
 gGeoManager->GetVolume("volGaseousArgon")->SetVisibility(showGasAr);
 gGeoManager->GetVolume("volGaseousArgon")->SetTransparency(85);


 gGeoManager->FindVolumeFast("volTPCActiveSmallestLongDrift")->SetVisibility(showLongActive);
 gGeoManager->FindVolumeFast("volTPCActiveSmallestLongDrift")->SetTransparency(80);
 gGeoManager->FindVolumeFast("volTPCActiveSmallestLongDrift")->SetLineColor(3);
 gGeoManager->FindVolumeFast("volTPCActiveMidLongDrift")->SetVisibility(showLongActive);
 gGeoManager->FindVolumeFast("volTPCActiveMidLongDrift")->SetTransparency(80);
 gGeoManager->FindVolumeFast("volTPCActiveMidLongDrift")->SetLineColor(3);
 gGeoManager->FindVolumeFast("volTPCActiveLargestLongDrift")->SetVisibility(showLongActive);
 gGeoManager->FindVolumeFast("volTPCActiveLargestLongDrift")->SetTransparency(80);
 gGeoManager->FindVolumeFast("volTPCActiveLargestLongDrift")->SetLineColor(3);
 gGeoManager->FindVolumeFast("volTPCActiveSmallestShortDrift")->SetVisibility(showShortActive);
 gGeoManager->FindVolumeFast("volTPCActiveSmallestShortDrift")->SetTransparency(80);
 gGeoManager->FindVolumeFast("volTPCActiveSmallestShortDrift")->SetLineColor(3);
 gGeoManager->FindVolumeFast("volTPCActiveMidShortDrift")->SetVisibility(showShortActive);
 gGeoManager->FindVolumeFast("volTPCActiveMidShortDrift")->SetTransparency(80);
 gGeoManager->FindVolumeFast("volTPCActiveMidShortDrift")->SetLineColor(3);
 gGeoManager->FindVolumeFast("volTPCActiveLargestShortDrift")->SetVisibility(showShortActive);
 gGeoManager->FindVolumeFast("volTPCActiveLargestShortDrift")->SetTransparency(80);
 gGeoManager->FindVolumeFast("volTPCActiveLargestShortDrift")->SetLineColor(3);

 gGeoManager->GetVolume("volFoamSouth")->SetLineColor(46);
 gGeoManager->GetVolume("volFoamSouth")->SetVisibility(showCryoShell);
 gGeoManager->GetVolume("volFoamSouth")->SetTransparency(85);
 gGeoManager->GetVolume("volFoamNorth")->SetLineColor(46);
 gGeoManager->GetVolume("volFoamNorth")->SetVisibility(showCryoShell);
 gGeoManager->GetVolume("volFoamNorth")->SetTransparency(85);
 gGeoManager->GetVolume("volFoamNorthNeck")->SetLineColor(46);
 gGeoManager->GetVolume("volFoamNorthNeck")->SetVisibility(showCryoShell);
 gGeoManager->GetVolume("volFoamNorthNeck")->SetTransparency(85);
 gGeoManager->GetVolume("volFoamEastWest")->SetLineColor(46);
 gGeoManager->GetVolume("volFoamEastWest")->SetVisibility(showCryoShell);
 gGeoManager->GetVolume("volFoamEastWest")->SetTransparency(85);
 gGeoManager->GetVolume("volFoamEastWestNeck")->SetLineColor(46);
 gGeoManager->GetVolume("volFoamEastWestNeck")->SetVisibility(showCryoShell);
 gGeoManager->GetVolume("volFoamEastWestNeck")->SetTransparency(85);
 gGeoManager->GetVolume("volFoamBottom")->SetLineColor(46);
 gGeoManager->GetVolume("volFoamBottom")->SetVisibility(showCryoShell);
 gGeoManager->GetVolume("volFoamBottom")->SetTransparency(85);
 gGeoManager->GetVolume("volFoamTop")->SetLineColor(46);
 gGeoManager->GetVolume("volFoamTop")->SetVisibility(showCryoShell);
 gGeoManager->GetVolume("volFoamTop")->SetTransparency(85);

 gGeoManager->GetVolume("volBottomConcreteShell")->SetLineColor(19);
 gGeoManager->GetVolume("volBottomConcreteShell")->SetVisibility(showCryoShell);
 gGeoManager->GetVolume("volBottomConcreteShell")->SetTransparency(90);
 gGeoManager->GetVolume("volNeckConcreteShell")->SetLineColor(19);
 gGeoManager->GetVolume("volNeckConcreteShell")->SetVisibility(showCryoShell);
 gGeoManager->GetVolume("volNeckConcreteShell")->SetTransparency(90);



 gGeoManager->GetTopNode();
 gGeoManager->CheckOverlaps(1e-5,"d");
 gGeoManager->PrintOverlaps();
 //gGeoManager->SetMaxVisNodes(70000);


  //gGeoManager->GetTopVolume()->Draw("ogl");
  //gGeoManager->FindVolumeFast("volDetEnclosure")->Draw("ogl");
  gGeoManager->FindVolumeFast("volWorld")->Draw("ogl");
  //gGeoManager->FindVolumeFast("volWorld")->Draw("");
  //gGeoManager->FindVolumeFast("volTPCLargestShortDrift")->Draw("ogl");



  TFile *tf = new TFile("lbne.root", "RECREATE");
 
  gGeoManager->Write();

  tf->Close();
}
