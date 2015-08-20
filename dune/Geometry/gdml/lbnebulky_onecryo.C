typedef struct _drawopt 
{
  const char* volume;
  int         color;
} drawopt;

dunebulky_onecryo(TString volName="")
{
  gSystem->Load("libGeom");
  gSystem->Load("libGdml");

  TGeoManager::Import("dunebulky_onecryo.gdml");

  drawopt optuboone[] = {
    {"volRockTop",        kOrange-7},
    {"volRockBottom",        kOrange-7},
    {"volUpperRockWithCavern",        kOrange-7},
    {"volLowerRockWithCavern",        kOrange-7},
    {0, 0}
  };

  for (int i=0;; ++i) {
    if (optuboone[i].volume==0) break;
      gGeoManager->FindVolumeFast(optuboone[i].volume)->SetLineColor(optuboone[i].color);
  }
  TList* mat = gGeoManager->GetListOfMaterials();
  TIter next(mat);
  TObject *obj;
  while (obj = next()) {
    obj->Print();
  }

  gGeoManager->GetTopNode();
  gGeoManager->CheckOverlaps(1e-16);
  gGeoManager->PrintOverlaps();
  gGeoManager->SetMaxVisNodes(70000);

  //gGeoManager->GetTopVolume()->Draw();
  //if ( ! volName.IsNull() ) gGeoManager->FindVolumeFast(volName)->Draw("ogl");
  gGeoManager->FindVolumeFast("volWorld")->Draw("ogl");


  TFile *tf = new TFile("dunebulky_onecryo.root", "RECREATE");
 
  gGeoManager->Write();

  tf->Close();
}
