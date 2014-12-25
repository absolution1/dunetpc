typedef struct _drawopt 
{
  const char* volume;
  int         color;
} drawopt;

lbne_geo(TString volName="")
{
  gSystem->Load("libGeom");
  gSystem->Load("libGdml");

  TGeoManager::Import("lbne35t4apa_v3_nowires.gdml");

  drawopt optuboone[] = {
// color in volumes later
    {"volCathode",	kOrange-6},
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
  gGeoManager->CheckOverlaps(1e-5);
  gGeoManager->PrintOverlaps();
  gGeoManager->SetMaxVisNodes(70000);

  gGeoManager->GetTopVolume()->Draw();
  //if ( ! volName.IsNull() ) gGeoManager->FindVolumeFast(volName)->Draw("ogl");
//  gGeoManager->FindVolumeFast("volCryostat")->Draw("X3D");


  TFile *tf = new TFile("lbne.root", "RECREATE");
 
  gGeoManager->Write();

  tf->Close();
}
