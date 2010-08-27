typedef struct _drawopt {
  const char* volume;
  int         color;
} drawopt;

void lbne_geo(TString volName="volWorld"){

gSystem->Load("libGeom");
gSystem->Load("libGdml");

TGeoManager::Import("lbne.gdml");

drawopt optLBNE[] = {
  {"volWorld",                 kWhite},
  {"volDetEnclosure",          kMagenta},
  {"volPit",                   kAzure},
  {"volCryostat",              kOrange},
  {"volInnerCryostat",         kYellow},
  {"volTPC",                   kOrange-5},
  {"volTPCCathode",            kRed},
  {"volTPCVertWall",           kCyan-5},
  {"volTPCHorizWall",          kOrange},
  {0, 0}
};

for (int i=0;; ++i) {
  if (optLBNE[i].volume==0) break;
    gGeoManager->FindVolumeFast(optLBNE[i].volume)->SetLineColor(optLBNE[i].color);
}

TList* mat = gGeoManager->GetListOfMaterials();
TIter next(mat);
TObject *obj;
while (obj = next()) {
 obj->Print();
}

 gGeoManager->CheckOverlaps(0.01);
 gGeoManager->PrintOverlaps();
 gGeoManager->SetMaxVisNodes(70000);

 //gGeoManager->GetTopVolume()->Draw();
 gGeoManager->FindVolumeFast(volName)->Draw();
}
