typedef struct _drawopt 
{
  const char* volume;
  int         color;
} drawopt;

void dune10kt_v1_geo(TString volName="")
{
  gSystem->Load("libGeom");
  gSystem->Load("libGdml");

//  TGeoManager::Import("dune10kt_v1_workspace.gdml");
  TGeoManager::Import("protodune_v4.gdml");
//  TGeoManager::Import("dune10kt_v1_3mmPitch.gdml");
//  TGeoManager::Import("dune10kt_v1_45deg.gdml");
//  TGeoManager::Import("dune10kt_v1_workspace_nowires.gdml");
//  TGeoManager::Import("dune10kt_v1_workspace.gdml");
//  TGeoManager::Import("dune10kt_v1_3mmPitch_workspace.gdml");
//  TGeoManager::Import("dune10kt_v1_45deg_workspace.gdml");
//  TGeoManager::Import("protoDune_v1.gdml");
  
  //TList* mat = gGeoManager->GetListOfMaterials();
  //TIter next(mat);

 TGeoIterator next(gGeoManager->GetTopVolume());
 TGeoNode *node = 0;

 gGeoManager->GetVolume("volSteelShell")->SetLineColor(19);
 gGeoManager->GetVolume("volSteelShell")->SetVisibility(1);
 gGeoManager->GetVolume("volSteelShell")->SetTransparency(90);

 gGeoManager->GetVolume("volGaseousArgon")->SetLineColor(kYellow-7);
 gGeoManager->GetVolume("volGaseousArgon")->SetVisibility(1);
 gGeoManager->GetVolume("volGaseousArgon")->SetTransparency(85);


 while ( (node=(TGeoNode*)next()) ){
   const char* nm = node->GetName();

   if( (strncmp(nm, "volCathode", 10) == 0) ){
     node->GetVolume()->SetLineColor(kOrange+3); 
     node->GetVolume()->SetVisibility(1);
     node->GetVolume()->SetTransparency(30);
   }
   if( (strncmp(nm, "volTPCActive", 12) == 0) ){
     node->GetVolume()->SetLineColor(kGreen-7); 
     node->GetVolume()->SetVisibility(1);
     node->GetVolume()->SetTransparency(80);
   }
   if( (strncmp(nm, "volAPAFrame", 11) == 0) ){
     node->GetVolume()->SetLineColor(kGray); 
     node->GetVolume()->SetVisibility(1);
     node->GetVolume()->SetTransparency(20);
   }
   if( (strncmp(nm, "volG10Board", 11) == 0) ){
     node->GetVolume()->SetLineColor(kMagenta-10); 
     node->GetVolume()->SetVisibility(1);
     node->GetVolume()->SetTransparency(40);
   }
   if( (strncmp(nm, "volTPCPlane", 11) == 0) ){
     node->GetVolume()->SetLineColor(kBlue-9); 
     node->GetVolume()->SetVisibility(1);
     node->GetVolume()->SetTransparency(50);
   }
   if( (strncmp(nm, "volTPCInner", 11) == 0) ){
     node->GetVolume()->SetLineColor(kWhite); 
     node->GetVolume()->SetVisibility(1);
     node->GetVolume()->SetTransparency(100);
   }
   if( (strncmp(nm, "volTPCOuter", 11) == 0) ){
     node->GetVolume()->SetLineColor(kWhite); 
     node->GetVolume()->SetVisibility(1);
     node->GetVolume()->SetTransparency(100);
   }
   if( (strncmp(nm, "volOpDetSensitive", 17) == 0) ){
     node->GetVolume()->SetLineColor(kRed-4); 
     node->GetVolume()->SetVisibility(1);
     node->GetVolume()->SetTransparency(10);
   }
   if( (strncmp(nm, "volWorld", 8) == 0) ){
     node->GetVolume()->SetLineColor(kOrange-7); 
     node->GetVolume()->SetVisibility(1);
     node->GetVolume()->SetTransparency(100);
   }
   if( (strncmp(nm, "volFoamPadding", 14) == 0) ){
     node->GetVolume()->SetLineColor(kCyan-10); 
     node->GetVolume()->SetVisibility(1);
     node->GetVolume()->SetTransparency(95);
   }
   if( (strncmp(nm, "volSteelSupport", 15) == 0) ){
     node->GetVolume()->SetLineColor(kGray); 
     node->GetVolume()->SetVisibility(1);
     node->GetVolume()->SetTransparency(95);
   }

 }

  gGeoManager->GetTopNode();
  //gGeoManager->CheckOverlaps(1e-5);
  //gGeoManager->PrintOverlaps();
  gGeoManager->SetMaxVisNodes(70000);

  //gGeoManager->GetTopVolume()->Draw("ogl");
  gGeoManager->FindVolumeFast("volWorld")->Draw("ogl");
  //gGeoManager->FindVolumeFast("volTPCPlaneUInner")->Draw("ogl");
  //if ( ! volName.IsNull() ) gGeoManager->FindVolumeFast(volName)->Draw("ogl");
//  gGeoManager->FindVolumeFast("volCryostat")->Draw("X3D");


  TFile *tf = new TFile("protodune_v4.root", "RECREATE");
 
  gGeoManager->Write();

  tf->Close();
}

