#include <vector>
#include <iostream>

#include "dune/RegCVN/func/RegCVNImageUtils.h"

cvn::RegCVNImageUtils::RegCVNImageUtils(){
  // Set a default image size
  SetPixelMapSize(2880,500);
  fViewReverse = {false,false,false};
  fNViews = 3;
}

cvn::RegCVNImageUtils::RegCVNImageUtils(unsigned int nWires, unsigned int nTDCs, unsigned int nViews){
  SetPixelMapSize(2880,500);
}

cvn::RegCVNImageUtils::~RegCVNImageUtils(){

}

void cvn::RegCVNImageUtils::SetViewReversal(bool reverseX, bool reverseY, bool reverseZ){
  fViewReverse = {reverseX,reverseY,reverseZ};
}

void cvn::RegCVNImageUtils::SetViewReversal(std::vector<bool> reverseViews){
  if(reverseViews.size() != 3){
    std::cout << "Expected three views for view reversals... using defaults." << std::endl;
  }
  else{
    SetViewReversal(reverseViews[0],reverseViews[1],reverseViews[2]);
  }
  return;
}

void cvn::RegCVNImageUtils::SetPixelMapSize(unsigned int nWires, unsigned int nTDCs){
  fPixelMapWires = nWires;
  fPixelMapTDCs = nTDCs;
}


float cvn::RegCVNImageUtils::ConvertToScaledCharge(float charge){
  return charge/500.;
}

void cvn::RegCVNImageUtils::ConvertPixelMapToImageVectorF(const cvn::RegPixelMap &pm, cvn::ImageVectorF &imageVec){

  SetPixelMapSize(pm.fNWire,pm.fNTdc);

  // Strip out the charge vectors and use these
  std::vector<float> v0pe = pm.fPEX;
  std::vector<float> v1pe = pm.fPEY;
  std::vector<float> v2pe = pm.fPEZ;

  ConvertChargeVectorsToImageVectorF(v0pe, v1pe, v2pe, imageVec);
}

void cvn::RegCVNImageUtils::ConvertChargeVectorsToImageVectorF(std::vector<float> &v0pe, std::vector<float> &v1pe,
                                                           std::vector<float> &v2pe, cvn::ImageVectorF &imageVec){

  cvn::ViewVectorF view0;
  cvn::ViewVectorF view1;
  cvn::ViewVectorF view2;

  ConvertChargeVectorsToViewVectors(v0pe, v1pe, v2pe, view0, view1, view2);

  cvn::ImageVectorF newImage = BuildImageVectorF(view0,view1,view2);

  imageVec = newImage;
}

void cvn::RegCVNImageUtils::ConvertChargeVectorsToViewVectors(std::vector<float> &v0pe, std::vector<float> &v1pe, std::vector<float> &v2pe,
                                       cvn::ViewVectorF& view0, cvn::ViewVectorF& view1, cvn::ViewVectorF& view2){

  // Reverse requested views
  if(fViewReverse[0]) ReverseView(v0pe);
  if(fViewReverse[1]) ReverseView(v1pe);
  if(fViewReverse[2]) ReverseView(v2pe);


  // Write the values to the three vectors
  for (unsigned int view = 0; view < fNViews; ++view){
    cvn::ViewVectorF viewChargeVec;
    //for (unsigned int wire = 0; wire < fPixelMapWires; ++wire){
    for (unsigned int wire = 0; wire < fPixelMapWires+1; ++wire){
      std::vector<float> wireTDCVec;
      //for (unsigned int time = 0; time < fPixelMapTDCs; ++time){
      for (unsigned int time = 0; time < fPixelMapTDCs+1; ++time){

        // Get the index for the pixel map
        unsigned int element = time + fPixelMapTDCs * wire;
	// FIXME: needs to add additional dimensions because of TF model
	if (time < fPixelMapTDCs && wire < fPixelMapWires){
	        // Scale Charges
	        float val = 0;
	        if(view == 0){ val = ConvertToScaledCharge(v0pe[element]); }
	        if(view == 1){ val = ConvertToScaledCharge(v1pe[element]); }
	        if(view == 2){ val = ConvertToScaledCharge(v2pe[element]); }
	        wireTDCVec.push_back(val);
	} else wireTDCVec.push_back(0.);
      }
      //wireTDCVec.push_back(0); // this is to increase 1 dim 280 -> 281
      viewChargeVec.push_back(wireTDCVec);
    }
    if(view == 0) view0 = viewChargeVec;
    if(view == 1) view1 = viewChargeVec;
    if(view == 2) view2 = viewChargeVec;
  }

  return;

}


void cvn::RegCVNImageUtils::ReverseView(std::vector<float> &peVec){

  std::vector<float> vecCopy(peVec.size(),0.); 

  for (unsigned int w = 0; w < fPixelMapWires; ++w)
  {
    // Get our new plane number
    unsigned int newPlane = fPixelMapWires - w - 1;

    for (unsigned int t = 0; t < fPixelMapTDCs; ++t)
    {
      float val = peVec[t + fPixelMapTDCs * w];
      vecCopy[t + fPixelMapTDCs * newPlane] = val;
    }
  }

  // Copy the values back into the original vector
  for(unsigned int e = 0; e < peVec.size(); ++e){
    float val = vecCopy[e];
    peVec[e] = val;
  }

}


cvn::ImageVectorF cvn::RegCVNImageUtils::BuildImageVectorF(cvn::ViewVectorF v0, cvn::ViewVectorF v1, cvn::ViewVectorF v2){

  // Tensorflow wants things in the arrangement <wires, TDCs, views>
  cvn::ImageVectorF image;
  for(unsigned int w = 0; w < v0.size(); ++w){
    std::vector<std::vector<float> > wireVec;
    for(unsigned int t = 0; t < v0[0].size(); ++t){
      std::vector<float> timeVec;
      timeVec.push_back(v0[w][t]);
      timeVec.push_back(v1[w][t]);
      timeVec.push_back(v2[w][t]);
      wireVec.push_back(timeVec);
    } // Loop over tdcs
    image.push_back(wireVec);
  } // Loop over wires
  
  

  return image;
}


