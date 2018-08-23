////////////////////////////////////////////////////////////////////////
// \file XYZCalibProtoDUNE.cxx
//
// \brief implementation of class for accessing (x,y,z) calibration data for ProtoDUNE
//
// \author jpaley@fnal.gov
// 
////////////////////////////////////////////////////////////////////////

// C++ language includes
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "math.h"
#include "stdio.h"

// LArSoft includes
#include "dune/Calib/XYZCalibProtoDUNE.h"

// nutools includes
#include "nutools/IFDatabase/Table.h"

// Framework includes
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//-----------------------------------------------
calib::XYZCalibProtoDUNE::XYZCalibProtoDUNE()
{
  fIsMC = true;
  fYZCorrLoaded = false;
  fXCorrLoaded = false;
  fNormCorrLoaded = false;
  fInterpolate = false;
  fCurrentTS = 0;
  fXCorrFileName="";
  fYZCorrFileName="";
  fNormCorrFileName="";
  fXCorrDBTag="";
  fYZCorrDBTag="";
  fNormCorrDBTag="";

}


//-----------------------------------------------
calib::XYZCalibProtoDUNE::XYZCalibProtoDUNE(
  fhicl::ParameterSet const& pset
)
{
  fIsMC = true;
  fYZCorrLoaded = false;
  fXCorrLoaded = false;
  fNormCorrLoaded = false;
  fInterpolate = false;
  fCurrentTS = 0;
  fXCorrFileName="";
  fYZCorrFileName="";
  fNormCorrFileName="";
  fXCorrDBTag="";
  fYZCorrDBTag="";
  fNormCorrDBTag="";
  Configure(pset);
}

//------------------------------------------------
bool calib::XYZCalibProtoDUNE::Configure(fhicl::ParameterSet const& pset)
{  
  fUseCondbXYZCorr = pset.get<bool>("UseCondbXYZCorr");
  fInterpolate     = pset.get<bool>("Interpolate");
  fXCorrFileName   = pset.get<std::string>("XCorrFileName");
  fYZCorrFileName  = pset.get<std::string>("YZCorrFileName");
  fNormCorrFileName= pset.get<std::string>("NormCorrFileName");
  fXCorrDBTag      = pset.get<std::string>("XCorrDBTag");
  fYZCorrDBTag     = pset.get<std::string>("YZCorrDBTag");
  fNormCorrDBTag   = pset.get<std::string>("NormCorrDBTag");
  return true;
}

//------------------------------------------------
bool calib::XYZCalibProtoDUNE::Update(uint64_t ts) 
{

  if (fYZCorrLoaded && ts != fCurrentTS) {
    fYZCorrHist.clear();
    fYZCorrLoaded = false;
  }

  if (fXCorrLoaded && ts != fCurrentTS) {
    fXCorrHist.clear();
    fXCorrLoaded = false;
  }

  if (fNormCorrLoaded && ts != fCurrentTS) {
    fNormCorr.clear();
    fNormCorrLoaded = false;
  }

  fCurrentTS = ts;
  // all done! 

  return true;
}

//------------------------------------------------
double calib::XYZCalibProtoDUNE::GetNormCorr(int chanId) 
{
  if (!fNormCorrLoaded) this->LoadNormCorr();

  if (fNormCorr.find(chanId) == fNormCorr.end()) {
    mf::LogError("XYZCalibProtoDUNE") << "Plane not found!";
    return 0.;
  }

  return fNormCorr[chanId].corr;
}

//------------------------------------------------
double calib::XYZCalibProtoDUNE::GetXCorr(int plane, double x) 
{
  if (!fXCorrLoaded) this->LoadXCorr();

  if (fXCorrHist.find(plane) == fXCorrHist.end()) {
    mf::LogError("XYZCalibProtoDUNE") << "Plane not found!";
    return 0.;
  }

  if (fInterpolate)
    return fXCorrHist[plane].Interpolate(x);
  else {
    int ix = fXCorrHist[plane].FindBin(x);
    return fXCorrHist[plane].GetBinContent(ix);
  }

}

//------------------------------------------------
double calib::XYZCalibProtoDUNE::GetYZCorr(int plane, int side, 
					   double y, double z) 
{
  if (!fYZCorrLoaded) this->LoadYZCorr();

  int chanId = plane*10+side;

  if (fYZCorrHist.find(chanId) == fYZCorrHist.end()) {
    mf::LogError("XYZCalibProtoDUNE") << "Plane not found!";
    return 0.;
  }

  if (fInterpolate)
    return fYZCorrHist[chanId].Interpolate(z,y);
  else {
    int iz = fYZCorrHist[chanId].GetXaxis()->FindBin(z);
    int iy = fYZCorrHist[chanId].GetYaxis()->FindBin(y);
    return fYZCorrHist[chanId].GetBinContent(iz,iy);
  }
  
}

//------------------------------------------------
bool calib::XYZCalibProtoDUNE::LoadNormCorr()
{
  if (!fUseCondbXYZCorr) return true;

  if (fNormCorrLoaded) return true;

  nutools::dbi::Table NormCorrTable;
  
  NormCorrTable.SetDetector("pdunesp");
  NormCorrTable.SetTableName("distcorrnorm");
  NormCorrTable.SetTableType(nutools::dbi::kConditionsTable);
  NormCorrTable.SetDataTypeMask(nutools::dbi::kDataOnly);
  if (fIsMC)
    NormCorrTable.SetDataTypeMask(nutools::dbi::kMCOnly);
  
  int normIdx = NormCorrTable.AddCol("norm","double");
  int normErrIdx = NormCorrTable.AddCol("norm_err","double");
  
  NormCorrTable.SetMinTSVld(fCurrentTS);
  NormCorrTable.SetMaxTSVld(fCurrentTS);
  NormCorrTable.SetTag(fNormCorrDBTag);

  NormCorrTable.SetVerbosity(100);

  bool readOk = false;
  if (!fNormCorrFileName.empty()) 
    readOk = NormCorrTable.LoadFromCSV(fNormCorrFileName);
  else
    readOk = NormCorrTable.Load();

  if (! readOk) {
    mf::LogError("XYZCalibProtoDUNE") << "Load from norm calib database table failed.";
    
    return false; //std::abort();

  }
  
  if (NormCorrTable.NRow() == 0) {
    mf::LogError("XYZCalibProtoDUNE") << "Number of rows in norm calib table is 0.  This should never be the case!";
    return false;
  }
  
  nutools::dbi::Row* row;
  uint64_t chan;
  for (int i=0; i<NormCorrTable.NRow(); ++i) {
    NormCorr_t norm;
    row = NormCorrTable.GetRow(i);      
    chan = row->Channel();
    row->Col(normIdx).Get(norm.corr);
    row->Col(normErrIdx).Get(norm.corr_err);
    fNormCorr[chan] = norm;
  }    

  fNormCorrLoaded = true;
  return true;
}

//------------------------------------------------
bool calib::XYZCalibProtoDUNE::LoadXCorr()
{
  if (!fUseCondbXYZCorr) return true;

  if (fXCorrLoaded) return true;

  nutools::dbi::Table XCorrTable;
  
  XCorrTable.SetDetector("pdunesp");
  XCorrTable.SetTableName("distcorrx");
  XCorrTable.SetTableType(nutools::dbi::kConditionsTable);
  XCorrTable.SetDataTypeMask(nutools::dbi::kDataOnly);
  if (fIsMC)
    XCorrTable.SetDataTypeMask(nutools::dbi::kMCOnly);
  
  int shapeIdx = XCorrTable.AddCol("shape","double");
  int shapeErrIdx = XCorrTable.AddCol("shape_err","double");
  int xIdx = XCorrTable.AddCol("x","double");
  int dxIdx = XCorrTable.AddCol("dx","double");
  
  XCorrTable.SetMinTSVld(fCurrentTS);
  XCorrTable.SetMaxTSVld(fCurrentTS);
  XCorrTable.SetTag(fXCorrDBTag);
  
  XCorrTable.SetVerbosity(100);

  bool readOk = false;
  if (!fXCorrFileName.empty()) 
    readOk = XCorrTable.LoadFromCSV(fXCorrFileName);
  else
    readOk = XCorrTable.Load();

  if (! readOk) {
    mf::LogError("XYZCalibProtoDUNE") << "Load from x calib database table failed.";
    return false; //std::abort();
  }
  
  if (XCorrTable.NRow() == 0) {
    mf::LogError("XYZCalibProtoDUNE") << "Number of rows in x calib table is 0.  This should never be the case!";
    return false;
  }
  
  nutools::dbi::Row* row;
  uint64_t chan;
  int plane;
  std::vector<int> planeVec;

  std::map<int,std::vector<XCorr_t> > fXCorr;

  for (int i=0; i<XCorrTable.NRow(); ++i) {
    row = XCorrTable.GetRow(i);      
    chan = row->Channel();
    plane = int(chan/10000);

    XCorr_t xcorr;
    row->Col(xIdx).Get(xcorr.x);
    row->Col(dxIdx).Get(xcorr.dx);
    row->Col(shapeIdx).Get(xcorr.corr);
    row->Col(shapeErrIdx).Get(xcorr.corr_err);
    
    if (fXCorr.find(plane) == fXCorr.end()) {
      planeVec.push_back(plane);
      std::vector<XCorr_t> xcorrVec;
      fXCorr[plane] = xcorrVec;      
    }

    fXCorr[plane].push_back(xcorr);

  }    

  // sort the x-corrections by x for easy look-up
  for (unsigned int i=0; i<planeVec.size(); ++i) {
    int ip = planeVec[i];
    std::sort(fXCorr[ip].begin(),fXCorr[ip].end());
    char hname[256];
    sprintf(hname,"xCorrHist_%d",ip);
    int nbins = int(fXCorr[ip].size());
    double xmin = fXCorr[ip][0].x;
    double xmax = fXCorr[ip][nbins-1].x;
    double bw = (xmax-xmin)/(nbins-1);
    fXCorrHist[ip] = TH1F(hname,"",nbins,xmin-bw,xmax+bw);
    for (unsigned int j=0; j<fXCorr[ip].size(); ++j) 
      fXCorrHist[ip].SetBinContent(j+1,fXCorr[ip][j].corr);
  }
    
  fXCorrLoaded = true;
  return true;

}

//------------------------------------------------
bool calib::XYZCalibProtoDUNE::LoadYZCorr()
{
  if (!fUseCondbXYZCorr) return true;

  if (fYZCorrLoaded) return true;

  nutools::dbi::Table YZCorrTable;
  
  YZCorrTable.SetDetector("pdunesp");
  YZCorrTable.SetTableName("distcorryz");
  YZCorrTable.SetTableType(nutools::dbi::kConditionsTable);
  YZCorrTable.SetDataTypeMask(nutools::dbi::kDataOnly);
  if (fIsMC)
    YZCorrTable.SetDataTypeMask(nutools::dbi::kMCOnly);
  
  int corrIdx = YZCorrTable.AddCol("corr","double");
  int corrErrIdx = YZCorrTable.AddCol("corr_err","double");
  int yIdx = YZCorrTable.AddCol("y","double");
  //  int dyIdx = YZCorrTable.AddCol("dy","double");
  int zIdx = YZCorrTable.AddCol("z","double");
  //  int dzIdx = YZCorrTable.AddCol("dz","double");
  
  YZCorrTable.SetMinTSVld(fCurrentTS);
  YZCorrTable.SetMaxTSVld(fCurrentTS);
  YZCorrTable.SetTag(fYZCorrDBTag);

  YZCorrTable.SetVerbosity(100);

  bool readOk = false;
  if (!fYZCorrFileName.empty()) 
    readOk = YZCorrTable.LoadFromCSV(fYZCorrFileName);
  else
    readOk = YZCorrTable.Load();

  if (! readOk) {
    mf::LogError("XYZCalibProtoDUNE") << "Load from yz calib database table failed.";
    return false; //std::abort();
  }
  
  if (YZCorrTable.NRow() == 0) {
    mf::LogError("XYZCalibProtoDUNE") << "Number of rows in yz calib table is 0.  This should never be the case!";
    return false;
  }
  
  nutools::dbi::Row* row;
  uint64_t chan;
  std::vector<int> planeVec;

  std::map<int,std::vector<YZCorr_t> > fYZCorr;

  for (int i=0; i<YZCorrTable.NRow(); ++i) {
    row = YZCorrTable.GetRow(i);      
    chan = row->Channel();
    int plane = chan/10000000;
    int side = (chan-plane*10000000)/1000000;
    int chanId = plane*10+side;

    YZCorr_t yzcorr;
    row->Col(yIdx).Get(yzcorr.y);
    //    row->Col(dyIdx).Get(yzcorr.dy);
    row->Col(zIdx).Get(yzcorr.z);
    //    row->Col(dzIdx).Get(yzcorr.dz);
    row->Col(corrIdx).Get(yzcorr.corr);
    row->Col(corrErrIdx).Get(yzcorr.corr_err);

    if (fYZCorr.find(chanId) == fYZCorr.end()) 
      planeVec.push_back(chanId);

    fYZCorr[chanId].push_back(yzcorr);
    
  }    

  // sort the (y,z)-corrections by y and then by z for easy look-up
  int nbinsy=0;
  int nbinsz=0;
  double ymin=0.;
  double ymax=0.;
  double zmin=0;
  double zmax=0.;
  double bwy=0.;
  double bwz=0.;

  for (unsigned int i=0; i<planeVec.size(); ++i) {
    int ip = planeVec[i];
    std::sort(fYZCorr[ip].begin(),fYZCorr[ip].end());
    if (i==0) {
      ymin = fYZCorr[ip][0].y;
      zmin = fYZCorr[ip][0].z;
      ymax = fYZCorr[ip][fYZCorr[ip].size()-1].y;
      zmax = fYZCorr[ip][fYZCorr[ip].size()-1].z;
      // now figure out how many z bins there are
      for (unsigned j=1; j<fYZCorr[ip].size(); ++j, ++nbinsz) 
	if (fYZCorr[ip][j].z < fYZCorr[ip][j-1].z) break;
      nbinsy = int(fYZCorr[ip].size())/++nbinsz;
      bwz = (zmax-zmin)/(nbinsz-1);
      bwy = (ymax-ymin)/(nbinsy-1);
    }
    char hname[256];
    sprintf(hname,"yzCorrHist_%d",i);
    fYZCorrHist[ip] = TH2F(hname,"",nbinsz,zmin-bwz,zmax+bwz,
			   nbinsy,ymin-bwy,ymax+bwy);
    for (unsigned int j=0; j<fYZCorr[ip].size(); ++j) 
      fYZCorrHist[ip].Fill(fYZCorr[ip][j].z,fYZCorr[ip][j].y,fYZCorr[ip][j].corr);
  }
  
  fYZCorrLoaded = true;
  return true;

}


