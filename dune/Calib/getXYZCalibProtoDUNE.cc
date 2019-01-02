#include "dune/Calib/XYZCalibProtoDUNE.h"

#include "nutools/IFDatabase/Table.h"
#include <getopt.h>
#include <iostream>

int gRun = -1;
std::string gDataType = "mc";
std::string gXCorrFile = "";
std::string gYZCorrFile = "";

//------------------------------------------------------------

void PrintUsage()
{
  std::cout << "Usage: getXYZCalibProtoDUNE -r|--run [run number] -d|--datatype [data|mc] -x|--xcorr [xcorr csv file] -y|--yzcorr [yzcorr csv file]" << std::endl;
}

//------------------------------------------------------------

bool ParseCLArgs(int argc, char* argv[])
{

  struct option long_options[] = {
    {"help",   0, 0, 'h'},
    {"run",   0, 0, 'r'},
    {"datatype",   0, 0, 'd'},
    {"xcorr",   0, 0, 'x'},
    {"yzcorr",   0, 0, 'y'},
    {0,0,0,0}
  };

  while (1) {
    int optindx;

    int c = getopt_long(argc,argv,"r:d:x:y:",long_options,&optindx);
        
    if (c==-1) break;
    
    switch(c) {
    case 'r':
      {
	int run = atoi(optarg);
	if (run < 0) {
	  std::cout << "Invalid run number." << std::endl;
	  exit(0);
	}
	gRun = run;
	break;
      }
    case 'd':
      {
	gDataType = optarg;
	break;
      }
    case 'x':
      {
	gXCorrFile = optarg;
	break;
      }
    case 'y':
      {
	gYZCorrFile = optarg;
	break;
      }
    default:
      break;
    }
  }

  if (gRun<0)
    return false;

  if ( gDataType != "mc" && gDataType != "data")
    return false;
  
  return true;
}

//------------------------------------------------------------

int main(int argc, char **argv)
{
  if (!ParseCLArgs(argc,argv)) {
    PrintUsage();
    return 1;
  }

  calib::XYZCalibProtoDUNE* xyzCalib = new calib::XYZCalibProtoDUNE();

  xyzCalib->SetIsMC((gDataType == "mc"));
  xyzCalib->SetUseCondb(true);
  if (! gXCorrFile.empty())
    xyzCalib->SetXCorrFileName(gXCorrFile);
  if (! gYZCorrFile.empty())
    xyzCalib->SetYZCorrFileName(gYZCorrFile);

  xyzCalib->Update(gRun);

  for (int ip=0; ip<3; ++ip) {
    std::cout << "Norm correction for plane " << ip << " = " 
	      << xyzCalib->GetNormCorr(ip) << std::endl;
    for (double x=-300.; x<=300.; x+=50.)
      std::cout << "xCorr(" << x << ") = " << xyzCalib->GetXCorr(ip,x) << std::endl;

    for (double y=50.; y<=600.; y+= 100.)
      for (double z=50.; z<=600.; z+= 100.)
	std::cout << "yzCorr(0," << y << "," << z << ") = " 
		  << xyzCalib->GetYZCorr(ip,0,y,z) << std::endl;
  }
  
  delete xyzCalib;
  
  return 0;
}

