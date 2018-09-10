#include "dune/Calib/LinCalibProtoDUNE.h"

#include "nutools/IFDatabase/Table.h"
#include <getopt.h>
#include <iostream>

int gRun = -1;
std::string gDataType = "data";
std::string gCSVFile = "";

//------------------------------------------------------------

void PrintUsage()
{
  std::cout << "Usage: getLinConstsProtoDUNE -r|--run [run number] -d|--datatype [data|mc] -f|--file [csv file] " << std::endl;
}

//------------------------------------------------------------

bool ParseCLArgs(int argc, char* argv[])
{

  struct option long_options[] = {
    {"help",   0, 0, 'h'},
    {"run",   0, 0, 'r'},
    {"datatype",   0, 0, 'd'},
    {"file",   0, 0, 'f'},
    {0,0,0,0}
  };

  while (1) {
    int optindx;

    int c = getopt_long(argc,argv,"hr:d:f:",long_options,&optindx);
        
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
    case 'f':
      {
	gCSVFile = optarg;
	break;
      }
    case 'h':
    default:
      {
	return false;
      }
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

  calib::LinCalibProtoDUNE* linCalib = new calib::LinCalibProtoDUNE();

  linCalib->SetIsMC((gDataType == "mc"));
  linCalib->SetUseCondb(true);
  if (! gCSVFile.empty())
    linCalib->SetCSVFileName(gCSVFile);

  linCalib->Update(gRun);

  calib::LinConsts_t lc = linCalib->GetLinConsts(100);
  std::cout << "Linearity constants correction for channel 100:" 
	    << std::endl; 
  std::cout << "\tGain = " << lc.gain
	    << "\n\tOffset = " << lc.offset << std::endl;

  delete linCalib;
  
  return 0;
}

