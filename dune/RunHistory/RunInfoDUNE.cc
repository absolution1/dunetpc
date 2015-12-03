#include "dune/RunHistory/RunHistoryDUNE.h"

#include "IFDatabase/Table.h"
#include <getopt.h>
#include <iostream>

int gRun = -1;
std::string gDetStr = "";
int gDetId = dune::RunHistoryDUNE::kUnknownDet;

//------------------------------------------------------------

void PrintUsage()
{
  std::cout << "Usage: RunInfo -d|--detector [detectorName, eg dune35t] -r|--run [run number]" << std::endl;

}

//------------------------------------------------------------

bool ParseCLArgs(int argc, char* argv[])
{

  struct option long_options[] = {
    {"help",   0, 0, 'h'},
    {"run",   0, 0, 'r'},
    {"detector",   0, 0, 'd'},
    {0,0,0,0}
  };

  while (1) {
    int optindx;

    int c = getopt_long(argc,argv,"r:d:h",long_options,&optindx);
        
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
	gDetStr = optarg;
	break;
      }
    case 'h':
      {
	PrintUsage();
	exit(0);
	break;
      }
    default:
      break;
    }
  }

  gDetId = dune::RunHistoryDUNE::DetNameToId(gDetStr);
  if ( gDetId != dune::RunHistoryDUNE::k35t)
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

  dune::RunHistoryDUNE* rh = new dune::RunHistoryDUNE(gDetId, gRun);
  
  std::cout << rh->RunNumber() << std::endl;

  delete rh;
  
  return 0;
}

