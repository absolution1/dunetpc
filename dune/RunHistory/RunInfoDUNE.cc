#include "dune/RunHistory/RunHistoryDUNE.h"

#include "nutools/IFDatabase/Table.h"
#include <getopt.h>
#include <iostream>

int gRun = -1;
std::string gDetStr = "";
int gDetId = dune::RunHistoryDUNE::kUnknownDet;
bool gDumpSCData = false;
bool gDumpASICSettings = false;
bool gPrintComponents = false;

//------------------------------------------------------------

void PrintUsage()
{
  std::cout << "Usage: RunInfo -d|--detector [detectorName, eg dune35t] -r|--run [run number] [options]" << std::endl;
  std::cout << "Options:" << std::endl;
  std::cout << "-S (--dumpSCData): dump slow controls data from database" << std::endl;
  std::cout << "-A (--dumpASICSettings): dump ASIC settings from database" << std::endl;
  std::cout << "-C (--components): print DAQ components" << std::endl;
}

//------------------------------------------------------------

bool ParseCLArgs(int argc, char* argv[])
{

  struct option long_options[] = {
    {"help",   0, 0, 'h'},
    {"run",   0, 0, 'r'},
    {"detector",   0, 0, 'd'},
    {"dumpSCData",   0, 0, 'S'},
    {"dumpASICSettings",   0, 0, 'A'},
    {"components",   0, 0, 'C'},
    {0,0,0,0}
  };

  while (1) {
    int optindx;

    int c = getopt_long(argc,argv,"r:d:hSCA",long_options,&optindx);
        
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
    case 'S':
      {
	gDumpSCData = true;
	break;
      }
    case 'A':
      {
	gDumpASICSettings = true;
	break;
      }
    case 'C':
      {
	gPrintComponents = true;
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
  
  std::cout << "Run " << rh->RunNumber() << ":" << std::endl;
  std::cout << "Cfg: " << rh->CfgLabel() << std::endl;
    
  std::cout << "RunType: " << rh->RunTypeAsString() << std::endl;
  std::cout << "Start time: " << rh->TStartAsString() << std::endl;
  std::cout << "Stop time: " << rh->TStopAsString() << std::endl;
  std::cout << "Duration: ";
  int duration = rh->Duration();
  int hours = duration/3600;
  int minutes = (duration-3600*hours)/60;
  int seconds = (duration-3600*hours-60*minutes);
  std::cout << hours << " hours, " << minutes << " minutes, " << seconds << " seconds" << std::endl;
  
  if (gPrintComponents) {
    std::cout << "Components: ";
    std::vector<std::string> comp = rh->Components();
    std::cout << comp[0];
    for (size_t i=1; i<comp.size(); ++i)
      std::cout << ", " << comp[i];
    std::cout << std::endl;
  }
  
  if (gDumpSCData) {
    rh->DumpSCData();
  }

  if (gDumpASICSettings) {
    rh->DumpASICSettings();
  }
  
  delete rh;
  
  return 0;
}

