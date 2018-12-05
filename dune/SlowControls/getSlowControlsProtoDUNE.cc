#include "dune/SlowControls/SlowControlsProtoDUNE.h"

#include "nutools/IFDatabase/Table.h"
#include <getopt.h>
#include <iostream>

uint32_t gTime = 0;
std::string gCSVFileName = "";

//------------------------------------------------------------

void PrintUsage()
{
  std::cout << "Usage: getSlowControlsProtoDUNE -t|--time [unix time] -f|--csvfile [/path/to/csv/file]" << std::endl;
}

//------------------------------------------------------------

bool ParseCLArgs(int argc, char* argv[])
{

  struct option long_options[] = {
    {"help",   0, 0, 'h'},
    {"time",   0, 0, 't'},
    {"csvfile", 0, 0, 'f'},
    {0,0,0,0}
  };

  while (1) {
    int optindx;

    int c = getopt_long(argc,argv,"t:f:",long_options,&optindx);
        
    if (c==-1) break;
    
    switch(c) {
    case 't':
      {
	int t1 = atoi(optarg);
	if (t1 < 0) {
	  std::cout << "Invalid time." << std::endl;
	  exit(0);
	}
	gTime = t1;
	break;
      }
    case 'f':
      {
	gCSVFileName = optarg;
	break;
      }
    default:
      break;
    }
  }
  
  if (gTime<0)
    return false;

  if (gCSVFileName == "")
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

  slowctrls::SlowControlsProtoDUNE* sc = new slowctrls::SlowControlsProtoDUNE();
  sc->SetUseCondb(true);
  sc->SetVerbosity(0);
  std::vector<std::string> chanList = {"NP04_DCS_01:Heinz_I",
				       "NP04_DCS_01:Heinz_V"};

  for (unsigned int i=0; i<chanList.size(); ++i) 
    sc->AddChannel(chanList[i],i+1);

  sc->SetTimeWindow(60);
  sc->SetSlowCtrlFileName(gCSVFileName);

  for (unsigned int i=0; i<chanList.size(); ++i) {
    double value = sc->GetValue(chanList[i],gTime);
    std::cout << chanList[i] << " : " << value << std::endl;
  }

  delete sc;
  
  return 0;
}

