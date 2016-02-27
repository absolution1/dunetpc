#include "dune/RunHistory/DetPedestalDUNE.h"

#include <getopt.h>
#include <iostream>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

bool gUseDB = false;
bool gPrintPedRun = false;
int gRun = -1;
std::string gDetName = "";
std::string gCfgFile = "";
std::vector<uint64_t> gChannelList{};

//------------------------------------------------------------

void PrintUsage()
{
  std::cout << "Usage: RunInfo [options] -r|--run [run number] -d|--detector [detector name, eg, dune35t]" << std::endl;
  std::cout << "Options:" << std::endl;
  std::cout << "-f (--file) [fhicl file] : use fhicl config file" << std::endl;
  std::cout << "-D (--database)          : get peds from database" << std::endl;
  std::cout << "-P (--pedestalRun)      : print pedestal run" << std::endl;
  std::cout << "-C [chan1,chan2,...]     : print only specific channels" << std::endl;
}

//------------------------------------------------------------

bool ParseCLArgs(int argc, char* argv[])
{

  if (argc < 2) return false;
  
  struct option long_options[] = {
    {"help",   0, 0, 'h'},
    {"run",   0, 0, 'r'},
    {"detector",   0, 0, 'd'},
    {"file",   0, 0, 'f'},
    {"database",   0, 0, 'D'},
    {"pedestalRun",   0, 0, 'P'},
    {0,0,0,0}
  };

  while (1) {
    int optindx;

    int c = getopt_long(argc,argv,"r:d:f:C:DPh",long_options,&optindx);
        
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
	gDetName = optarg;
	break;
      }
    case 'D':
      {
	gUseDB = true;
	break;
      }
    case 'P':
      {
	gPrintPedRun = true;
	break;
      }
    case 'f':
      {
	gCfgFile = optarg;
	break;
      }
    case 'h':
      {
	PrintUsage();
	exit(0);
	break;
      }
    case 'C':
      {
	std::vector<std::string> cList;
	std::string channels = optarg;
	boost::tokenizer< boost::escaped_list_separator<char> > tok(channels);
	cList.assign(tok.begin(),tok.end());
	for (size_t i=0; i<cList.size(); ++i) {
	  try {
	    gChannelList.push_back(boost::lexical_cast<uint64_t>(cList[i]));
	  }
	  catch (const boost::bad_lexical_cast &) {
	    std::cout << cList[i] << " does not appear to be a channel id, skipping"
		      << std::endl;
	    continue;
	  }
	}
	break;
      }
    default:
      break;
    }
  }

  return true;
}

//------------------------------------------------------------

int main(int argc, char **argv)
{
  if (!ParseCLArgs(argc,argv)) {
    PrintUsage();
    return 1;
  }

  dune::DetPedestalDUNE* ped = new dune::DetPedestalDUNE(gDetName);

  /*
  if (gCfgFile != "")
    ped->Configure(gCfgFile);
  */

  if (gUseDB) {
    ped->SetUseDB(gUseDB);
    ped->Update(gRun);
  }

  if (gPrintPedRun) {
    std::cout << "Pedestal run: " << ped->VldTimeUsed() << std::endl;
  }
  
  if (gChannelList.empty()) 
    ped->PrintAllValues();
  else {
    for (size_t i=0; i<gChannelList.size(); ++i)
      std::cout << "Channel: " << gChannelList[i]
		<< ", Mean = " << ped->PedMean(gChannelList[i])
		<< ", RMS = " << ped->PedRms(gChannelList[i])
		<< std::endl;
  }
  
  delete ped;
  
  return 0;
}

