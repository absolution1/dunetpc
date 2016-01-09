#ifndef _NLPlotMkr
#define _NLPlotMkr

#include <fstream>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TPaveText.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TDatime.h>
#include <TMultiGraph.h>
#include <TGaxis.h>
#include <TPaveStats.h>

// REMOVE includes below to be removed when I add partition to the OnMon header
#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TSQLStatement.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cstring>

// NOTE:  Go through ALL of the above and remove what isn't necessary!!!


using namespace std;

//
// Define global variables to be used in NearlinePlotMaker.C
//

ifstream inFile;

TDatime *TNowGMT; // current GMT time
TDatime *TNow;    // current local time
TDatime *SRtime;
Int_t GMToffset;

unsigned int run;
unsigned int subrun;

UInt_t LastRun;
UInt_t LastSR;




#endif
