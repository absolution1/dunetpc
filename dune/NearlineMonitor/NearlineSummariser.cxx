#include <iostream>
#include <string>
#include <vector>
#include <memory>

//ROOT

#include "TFile.h"
#include "TTree.h"
#include "TList.h"
#include "TKey.h"
#include "TH1.h"

const int NearlineMinorVersion=1;
const int NearlineMajorVersion=0;

const std::string NearlineInstanceLabel = "nearlineana";

bool checkFileHasDir(TFile* fp, std::string dirName);
TH1* getHistogram(TFile* fp, std::string histName, std::string newHistName);

struct NearlineFileInfo{

  bool HasHeader;
  unsigned int Run;
  unsigned int Subrun;
  int          FirstEvent;
  int          LastEvent;
  int          Nevents;
  unsigned int StartYear;
  unsigned int EndYear;
  unsigned int StartMonth;
  unsigned int EndMonth;
  unsigned int StartDay;
  unsigned int EndDay;
  double       StartHour;
  double       EndHour;
  unsigned long long int StartTime;
  unsigned long long int EndTime;
  
  //NearlineVersionNumbers
  int ThisNearlineMinorVersion;
  int ThisNearlineMajorVersion;

  NearlineFileInfo(TFile *fp, std::string headerTreeName = NearlineInstanceLabel + "/Header", std::string versionHistName = NearlineInstanceLabel + "/hist_nearline_version"):
    HasHeader(false),
    Run(0),
    Subrun(0),
    FirstEvent(1e9),
    LastEvent(-1),
    Nevents(0),
    StartYear(0),
    EndYear(0),
    StartMonth(0),
    EndMonth(0),
    StartDay(0),
    EndDay(0),
    StartHour(0.0),
    EndHour(0.0),
    StartTime(-1), // this is an unsigned int so it will default to a huge number
    EndTime(0),
    ThisNearlineMinorVersion(0),
    ThisNearlineMajorVersion(0)
  {

    getHeaderInfo(fp, headerTreeName);
    getNearlineVersion(fp, versionHistName);

  }//constructor

  void getHeaderInfo(TFile *fp, std::string headerTreeName){
    TTree *Header = (TTree*) fp->Get(headerTreeName.c_str());
    if(!Header){
      std::cout << "NearlineSummariser: ERROR: " << fp->GetName() << " doesn't contain header tree: " << headerTreeName << std::endl;
      return;
    }
    
    Header->SetBranchAddress("Run",&Run);
    Header->SetBranchAddress("Subrun",&Subrun);
    Header->SetBranchAddress("FirstEvent",&FirstEvent);
    Header->SetBranchAddress("LastEvent",&LastEvent);
    Header->SetBranchAddress("Nevents",&Nevents);
    Header->SetBranchAddress("StartYear",&StartYear);
    Header->SetBranchAddress("StartMonth",&StartMonth);
    Header->SetBranchAddress("StartDay",&StartDay);
    Header->SetBranchAddress("StartHour",&StartHour);
    Header->SetBranchAddress("EndYear",&EndYear);
    Header->SetBranchAddress("EndMonth",&EndMonth);
    Header->SetBranchAddress("EndDay",&EndDay);
    Header->SetBranchAddress("EndHour",&EndHour);
    
    Header->GetEntry(0);
    
    HasHeader = true;

  }//getHeaderInfo

  void getNearlineVersion(TFile* fp, std::string versionHistName){

    TH1* hist = (TH1*) fp->Get(versionHistName.c_str());

    if(!hist){
      std::cout << "NearlineSummariser: ERROR: " << fp->GetName() << " doesn't contain version number hist: " << versionHistName << std::endl;
      return;
    }

    ThisNearlineMinorVersion = hist->GetBinContent(1);
    ThisNearlineMajorVersion = hist->GetBinContent(2);

  }//getNearlineVersion

  friend std::ostream & operator << (std::ostream &os, NearlineFileInfo &rhs){
    os << "Run " << rhs.Run
       << " StartYear " << rhs.StartYear
       << " StartMonth " << rhs.StartMonth
       << " StartDay " << rhs.StartDay
       << " StartHour " << rhs.StartHour
       << " StartTime " << rhs.StartTime;
    return os;
  }//operator <<

};//NearlineFileInfo


////////////////////////////////////////////////////////////////////////////////

struct NearlineRunPlotSet{

  std::map<std::string,std::shared_ptr<TH1>> MapNameHists;
  NearlineFileInfo FileInfo;
  TFile *File;
  NearlineRunPlotSet(TFile *fp):
    FileInfo(fp),
    File(fp)
  {
    
  }//NearlinePlots
  void loadHistogram(std::string dirName, std::string histName){
    std::string newHistName = histName + "_run" + std::to_string(FileInfo.Run);
    std::string histFullName;
    if(dirName=="") histFullName = histName;
    else histFullName = dirName + "/" + histName;
    std::shared_ptr<TH1> thisHisto(getHistogram(File, histFullName, newHistName));
    MapNameHists.insert(std::pair<std::string,std::shared_ptr<TH1>>(histName, thisHisto));
  }
};//NearlineRunPlotSet

////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv){
  
  std::vector<std::string> vecFileNames;
  for(int i=1;i<argc;i++){
    std::string thisFileName(argv[i]);
    vecFileNames.push_back(thisFileName);
  }

  std::vector<NearlineRunPlotSet> VecRunPlotSet;
  std::vector<std::string> VecPlotNames;
  VecPlotNames.push_back("hped_per_event_chan_0");
  VecPlotNames.push_back("hped_per_event_chan_128");
  VecPlotNames.push_back("hped_per_event_chan_256");
  VecPlotNames.push_back("hped_per_event_chan_384");
  VecPlotNames.push_back("hped_per_event_chan_512");
  VecPlotNames.push_back("hped_per_event_chan_640");
  VecPlotNames.push_back("hped_per_event_chan_768");
  VecPlotNames.push_back("hped_per_event_chan_896");
  VecPlotNames.push_back("hped_per_event_chan_1024");
  VecPlotNames.push_back("hped_per_event_chan_1152");
  VecPlotNames.push_back("hped_per_event_chan_1280");
  VecPlotNames.push_back("hped_per_event_chan_1408");
  VecPlotNames.push_back("hped_per_event_chan_1536");
  VecPlotNames.push_back("hped_per_event_chan_1664");
  VecPlotNames.push_back("hped_per_event_chan_1792");
  VecPlotNames.push_back("hped_per_event_chan_1920");



  for(auto thisFileName: vecFileNames){
    TFile *fp = TFile::Open(thisFileName.c_str());
    if(fp){
      if(!checkFileHasDir(fp, NearlineInstanceLabel)) continue;
    }

    NearlineRunPlotSet thisRunPlotSet(fp);

    //DEBUG
    std::cout << "NearlineSummariser: INFO: "
              << thisRunPlotSet.FileInfo
              // << thisRunPlotSet.FileInfo.StartYear 
              // << " " << thisRunPlotSet.FileInfo.StartMonth
              // << " " << thisRunPlotSet.FileInfo.StartDay
              // << " " << thisRunPlotSet.FileInfo.StartHour
              // << " " << thisRunPlotSet.FileInfo.StartTime
              << "\n";

    std::cout << "NearlineSummariser: INFO: ThisNearlineVersion: "
              << thisRunPlotSet.FileInfo.ThisNearlineMajorVersion 
              << "." << thisRunPlotSet.FileInfo.ThisNearlineMinorVersion
              << " NearlineVersion (Summariser): "
              << NearlineMajorVersion 
              << "." << NearlineMinorVersion
              << std::endl;

    for(auto thisPlotName: VecPlotNames){
      thisRunPlotSet.loadHistogram(NearlineInstanceLabel, thisPlotName);
    }//VecPlotNames

    // thisRunPlotSet.loadHistogram(NearlineInstanceLabel , "hped_per_event_chan_0");
    // thisRunPlotSet.loadHistogram(NearlineInstanceLabel , "hped_per_event_chan_128");
    // thisRunPlotSet.loadHistogram(NearlineInstanceLabel , "hped_per_event_chan_256");
    // thisRunPlotSet.loadHistogram(NearlineInstanceLabel , "hped_per_event_chan_384");
    // thisRunPlotSet.loadHistogram(NearlineInstanceLabel , "hped_per_event_chan_512");
    // thisRunPlotSet.loadHistogram(NearlineInstanceLabel , "hped_per_event_chan_640");
    // thisRunPlotSet.loadHistogram(NearlineInstanceLabel , "hped_per_event_chan_768");
    // thisRunPlotSet.loadHistogram(NearlineInstanceLabel , "hped_per_event_chan_896");
    // thisRunPlotSet.loadHistogram(NearlineInstanceLabel , "hped_per_event_chan_1024");
    // thisRunPlotSet.loadHistogram(NearlineInstanceLabel , "hped_per_event_chan_1152");
    // thisRunPlotSet.loadHistogram(NearlineInstanceLabel , "hped_per_event_chan_1280");
    // thisRunPlotSet.loadHistogram(NearlineInstanceLabel , "hped_per_event_chan_1408");
    // thisRunPlotSet.loadHistogram(NearlineInstanceLabel , "hped_per_event_chan_1536");
    // thisRunPlotSet.loadHistogram(NearlineInstanceLabel , "hped_per_event_chan_1664");
    // thisRunPlotSet.loadHistogram(NearlineInstanceLabel , "hped_per_event_chan_1792");
    // thisRunPlotSet.loadHistogram(NearlineInstanceLabel , "hped_per_event_chan_1920");

    VecRunPlotSet.push_back(thisRunPlotSet);

    //    if(fp) checkFileHasDir(fp, "foo");
    //    listHistogramsInFile(thisFileName);

  }//vecFileNames

  for(auto thisRunPlotSet: VecRunPlotSet){
    std::cout << "NearlineSummariser: INFO: " << thisRunPlotSet.FileInfo << std::endl;
    
    // for(auto thisPlotName: VecPlotNames){
    //   auto thisHist =thisRunPlotSet.MapNameHists.at(thisPlotName);
    //   double mean = (thisHist)->GetMean();
    //   double rms = (thisHist)->GetRMS();
    //   std::cout << "NearlineSummariser: INFO: Name " << thisPlotName << " count " << thisHist.use_count() << " mean " << mean << " rms " << rms << std::endl;
    // }//thisPlotName

    // for(auto thisNameHistPair: thisRunPlotSet.MapNameHists){
    //   double mean = (thisNameHistPair.second)->GetMean();
    //   double rms = (thisNameHistPair.second)->GetRMS();
    //   std::cout << "NearlineSummariser: INFO: Name " << thisNameHistPair.first << " mean " << mean << " rms " << rms << std::endl;
    // }//MapNameHists
    
  }//VecRunPlotSet


  return 0;
}

////////////////////////////////////////////////////////////////////////////////

bool checkFileHasDir(TFile* fp, std::string dirName){

  TList* listOfKeys = fp->GetListOfKeys();
  TIter next(listOfKeys);
  TKey* key;
  while(( key = (TKey*)next() )){
    std::string thisKeyName(key->GetName());
    if(thisKeyName == dirName){
      //      std::cout << "NearlineSummariser: INFO: Got dirName: " << dirName << std::endl;
      return true;
    }
  }
  std::cout << "NearlineSummariser: ERROR: " << fp->GetName() << " doesn't contain a dir: " << dirName << std::endl;
  return false;
}


////////////////////////////////////////////////////////////////////////////////

TH1* getHistogram(TFile* fp, std::string histName, std::string newHistName){
  
  TH1* hist;
  TH1* outhist;

  //  std::cout << "NearlineSummariser: INFO: Getting: " << histName << std::endl;
  
  fp->GetObject(histName.c_str(), hist);

  if(!hist){
    std::cout << "NearlineSummariser: ERROR: Failed to get histogram: " << histName << "\n";
    return NULL;
  }

  //  std::cout << "NearlineSummariser: INFO: Got: " << hist->GetName() << std::endl;
  
  outhist = (TH1*) hist->Clone(newHistName.c_str());
  outhist->SetDirectory(0);

  //  std::cout << "NearlineSummariser: INFO: Returning: " << outhist->GetName() << std::endl;

  return outhist;
}
////////////////////////////////////////////////////////////////////////////////
