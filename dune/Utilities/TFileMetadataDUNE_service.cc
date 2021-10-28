////////////////////////////////////////////////////////////////////////
// Name:  TFileMetadataDUNE_service.cc.  
//
// Purpose:  generate DUNE-specific sam metadata for root Tfiles (histogram or ntuple files).
//
// FCL parameters: GenerateTFileMetadata: This needs to be set to "true" in the fcl file
//				    to generate metadata (default value: false)
//	     	   dataTier: Currrently this needs to be parsed by the user
//		     	     for ntuples, dataTier = root-tuple; 
//		             for histos, dataTier = root-histogram
//		             (default value: root-tuple)
//	           fileFormat: This is currently specified by the user,
//			       the fileFormat for Tfiles is "root" (default value: root)	
//
// Other notes: 1. This service uses the ART's standard file_catalog_metadata service
//		to extract some of the common (common to both ART and TFile outputs)
//	        job-specific metadata parameters, so, it is important to call this
//  		service in your fcl file
//		stick this line in your "services" section of fcl file:
//		FileCatalogMetadata:  @local::art_file_catalog_mc
//	
//              2. When you call FileCatalogMetadata service in your fcl file, and if 
//		you have (art) root Output section in your fcl file, and if you do not  
//		have "dataTier" specified in that section, then this service will throw
//		an exception. To avoid this, either remove the entire root Output section
//		in your fcl file (and remove art stream output from your end_paths) or 
//		include appropriate dataTier information in the section.If you are only
//		running analysis job, best way is to not include any art root Output section.
//
//	        3. This service is exclusively written to work with production (in other
//		words for jobs submitted through grid). Some of the metadata parameters
//		(output TFileName, filesize, Project related details) are captured/updated
//		during and/or after the workflow. 
//	
//
// Created:  1-Nov-2017,  T. Junk
//  based on the MicroBooNE example by S. Gollapinni
//
////////////////////////////////////////////////////////////////////////

#include <string>
#include <sstream>
#include <iomanip>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include "dune/Utilities/TFileMetadataDUNE.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include "dune/Utilities/FileCatalogMetadataDUNE.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/System/FileCatalogMetadata.h"
#include "art/Utilities/OutputFileInfo.h"
#include "art_root_io/RootDB/SQLite3Wrapper.h"
#include "art_root_io/RootDB/SQLErrMsg.h"
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTimeStamp.h"
#include <ctime>
#include <stdio.h>
#include <time.h>

using namespace std;


//--------------------------------------------------------------------

// Constructor.
util::TFileMetadataDUNE::TFileMetadataDUNE(fhicl::ParameterSet const& pset, 
					   art::ActivityRegistry& reg):
  fGenerateTFileMetadata(false)							
{
  reconfigure(pset);

  // Register for callbacks.
 
  reg.sPostBeginJob.watch(this, &TFileMetadataDUNE::postBeginJob);
  reg.sPostOpenFile.watch(this, &TFileMetadataDUNE::postOpenFile);
  reg.sPostEndJob.watch(this, &TFileMetadataDUNE::postEndJob);
  reg.sPostProcessEvent.watch(this, &TFileMetadataDUNE::postEvent);
  reg.sPostBeginSubRun.watch(this, &TFileMetadataDUNE::postBeginSubRun);
}

//--------------------------------------------------------------------
// Destructor.
util::TFileMetadataDUNE::~TFileMetadataDUNE()
{
}

//--------------------------------------------------------------------
// Set service paramters
void util::TFileMetadataDUNE::reconfigure(fhicl::ParameterSet const& pset)
{    
  // Get parameters.
  fGenerateTFileMetadata = pset.get<bool>("GenerateTFileMetadata", false);
  fJSONFileName = pset.get<std::string>("JSONFileName");
    
  if (!fGenerateTFileMetadata) return;

  md.fdata_tier  	   = pset.get<std::string>("dataTier","root-tuple");	   
  md.ffile_format	   = pset.get<std::string>("fileFormat","root");              
}

//--------------------------------------------------------------------
// PostBeginJob callback.
// Insert per-job metadata via TFileMetadata service.
void util::TFileMetadataDUNE::postBeginJob()
{ 
  // only generate metadata when this is true
  if (!fGenerateTFileMetadata) return;
    
  // get the start time  
  md.fstart_time = time(0); 
  
  // Get art metadata service and extract paramters from there
  art::ServiceHandle<art::FileCatalogMetadata> artmds;
  
  art::FileCatalogMetadata::collection_type artmd;
  artmds->getMetadata(artmd);
  
  for(auto const & d : artmd)
    mdmap[d.first] = d.second;
    
  std::map<std::string,std::string>::iterator it;
  
  // if a certain paramter/key is not found, assign an empty string value to it
  
  if ((it=mdmap.find("applicationFamily"))!=mdmap.end()) std::get<0>(md.fapplication) = it->second;
  else std::get<0>(md.fapplication) = "\" \"";   

  if ((it=mdmap.find("application.family"))!=mdmap.end()) std::get<0>(md.fapplication) = it->second;
  else std::get<0>(md.fapplication) = "\" \"";   
   
  if ((it=mdmap.find("process_name"))!=mdmap.end()) std::get<1>(md.fapplication) = it->second;
  else std::get<1>(md.fapplication) = "\" \"";  

  if ((it=mdmap.find("art.process_name"))!=mdmap.end()) std::get<1>(md.fapplication) = it->second;
  else std::get<1>(md.fapplication) = "\" \"";  
  
  if ((it=mdmap.find("applicationVersion"))!=mdmap.end()) std::get<2>(md.fapplication) = it->second;
  else  std::get<2>(md.fapplication) = "\" \"";   

  if ((it=mdmap.find("application.version"))!=mdmap.end()) std::get<2>(md.fapplication) = it->second;
  else  std::get<2>(md.fapplication) = "\" \"";   
  
  if ((it=mdmap.find("group"))!=mdmap.end()) md.fgroup = it->second;
  else md.fgroup = "\" \"";   
    
  if ((it=mdmap.find("file_type"))!=mdmap.end()) md.ffile_type = it->second;
  else  md.ffile_type = "\" \"";  
    
  if ((it=mdmap.find("run_type"))!=mdmap.end()) frunType = it->second;
  else frunType = "\" \"";         	     	

  if ((it=mdmap.find("art.run_type"))!=mdmap.end()) frunType = it->second;
  else frunType = "\" \"";         	     	
}


//--------------------------------------------------------------------        
// PostOpenFile callback.
void util::TFileMetadataDUNE::postOpenFile(std::string const& fn)
{
  if (!fGenerateTFileMetadata) return;
  
  // save parent input files here
  md.fParents.insert(fn);
  
}

//--------------------------------------------------------------------  	
// PostEvent callback.
void util::TFileMetadataDUNE::postEvent(art::Event const& evt, art::ScheduleContext)
{
 
  if(!fGenerateTFileMetadata) return;	
  
  art::RunNumber_t run = evt.run();
  art::SubRunNumber_t subrun = evt.subRun();
  art::EventNumber_t event = evt.event();
  art::SubRunID srid = evt.id().subRunID();
      
  // save run, subrun and runType information once every subrun    
  if (fSubRunNumbers.count(srid) == 0){
    fSubRunNumbers.insert(srid);
    md.fruns.push_back(make_tuple(run, subrun, frunType));   
  }
  
  // save the first event
  if (md.fevent_count == 0) md.ffirst_event = event;
  md.flast_event = event;
  // event counter
  ++md.fevent_count;
    
}

//--------------------------------------------------------------------  	
// PostSubRun callback.
void util::TFileMetadataDUNE::postBeginSubRun(art::SubRun const& sr)
{

  if(!fGenerateTFileMetadata) return;

  art::RunNumber_t run = sr.run();
  art::SubRunNumber_t subrun = sr.subRun();
  art::SubRunID srid = sr.id();

  // save run, subrun and runType information once every subrun
  if (fSubRunNumbers.count(srid) == 0){
    fSubRunNumbers.insert(srid);
    md.fruns.push_back(make_tuple(run, subrun, frunType));
  }
}

//--------------------------------------------------------------------
// PostCloseFile callback.
void util::TFileMetadataDUNE::postEndJob()
{
	
  // Do nothing if generating TFile metadata is disabled.	
  if(!fGenerateTFileMetadata) return;	
  
  // get metadata from the FileCatalogMetadataDUNE service, which is filled on its construction
        
  art::ServiceHandle<util::FileCatalogMetadataDUNE> paramhandle; 
  md.fMCGenerators =		paramhandle->MCGenerators();			  
  md.fMCOscillationP =          paramhandle->MCOscillationP();         
  md.fMCTriggerListVersion =	paramhandle->MCTriggerListVersion();		  
  md.fMCBeamEnergy =		paramhandle->MCBeamEnergy();			  
  md.fMCBeamFluxID =		paramhandle->MCBeamFluxID();			  
  md.fMCName =	                paramhandle->MCName();	                
  md.fMCDetectorType =          paramhandle->MCDetectorType();         
  md.fMCNeutrinoFlavors =	paramhandle->MCNeutrinoFlavors();	
  md.fMCMassHierarchy =         paramhandle->MCMassHierarchy();        
  md.fMCMiscellaneous =         paramhandle->MCMiscellaneous();        
  md.fMCGeometryVersion =	paramhandle->MCGeometryVersion();	
  md.fMCOverlay =		paramhandle->MCOverlay();		
  md.fDataRunMode =		paramhandle->DataRunMode();			  
  md.fDataDetectorType =        paramhandle->DataDetectorType();	
  md.fDataName =		paramhandle->DataName();		
  md.fStageName =               paramhandle->StageName();              

  //update end time
  md.fend_time = time(0);

  // convert start and end times into time format: Year-Month-DayTHours:Minutes:Seconds
  char endbuf[80], startbuf[80];
  struct tm tstruct;
  tstruct = *localtime(&md.fend_time);
  strftime(endbuf,sizeof(endbuf),"%Y-%m-%dT%H:%M:%S",&tstruct);
  tstruct = *localtime(&md.fstart_time);
  strftime(startbuf,sizeof(startbuf),"%Y-%m-%dT%H:%M:%S",&tstruct);
  
  // open a json file and write everything from the struct md complying to the 
  // samweb json format. This json file holds the below information temporarily. 
  // If you submitted a grid job invoking this service, the information from 
  // this file is appended to a final json file and this file will be removed
  
  std::ofstream jsonfile;
  jsonfile.open(fJSONFileName);
  jsonfile<<"{\n  \"application\": {\n    \"family\": "<<std::get<0>(md.fapplication)<<",\n    \"name\": ";
  jsonfile<<std::get<1>(md.fapplication)<<",\n    \"version\": "<<std::get<2>(md.fapplication)<<"\n  },\n  ";
  jsonfile<<"\"data_tier\": \""<<md.fdata_tier<<"\",\n  ";
  jsonfile<<"\"event_count\": "<<md.fevent_count<<",\n  ";
  //jsonfile<<"\"fcl.name\": \""<<md.ffcl_name<<"\",\n  ";
  //jsonfile<<"\"fcl.version\":  \""<<md.ffcl_version<<"\",\n  ";
  jsonfile<<"\"file_format\": \""<<md.ffile_format<<"\",\n  ";
  jsonfile<<"\"file_type\": "<<md.ffile_type<<",\n  ";
  jsonfile<<"\"first_event\": "<<md.ffirst_event<<",\n  ";
  jsonfile<<"\"group\": "<<md.fgroup<<",\n  ";
  jsonfile<<"\"last_event\": "<<md.flast_event<<",\n  ";
  //if (md.fdataTier != "generated"){
  unsigned int c=0;
  jsonfile<<"\"parents\": [\n";
  for(auto parent : md.fParents) {
    c++;
    size_t n = parent.find_last_of('/');
    size_t f1 = (n == std::string::npos ? 0 : n+1);
    jsonfile<<"    {\n     \"file_name\": \""<<parent.substr(f1)<<"\"\n    }";
    if (md.fParents.size()==1 || c==md.fParents.size()) jsonfile<<"\n";
    else jsonfile<<",\n"; 
  }      
  jsonfile<<"  ],\n  "; 
  //}   
  c=0;
  jsonfile<<"\"runs\": [\n";
  for(auto &t : md.fruns){
    c++;
    jsonfile<<"    [\n     "<<std::get<0>(t)<<",\n     "<<std::get<1>(t)<<",\n     "<<std::get<2>(t)<<"\n    ]";
    if (md.fruns.size()==1 || c==md.fruns.size()) jsonfile<<"\n";
    else jsonfile<<",\n"; 
  }
  jsonfile<<"  ],\n";          

  if (md.fMCGenerators!="") jsonfile << "\"lbne_MC.generators\": \"" << md.fMCGenerators << "\",\n";
  if (md.fMCOscillationP!="") jsonfile << "\"lbne_MC.oscillationP\": \"" << md.fMCOscillationP << "\",\n";
  if (md.fMCTriggerListVersion!="") jsonfile << "\"lbne_MC.trigger-list-version\": \"" << md.fMCTriggerListVersion << "\",\n";
  if (md.fMCBeamEnergy!="") jsonfile << "\"lbne_MC.beam_energy\": \"" << md.fMCBeamEnergy << "\",\n";
  if (md.fMCBeamFluxID!="") jsonfile << "\"lbne_MC.beam_flux_ID\": \"" << md.fMCBeamFluxID << "\",\n";
  if (md.fMCName!="") jsonfile << "\"lbne_MC.name\": \"" << md.fMCName << "\",\n";
  if (md.fMCDetectorType!="") jsonfile << "\"lbne_MC.detector_type\": \"" << md.fMCDetectorType << "\",\n";
  if (md.fMCNeutrinoFlavors!="") jsonfile << "\"lbne_MC.neutrino_flavors\": \"" << md.fMCNeutrinoFlavors << "\",\n";
  if (md.fMCMassHierarchy!="") jsonfile << "\"lbne_MC.mass_hierarchy\": \"" << md.fMCMassHierarchy << "\",\n";
  if (md.fMCMiscellaneous!="") jsonfile << "\"lbne_MC.miscellaneous\": \"" << md.fMCMiscellaneous << "\",\n";
  if (md.fMCGeometryVersion!="") jsonfile << "\"lbne_MC.geometry_version\": \"" << md.fMCGeometryVersion << "\",\n";
  if (md.fMCOverlay!="") jsonfile << "\"lbne_MC.overlay\": \"" << md.fMCOverlay << "\",\n";
  if (md.fDataRunMode!="") jsonfile << "\"lbne_data.run_mode\": \"" << md.fDataRunMode << "\",\n";
  if (md.fDataDetectorType!="") jsonfile << "\"lbne_data.detector_type\": \"" << md.fDataDetectorType << "\",\n";
  if (md.fDataName!="") jsonfile << "\"lbne_data.name\": \"" << md.fDataName << "\",\n";
  // fStageName appears not to be in our metadata spec

  // put these at the end because we know they'll be there and the last one needs to not have a comma

  jsonfile<<"\"start_time\": \""<<startbuf<<"\",\n";
  jsonfile<<"\"end_time\": \""<<endbuf<<"\"\n";
  
  //jsonfile<<"  \"ub_project.name\": \""<<md.fproject_name<<"\",\n  ";
  //jsonfile<<"\"ub_project.stage\": \""<<md.fproject_stage;
  //jsonfile<<"\",\n  \"ub_project.version\": \""<<md.fproject_version<<"\"\n";
  
  jsonfile<<"}\n";
  jsonfile.close();  
}



namespace util{

  DEFINE_ART_SERVICE(util::TFileMetadataDUNE)
}//namespace util
