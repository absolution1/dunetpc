////////////////////////////////////////////////////////////////////////
//Adapted for ProtoDUNE by Owen Goodwin (ogoodwin@fnal.gov) from
//  \file DetectorProperties_service.cc
//
////////////////////////////////////////////////////////////////////////
// Framework includes
// LArSoft includes
#include "dune/Protodune/singlephase/DetectorServices/Services/DetectorPropertiesServiceProtoDUNEsp.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include "lardataalg/DetectorInfo/LArProperties.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "lardata/DetectorInfoServices/ServicePack.h" // lar::extractProviders()
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
// Art includes
#include "art_root_io/RootDB/SQLite3Wrapper.h"
#include "fhiclcpp/make_ParameterSet.h"


#include "IFDH_service.h"

#include "TFile.h"

namespace spdp{
  //--------------------------------------------------------------------
  DetectorPropertiesServiceProtoDUNEsp::DetectorPropertiesServiceProtoDUNEsp
    (fhicl::ParameterSet const& pset, art::ActivityRegistry &reg)
    : fInheritNumberTimeSamples(pset.get<bool>("InheritNumberTimeSamples", false))
  {
    // Dummy handle to ensure that the DetectorClocksService callbacks are invoked first
    art::ServiceHandle<detinfo::DetectorClocksService const>{};

    // Register for callbacks.
    reg.sPostOpenFile.watch   (this, &DetectorPropertiesServiceProtoDUNEsp::postOpenFile);
    reg.sPreBeginRun.watch (this, &DetectorPropertiesServiceProtoDUNEsp::preBeginRun);


    isNewRun=true;
/*
    // obtain the required dependency service providers and create our own
    const geo::GeometryCore* geo = lar::providerFrom<geo::Geometry>();

    const detinfo::LArProperties* lp = lar::providerFrom<detinfo::LArPropertiesService>();

    const detinfo::DetectorClocks* clks = art::ServiceHandle<detinfo::DetectorClocksService const>()->;

    fProp = std::make_unique<detinfo::DetectorPropertiesStandard>(pset,geo,lp,clks);
    */
    fProp = std::make_unique<spdp::DetectorPropertiesProtoDUNEsp>(pset,
      lar::extractProviders<
        geo::Geometry,
        detinfo::LArPropertiesService
        >(),
        std::set<std::string>({ "InheritNumberTimeSamples" })
      );

    // at this point we need and expect the provider to be fully configured
    fProp->CheckIfConfigured();

    // Save the parameter set.
    fPS = pset;


  }
  //--------------------------------------------------------------------
  void DetectorPropertiesServiceProtoDUNEsp::reconfigure(fhicl::ParameterSet const& p)
  {
    fProp->ValidateAndConfigure(p, { "InheritNumberTimeSamples" });

    // Save the parameter set.
    fPS = p;
    return;
  }

    void DetectorPropertiesServiceProtoDUNEsp::preBeginRun(const art::Run& run)
  {
      isNewRun=true;
      std::cout<<"New run, Num is: "<<run.run()<<" Updating DetectorProperties."<<std::endl;
      // make it into a TTimeStamp

      //Get run dependent temperature  correction

      //  Before 11/17/18, temperature = 87.68 and is stable
      // Between 11/17/18 and 3/1/19, temperature = 87.36 and has a rather large fluctuation. first run in this period 5903
      // After 3/1/19, temperature = 87.65 and is stable.  first run 6930
      fProp->UpdateTemp(run.run());
    }


  //--------------------------------------------------------------------
  //  Callback called after input file is opened.
  void DetectorPropertiesServiceProtoDUNEsp::postOpenFile(const std::string& filename)
  {








    // std::cout<<"new run?"<<isNewRun<<std::endl;
    if(isNewRun){

      auto start = filename.rfind("/"); //finds the final "/"
      if (start == std::string::npos)
        {
          start = 0;
        }
      else
        {
          start += 1;
        }

      int end = filename.length(); //last postion
      std::string filename_s=(filename.substr(start, end-start)); //creates string of just file name, not path
      // std::cout<<"filename:"<<filename_s<<std::endl;
      art::ServiceHandle<ifdh_ns::IFDH> ifdh;
      fProp->UpdateHV(filename_s); //pass file name to be able update HV value from MetaData (if requested)
      auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
      fProp->UpdateReadoutWindowSize(clockData, filename_s); //update Readout window value from filename
      isNewRun=false;

    }


    // Use this method to figure out whether to inherit configuration
    // parameters from previous jobs.
    // There is no way currently to correlate parameter sets saved in
    // sqlite RootFileDB with process history (from MetaData tree).
    // Therefore, we use the approach of scanning every historical
    // parameter set in RootFileDB, and finding all parameter sets
    // that appear to be DetectorPropertiesService configurations.  If all
    // historical parameter sets are in agreement about the value of
    // an inherited parameter, then we accept the historical value,
    // print a message, and override the configuration parameter.  In
    // cases where the historical configurations are not in agreement
    // about the value of an inherited parameter, we ignore any
    // historical parameter values that are the same as the current
    // configured value of the parameter (that is, we resolve the
    // conflict in favor of parameters values that are different than
    // the current configuration).  If two or more historical values
    // differ from the current configuration, throw an exception.
    // Note that it is possible to give precendence to the current
    // configuration by disabling inheritance for that configuration
    // parameter.
    // Don't do anything if no parameters are supposed to be inherited.
    if(!fInheritNumberTimeSamples) return;
    // The only way to access art service metadata from the input file
    // is to open it as a separate TFile object.  Do that now.
    if(filename.size() != 0) {
      TFile* file = TFile::Open(filename.c_str(), "READ");
      if(file != 0 && !file->IsZombie() && file->IsOpen()) {
        // Open the sqlite datatabase.
        art::SQLite3Wrapper sqliteDB(file, "RootFileDB");
        // Loop over all stored ParameterSets.
        unsigned int iNumberTimeSamples = 0;  // Combined value of NumberTimeSamples.
        unsigned int nNumberTimeSamples = 0;  // Number of NumberTimeSamples parameters seen.
        sqlite3_stmt * stmt = 0;
        sqlite3_prepare_v2(sqliteDB, "SELECT PSetBlob from ParameterSets;", -1, &stmt, NULL);
        while (sqlite3_step(stmt) == SQLITE_ROW) {
          auto ps = fhicl::ParameterSet::make(reinterpret_cast<char const *>(sqlite3_column_text(stmt, 0)));
          // Is this a DetectorPropertiesService parameter set?
          bool psok = isDetectorPropertiesServiceProtoDUNEsp(ps);
          if(psok) {
            // Check NumberTimeSamples
            //            if(fInheritNumberTimeSamples) {
            unsigned int newNumberTimeSamples = ps.get<unsigned int>("NumberTimeSamples");

            // Ignore parameter values that match the current configuration.

            if(newNumberTimeSamples != fPS.get<unsigned int>("NumberTimeSamples")) {
              if(nNumberTimeSamples == 0)
                iNumberTimeSamples = newNumberTimeSamples;
              else if(newNumberTimeSamples != iNumberTimeSamples) {
                throw cet::exception(__FUNCTION__)
                  << "Historical values of NumberTimeSamples do not agree: "
                  << iNumberTimeSamples << " " << newNumberTimeSamples << "\n" ;
              }
              ++nNumberTimeSamples;
              //            }
            }
          }
        }
        // Done looping over parameter sets.
        // Now decide which parameters we will actually override.
        if(// fInheritNumberTimeSamples &&
           nNumberTimeSamples != 0 &&
           iNumberTimeSamples != fProp->NumberTimeSamples()) {
          mf::LogInfo("DetectorPropertiesServiceProtoDUNEsp")
            << "Overriding configuration parameter NumberTimeSamples using historical value.\n"
            << "  Configured value:        " << fProp->NumberTimeSamples() << "\n"
            << "  Historical (used) value: " << iNumberTimeSamples << "\n";
          fProp->SetNumberTimeSamples(iNumberTimeSamples);
        }
      }
      // Close file.
      if(file != 0) {
        if(file->IsOpen())
          file->Close();
        delete file;
      }
    }

  }
  //--------------------------------------------------------------------
  //  Determine whether a parameter set is a DetectorPropertiesService configuration.

  bool DetectorPropertiesServiceProtoDUNEsp::isDetectorPropertiesServiceProtoDUNEsp
    (const fhicl::ParameterSet& ps) const
  {
    // This method uses heuristics to determine whether the parameter
    // set passed as argument is a DetectorPropertiesService configuration
    // parameter set.

    return
         (ps.get<std::string>("service_type", "") == "DetectorPropertiesService")
      && (ps.get<std::string>("service_provider", "") == "DetectorPropertiesServiceProtoDUNEsp")
      ;
#if 0
    // old heuristics here:
    std::string s;
    double d;
    int i;
    unsigned int u;

    bool result = !ps.get_if_present("module_label", s);
    result = result && ps.get_if_present("TriggerOffset", i);
    result = result && ps.get_if_present("SamplingRate", d);
    result = result && ps.get_if_present("NumberTimeSamples", u);
    result = result && ps.get_if_present("ReadOutWindowSize", u);

    return result;
#endif // 0
  }
} // namespace spdf
DEFINE_ART_SERVICE_INTERFACE_IMPL(spdp::DetectorPropertiesServiceProtoDUNEsp, detinfo::DetectorPropertiesService)
