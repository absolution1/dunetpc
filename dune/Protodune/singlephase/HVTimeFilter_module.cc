////////////////////////////////////////////////////////////////////////
//
// HVTimeFilter class
// 
// author: Owen Goodwin
// email: owen.goodwin@manchester.ac.uk
//
// - A filter to select events between given start and end dates and times.
//    - Dates should be passed in the form YYYYMMDD
//    - Times should be passes in the form HHMMSS
//    - All dates and times should be in UTC
//
////////////////////////////////////////////////////////////////////////
// C++
extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}
#include <math.h>
#include <algorithm>
#include <iostream>
#include <fstream>
// ROOT
#include "TMath.h"
#include "TTimeStamp.h"
/// Framework 
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 

#include "lardataobj/RawData/RDTimeStamp.h"
#include "dune/DuneObj/ProtoDUNETimeStamp.h"
///filters for events, etc
namespace filter {
        class HVTimeFilter : public art::EDFilter  {
        public:
                explicit HVTimeFilter(fhicl::ParameterSet const& ); 
                virtual ~HVTimeFilter();      
                uint64_t GetRawDecoderInfo(art::Event & e);
                bool filter(art::Event& evt);
                void beginJob();
                void endJob();
                void reconfigure(fhicl::ParameterSet const& p);
        private:
                UInt_t fDateRangeLow; 
                UInt_t fTimeRangeLow; 
                UInt_t fDateRangeHigh; 
                UInt_t fTimeRangeHigh;

                std::vector<std::pair<UInt_t,UInt_t>> fTimeRanges; 

                long long RDTSTime;
                double RDTSTimeSec;
                double PrevRDTSTimeSec;  
                double RDTSTimeNano; 
                int RDTSTrigger;
                art::Handle< std::vector<raw::RDTimeStamp> > RDTimeStampHandle;

        }; //class HVTimeFilter
}
void filter::HVTimeFilter::beginJob() { }
void filter::HVTimeFilter::endJob() { }
void filter::HVTimeFilter::reconfigure(fhicl::ParameterSet const& p) {
        fDateRangeLow   = p.get<UInt_t>("DateRangeLow", 0);   // YYYYMMDD
        fTimeRangeLow   = p.get<UInt_t>("TimeRangeLow", 0);   // HHMMSS
        fDateRangeHigh  = p.get<UInt_t>("DateRangeHigh", 0);  // YYYYMMDD
        fTimeRangeHigh  = p.get<UInt_t>("TimeRangeHigh", 0);  // HHMMSS
        fTimeRanges  = p.get<std::vector <std::pair<UInt_t,UInt_t >>>("TimeRanges");

 }  


filter::HVTimeFilter::HVTimeFilter(fhicl::ParameterSet const& pset) {
        this->reconfigure(pset);
}
filter::HVTimeFilter::~HVTimeFilter() { }


uint64_t filter::HVTimeFilter::GetRawDecoderInfo(art::Event & e){
    LOG_INFO("BeamEvent") << "\n";
    LOG_INFO("BeamEvent") << "Getting Raw Decoder Info" << "\n";
    e.getByLabel("timingrawdecoder","daq",RDTimeStampHandle);
    LOG_INFO("BeamEvent") << "RDTS valid? " << RDTimeStampHandle.isValid() << "\n";
    for (auto const & RDTS : *RDTimeStampHandle){
        LOG_INFO("BeamEvent") << "High: " << RDTS.GetTimeStamp_High() << "\n";
        LOG_INFO("BeamEvent") << "Low: " << RDTS.GetTimeStamp_Low() << "\n";

        uint64_t high = RDTS.GetTimeStamp_High();
        uint64_t low  = RDTS.GetTimeStamp_Low();

        high = high << 32;
        uint64_t joined = (high | low);

        LOG_INFO("BeamEvent") << "Raw Decoder Timestamp: " << joined << "\n";

        RDTSTime = joined;
        RDTSTrigger = RDTS.GetFlags();
        LOG_INFO("BeamEvent") << "Trigger: " << RDTSTrigger << "\n";

        //Separates seconds portion of the ticks 
        //From the nanoseconds
        long long RDTSTickSec = (RDTSTime * 2) / (int)(TMath::Power(10,8));
        RDTSTickSec = RDTSTickSec * (int)(TMath::Power(10,8)) / 2;
        long long RDTSTickNano = RDTSTime - RDTSTickSec;

        //Units are 20 nanoseconds ticks
        RDTSTimeSec  = 20.e-9 * RDTSTickSec;
        RDTSTimeNano = 20.    * RDTSTickNano;


  }
  return RDTSTimeSec;
}

bool filter::HVTimeFilter::filter(art::Event &evt) {   

        TTimeStamp * evtTTS;
        evtTTS = new TTimeStamp(GetRawDecoderInfo(evt));
        // if (evtTime.timeHigh() == 0) { evtTTS = new TTimeStamp(evtTime.timeLow()); }
        // else { evtTTS = new TTimeStamp(evtTime.timeHigh(), evtTime.timeLow()); }
        std::cout << "Event time:  " << evtTTS -> AsString() << std::endl;
        // Requested time range lower end

        for (auto TimeRange : fTimeRanges){ //loop through beam side APAs

            fTimeRangeLow=TimeRange.first;
            fTimeRangeHigh=TimeRange.second;


                        // Check that input date is in correct format
            
            if (fDateRangeHigh > 99999999 || fDateRangeLow > 99999999) {
                    std::cout << "Warning: please provide date in format YYYYMMDD, event time "
                              << "filter returning false." << std::endl; 
                    return false;
            }
            if (fDateRangeHigh > 0 && fDateRangeHigh < 10000000) {
                    std::cout << "Warning: please provide date in format YYYYMMDD, event time "
                              << "filter returning false." << std::endl; 
                    return false;
            }
            if (fDateRangeLow > 0 && fDateRangeLow < 10000000) {
                    std::cout << "Warning: please provide date in format YYYYMMDD, event time "
                              << "filter returning false." << std::endl; 
                    return false;
            }
            // Check that input times are in correct format
            if (fTimeRangeHigh > 999999 || fTimeRangeLow > 999999) {
                    std::cout << "Warning: please provide time in format HHMMSS, event time "
                              << "filter returning false.1" << std::endl; 
                    return false;

            }
            //no idea what this is supposed to be doing!! will reject any times before 10am
            // 
            // if (fTimeRangeHigh > 0 && fTimeRangeHigh < 100000) {
            //         std::cout << "Warning: please provide time in format HHMMSS, event time "
            //                   << "filter returning false.2" << std::endl; 
            //         return false;
            // }
            // if (fTimeRangeLow > 0 && fTimeRangeLow < 100000) {
            //         std::cout << "Warning: please provide time in format HHMMSS, event time "
            //                   << "filter returning false.3" << std::endl; 
            //         return false;
            // }
            // Event time
            //art::Timestamp evtTime = evt.time();


                // all the checking he dies then ask if its in the correct range

                        TTimeStamp * ttsLow(nullptr); 
            if (fDateRangeLow != 0) {
                    if (fTimeRangeLow != 0) { 
                            ttsLow = new TTimeStamp(fDateRangeLow, fTimeRangeLow, 0u); 
                    }
                    else { 
                            ttsLow = new TTimeStamp(fDateRangeLow, 0u, 0u); 
                            std::cout << "Warning: No start time given for event time filter, "
                                      << "assuming 00:00:00" << std::endl;
                    }
            }
            // Requested time range higher end
            TTimeStamp * ttsHigh(nullptr);
            if (fDateRangeHigh != 0) {
                    if (fTimeRangeHigh != 0) { 
                            ttsHigh = new TTimeStamp(fDateRangeHigh, fTimeRangeHigh, 0u); 
                    }
                    else { 
                            std::cout << "Warning: No end time given for event time filter, assuming "
                                      << "23:59:59" << std::endl;
                            ttsHigh = new TTimeStamp(fDateRangeHigh, 235959u, 0u); 
                    }
            }
            // Filter decision
            std::cout << "Lower Limit:  " << ttsLow -> AsString() << std::endl;
            std::cout << "Upper Limit:  " << ttsHigh -> AsString() << std::endl;
    
                
            if (evtTTS -> GetSec() > ttsLow -> GetSec() && 
                evtTTS -> GetSec() < ttsHigh -> GetSec()) { return false; }
                    

    }
    return true;    //need this to stop the get to end error make sure to check this is working as intended
}
 DEFINE_ART_MODULE(filter::HVTimeFilter) 
