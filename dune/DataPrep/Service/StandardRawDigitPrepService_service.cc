// StandardRawDigitPrepService_service.cc

#include "StandardRawDigitPrepService.h"
#include <iostream>
#include <sstream>
#include <set>
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "lardataobj/RawData/RawDigit.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
#include "dune/DuneInterface/AdcChannelData.h"
#include "dune/DuneInterface/ChannelMappingService.h"
#include "dune/DuneInterface/RawDigitExtractService.h"
#include "dune/DuneInterface/AdcMitigationService.h"
#include "dune/DuneInterface/AdcSignalFindingService.h"
#include "dune/DuneInterface/AdcNoiseRemovalService.h"
#include "dune/DuneInterface/PedestalEvaluationService.h"
#include "dune/DuneInterface/AdcDeconvolutionService.h"
#include "dune/DuneInterface/AdcRoiBuildingService.h"
#include "dune/DuneInterface/AdcWireBuildingService.h"
#include "dune/DuneInterface/AdcChannelDataCopyService.h"
#include "dune/ArtSupport/DuneToolManager.h"
#include "dune/DuneInterface/Tool/AdcChannelTool.h"

using std::string;
using std::cout;
using std::endl;
using std::vector;
using std::set;
using std::ostringstream;
using raw::RawDigit;

//**********************************************************************

StandardRawDigitPrepService::
StandardRawDigitPrepService(fhicl::ParameterSet const& pset, art::ActivityRegistry&)
: m_LogLevel(1),
  m_ChannelStatusOnline(false),
  m_WiresWithoutROIFlag(2),
  m_DoDump(false), m_DumpChannel(0), m_DumpTick(0),
  m_pChannelMappingService(0),
  m_pChannelStatusProvider(nullptr),
  m_pExtractSvc(nullptr),
  m_pmitigateSvc(nullptr),
  m_pAdcSignalFindingService(nullptr),
  m_pNoiseRemoval(nullptr),
  m_pPedestalEvaluation(nullptr),
  m_pDeconvolutionService(nullptr),
  m_pRoiBuildingService(nullptr),
  m_pWireBuildingService(nullptr),
  m_pAdcChannelDataCopyService(nullptr) {
  const string myname = "StandardRawDigitPrepService::ctor: ";
  pset.get_if_present<int>("LogLevel", m_LogLevel);
  m_SkipBad        = pset.get<bool>("SkipBad");
  m_SkipNoisy      = pset.get<bool>("SkipNoisy");
  pset.get_if_present<bool>("ChannelStatusOnline", m_ChannelStatusOnline);
  m_DoMitigation = pset.get<bool>("DoMitigation");
  m_DoEarlySignalFinding = pset.get<bool>("DoEarlySignalFinding");
  m_DoNoiseRemoval       = pset.get<bool>("DoNoiseRemoval");
  m_DoPedestalAdjustment = pset.get<bool>("DoPedestalAdjustment");
  m_DoDeconvolution      = pset.get<bool>("DoDeconvolution");
  m_DoROI                = pset.get<bool>("DoROI");
  m_DoWires              = pset.get<bool>("DoWires");
  m_DoIntermediateStates = pset.get<bool>("DoIntermediateStates");
  pset.get_if_present<bool>("DoDump", m_DoDump);
  pset.get_if_present<bool>("DoDump", m_DoDump);
  pset.get_if_present<unsigned int>("WiresWithoutROIFlag", m_WiresWithoutROIFlag);
  pset.get_if_present<unsigned int>("DumpTick", m_DumpTick);
  pset.get_if_present<vector<string>>("DisplayTools", m_DisplayTools);
  if ( m_LogLevel ) cout << myname << "Fetching extract service." << endl;
  m_pExtractSvc = &*art::ServiceHandle<RawDigitExtractService>();
  if ( m_SkipBad || m_SkipNoisy ) {
    if ( m_ChannelStatusOnline ) {
      if ( m_LogLevel ) cout << myname << "Fetching channel mapping service." << endl;
      m_pChannelMappingService = &*art::ServiceHandle<ChannelMappingService>();
      if ( m_LogLevel ) cout << myname << "  Channel mapping service: @"
                             << m_pChannelMappingService << endl;
    }
    if ( m_LogLevel ) cout << myname << "Fetching channel status provider." << endl;
    m_pChannelStatusProvider = &art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider();
    if ( m_LogLevel ) cout << myname << "  Channel status provider: @"
                           <<  m_pChannelStatusProvider << endl;
  }
  if ( m_LogLevel ) cout << myname << "  Extract service: @" <<  m_pExtractSvc << endl;
  if ( m_DoMitigation ) {
    if ( m_LogLevel ) cout << myname << "Fetching mitigation service." << endl;
    m_pmitigateSvc = &*art::ServiceHandle<AdcMitigationService>();
    if ( m_LogLevel ) cout << myname << "  Mitigation service: @" <<  m_pmitigateSvc << endl;
  }
  if ( m_DoEarlySignalFinding ) {
    if ( m_LogLevel ) cout << myname << "Fetching signal finding service." << endl;
    m_pAdcSignalFindingService = &*art::ServiceHandle<AdcSignalFindingService>();
    if ( m_LogLevel ) cout << myname << "  Signal finding service: @" <<  m_pAdcSignalFindingService << endl;
  }
  if ( m_DoNoiseRemoval ) {
    if ( m_LogLevel ) cout << myname << "Fetching noise removal service." << endl;
    m_pNoiseRemoval = &*art::ServiceHandle<AdcNoiseRemovalService>();
    if ( m_LogLevel ) cout << myname << "  Noise removal service: @" <<  m_pNoiseRemoval << endl;
  }
  if ( m_DoPedestalAdjustment ) {
    if ( m_LogLevel ) cout << myname << "Fetching pedestal evaluation service." << endl;
    m_pPedestalEvaluation = &*art::ServiceHandle<PedestalEvaluationService>();
    if ( m_LogLevel ) cout << myname << "  Pedestal evalution service: @" <<  m_pPedestalEvaluation << endl;
  }
  if ( m_DoDeconvolution ) {
    if ( m_LogLevel ) cout << myname << "Fetching deconvolution service." << endl;
    m_pDeconvolutionService = &*art::ServiceHandle<AdcDeconvolutionService>();
    if ( m_LogLevel ) cout << myname << "  Deconvolution service: @" <<  m_pDeconvolutionService << endl;
  }
  if ( m_DoROI ) {
    if ( m_LogLevel ) cout << myname << "Fetching ROI building service." << endl;
    m_pRoiBuildingService = &*art::ServiceHandle<AdcRoiBuildingService>();
    if ( m_LogLevel ) cout << myname << "  ROI building service: @" <<  m_pRoiBuildingService << endl;
  }
  if ( m_DoWires ) {
    if ( m_LogLevel ) cout << myname << "Fetching wire building service." << endl;
    m_pWireBuildingService = &*art::ServiceHandle<AdcWireBuildingService>();
    if ( m_LogLevel ) cout << myname << "  Wire building service: @" <<  m_pWireBuildingService << endl;
  }
  if ( m_DoIntermediateStates > 0 ) {
    if ( m_LogLevel ) cout << myname << "Fetching intermediate state copying building service." << endl;
    m_pAdcChannelDataCopyService = &*art::ServiceHandle<AdcChannelDataCopyService>();
    if ( m_LogLevel ) cout << myname << "  Intermediate state copying service: @" <<  m_pAdcChannelDataCopyService << endl;
  }
  if ( m_DoWires && !m_DoROI ) {
    if ( m_WiresWithoutROIFlag >= 1 ) cout << myname << "WARNING: Wire building requested without ROI building." << endl;
    if ( m_WiresWithoutROIFlag >= 3 ) abort();
  }
  if ( m_DisplayTools.size() ) {
    DuneToolManager* ptm = DuneToolManager::instance("");
    if ( ptm == nullptr ) {
      cout << myname << "ERROR: Unable to retrieve tool manaager." << endl;
    } else {
      cout << myname << "  Fetching display tools: " << endl;
      for ( string tname : m_DisplayTools ) {
        if ( m_LogLevel ) cout << myname << "     Fetching " << tname << endl;
        auto padv = ptm->getPrivate<AdcChannelTool>(tname);
        if ( padv ) {
          if ( m_LogLevel ) cout << myname << "    Display tool " << tname << ": @" << padv.get() << endl;
          m_DisplayToolPtrs.push_back(std::move(padv));
          m_FoundDisplayToolNames.push_back(tname);
        } else {
          cout << myname << "ERROR: Unable to retrieve display tool " << tname << endl;
        }
      }
    }
  }
  if ( m_LogLevel >=1 ) print(cout, myname);
}

//**********************************************************************

int StandardRawDigitPrepService::
prepare(AdcChannelDataMap& datamap,
        std::vector<recob::Wire>* pwires, WiredAdcChannelDataMap* pintStates) const {
  const string myname = "StandardRawDigitPrepService:prepare: ";
  // Extract digits.
  if ( m_LogLevel >= 2 ) {
    cout << myname << "Processing digits..." << endl;
    cout << myname << "  # input digits: " << datamap.size() << endl;
  }
  int nbad = 0;
  unsigned int ichan = m_DumpChannel;
  unsigned int isig = m_DumpTick;
  set<string> snames;
  vector<AdcChannel> skipChannels;   // Channels to drop from the data map.
  for ( AdcChannelDataMap::value_type& iacd : datamap ) {
    AdcChannel chan = iacd.first;
    AdcChannelData& data = iacd.second;
    if ( m_LogLevel >= 3 ) cout << myname << "Processing digit for channel " << chan << endl;
    if ( data.digit == nullptr ) {
      if ( m_LogLevel >= 2 ) cout << myname << "Skipping null digit." << endl;
      skipChannels.push_back(chan);
      continue;
    }
    m_pExtractSvc->extract(data);
    if ( chan != data.digit->Channel() || chan != data.channel ) {
      cout << myname << "ERROR: Inconsistent channel number!" << endl;
      skipChannels.push_back(chan);
      continue;
    }
    if ( m_SkipBad || m_SkipNoisy ) {
      unsigned int chanstat = chan;
      if ( m_ChannelStatusOnline ) {
        unsigned int chanon = m_pChannelMappingService->online(chan);
        chanstat = chanon;
      }
      if ( m_SkipBad && m_pChannelStatusProvider->IsBad(chanstat) ) {
        if ( m_LogLevel >= 3 ) cout << myname << "Skipping bad channel " << chanstat << endl;
        skipChannels.push_back(chan);
        continue;
      }
      if ( m_SkipNoisy && m_pChannelStatusProvider->IsNoisy(chanstat) ) {
        if ( m_LogLevel >= 3 ) cout << myname << "Skipping noisy channel " << chanstat << endl;
        skipChannels.push_back(chan);
        continue;
      }
    }
    string state = "extracted";
    if ( pintStates != nullptr && pintStates->hasData(state) && m_pAdcChannelDataCopyService != nullptr ) {
      if ( m_LogLevel >= 3 ) cout << myname << "Saving intermediate state " << state << "." << endl;
      m_pAdcChannelDataCopyService->copy(data, pintStates->dataMaps[state][chan]);
      snames.insert(state);
    }
    if ( m_DoMitigation ) {
      m_pmitigateSvc->update(data);
      string state = "mitigated";
      if ( pintStates != nullptr && pintStates->hasData(state) && m_pAdcChannelDataCopyService != nullptr ) {
        if ( m_LogLevel >= 3 ) cout << myname << "Saving intermediate state " << state << "." << endl;
        m_pAdcChannelDataCopyService->copy(data, pintStates->dataMaps[state][chan]);
        snames.insert(state);
      }
    }
    if ( m_DoEarlySignalFinding ) {
      m_pAdcSignalFindingService->find(data);
    }
  }
  // Remove the skipped channels from the data map.
  for ( AdcChannel chan : skipChannels ) {
    datamap.erase(chan);
  }
  if ( m_DoDump ) {
    cout << myname << "Dumping channel " << m_DumpChannel << ", Tick " << isig << endl;
    if ( datamap.find(m_DumpChannel) == datamap.end() ) {
      if ( datamap.size() == 0 ) {
        cout << "Prepared data is empty." << endl;
      } else {
        cout << "Prepared data does not include channel " << m_DumpChannel << endl;
        cout << "  First channel is " << datamap.begin()->first << endl;
        cout << "   Last channel is " << datamap.end()->first << endl;
      }
    } else {
      cout << myname << "    Pedestal: " << datamap[ichan].pedestal << endl;
      if ( isig >= datamap[ichan].raw.size() ) {
        cout << myname << "Raw data does not include tick. Size is " << datamap[ichan].raw.size() << "." << endl;
      } else {
        cout << myname << "         raw: " << datamap[ichan].raw[isig] << endl;
        cout << myname << "        flag: " << datamap[ichan].flags[isig] << endl;
        cout << myname << "   After ext: " << datamap[ichan].samples[isig] << endl;
      }
    }
  }
  if ( m_DoNoiseRemoval ) {
    m_pNoiseRemoval->update(datamap);
    string state = "noiseRemoved";
    if ( pintStates != nullptr && pintStates->hasData(state) && m_pAdcChannelDataCopyService != nullptr ) {
      if ( m_LogLevel >= 3 ) cout << myname << "Saving intermediate state " << state << "." << endl;
      for ( const auto& idat : datamap ) {
        AdcChannel chan = idat.first;
        const AdcChannelData& data = idat.second;
        m_pAdcChannelDataCopyService->copy(data, pintStates->dataMaps[state][chan]);
        snames.insert(state);
      }
    }
  }
  if ( m_DoDeconvolution ) {
    for ( AdcChannelDataMap::value_type& chdata : datamap ) {
      m_pDeconvolutionService->update(chdata.second);
    }
  }
  if ( m_DoDump ) cout << myname << "   After dco: " << datamap[ichan].samples[isig] << endl;
  if ( m_DoPedestalAdjustment ) {
    for ( auto& chdata : datamap ) {
      AdcChannelData& data = chdata.second;
      AdcSignal ped = 0.0;
      m_pPedestalEvaluation->evaluate(data, &ped);
      for ( AdcSignal& sig : data.samples ) sig -= ped;
    }
  }
  if ( m_DoROI ) {
    for ( auto& chdata : datamap ) {
      AdcChannelData& acd = chdata.second;
      m_pRoiBuildingService->build(acd);
      if ( m_DoDump && chdata.first == ichan ) {
        vector<AdcRoi> rois;
        for ( AdcRoi roi : acd.rois ) if ( isig >= roi.first && isig <= roi.second ) rois.push_back(roi);
        cout << myname << "   After roi: " << rois.size() << " of " << acd.rois.size() << " match: ";
        for ( AdcRoi roi : rois ) cout << " (" << roi.first << " , " << roi.second << ")";
        cout << endl;
      }
    }
  }
  if ( m_DoWires ) {
    if ( ! m_DoROI ) {
      if ( m_WiresWithoutROIFlag >= 2 ) cout << myname << "WARNING: Wire building requested without ROI building." << endl;
      if ( m_WiresWithoutROIFlag >= 3 ) abort();
    }
    for ( auto& chdata : datamap ) {
      AdcChannelData& acd = chdata.second;
      m_pWireBuildingService->build(acd, pwires);
      if ( m_DoDump && acd.channel==ichan ) {
        cout << myname << "        Wire: " << pwires->back().Signal().at(isig) << endl;
      }
    }
    for ( string sname : snames ) {
      WiredAdcChannelDataMap& intStates = *pintStates;
      auto inamedacdmap = intStates.dataMaps.find(sname);
      if ( inamedacdmap == intStates.dataMaps.end() ) {
        cout << myname << "WARNING: State " << sname << " does not have data." << endl;
        continue;
      }
      auto inamedwires = intStates.wires.find(sname);
      if ( inamedwires == intStates.wires.end() ) {
        cout << myname << "WARNING: State " << sname << " does not have a wire container." << endl;
        continue;
      }
      AdcChannelDataMap& acdmapState = inamedacdmap->second;
      std::vector<recob::Wire>* pwiresState = inamedwires->second;
      for ( auto& chdata : acdmapState ) {
        AdcChannelData& acd = chdata.second;
        // Create a single ROI.
        acd.signal.clear();
        acd.signal.resize(acd.samples.size(), true);
        acd.roisFromSignal();
        // Build wires.
        m_pWireBuildingService->build(acd, pwiresState);
        if ( m_DoDump && acd.channel==ichan ) {
          cout << myname << "        State " << sname << " wire: " << pwires->back().Signal().at(isig) << endl;
        }
      }
    }
  }
  if ( m_DisplayToolPtrs.size() ) {
    unsigned int idto = 0;
    if ( m_LogLevel >= 2 ) cout << myname << "Running display tools." << endl;
    for ( const AdcChannelToolPtr& padv : m_DisplayToolPtrs ) {
      if ( m_LogLevel >= 2 ) {
        cout << myname << "Running tool " << m_FoundDisplayToolNames[idto] << endl;
      }
      padv->viewMap(datamap);
      ++idto;
    }
  }

  return nbad;
}

//**********************************************************************

std::ostream& StandardRawDigitPrepService::
print(std::ostream& out, std::string prefix) const {
  out << prefix << "StandardRawDigitPrepService:"                      << endl;
  out << prefix << "             LogLevel: " << m_LogLevel             << endl;
  out << prefix << "              SkipBad: " << m_SkipBad              << endl;
  out << prefix << "            SkipNoisy: " << m_SkipNoisy            << endl;
  out << prefix << "  ChannelStatusOnline: " << m_ChannelStatusOnline  << endl;
  out << prefix << "         DoMitigation: " << m_DoMitigation         << endl;
  out << prefix << " DoEarlySignalFinding: " << m_DoEarlySignalFinding << endl;
  out << prefix << "       DoNoiseRemoval: " << m_DoNoiseRemoval       << endl;
  out << prefix << "      DoDeconvolution: " << m_DoDeconvolution      << endl;
  out << prefix << " DoPedestalAdjustment: " << m_DoPedestalAdjustment << endl;
  out << prefix << "                DoROI: " << m_DoROI                << endl;
  out << prefix << "              DoWires: " << m_DoWires              << endl;
  out << prefix << "               DoDump: " << m_DoDump               << endl;
  out << prefix << " DoIntermediateStates: " << m_DoIntermediateStates << endl;
  if ( m_DoDump ) {
    out << prefix << "          DumpChannel: " << m_DumpChannel          << endl;
    out << prefix << "             DumpTick: " << m_DumpTick             << endl;
  }
  if ( m_DisplayTools.size() ) {
    cout << prefix << "        Display tools:";
    for ( string tname : m_DisplayTools ) {
       cout << " " << tname;
    }
    cout << endl;
  } else {
    cout << prefix << "No display tools." << endl;
  }
  return out;
}

//**********************************************************************

DEFINE_ART_SERVICE_INTERFACE_IMPL(StandardRawDigitPrepService, RawDigitPrepService)

//**********************************************************************

