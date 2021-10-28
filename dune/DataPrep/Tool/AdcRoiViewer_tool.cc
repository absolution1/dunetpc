// AdcRoiViewer_tool.cc

#include "AdcRoiViewer.h"
#include "dune/DuneInterface/Tool/AdcChannelStringTool.h"
#include "dune/ArtSupport/DuneToolManager.h"
#include "dune/DuneInterface/Tool/RunDataTool.h"
#include "dune/DuneInterface/Tool/TimeOffsetTool.h"
#include "dune/DuneInterface/Tool/IndexRangeTool.h"
#include "dune/DuneCommon/Utility/gausTF1.h"
#include "dune/DuneCommon/Utility/coldelecResponse.h"
#include "dune/DuneCommon/Utility/quietHistFit.h"
#include "dune/DuneCommon/Utility/shiftHistFit.h"
#include "dune/DuneCommon/Utility/StringManipulator.h"
#include "dune/DuneCommon/Utility/LineColors.h"
#include "dune/DuneCommon/Utility/GausStepFitter.h"
#include "dune/DuneCommon/Utility/GausRmsFitter.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include "TH1F.h"
#include "TH2F.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TF1.h"
#include "TTimeStamp.h"

using std::string;
using std::to_string;
using std::cout;
using std::endl;
using std::ostringstream;
using std::istringstream;
using fhicl::ParameterSet;

using Index = AdcRoiViewer::Index;
using Name = AdcRoiViewer::Name;
using NameVector = std::vector<Name>;
using ParameterSetVector = std::vector<ParameterSet>;
using FloatVector = std::vector<float>;
using IntVector = std::vector<int>;
using ManVector = std::vector<TPadManipulator*>;
using ManVectorMap = std::vector<Name, ManVector>;

namespace {
string timeString(int itim) {
  TTimeStamp ts(itim);
  string stim = ts.AsString("s");
  stim += " UTC";
  return stim;
}
}  // end unnamed namespace

//**********************************************************************
// Subclass methods.
//**********************************************************************

AdcRoiViewer::State::~State() {
  for ( HistInfoMap::value_type ihin : sumHistTemplates ) {
    TH1*& ph = ihin.second.ph;
    delete ph;
    ph = nullptr;
  }
  for ( HistMap::value_type ihst : sumHists ) {
    TH1*& ph = ihst.second;
    delete ph;
    ph = nullptr;
  }
  for ( HistMap::value_type ihst : chanSumHists ) {
    TH1*& ph = ihst.second;
    delete ph;
    ph = nullptr;
  }
}

//**********************************************************************

TH1* AdcRoiViewer::State::getSumHist(Name hname) {
  HistMap::iterator ihst = sumHists.find(hname);
  if ( ihst == sumHists.end() ) return nullptr;
  return ihst->second;
}

//**********************************************************************

Name AdcRoiViewer::State::getSumFitName(Name hnam) const {
  NameMap::const_iterator ifit = sumFitNames.find(hnam);
  if ( ifit == sumFitNames.end() ) return "";
  return ifit->second;
}

//**********************************************************************

Name AdcRoiViewer::State::getSumPlotName(Name hnam) const {
  NameMap::const_iterator iplt = sumPlotNames.find(hnam);
  if ( iplt == sumPlotNames.end() ) return "";
  return iplt->second;
}

//**********************************************************************

float AdcRoiViewer::State::getSumPlotWidth(Name hnam) const {
  FloatMap::const_iterator iplt = sumPlotWidths.find(hnam);
  if ( iplt == sumPlotWidths.end() ) return 0.0;
  return iplt->second;
}

//**********************************************************************

Name AdcRoiViewer::State::getChanSumHistTemplateName(Name hnam) const {
  NameMap::const_iterator ihst = chanSumHistTemplateNames.find(hnam);
  if ( ihst == chanSumHistTemplateNames.end() ) return nullptr;
  return ihst->second;
}

//**********************************************************************

Name AdcRoiViewer::State::getChanSumHistVariableType(Name hnam) const {
  NameMap::const_iterator ihst = chanSumHistVariableTypes.find(hnam);
  if ( ihst == chanSumHistVariableTypes.end() ) return nullptr;
  return ihst->second;
}

//**********************************************************************

Name AdcRoiViewer::State::getChanSumHistErrorType(Name hnam) const {
  NameMap::const_iterator ihst = chanSumHistErrorTypes.find(hnam);
  if ( ihst == chanSumHistErrorTypes.end() ) return nullptr;
  return ihst->second;
}

//**********************************************************************

Name AdcRoiViewer::State::getChanSumPlotName(Name hnam) const {
  NameMap::const_iterator ihst = chanSumPlotNames.find(hnam);
  if ( ihst == chanSumPlotNames.end() ) return nullptr;
  return ihst->second;
}

//**********************************************************************

Index AdcRoiViewer::State::getChannelStatus(Index icha) const {
  const Name myname = "AdcRoiViewer::State::getChannelStatus: ";
  IndexByIndexMap::const_iterator ichs = channelStatuses.find(icha); 
  if ( ichs == channelStatuses.end() ) {
    cout << myname << "WARNING: Status not found for channel " << icha << endl;
    return 0;
  }
  return ichs->second;
}

//**********************************************************************

Index AdcRoiViewer::State::getChannelStatus(Name hnam) const {
  const Name myname = "AdcRoiViewer::State::getChannelStatus: ";
  IndexByNameMap::const_iterator ichc = sumHistChannels.find(hnam); 
  if ( ichc == sumHistChannels.end() ) {
    cout << myname << "WARNING: Channel not found for chansum histogram " << hnam << endl;
    return 0;
  }
  return getChannelStatus(ichc->second);
}

//**********************************************************************
// Class methods.
//**********************************************************************

AdcRoiViewer::AdcRoiViewer(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")),
  m_SigThresh(ps.get<float>("SigThresh")),
  m_TickBorder(ps.get<Index>("TickBorder")),
  m_RoiHistOpt(ps.get<int>("RoiHistOpt")),
  m_FitOpt(ps.get<int>("FitOpt")),
  m_RoiPlotOpt(ps.get<Index>("RoiPlotOpt")),
  m_MaxRoiPlots(ps.get<int>("MaxRoiPlots")),
  m_RoiPlotPadX(ps.get<Index>("RoiPlotPadX")),
  m_RoiPlotPadY(ps.get<Index>("RoiPlotPadY")),
  m_StartTime(ps.get<time_t>("StartTime")),
  m_PulserStepCharge(ps.get<float>("PulserStepCharge")),
  m_PulserDacOffset(ps.get<float>("PulserDacOffset")),
  m_PulserChargeUnit(ps.get<string>("PulserChargeUnit")),
  m_SumNegate(ps.get<bool>("SumNegate")),
  m_SumPlotPadX(ps.get<Index>("SumPlotPadX")),
  m_SumPlotPadY(ps.get<Index>("SumPlotPadY")),
  m_ChannelLineModulus(ps.get<Index>("ChannelLineModulus")),
  m_ChannelLinePattern(ps.get<IndexVector>("ChannelLinePattern")),
  m_RunDataTool(ps.get<string>("RunDataTool")),
  m_TickOffsetTool(ps.get<string>("TickOffsetTool")),
  m_RoiRootFileName(ps.get<string>("RoiRootFileName")),
  m_SumRootFileName(ps.get<string>("SumRootFileName")),
  m_ChanSumRootFileName(ps.get<string>("ChanSumRootFileName")),
  m_ChannelRanges(ps.get<NameVector>("ChannelRanges")),
  m_PlotLabels(ps.get<NameVector>("PlotLabels")),
  m_state(new AdcRoiViewer::State)
{
  const string myname = "AdcRoiViewer::ctor: ";
  if ( m_LogLevel >=2 ) cout << myname << "Begin constructing tool." << endl;
  string stringBuilder = "adcStringBuilder";
  DuneToolManager* ptm = DuneToolManager::instance();
  m_adcStringBuilder = ptm->getShared<AdcChannelStringTool>(stringBuilder);
  if ( m_adcStringBuilder == nullptr ) {
    cout << myname << "WARNING: AdcChannelStringTool not found: " << stringBuilder << endl;
  }
  if ( m_RunDataTool.size() ) {
    m_pRunDataTool = ptm->getShared<RunDataTool>(m_RunDataTool);
    if ( m_pRunDataTool == nullptr ) {
      cout << myname << "WARNING: RunDataTool not found: " << m_RunDataTool << endl;
    } else {
      cout << myname << "Found run data tool." << endl;
    }
  }
  if ( m_TickOffsetTool.size() ) {
    m_pTickOffsetTool = ptm->getShared<TimeOffsetTool>(m_TickOffsetTool);
    if ( m_pTickOffsetTool == nullptr ) {
      cout << myname << "WARNING: Tick offset tool not found: " << m_TickOffsetTool << endl;
    } else {
      cout << myname << "Found tick offset tool." << endl;
    }
  }
  if ( m_ChannelRangeTool.size() ) {
    m_pChannelRangeTool = ptm->getShared<IndexRangeTool>(m_ChannelRangeTool);
    if ( m_pChannelRangeTool == nullptr ) {
      cout << myname << "WARNING: Index range tool not found: " << m_ChannelRangeTool << endl;
    } else {
      cout << myname << "Found channel range tool." << endl;
    }
  }
  // Build the label substitutions.
  for ( Index ilab=0; ilab<m_PlotLabels.size(); ++ilab ) {
    Name stxt = "%LAB" + to_string(ilab) + "%";
    m_plotLabelSubs[stxt] = m_PlotLabels[ilab];
  }
  // Build the summary template histograms.
  // The summary histogram for each channel is created the first time it is encountered in th data.
  ParameterSetVector pshists = ps.get<ParameterSetVector>("SumHists");
  string stimepre = "Time since " + timeString(m_StartTime) + " ";
  for ( const ParameterSet& psh : pshists ) {
    Name hvarx = psh.get<Name>("var");
    Name hvary;
    if ( hvarx == "timingPhase_fitToffPulserMod10" ) {
      hvarx = "fitToffPulserMod10";
      hvary = "timingPhase";
    }
    if ( hvarx == "event_fitToffPulser" ) {
      hvarx = "fitToffPulser";
      hvary = "event";
    }
    Name hnam  = psh.get<Name>("name");
    setPlotLabels(hnam);
    Name httl  = psh.get<Name>("title");
    setPlotLabels(httl);
    if ( getState().sumHistTemplates.find(hnam) != getState().sumHistTemplates.end() ) {
      cout << myname << "ERROR: Duplicate summary template name: " << hnam << endl;
      continue;
    }
    int nbin   = psh.get<int>("nbin");
    float xmin = psh.get<float>("xmin");
    float xmax = psh.get<float>("xmax");
    Name sfit;
    psh.get_if_present("fit", sfit);
    Name plotName;
    psh.get_if_present("plot", plotName);
    float plotWidth = 0.0;
    psh.get_if_present("pwid", plotWidth);
    Name xlab = hvarx;
    if      ( hvarx == "fitHeight"    ) xlab = "Fit height% [SUNIT]%";
    else if ( hvarx == "fitHeightNeg" ) xlab = "-(Fit height)% [SUNIT]%";
    else if ( hvarx == "fitHeightGain" ) {
      xlab = "Fit height gain% [SUNIT]%";
      Name sden = m_PulserChargeUnit;
      if ( sden.size() ) {
        if ( sden.find(" ") != string::npos ) sden = "(" + sden + ")";
        xlab = "Fit height gain [%((SUNIT))%/" + sden + "]";
      }
    }
    else if ( hvarx == "fitWidth"     ) xlab = "Fit width [Tick]";
    else if ( hvarx == "fitPos"       ) xlab = "Fit position [Tick]";
    else if ( hvarx == "fitPosRem"    ) xlab = "Fit position tick remainder [Tick]";
    else if ( hvarx == "fitPosPulser" )
      xlab = "Fit position wrt pulser [Tick]";
    else if ( hvarx == "fitToffPulser" )
      xlab = "Offset fit position wrt pulser [Tick]";
    else if ( hvarx == "fitToffPulserMod10" )
      xlab = "mod_{10}(offset fit position wrt pulser) [Tick]";
    else if ( hvarx == "fitChiSquare" ) xlab = "Fit #chi^{2}";
    else if ( hvarx == "fitChiSquareDof" ) xlab = "Fit #chi^{2}/DOF";
    else if ( hvarx == "fitCSNorm" ) xlab = "Normalized fit #chi^{2}";
    else if ( hvarx == "fitCSNormDof" ) xlab = "Normalized fit #chi^{2}/DOF";
    else if ( hvarx == "sigArea" ) xlab = "Area [%ASUNIT%]";
    else if ( hvarx == "sigAreaNeg" ) xlab = "-Area [%ASUNIT%]";
    else if ( hvarx == "sigWidth" ) xlab = "Width [Tick]";
    else if ( hvarx == "timeSec"   || hvarx == "procTimeSec" ) xlab = stimepre + " [sec]";
    else if ( hvarx == "timeHour"  || hvarx == "procTimeHour") xlab = stimepre + " [hour]";
    else if ( hvarx == "timeDay"   || hvarx == "procTimeDay" ) xlab = stimepre + " [day]";
    else {
      cout << myname << "WARNING: Unknown summary variable: " << hvarx << endl;
    }
    Name ylab;
    if ( hvary.size() ) {
      if ( hvary == "timingPhase" ) ylab = "Timing phase [Tick]";
      if ( hvary == "event" ) ylab = "Event";
    }
    TH1* ph = nullptr;
    if ( hvary == "" ) {
      ph = new TH1F(hnam.c_str(), httl.c_str(), nbin, xmin, xmax);
      ph->GetYaxis()->SetTitle("# ROI");
      ph->Sumw2();  // Needed for likelihood fit
    } else {
      int nbiny  = psh.get<int>("nbiny");
      float ymin = psh.get<float>("ymin");
      float ymax = psh.get<float>("ymax");
      ph = new TH2F(hnam.c_str(), httl.c_str(), nbin, xmin, xmax, nbiny, ymin, ymax);
      ph->GetYaxis()->SetTitle(ylab.c_str());
    }
    ph->SetDirectory(nullptr);
    ph->SetLineWidth(2);
    ph->GetXaxis()->SetTitle(xlab.c_str());
    // Add fit to template so it will be used for each child histogram.
    //if ( sfit.size() ) {
    //  if ( m_LogLevel >= 1 ) cout << myname << "Adding fitter " << sfit
    //                              << " to hist template " << hnam << endl;
    //  TF1* pf = new TF1(sfit.c_str(), sfit.c_str());
    //  ph->GetListOfFunctions()->AddLast(pf);
    //  ph->GetListOfFunctions()->SetOwner(kTRUE);
    //}
    HistInfo& hin = getState().sumHistTemplates[hnam];
    hin.ph = ph;
    hin.varx = hvarx;
    hin.vary = hvary;
    hin.plotName = plotName;
    hin.plotWidth = plotWidth;
    hin.fitName = sfit;
  }
  // Build the channel summary histograms.
  ParameterSetVector pcshists = ps.get<ParameterSetVector>("ChanSumHists");
  for ( const ParameterSet& psh : pcshists ) {
    Name hnam0   = psh.get<Name>("name");    // Name for this histogram
    Name httl0   = psh.get<Name>("title");   // Title for this histogram
    Name vhnam   = psh.get<Name>("valHist"); // Name of the template for the histogram used to fill
    Name valType = psh.get<Name>("valType"); // Type of variable extracted from histogram
    Name etype   = psh.get<Name>("errType"); // Type of variable extracted from histogram
    Name crname0  = psh.get<Name>("cr");     // Name of the channel range for this histogram
    Name plname  = psh.get<Name>("plot");    // Name of the plot file for this histogram
    Name spran   = psh.get<Name>("pran");    // Plot range with format "ymin:ymax"
    Index nbins  = psh.get<Index>("nbins");   // # bins for #ROI vs var plot. 0 means var vs chan plot.;
    if ( hnam0.size() == 0 ) {
      cout << myname << "ERROR: Channel summary histogram name is missing." << endl;
      continue;
    }
    if ( m_pChannelRangeTool == nullptr ) {
      cout << myname << "ERROR: Channel range tool not found." << endl;
      continue;
    }
    HistInfoMap::const_iterator ivh = getState().sumHistTemplates.find(vhnam);
    if ( ivh == getState().sumHistTemplates.end() || ivh->second.ph == nullptr ) {
      cout << myname << "ERROR: Channel summary histogram value histogram not found: " << vhnam << endl;
      continue;
    }
    const NameVector valTypes = {"entries", "count", "mean", "peak", "rms", "sum", "fitMean", "fitSigma"};
    if ( std::find(valTypes.begin(), valTypes.end(), valType) == valTypes.end() ) {
      cout << myname << "ERROR: Channel summary histogram has invalid variable type: " << valType << endl;
      continue;
    }
    const NameVector errTypes = {"none", "zero", "rms", "meanError", "rmsError", "fitSigma"};
    if ( std::find(errTypes.begin(), errTypes.end(), etype) == errTypes.end() ) {
      cout << myname << "ERROR: Channel summary histogram has invalid error type: \"" << etype << '"' << endl;
      continue;
    }
    TH1* phval = ivh->second.ph;
    Name valLabel = phval->GetXaxis()->GetTitle();
    Name yttl = "Unknown";
    if ( valType == "mean" ) {
      yttl = "Mean of " + valLabel;
    } else if ( valType == "peak" ) {
      yttl = "Peak of " + valLabel;
    } else if ( valType == "rms" ) {
      yttl = "Sum of " + valLabel;
    } else if ( valType == "rms" ) {
      yttl = "RMS of " + valLabel;
    } else if ( valType == "fitMean" ) {
      yttl = "Fit mean of " + valLabel;
    } else if ( valType == "fitSigma" ) {
      yttl = "Fit sigma of " + valLabel;
    } else if ( valType == "entries" ) {
      yttl = "# ROI";
    } else if ( valType == "count" ) {
      yttl = "# ROI";
    }
    // Find y-range for plot.
    bool havePlotYMin = false;
    bool havePlotYMax = false;
    float plotYMin = 0;
    float plotYMax = 0;
    Name plotYOpt;
    if ( spran.size() ) {
      Name::size_type ipos = spran.find(":");
      if ( ipos == Name::npos ) {
        cout << myname << "WARNING: Channel summary range specifcation must include \":\"" << endl;
      } else {
        Name::size_type jpos = spran.find(":", ipos+1);
        Name spmin = spran.substr(0, ipos);
        if ( spmin.size() ) {
          istringstream ssin(spmin);
          ssin >> plotYMin;
          havePlotYMin = true;
        } 
        Name spmax = spran.substr(ipos+1, jpos-ipos);
        if ( spmax.size() ) {
          istringstream ssin(spmax);
          ssin >> plotYMax;
          havePlotYMax = true;
        } 
        if ( jpos != Name::npos ) {
          plotYOpt = spran.substr(jpos+1);
        }
      }
    }
    // Loop over channel ranges. Value "list" means all; otherwise just the one given.
    NameVector crns;
    if ( crname0 == "list" ) crns = m_ChannelRanges;
    else crns.push_back(crname0);
    if ( crns.size() == 0 ) crns.push_back("all");
    for ( Name crname : crns ) {
      if ( m_LogLevel >= 2 ) cout << myname << "Creating channel summary histograms for channel range "
                                  << crname << endl;
      IndexRange cr = m_pChannelRangeTool->get(crname);
      if ( ! cr.isValid() ) {
        cout << myname << "ERROR: Channel range " << crname << " not found." << endl;
        continue;
      }
      StringManipulator smhnam(hnam0, false);
      smhnam.replace("%CRNAME%", cr.name);
      smhnam.replace("%CRLABEL%", cr.label());
      smhnam.replace("%CRLABEL1%", cr.label(1));
      smhnam.replace("%CRLABEL2%", cr.label(2));
      Name hnam = smhnam.str();
      if ( getState().chanSumHists.find(hnam) != getState().chanSumHists.end() ) {
        cout << myname << "ERROR: Duplicate channel summary histogram name: " << hnam << endl;
        continue;
      }
      setPlotLabels(hnam);
      StringManipulator smttl(httl0, false);
      smttl.replace("%CRNAME%", cr.name);
      smttl.replace("%CRLABEL%", cr.label());
      smttl.replace("%CRLABEL1%", cr.label(1));
      smttl.replace("%CRLABEL2%", cr.label(2));
      Name httl = smttl.str();
      setPlotLabels(httl);
      TH1* phf = nullptr;
      if ( nbins == 0 ) {
        phf = new TH1F(hnam.c_str(), httl.c_str(), cr.size(), cr.begin, cr.end);
        phf->GetXaxis()->SetTitle("Channel");
        phf->GetYaxis()->SetTitle(yttl.c_str());
      } else {
        phf = new TH1F(hnam.c_str(), httl.c_str(), nbins, plotYMin, plotYMax);
        phf->GetXaxis()->SetTitle(yttl.c_str());
        phf->GetYaxis()->SetTitle("# channels");
        phf->SetLineWidth(2);
      }
      phf->SetDirectory(nullptr);
      phf->SetStats(0);
      if ( cr.size() < 400 ) phf->SetLineWidth(2);
      if ( etype == "none" ) phf->SetMarkerStyle(2);
      else phf->SetMarkerStyle(0);  // Draw error bars instead of markers
      StringManipulator smplt(plname, false);
      smplt.replace("%HNAME%", hnam0);
      smplt.replace("%CRNAME%", cr.name);
      smplt.replace("%CRLABEL%", cr.label());
      smplt.replace("%CRLABEL1%", cr.label(1));
      smplt.replace("%CRLABEL2%", cr.label(2));
      plname = smplt.str();
      setPlotLabels(plname);
      getState().chanSumHists[hnam] = phf;
      getState().chanSumHistTemplateNames[hnam] = vhnam;
      getState().chanSumHistTypes[hnam] = nbins > 0;
      getState().chanSumHistVariableTypes[hnam] = valType;
      getState().chanSumHistErrorTypes[hnam] = etype;
      getState().chanSumPlotNames[hnam] = plname;
      if ( havePlotYMin ) getState().chanSumPlotYMins[hnam] = plotYMin;
      if ( havePlotYMax ) getState().chanSumPlotYMaxs[hnam] = plotYMax;
      if ( havePlotYMin || havePlotYMax ) getState().chanSumPlotYOpts[hnam] = plotYOpt;
      getState().chanSumChaBegin[hnam] = cr.begin;
      getState().chanSumChaEnd[hnam] = cr.end;
      if ( m_LogLevel >= 3 ) {
        cout << myname << "  Histogram name: " << hnam << endl;
        cout << myname << "      Value type: " << valType << endl;
        cout << myname << "      Error type: " << etype << endl;
        cout << myname << "  Histogram name: " << hnam << endl;
        cout << myname << "  Histogram name: " << hnam << endl;
        cout << myname << "       Plot name: " << plname << endl;
        cout << myname << "       Plot ymin: ";
        if ( havePlotYMin ) cout << plotYMin;
        cout << endl;
        cout << myname << "       Plot ymax: ";
        if ( havePlotYMax ) cout << plotYMax;
        cout << endl;
        cout << "       Plot yopt: " << plotYOpt << endl;
      }
    }  // End loop over channel ranges
  }  // End loop over channel summmary histogram configurations
  // Display the configuration.
  if ( m_LogLevel>= 1 ) {
    cout << myname << "             LogLevel: " << m_LogLevel << endl;
    cout << myname << "           RoiHistOpt: " << m_RoiHistOpt << endl;
    cout << myname << "            SigThresh: " << m_SigThresh << endl;
    cout << myname << "           TickBorder: " << m_TickBorder << endl;
    cout << myname << "               FitOpt: " << m_FitOpt << endl;
    cout << myname << "            StartTime: " << m_StartTime
                   << " (" << timeString(m_StartTime) << ")" << endl;
    cout << myname << "     PulserStepCharge: " << m_PulserStepCharge << endl;
    cout << myname << "      PulserDacOffset: " << m_PulserDacOffset << endl;
    cout << myname << "     PulserChargeUnit: " << m_PulserChargeUnit << endl;
    cout << myname << "           RoiPlotOpt: " << m_RoiPlotOpt << endl;
    cout << myname << "          MaxRoiPlots: " << m_MaxRoiPlots << endl;
    cout << myname << "          RoiPlotPadX: " << m_RoiPlotPadX << endl;
    cout << myname << "          RoiPlotPadY: " << m_RoiPlotPadY << endl;
    cout << myname << "            SumNegate: " << (m_SumNegate ? "true" : "false") << endl;
    cout << myname << "          SumPlotPadX: " << m_SumPlotPadX << endl;
    cout << myname << "          SumPlotPadY: " << m_SumPlotPadY << endl;
    cout << myname << "        ChannelRanges: [";
    bool first = true;
    for ( string crn : m_ChannelRanges ) {
      if ( first ) first = false;
      else cout << ", ";
      cout << crn;
    }
    cout << "]" << endl;
    cout << myname << "  ChannelLineModulus: " << m_ChannelLineModulus << endl;
    cout << myname << "  ChannelLinePattern: {";
    first = true;
    for ( Index icha : m_ChannelLinePattern ) {
      if ( ! first ) cout << ", ";
      first = false;
      cout << icha;
    }
    cout << "}" << endl;
    cout << myname << "      RoiRootFileName: " << m_RoiRootFileName << endl;
    cout << myname << "      SumRootFileName: " << m_SumRootFileName << endl;
    cout << myname << "  ChanSumRootFileName: " << m_ChanSumRootFileName << endl;
    if ( getState().sumHistTemplates.size() == 0 ) {
      cout << myname << "  No summary histograms" << endl;
    } else {
      cout << myname << "         SumHists:" << endl;
      for ( const HistInfoMap::value_type& ish : getState().sumHistTemplates ) {
        const HistInfo& hin = ish.second;
        cout << myname << "                   ";
        cout << hin.ph->GetName() << "(" << hin.varx;
        if ( hin.vary.size() ) cout << "," << hin.vary;
        cout << ")";
        if ( hin.fitName.size() ) cout << " fit=" << hin.fitName;
        if ( hin.plotName.size() ) cout << " plot=" << hin.plotName;
        cout << endl;
      }
    }
    if ( getState().chanSumHists.size() == 0 ) {
      cout << myname << "  No channel summary histograms" << endl;
    } else {
      cout << myname << "     ChanSumHists:" << endl;
      for ( HistMap::value_type ihst : getState().chanSumHists ) {
        TH1* ph = ihst.second;
        cout << myname << "                 " << ph->GetName() << endl;
      }
    }
    cout << myname << "        PlotLabels: [";
    first = true;
    for ( Name slab : m_PlotLabels ) {
      if ( first ) first = false;
      else cout << ", ";
      cout << slab;
    }
    cout << "]" << endl;
    cout << myname << "      RunDataTool: \"" << m_RunDataTool << "\" @ "
         << m_pRunDataTool << endl;
    cout << myname << "   TickOffsetTool: \"" << m_TickOffsetTool << "\" @ "
         << m_pTickOffsetTool << endl;
  }
  if ( m_LogLevel >=2 ) cout << myname << "End constructing tool." << endl;
}

//**********************************************************************

AdcRoiViewer::~AdcRoiViewer() {
  const string myname = "AdcRoiViewer::dtor: ";
  close(nullptr);
}

//**********************************************************************

DataMap AdcRoiViewer::view(const AdcChannelData& acd) const {
  DataMap res;
  doView(acd, m_LogLevel, res);
  if ( m_RoiRootFileName.size() ) {
    writeRoiHists(res, m_LogLevel);
  }
  return res;
}

//**********************************************************************

DataMap AdcRoiViewer::viewMap(const AdcChannelDataMap& acds) const {
  const string myname = "AdcRoiViewer::viewMap: ";
  DataMap ret;
  Index ncha = 0;
  Index nroi = 0;
  Index nfail = 0;
  DataMap::IntVector failedChannels;
  bool save = m_RoiRootFileName.size();
  Index ndm = save ? acds.size() : 1;
  int dbg = m_LogLevel > 3 ? m_LogLevel - 2 : 0;
  Index nacd = acds.size();
  DataMapVector dms;  // Cache to hold results from doView
  dms.reserve(ndm);
  Index nroiLimit = save ? 1000 : 1000; // Clear cache after we get this many ROIs.
  Index nroiCached = 0;
  for ( const AdcChannelDataMap::value_type& iacd : acds ) {
    const AdcChannelData& acd = iacd.second;
    if ( m_LogLevel >= 3 ) {
      cout << myname << "Processing channel " << acd.channel()
           << " (" << ncha << "/" << nacd << ")" << endl;
    }
    dms.emplace_back();
    DataMap& dm = dms.back();
    doView(acd, dbg, dm);
    if ( dm.status() ) {
      ++nfail;
      failedChannels.push_back(acd.channel());
    }
    ++ncha;
    nroi += dm.getInt("roiCount");
    nroiCached += dm.getInt("roiCount");
    if ( nroiCached > nroiLimit ) {
      if ( m_LogLevel >= 3 ) cout << myname << "  Clearing result cache." << endl;
      if ( save && dms.size() ) writeRoiHists(dms, dbg);
      dms.clear();
      nroiCached = 0;
    }
  }
  if ( save && dms.size() ) writeRoiHists(dms, dbg);
  ret.setInt("roiChannelCount", ncha);
  ret.setInt("roiFailedChannelCount", nfail);
  ret.setIntVector("roiFailedChannels", failedChannels);
  ret.setInt("roiCount", nroi);
  return ret;
}

//**********************************************************************

DataMap AdcRoiViewer::close(const DataMap* pdmin) {
  const string myname = "AdcRoiViewer::close: ";
  DataMap ret;
  if ( getState().closeCount ) return ret.setStatus(1);
  ++getState().closeCount;
  if ( m_LogLevel >= 1 ) {
    cout << myname << "Closing." << endl;
    cout << myname << "  Event count: " << getState().eventCallCount.size() << endl;
    cout << myname << "   Call count: " << getState().callCount << endl;
    cout << myname << "  Sample unit: " << getState().cachedSampleUnit << endl;
  }
  if ( m_RoiPlotOpt == 2 ) {
    Index ntpm = getState().roiPads.size();
    if (  m_LogLevel >= 2 ) cout << myname << "Printing " << ntpm << " ROI pad"
                                 << (ntpm == 1 ? "" : "s") << endl;
    for ( TpmMap::value_type& itpm : getState().roiPads ) {
      Index icha = itpm.first;
      TpmPtr& pmantop = itpm.second;
      if ( pmantop ) {
        Name plotFileName = getState().roiPadNames[icha];
        if ( m_LogLevel >=3 ) cout << myname << "Writing " << plotFileName << endl;
        pmantop->print(plotFileName);
        pmantop.reset(nullptr);
      }
    }
  }
  if ( getState().sumHists.size() ) {
    fitSumHists();
    writeSumHists();
    writeSumPlots(pdmin);
  }
  if ( getState().chanSumHists.size() ) {
    fillChanSumHists();
    writeChanSumHists();
    writeChanSumPlots();
  }
  return ret;
}

//**********************************************************************

int AdcRoiViewer::doView(const AdcChannelData& acd, int dbg, DataMap& res) const {
  const string myname = "AdcRoiViewer::doView: ";
  unsigned int ievt = acd.event();
  ++getState().callCount;
  if ( getState().eventCallCount.find(ievt) == getState().eventCallCount.end() ) {
    getState().eventCallCount[ievt] = 0;
  } else {
    ++getState().eventCallCount[ievt];
  }
  unsigned int nraw = acd.raw.size();
  unsigned int nsam = acd.samples.size();
  unsigned int ntickChannel = nsam > nraw ? nsam : nraw;
  unsigned int nroiRaw = acd.rois.size();
  bool doHist = m_RoiHistOpt != 0;
  bool histRelativeTick = false;
  int histType = 0;
  if ( doHist ) {
    if        ( m_RoiHistOpt ==  1 ) {
      histType = 1;
    } else if ( m_RoiHistOpt ==  2 ) {
      histType = 2;
    } else if ( m_RoiHistOpt == 11 ) {
      histType = 1;
      histRelativeTick = true;
    } else if ( m_RoiHistOpt == 12 ) {
      histType = 2;
      histRelativeTick = true;
    } else {
      cout << myname << "Invalid value for RoiHistOpt: " << m_RoiHistOpt << endl;
      return res.setStatus(1).status();
    }
  }
  if ( dbg >=2 ) cout << myname << "Processing channel " << acd.channel() << "."
                      << " Input ROI count is " << nroiRaw << endl;
  DataMap::HistVector roiHists;
  DataMap::FloatVector roiSigMins;
  DataMap::FloatVector roiSigMaxs;
  DataMap::FloatVector roiSigAreas;
  DataMap::FloatVector roiFitHeights;
  DataMap::FloatVector roiFitWidths;
  DataMap::FloatVector roiFitPositions;
  DataMap::FloatVector roiFitChiSquares;
  DataMap::FloatVector roiFitChiSquareDofs;
  DataMap::IntVector roiTickMins;
  DataMap::IntVector roiTickMaxs;
  DataMap::IntVector roiTimes;
  DataMap::IntVector roiNUnderflows;
  DataMap::IntVector roiNOverflows;
  DataMap::IntVector tick1;
  DataMap::IntVector ntick;
  DataMap::IntVector roiFitStats;
  Index nroi = 0;
  for ( unsigned int iroiRaw=0; iroiRaw<nroiRaw; ++iroiRaw ) {
    AdcRoi roi = acd.rois[iroiRaw];
    if ( dbg >=3 ) cout << myname << "  ROI " << nroi << "(raw " << iroiRaw << "): ["
                        << roi.first << ", " << roi.second << "]" << endl;
    ostringstream sshnam;
    sshnam << "hroi_run_%0RUN%_evt%0EVENT%_chan%0CHAN%_roi";
    if ( nroi < 100 ) sshnam << "0";
    if ( nroi <  10 ) sshnam << "0";
    sshnam << nroi;
    string hnam = AdcChannelStringTool::AdcChannelStringTool::build(m_adcStringBuilder, acd, sshnam.str());
    ostringstream sshttl;
    sshttl << "Run %RUN% event %EVENT% channel %CHAN% ROI " << nroi;
    sshttl << " ;Tick ;";
    if ( histType == 1 ) sshttl << "Signal% [SUNIT]%";
    if ( histType == 2 ) sshttl << "ADC count";
    string httl = AdcChannelStringTool::build(m_adcStringBuilder, acd, sshttl.str());
    unsigned int isam1 = roi.first;
    unsigned int isam2 = roi.second + 1;
    // Check position if this a ROI to keep.
    if ( m_TickBorder > 0 ) {
      if ( isam1 < m_TickBorder ) continue;
      if ( isam2 + m_TickBorder > nsam ) continue;
    }
    float x1 = histRelativeTick ? 0.0 : isam1;
    float x2 = histRelativeTick ? isam2 - isam1 : isam2;
    TH1* ph = new TH1F(hnam.c_str(), httl.c_str(), isam2-isam1, x1, x2);
    ph->SetDirectory(nullptr);
    //ph->Sumw2();  // Likelihood fit needs weights.
    ph->SetStats(0);
    ph->SetLineWidth(2);
    unsigned int ibin = 0;
    float sigmin = 0.0;
    float sigmax = 0.0;
    float sigarea = 0.0;
    int roiTickMin = 0;
    int roiTickMax = 0;
    Index nunder = 0;
    Index nover = 0;
    for ( unsigned int isam=isam1; isam<isam2; ++isam ) {
      float sig = 0.0;
      if ( histType == 1 && isam<nsam ) sig = acd.samples[isam];
      if ( histType == 2 && isam<nraw ) sig = isam<nraw ? acd.raw[isam] : 0.0;
      if ( ibin == 0 ) {
        sigmin = sig;
        sigmax = sig;
      } else {
        if ( sig < sigmin ) {
          sigmin = sig;
          roiTickMin = ibin;
        }
        if ( sig > sigmax ) {
          sigmax = sig;
          roiTickMax = ibin;
        }
      }
      sigarea += sig;
      ph->SetBinContent(++ibin, sig);
      AdcFlag flag = acd.flags.size() > isam ? acd.flags[isam] : 0;
      if ( flag == AdcUnderflow ) ++nunder;
      if ( flag == AdcOverflow ) ++nover;
    }
    // Check height if this a ROI to keep.
    if ( m_SigThresh < 0.0 && sigmin > m_SigThresh ) continue;
    if ( m_SigThresh > 0.0 && sigmax < m_SigThresh ) continue;
    ++nroi;
    roiHists.push_back(ph);
    roiTickMins.push_back(roiTickMin);
    roiNUnderflows.push_back(nunder);
    roiNOverflows.push_back(nover);
    roiTickMaxs.push_back(roiTickMax);
    roiTimes.push_back(int(acd.time()) - int(m_StartTime));
    roiSigMins.push_back(sigmin);
    roiSigMaxs.push_back(sigmax);
    roiSigAreas.push_back(sigarea);
    tick1.push_back(isam1);
    ntick.push_back(isam2 - isam1);
    if ( m_FitOpt == 1 ) {
      if ( dbg >= 3 ) cout << "  Fitting with coldelecResponse" << endl;
      bool isNeg = fabs(sigmin) > sigmax;
      double h = isNeg ? sigmin : sigmax;
      //double shap = 2.5*ph->GetRMS();  // No! Negative entries break RMS calculation.
      double shap = 0.8*fabs(sigarea)/fabs(h);
      double t0 = x1 + (isNeg ? roiTickMin : roiTickMax) - shap;
      TF1* pf = coldelecResponseTF1(h, shap, t0, "coldelec");
      TF1* pfinit = coldelecResponseTF1(h, shap, t0, "coldelec");
      // The following was very slow when function was written to a file.
      //TF1* pfinit = dynamic_cast<TF1*>(pf->Clone("coldelec0"));
      pfinit->SetLineColor(3);
      pfinit->SetLineStyle(2);
      string fopt = "0";
      fopt = "WWB";
      //fopt = "LWB";  // Use likelihood fit to include empty bins. Do we want this here?
      if ( dbg < 3 ) fopt += "Q";
      // Do shifted fit to avoid numerical issues if t0 is large.
      double xshift = 0.0;
      if ( t0 > 500.0 ) {
        xshift = 500*(long(t0)/500);
      }
      int fstat = shiftHistFit(ph, pf, fopt.c_str(), 2, xshift);
      ph->GetListOfFunctions()->AddLast(pfinit, "0");
      ph->GetListOfFunctions()->Last()->SetBit(TF1::kNotDraw, true);
      ph->GetListOfFunctions()->SetOwner(kTRUE);  // So the histogram owns pfinit
      roiFitHeights.push_back(pf->GetParameter(0));
      roiFitWidths.push_back(pf->GetParameter(1));
      roiFitPositions.push_back(pf->GetParameter(2));
      roiFitStats.push_back(fstat);
      float cs = pf->GetChisquare();
      int ndf = pf->GetNDF();
      float csn = ndf > 0 ? cs/float(ndf) : -1.0;
      roiFitChiSquares.push_back(cs);
      roiFitChiSquareDofs.push_back(csn);
      delete pf;
      //delete pfinit;  This give error: list accessing deleted object
    }
  }
  writeRoiPlots(roiHists, acd);
  res.setInt("roiEvent",   acd.event());
  res.setInt("roiRun",     acd.run());
  res.setInt("roiSubRun",  acd.subRun());
  res.setInt("roiChannel", acd.channel());
  res.setInt("roiCount", nroi);
  res.setInt("roiRawCount", nroiRaw);
  res.setInt("roiNTickChannel", ntickChannel);
  res.setIntVector("roiTick0s", tick1);
  res.setIntVector("roiNTicks", ntick);
  res.setIntVector("roiNUnderflows", roiNUnderflows);
  res.setIntVector("roiNOverflows", roiNOverflows);
  res.setIntVector("roiTickMins", roiTickMins);
  res.setIntVector("roiTickMaxs", roiTickMaxs);
  res.setIntVector("roiTimes", roiTimes);
  res.setFloatVector("roiSigMins", roiSigMins);
  res.setFloatVector("roiSigMaxs", roiSigMaxs);
  res.setFloatVector("roiSigAreas", roiSigAreas);
  res.setHistVector("roiHists", roiHists, true);
  if ( roiFitHeights.size() ) {
    res.setFloatVector("roiFitHeights", roiFitHeights);
    res.setFloatVector("roiFitWidths", roiFitWidths);
    res.setFloatVector("roiFitPositions", roiFitPositions);
    res.setIntVector("roiFitStats", roiFitStats);
    res.setFloatVector("roiFitChiSquares", roiFitChiSquares);
    res.setFloatVector("roiFitChiSquareDofs", roiFitChiSquareDofs);
  }
  fillSumHists(acd, res);
  if ( acd.run() != AdcChannelData::badIndex() ) {
    if ( getState().cachedRunCount == 0 ) {
      getState().cachedRun = acd.run();
      getState().cachedRunCount = 1;
    } else {
      if ( acd.run() != getState().cachedRun ) {
        getState().cachedRun = acd.run();
        ++getState().cachedRunCount;
      }
    }
  }
  if ( getState().cachedSampleUnit.size() == 0 ) {
    getState().cachedSampleUnit = acd.sampleUnit;
  }
  getState().channelStatuses[acd.channel()] = acd.channelStatus();
  return res.status();
}

//**********************************************************************

void AdcRoiViewer::writeRoiHists(const DataMap& dm, int dbg) const {
  DataMapVector dms(1, dm);
  writeRoiHists(dms, dbg);
}

//**********************************************************************

void AdcRoiViewer::writeRoiHists(const DataMapVector& dms, int dbg) const {
  const string myname = "AdcRoiViewer::writeRoiHists: ";
  if ( m_RoiRootFileName.size() == 0 ) return;
  TDirectory* savdir = gDirectory;
  string ofrnameOld = "";
  TFile* pfile = nullptr;
  for ( const DataMap& dm : dms ) {
    AdcChannelData acd;
    acd.setEventInfo(
      dm.getInt("roiRun"),
      dm.getInt("roiEvent"),
      dm.getInt("roiSubRun")
    );
    acd.setChannelInfo(dm.getInt("roiChannel"));
    string ofrname = AdcChannelStringTool::build(m_adcStringBuilder, acd, m_RoiRootFileName);
    if ( ofrname != ofrnameOld ) {
      if ( pfile != nullptr ) pfile->Close();
      delete pfile;
      if ( m_LogLevel >= 2 ) cout << myname << "Writing histograms to " << ofrname << endl;
      pfile = TFile::Open(ofrname.c_str(), "UPDATE");
      ofrnameOld = ofrname;
    }
    const DataMap::HistVector& roiHists = dm.getHistVector("roiHists");
    for ( TH1* ph : roiHists ) {
      TH1* phnew = dynamic_cast<TH1*>(ph->Clone());
      phnew->GetListOfFunctions()->SetOwner(kTRUE);  // So the histogram owns pfinit
      phnew->Write();
      if ( dbg >= 3 ) cout << myname << "  Wrote " << phnew->GetName() << endl;
    }
  }
  if ( pfile != nullptr ) pfile->Close();
  delete pfile;
  savdir->cd();
}

//**********************************************************************

void AdcRoiViewer::writeRoiPlots(const HistVector& hsts, const AdcChannelData& acd) const {
  const string myname = "AdcRoiViewer::writeRoiPlots: ";
  if ( m_MaxRoiPlots >=0 && getState().nRoiPlot >= Index(m_MaxRoiPlots) ) return;
  Index npadx = m_RoiPlotPadX;
  Index npady = m_RoiPlotPadY;
  Index npad = npadx*npady;
  if ( npad == 0 ) return;
  Index wpadx = 1400;
  Index wpady = 1000;
  TpmPtr pmantopLocal;
  TpmPtr& pmantop = m_RoiPlotOpt == 2 ? getState().roiPads[acd.channel()] : pmantopLocal;
  Name plotFileName;
  Index ipad = 0;
  if ( m_RoiPlotOpt == 2 ) {
    // Fetch the print name.
    plotFileName = getState().roiPadNames[acd.channel()];
    // Find the first empty sub pad.
    if ( pmantop && npad > 1 ) {
      for ( ipad=0; ipad<npad; ++ipad ) {
        if ( pmantop->man(ipad)->hist() == nullptr ) break;
      }
      if ( ipad >= npad ) {
        cout << myname << "FATAL: ROI pad is full." << endl;
        abort();
      }
    }
  }
  Index ihst = 0;
  for ( TH1* ph : hsts ) {
    if ( ph == nullptr ) continue;
    Name hnam = ph->GetName();
    if (  ! pmantop ) {
      plotFileName = hnam.substr(1) + ".png";  // Strip leading h from histogram name.
      if ( m_RoiPlotOpt == 2 ) {
        ostringstream sscha;
        sscha << acd.channel();
        string scha = sscha.str();
        while ( scha.size() < 6 ) scha = "0" + scha;
        ostringstream sspag;
        sspag << getState().roiPadCounts[acd.channel()];
        string spag = sspag.str();
        while ( spag.size() < 3 ) spag = "0" + spag;
        plotFileName = "roi_chan" + scha + "_" + spag + ".png";
        getState().roiPadNames[acd.channel()] = plotFileName;
      }
      ipad = 0;
      pmantop.reset(new TPadManipulator);
      pmantop->setCanvasSize(wpadx, wpady);
      if ( npad > 1 ) pmantop->split(npadx, npady);
      if (  m_LogLevel >= 3 ) cout << myname << "  Creating plots for " << plotFileName << endl;
      if (  m_LogLevel >= 4 ) cout << myname << "    Plotting " << ph->GetName() << endl;
    }
    TPadManipulator* pman = pmantop->man(ipad);
    pman->add(ph, "hist", false);
    pman->addHistFun(0);
    TF1* pfit = ph->GetFunction("coldelec");
    if ( true ) {
      NameVector labs;
      double area = ph->Integral();
      ostringstream ssout;
      ssout.precision(3);
      ssout.setf(std::ios_base::fixed);
      ssout << "Area: " << area;
      labs.push_back(ssout.str());
      if ( pfit != nullptr ) {
        double height = pfit->GetParameter("Height");
        double shaping = pfit->GetParameter("Shaping");
        double t0 = pfit->GetParameter("T0");
        ssout.str("");
        ssout << "Height: " << height;
        labs.push_back(ssout.str());
        ssout.str("");
        ssout << "Shaping: " << shaping << " tick";
        labs.push_back(ssout.str());
        ssout.str("");
        ssout.precision(2);
        ssout << "Position: " << t0 << " tick";
        labs.push_back(ssout.str());
        ssout.str("");
        ssout.precision(1);
        ssout << "#chi^{2}: " << pfit->GetChisquare();
        labs.push_back(ssout.str());
        Index chanStat = acd.channelStatus();
        if ( chanStat == AdcChannelStatusBad ) labs.push_back("Bad channel");
        if ( chanStat == AdcChannelStatusNoisy ) labs.push_back("Noisy channel");
      }
      double xlab = 0.70;
      double ylab = 0.80;
      double dylab = 0.04;
      for ( Name lab : labs ) {
        TLatex* pptl = nullptr;
        pptl = new TLatex(xlab, ylab, lab.c_str());
        pptl->SetNDC();
        pptl->SetTextFont(42);
        pptl->SetTextSize(dylab);
        pman->add(pptl);
        ylab -= 1.2*dylab;
      }
    }
    pman->addAxis();
    pman->addHorizontalLine(0.0);
    pman->showUnderflow();
    pman->showOverflow();
    ++ipad;
    if ( ipad >= npad || (++ihst >= hsts.size() && m_RoiPlotOpt != 2) ) {
      if (  m_LogLevel >= 3 ) cout << myname << "  Writing " << plotFileName << endl;
      pmantop->print(plotFileName);
      pmantop.reset(nullptr);
      ++getState().roiPadCounts[acd.channel()];
      ipad = 0;
      ++getState().nRoiPlot;
      if ( m_MaxRoiPlots >=0 && getState().nRoiPlot >= Index(m_MaxRoiPlots) ) return;
    }
  }
}

//**********************************************************************

void AdcRoiViewer::fillSumHists(const AdcChannelData& acd, const DataMap& dm) const {
  const string myname = "AdcRoiViewer::fillSumHists: ";
  // Fetch the run data.
  RunData rdat;
  if ( m_pRunDataTool != nullptr ) {
    rdat = m_pRunDataTool->runData(acd.run(), acd.subRun());
    RunData& rdatOld = getState().runData;
    if ( rdat.isValid() && ! rdatOld.isValid() ) {
      if ( m_LogLevel >= 2 ) cout << myname << "Setting run data." << endl;
      rdatOld = rdat;
    } else if ( rdat.isValid() && rdatOld.isValid() ) {
      if ( rdat.run() != rdatOld.run() ) {
        cout << myname << "Ignoring unexpected change in run number: " << rdatOld.run()
             << " --> " << rdat.run();
      }
    } else if ( ! rdat.isValid() ) {
      if ( m_LogLevel >= 3 ) cout << myname << "Run data not found." << endl;
    }
  }
  float pulserQin = 0.0;
  bool havePulserAmplitude = rdat.havePulserAmplitude() && rdat.havePulserSource();
  bool havePulserPeriod = rdat.havePulserPeriod();
  bool haveQin = false;
  if ( havePulserAmplitude ) {
    int qfac = rdat.pulserAmplitude();
    //if ( rdat.pulserSource() == 2 && qfac > 0 ) --qfac;     // Should we do this??
    pulserQin = (qfac - m_PulserDacOffset)*m_PulserStepCharge;
    haveQin = pulserQin != 0.0;
    if ( ! haveQin ) {
      cout << myname << "WARNING: Pulser charge evaluates to zero." << endl;
    }
  }
  Index pulserPeriod = 0;
  if ( havePulserPeriod ) {
    pulserPeriod = rdat.pulserPeriod();
    if ( pulserPeriod == 0 ) {
      havePulserPeriod = false;
      cout << myname << "WARNING: Pulser period is zero." << endl;
    }
  }
  // Fetch the tick offset.
  TimeOffsetTool::Data tdat;
  bool haveTickOffset = false;
  long tickOffset = 0;
  bool haveTickOffsetPulserMod = false;
  Index tickOffsetPulserMod = 0;   // Tick offset modulus the pulser period [0, pulserPeriod).
  double timingPhase = 0.0;   // Phase  of the timing clock (0,1]
  if ( m_pTickOffsetTool != nullptr ) {
    tdat.run = acd.run();
    tdat.subrun = acd.subRun();
    tdat.event = acd.event();
    tdat.channel = acd.channel();
    tdat.fembID = acd.fembID();
    tdat.triggerClock = acd.triggerClock();
    TimeOffsetTool::Offset off = m_pTickOffsetTool->offset(tdat);
    if ( off.isValid() ) {
      haveTickOffset = true;
      tickOffset = off.value;
      timingPhase = off.rem;
    } else {
      cout << myname << "Unable to retrieve tick offset for run " << tdat.run << "-" << tdat.subrun
           << " event " << tdat.event << " channel " << tdat.channel << endl;
    }
    if ( haveTickOffset && havePulserPeriod ) {
      long toff = tickOffset;
      long period = pulserPeriod;
      toff = toff % period;
      if ( toff < 0 ) toff += period;
      tickOffsetPulserMod = toff;
      haveTickOffsetPulserMod = true;
    }
  }
  // Loop over summary histogram templates.
  Index nhst = 0;
  Index nhstGood = 0;
  for ( const HistInfoMap::value_type ish : getState().sumHistTemplates ) {
    const HistInfo& hin0 = ish.second;
    ++nhst;
    Name varx = hin0.varx;
    if ( m_SumNegate ) {
      if ( varx == "fitHeight" ) varx = "fitHeightNeg";
      if ( varx == "sigArea" ) varx = "sigAreaNeg";
    }
    Name vary = hin0.vary;
    TH1* ph0 = hin0.ph;
    Name fitName = hin0.fitName;
    Name plotNameTemplate = hin0.plotName;
    FloatVector vals;
    IntVector ivals;
    if      ( varx == "sigArea" )            vals = dm.getFloatVector("roiSigAreas");
    else if ( varx == "sigAreaNeg" )         vals = dm.getFloatVector("roiSigAreas");
    else if ( varx == "sigWidth" )          ivals = dm.getIntVector("roiNTicks");
    else if ( varx == "sigTick0" )          ivals = dm.getIntVector("roiTick0s");
    else if ( varx == "sigTick0Pulser" )    ivals = dm.getIntVector("roiTick0s");
    else if ( varx == "sigTickMax" )        ivals = dm.getIntVector("roiTickMaxs");
    else if ( varx == "sigTickMaxPulser" )  ivals = dm.getIntVector("roiTickMaxs");
    else if ( varx == "fitHeight"    )       vals = dm.getFloatVector("roiFitHeights");
    else if ( varx == "fitHeightNeg" )       vals = dm.getFloatVector("roiFitHeights");
    else if ( varx == "fitHeightGain" )      vals = dm.getFloatVector("roiFitHeights");
    else if ( varx == "fitWidth"     )       vals = dm.getFloatVector("roiFitWidths");
    else if ( varx == "fitPos"  )            vals = dm.getFloatVector("roiFitPositions");
    else if ( varx == "fitPosRem"   )        vals = dm.getFloatVector("roiFitPositions");
    else if ( varx == "fitPosPulser" )       vals = dm.getFloatVector("roiFitPositions");
    else if ( varx == "fitToffPulser" )      vals = dm.getFloatVector("roiFitPositions");
    else if ( varx == "fitToffPulserMod10" ) vals = dm.getFloatVector("roiFitPositions");
    else if ( varx == "fitStat" )           ivals = dm.getIntVector("roiFitStats");
    else if ( varx == "fitChiSquare" )       vals = dm.getFloatVector("roiFitChiSquares");
    else if ( varx == "fitChiSquareDof" )    vals = dm.getFloatVector("roiFitChiSquareDofs");
    else if ( varx == "fitCSNorm" )          vals = dm.getFloatVector("roiFitChiSquares");
    else if ( varx == "fitCSNormDof" )       vals = dm.getFloatVector("roiFitChiSquareDofs");
    else if ( varx == "timeSec" )           ivals = dm.getIntVector("roiTimes");
    else if ( varx == "timeHour" )          ivals = dm.getIntVector("roiTimes");
    else if ( varx == "timeDay" )           ivals = dm.getIntVector("roiTimes");
    else if ( varx == "procEvent" )         ivals.push_back(acd.event());
    else if ( varx == "procTimeSec" )       ivals.push_back(int(acd.time()) - int(m_StartTime));
    else if ( varx == "procTimeHour" )      ivals.push_back(int(acd.time()) - int(m_StartTime));
    else if ( varx == "procTimeDay" )       ivals.push_back(int(acd.time()) - int(m_StartTime));
    else {
      if ( m_LogLevel >= 2 ) {
        cout << myname << "ERROR: Invalid variable name: " << varx << endl;
        continue;
      }
    }
    Name hnam0 = ph0->GetName();
    Name hnam = AdcChannelStringTool::build(m_adcStringBuilder, acd, hnam0);
    Name xlab = AdcChannelStringTool::build(m_adcStringBuilder, acd, ph0->GetXaxis()->GetTitle());
    Name ylab = ph0->GetYaxis()->GetTitle();
    TH1* ph = getState().getSumHist(hnam);
    if ( ivals.size() && !vals.size() ) for ( int ival : ivals ) vals.push_back(ival);
    if ( varx == "fitPosRem" ) for ( float& val : vals ) val = std::remainder(val,1);
    if ( varx == "fitPosPulser" ) {
      if ( ! havePulserPeriod ) {
        cout << myname << "WARNING: Cannot evaluate " << varx << " without pulser period" << endl;
        continue;
      }
      for ( float& val : vals ) val = fmod(val, pulserPeriod);
    }
    if ( varx == "fitToffPulser" || varx == "fitToffPulserMod10" ||
         varx == "sigTick0Pulser" || varx == "sigTickMaxPulser" ) {
      if ( ! haveTickOffsetPulserMod ) {
        cout << myname << "WARNING: Cannot evaluate " << varx << " without timing offset and pulser period" << endl;
        continue;
      }
      for ( float& val : vals ) val = fmod(val + pulserPeriod + tickOffsetPulserMod, pulserPeriod);
      if ( varx == "fitToffPulserMod10" ) {
        for ( float& val : vals ) val = fmod(val, 10.0);
      }
    }
    float varfac = 1.0;
    if ( varx == "fitHeightNeg" ) varfac = -1.0;
    if ( varx == "sigAreaNeg" ) varfac = -1.0;
    if ( varx == "fitCSNorm" ||  varx == "fitCSNormDof" ) {
      float pedrms = acd.pedestalRms;
      if ( pedrms > 0.0 ) varfac = 1.0/(pedrms*pedrms);
      else varfac = 0.0;
    }
    if ( varx == "fitHeightGain" ) {
      if ( ! haveQin ) {
        cout << myname << "WARNING: Cannot evaluate " << varx << " without Qin" << endl;
        continue;
      }
      varfac = 1.0/pulserQin;
    }
    if ( varx == "timeHour" || varx == "procTimeHour" ) varfac = 1/3600.0;
    if ( varx == "timeDay" || varx == "procTimeDay" ) varfac = 1/(24*3600.0);
    if ( varfac != 1.0 ) for ( float& val : vals ) val *= varfac;
    // Create histogram if it does not yet exist or is empty and we have data here.
    //if ( ph == nullptr && vals.size() ) {
    //if ( ph == nullptr ) {
    if ( (ph == nullptr) || (ph->GetEntries() == 0 && vals.size()) ) {
      bool replacingHistogram = ph != nullptr;
      // Fetch the vector that lists the summary histograms to be plotted for the current template.
      // Find or make a place to record this histogram.
      HistVector empty;
      bool havePlot = plotNameTemplate.size();
      HistVector& plotHists = havePlot ? getState().sumPlotHists[plotNameTemplate] : empty;
      HistVector::iterator iplotHist = plotHists.end();
      if ( replacingHistogram ) {
        if ( havePlot ) {
          iplotHist = find(plotHists.begin(), plotHists.end(), ph);
          if ( iplotHist == plotHists.end() ) {
            cout << myname << "ERROR: Unable to find histogram in plot name list." << endl;
            plotHists.clear();
            iplotHist = plotHists.end();
          }
          *iplotHist = nullptr;
        }
        delete ph;
        ph = nullptr;
        if ( m_LogLevel >= 2 ) cout << myname << "Replacing histogram " << hnam << endl;
      } else {
        if ( havePlot ) {
          plotHists.push_back(nullptr);
          iplotHist = plotHists.end();
          --iplotHist;
        }
        if ( m_LogLevel >= 2 ) cout << myname << "Creating histogram " << hnam << endl;
      }
      Name httl0 = ph0->GetTitle();
      Name httl = AdcChannelStringTool::build(m_adcStringBuilder, acd, httl0);
      int nbin = ph0->GetNbinsX();
      float xmin = ph0->GetXaxis()->GetXmin();
      float xmax = ph0->GetXaxis()->GetXmax();
      // If xmin > xmax and xmin > 0, then we center histogram on median and use width = xmin.
      // If also xmax >0, then we round the first bin edge to that value.
      if ( vals.size() && xmin > xmax && xmin > 0.0 ) {
        FloatVector tmpvals = vals;
        std::sort(tmpvals.begin(), tmpvals.end());
        Index nval = tmpvals.size();
        if ( m_LogLevel >= 3 ) cout << myname << "  Centering histogram on median of "
                                    << nval << " value" << (nval==1 ? "" : "s") << endl;
        float xmed = 0.5*(tmpvals[(nval-1)/2] + tmpvals[nval/2]);
        float width = xmin;
        xmin = xmed - 0.5*width;
        bool roundXmin = xmax > 0.0;
        if ( roundXmin ) {
          //float rexp = log10(width/50.0);
          //rexp = rexp > 0 ? int(rexp) : int(rexp-1);
          //float rfac = pow(10, rexp);
          float rfac = xmax;
          xmin = rfac*int(xmin/rfac + (xmin > 0.0 ? 0.5 : -0.5));
        }
        xmax = xmin + width;
      // Make sure we have xmax > xmin so Root won't complain when we plot.
      } else if ( xmin >= xmax ) {
        if ( xmin > 0.0 ) {
          xmax = xmin;
          xmin = 0;
        } else if ( xmin < 0.0 ) {
          xmax = 0.0;
        } else {
          xmin = 0.0;
          xmax = 1.0;
        }
      }
      if ( m_LogLevel >= 3 ) {
        cout << myname << "   Name: " << hnam << endl;
        cout << myname << "  Title: " << httl << endl;
        cout << myname << "   nbin: " << nbin << endl;
        cout << myname << "   xmin: " << xmin << endl;
        cout << myname << "   xmax: " << xmax << endl;
        cout << myname << "    fit: " << fitName << endl;
      }
      bool isTH2 = dynamic_cast<TH2*>(ph0);
      if ( ! isTH2 ) {
        ph = new TH1F(hnam.c_str(), httl.c_str(), nbin, xmin, xmax);
      } else {
        int nbiny = ph0->GetNbinsY();
        float ymin = ph0->GetYaxis()->GetXmin();
        float ymax = ph0->GetYaxis()->GetXmax();
        ph = new TH2F(hnam.c_str(), httl.c_str(), nbin, xmin, xmax, nbiny, ymin, ymax);
      }
      ph->SetDirectory(nullptr);
      ph->SetStats(0);
      ph->Sumw2();  // Needed for likelihood fit.
      ph->SetLineWidth(2);
      ph->GetXaxis()->SetTitle(xlab.c_str());
      ph->GetYaxis()->SetTitle(ylab.c_str());
      if ( ph0->GetListOfFunctions()->GetEntries() ) {
        TF1* pf = dynamic_cast<TF1*>(ph0->GetListOfFunctions()->At(0)->Clone());
        ph->GetListOfFunctions()->AddLast(pf);
        ph->GetListOfFunctions()->SetOwner(kTRUE);
      }
      getState().sumHistChannels[hnam] = acd.channel();
      getState().sumHists[hnam] = ph;
      getState().sumFitNames[hnam] = fitName;
      if ( plotNameTemplate.size() ) {
        Name plotNameHist = AdcChannelStringTool::build(m_adcStringBuilder, acd, plotNameTemplate);
        if ( havePlot ) *iplotHist = ph;
        if ( ! replacingHistogram ) {
          getState().sumPlotNames[hnam] = plotNameHist;
          getState().sumPlotWidths[hnam] = hin0.plotWidth;
        }
      }
    }
    if ( m_LogLevel >= 3 ) cout << myname << "Filling summary histogram " << hnam
                                << ". Fill count is " << vals.size() << endl;
    FloatVector csds = dm.getFloatVector("roiFitChiSquareDofs");
    IntVector fstats = dm.getIntVector("roiFitStats");
    bool checkFit = varx.substr(0,3) == "fit";
    if ( csds.size() == 0 || fstats.size() == 0 ) checkFit = false;
    double chiSquareDofMax = 0.0;   // was 1000; should be config param?
    if ( checkFit ) {
      if ( csds.size() != vals.size() ) {
      cout << "ERROR: Variable and chi-square/DF vectors have different sizes." << endl;
        checkFit = false;
      }
      if ( fstats.size() != vals.size() ) {
        cout << "ERROR: Variable and fit status vectors have different sizes: "
             << vals.size() << " != " << fstats.size()  << endl;
        checkFit = false;
      }
    }
    Index nval = 0;
    Index nvalSkip = 0;
    bool logerr = m_LogLevel >= 3 && nhstGood == 0;
    for ( Index ival=0; ival<vals.size(); ++ival ) {
      ++nval;
      if ( checkFit ) {
        int fstat = fstats[ival];
        float csd = csds[ival];
        if ( fstat ) {
          if ( logerr) cout << myname << "WARNING: Skipping entry with fit status " << fstat
                            << " (chi-square/DOF = " << csd << ")" << endl;
          ++nvalSkip;
          continue;
        } else if ( chiSquareDofMax > 0.0 && csd > chiSquareDofMax ) {
          if ( logerr) cout << myname << "WARNING: Skipping entry with chi-square/DOF = " << csd << endl;
          ++nvalSkip;
          continue;
        }
      }
      float val = vals[ival];
      double valy = 0.0;
      if ( vary == "timingPhase" ) valy = timingPhase;
      if ( vary == "event" ) valy = acd.event();
      if ( vary == "" ) {
        ph->Fill(val);
      } else {
        ph->Fill(val, valy);
      }
    }
    // Show skips for the first histogram only.
    if ( nvalSkip ) {
      cout << myname << "WARNING: Skipped " << nvalSkip << " of " << nval
           << " entries due to bad fit for histogram " << hnam << "." << endl;
    }
    ++nhstGood;
  }
  if ( nhstGood != nhst ) {
    cout << myname << "WARNING: Only filled " << nhstGood << " of " << nhst << " histograms." << endl;
  }
}

//**********************************************************************

void AdcRoiViewer::fitSumHists() const {
  const string myname = "AdcRoiViewer::fitSumHists: ";
  if ( getState().sumHists.size() == 0 ) {
    cout << myname << "No summary histograms found." << endl;
    return;
  }
  if ( m_LogLevel >= 1 ) cout << myname << "Fitting summary histograms. Count is "
                              << getState().sumHists.size() << "." << endl;
  for ( HistMap::value_type ihst : getState().sumHists ) {
    TH1* ph = ihst.second;
    string hnam = ph->GetName();
    string fitName = getState().getSumFitName(hnam);
    bool fitDone = false;
    if ( m_LogLevel >= 3 ) cout << myname << "Fitting hist " << ph->GetName() << " with " << fitName << endl;
    if ( fitName.size() ) {
      TF1* pf = nullptr;
      int binMax = ph->GetMaximumBin();
      double mean0 = ph->GetBinLowEdge(binMax);
      double sigma0 = ph->GetRMS();
      double height0 = ph->GetMaximum();
      // Use gaus step fit.
      if ( fitName.substr(0,5) == "sgaus" ) {
        if ( m_LogLevel >= 4 ) cout << myname << "  Doing gaus step fit." << endl;
        if ( fitName.size() > 5 ) {
          istringstream ssin(fitName.substr(5));
          ssin >> sigma0;
        }
        GausStepFitter gsf(mean0, sigma0, height0, fitName, "WWS");
        fitDone = gsf.fit(ph) == 0;
      // Use gaus from RMS.
      } else if ( fitName.substr(0,5) == "rgaus" ) {
        if ( m_LogLevel >= 4 ) cout << myname << "  Doing fixed rms fit." << endl;
        double sigma0 = 0.0;
        double nsigma = 4.0;
        Name spar1 = fitName.substr(5);
        Name spar2;
        if ( spar1.size() ) {
          if ( spar1[0] == '_' ) spar1 = spar1.substr(1);
          string::size_type ipos = spar1.find("_");
          if ( ipos != string::npos ) {
            spar2 = spar1.substr(ipos+1);
            spar1 = spar1.substr(0, ipos);
            istringstream ssin2(spar2);
            ssin2 >> nsigma;
          }
          istringstream ssin1(spar1);
          ssin1 >> sigma0;
        }
        GausRmsFitter grf(sigma0, nsigma, fitName);
        if ( m_LogLevel >=4 ) grf.setLogLevel(m_LogLevel - 3);
        if ( grf.fit(ph, mean0) == 0 ) {
          fitDone = true;
        } else {
          pf = new TF1("mygaus", "gaus");
        }
      } else {
        pf = new TF1(fitName.c_str(), fitName.c_str());
      }
      if ( m_LogLevel >= 4 && pf != nullptr ) {
        cout << myname << "  Created function " << pf->GetName() << " at " << std::hex << pf << endl;
      }
      if ( ! fitDone ) {
        if ( m_LogLevel >= 4 ) cout << myname << "  Doing unconstrained fit" << endl;
        string fopt = "LWWS";
        if ( m_LogLevel < 4 ) fopt += "Q"; 
        int fstat = quietHistFit(ph, pf, fopt.c_str());
        if ( fstat != 0 ) {
          cout << myname << "  WARNING: Fit " << pf->GetName() << " of " << ph->GetName() << " returned " << fstat << endl;
          ph->GetListOfFunctions()->Clear();   // Otherwise we may get a crash when we try to view saved copy of histo
        } else {
          if ( m_LogLevel >=4 ) cout << myname << "  Fit succeeded." << endl;
        }
        ph->GetListOfFunctions()->SetOwner(kTRUE);  // So the histogram owns pf
      }
      delete pf;
    }
  }
}

//**********************************************************************

void AdcRoiViewer::writeSumHists() const {
  const string myname = "AdcRoiViewer::writeSumHists: ";
  bool saveHist = m_SumRootFileName.size();
  if ( ! saveHist ) return;
  if ( getState().sumHists.size() == 0 ) {
    cout << myname << "No summary histograms found." << endl;
    return;
  }
  TDirectory* savdir = gDirectory;
  Name ofrname = m_SumRootFileName;
  TFile* pfile = TFile::Open(ofrname.c_str(), "UPDATE");
  saveHist = pfile != nullptr && pfile->IsOpen();
  if ( ! saveHist ) {
    cout << myname << "ERROR: Unable to open output file " << ofrname << endl;
    return;
  }
  if ( m_LogLevel >= 1 ) cout << myname << "Writing summary histograms. Count is "
                              << getState().sumHists.size() << "." << endl;
  for ( HistMap::value_type ihst : getState().sumHists ) {
    TH1* ph = ihst.second;
    TH1* phnew = dynamic_cast<TH1*>(ph->Clone());
    if ( saveHist ) phnew->Write();
    if ( m_LogLevel >= 2 ) cout << myname << "  Wrote " << phnew->GetName() << endl;
  }
  pfile->Close();
  if ( m_LogLevel >= 1 ) cout << myname << "Closed summary histogram file " << ofrname << endl;
  savdir->cd();
  delete pfile;
}

//**********************************************************************

void AdcRoiViewer::writeSumPlots(const DataMap* pdmin) const {
  const string myname = "AdcRoiViewer::writeSumPlots: ";
  int dbgin = pdmin == nullptr ? 0 : pdmin->getInt("dbg", 0);
  Index npad = 0;
  Index npadx = m_SumPlotPadX;
  Index npady = m_SumPlotPadY;
  Index wpadx = 1400;
  Index wpady = 1000;
  npad = npadx*npady;
  Index nvec = getState().sumPlotHists.size();
  if ( npad == 0 ) return;
  if (  m_LogLevel >= 1 ) cout << myname << "Plotting " << nvec << " set"
                               << (nvec == 1 ? "" : "s") << " of summary histograms " << endl;
  for ( const HistVectorMap::value_type ihv : getState().sumPlotHists ) {
    Name plotNameTemplate = ihv.first;
    const HistVector& hsts = ihv.second;
    TPadManipulator* pmantop = nullptr;
    Index ipad = 0;
    Name plotFileName;
    for ( Index ihst=0; ihst<hsts.size(); ++ihst ) {
      TH1* ph = hsts[ihst];
      if ( ph == nullptr ) {
        cout << myname << "WARNING: Histogram " << ihst << " not found for template " << plotNameTemplate << endl;
      } else {
        Name hnam = ph->GetName();
        bool haveExpectedValue = pdmin != nullptr && pdmin->haveFloat(hnam);
        float expValue = haveExpectedValue ? pdmin->getFloat(hnam) : 0.0;
        if ( dbgin ) {
          cout << myname << "  Processing histogram " << hnam << endl;
          if ( haveExpectedValue ) cout << myname << "    Expected value: " << expValue << endl;
          else cout << myname << "    No expected value." << endl;
        }
        if ( pmantop == nullptr ) {
          plotFileName = getState().getSumPlotName(hnam);
          if ( plotFileName.size() == 0 ) {
            cout << myname << "ERROR: Plot file name is not assigned for " << hnam << endl;
            break;
          }
          ipad = 0;
          pmantop = new TPadManipulator;
          if ( npadx && npady ) pmantop->setCanvasSize(wpadx, wpady);
          if ( npad > 1 ) pmantop->split(npadx, npady);
          if (  m_LogLevel >= 2 ) cout << myname << "  Creating plots for " << plotFileName << endl;
        }
        if (  m_LogLevel >= 3 ) cout << myname << "    Plotting " << ph->GetName() << endl;
        TPadManipulator* pman = pmantop->man(ipad);
        pman->add(ph, "hist", false);
        if ( ph->GetListOfFunctions()->GetEntries() ) {
          //dynamic_cast<TF1*>(pman->hist()->GetListOfFunctions()->At(0))->SetNpx(2000);
          pman->addHistFun(0);
        }
        pman->addAxis();
        pman->showUnderflow();
        pman->showOverflow();
        float plotWidth = getState().getSumPlotWidth(hnam);
        if ( plotWidth > 0.0 ) {
          int binMax = ph->GetMaximumBin();
          if ( binMax && binMax <= ph->GetNbinsX() ) {
            float xCen = ph->GetBinLowEdge(binMax);
            float xmin = xCen - 0.5*plotWidth;
            float xmax = xCen + 0.5*plotWidth;
            pman->setRangeX(xmin, xmax);
          }
        }
        NameVector labs;
        bool showMean = true;
        if ( showMean ) {
          double mean = ph->GetMean();
          double sigm = ph->GetRMS();
          float meanDen = haveExpectedValue ? expValue : mean;
          ostringstream ssout;
          ssout.precision(3);
          ssout.setf(std::ios_base::fixed);
          ssout.str("");
          ssout << "# ROI: " << int(ph->GetEntries() + 0.1);
          labs.push_back(ssout.str());
          ssout.str("");
          ssout << "Hist mean: " << mean;
          labs.push_back(ssout.str());
          ssout.str("");
          ssout << "Hist RMS: " << sigm;
          ssout.precision(1);
          if ( meanDen != 0.0 ) ssout << " [" << 100*sigm/meanDen << "%]";
          ssout.precision(3);
          labs.push_back(ssout.str());
          if ( haveExpectedValue ) {
            ssout.str("");
            float bias = mean - expValue;
            ssout << "Hist bias: " << bias;
            ssout.precision(1);
            if ( expValue != 0.0 ) ssout << " [" << 100*bias/expValue << "%]";
            ssout.precision(3);
            labs.push_back(ssout.str());
          }
        }
        TF1* pffit = dynamic_cast<TF1*>(ph->GetListOfFunctions()->Last());
        if ( pffit != nullptr ) {
          string fnam = pffit->GetName();
          double mean = pffit->GetParameter("Mean");
          double sigm = pffit->GetParameter("Sigma");
          float meanDen = haveExpectedValue ? expValue : mean;
          labs.push_back("Fitter: " + fnam);
          ostringstream ssout;
          ssout.precision(3);
          ssout.setf(std::ios_base::fixed);
          ssout << "Fit mean: " << mean;
          labs.push_back(ssout.str());
          ssout.str("");
          ssout << "Fit sigma: " << sigm;
          ssout.precision(1);
          if ( meanDen != 0.0 ) ssout << " [" << 100*sigm/meanDen << "%]";
          ssout.precision(3);
          labs.push_back(ssout.str());
          if ( haveExpectedValue ) {
            ssout.str("");
            float bias = mean - expValue;
            ssout << "Fit bias: " << bias;
            ssout.precision(1);
            if ( expValue != 0.0 ) ssout << " [" << 100*bias/expValue << "%]";
            ssout.precision(3);
            labs.push_back(ssout.str());
          }
        }
        {
          ostringstream ssout;
          ssout.str("");
          Index chanStat = getState().getChannelStatus(ph->GetName());
          if ( chanStat == AdcChannelStatusBad ) labs.push_back("Bad channel");
          if ( chanStat == AdcChannelStatusNoisy ) labs.push_back("Noisy channel");
        }
        double xlab = 0.65;
        double ylab = 0.80;
        double dylab = 0.04;
        for ( Name lab : labs ) {
          TLatex* pptl = nullptr;
          pptl = new TLatex(xlab, ylab, lab.c_str());
          pptl->SetNDC();
          pptl->SetTextFont(42);
          pptl->SetTextSize(dylab);
          pman->add(pptl);
          ylab -= 1.2*dylab;
        }
      }
      ++ipad;
      if ( ipad >= npad || ihst+1 >= hsts.size() ) {
        if ( pmantop == nullptr ) {
          if (  m_LogLevel >= 2 ) cout << myname << "  Not writing empty plot" << endl;
        } else {
          if (  m_LogLevel >= 2 ) cout << myname << "  Writing " << plotFileName << endl;
          pmantop->print(plotFileName);
          delete pmantop;
          pmantop = nullptr;
        }
        ipad = 0;
      }
    }
  }
}

//**********************************************************************

void AdcRoiViewer::fillChanSumHists() const {
  const string myname = "AdcRoiViewer::fillChanSumHists: ";
  if ( getState().chanSumHists.size() == 0 ) {
    cout << myname << "No channel summary histograms found." << endl;
    return;
  }
  Name ofrname = m_ChanSumRootFileName;
  if ( m_LogLevel >= 1 ) cout << myname << "Filling channel summary histograms. Count is "
                              << getState().chanSumHists.size() << "." << endl;
  for ( HistMap::value_type ihst : getState().chanSumHists ) {
    TH1* ph = ihst.second;
    Name hnam = ph->GetName();
    Name hnamTemplate = getState().getChanSumHistTemplateName(hnam);
    if ( hnamTemplate.size() == 0 ) {
      cout << myname << "ERROR: Summary template histogram name not found for channel summary " << hnam << endl;
      continue;
    }
    Name vartype = getState().getChanSumHistVariableType(hnam);
    Name errtype = getState().getChanSumHistErrorType(hnam);
    double errmin = 0.0;
    if ( errtype.substr(0,3) == "rms" && errtype.size() > 3 ) {
      string serrmin = errtype.substr(3);
      istringstream ssermin(serrmin);
      ssermin >> errmin;
      errtype = "rms";
    }
    bool isDist = getState().chanSumHistTypes[hnam];
    if ( vartype.size() == 0 ) {
      cout << myname << "ERROR: Variable type name not found for " << hnam << endl;
      continue;
    }
    AdcChannelData acd;
    acd.setEventInfo(getState().cachedRun, 0);
    acd.sampleUnit = getState().cachedSampleUnit;
    Index ncha = 0;
    Index nchaGood = 0;
    int logthresh = 3;
    Index icha1 = getState().chanSumChaBegin[hnam];
    Index icha2 = getState().chanSumChaEnd[hnam];
    //for ( int ibin=1; ibin<=ph->GetNbinsX(); ++ibin ) {
    //  Index icha = ph->GetBinCenter(ibin);
    for ( Index icha=icha1; icha<icha2; ++icha ) {
      Index ichaBin = icha + 1 - icha1;
      acd.setChannelInfo(icha);
      Name hnam = AdcChannelStringTool::build(m_adcStringBuilder, acd, hnamTemplate);
      //Index chanStat = getState().getChannelStatus(hnam);
      //if ( chanStat ) continue;
      ++ncha;
      TH1* phvar = getState().getSumHist(hnam);
      if ( phvar == nullptr ) {
        if ( m_LogLevel >= logthresh )
          cout << myname << "Unable to find sum hist " << hnam << endl;
        continue;
      }
      float val = 0.0;
      if ( vartype == "mean" ) {
        val = phvar->GetMean();
      } else if ( vartype == "peak" ) {
        int ibin = phvar->GetMaximumBin();
        val = phvar->GetBinCenter(ibin);
      } else if ( vartype == "rms" ) {
        val = phvar->GetRMS();
      } else if ( vartype == "entries" ) {
        val = phvar->GetEntries();
      } else if ( vartype == "count" ) {
        val = phvar->Integral();
      } else if ( vartype == "sum" ) {
        val = phvar->Integral()*phvar->GetMean();
      } else if ( vartype.substr(0,3) == "fit" ) {
        Index nfun = phvar->GetListOfFunctions()->GetEntries();
        TF1* pf = nfun ? dynamic_cast<TF1*>(phvar->GetListOfFunctions()->At(0)) : nullptr;
        if ( pf == nullptr ) {
          if ( m_LogLevel >= logthresh )
            cout << myname << "Unable to find find fit for sum hist " << hnam << endl;
          continue;
        }
        bool doRat = vartype.substr(3,3) == "rat";
        string::size_type ipos = doRat ? 6 : 3;
        Name spar = vartype.substr(ipos);
        int ipar = pf->GetParNumber(spar.c_str());
        if ( ipar < 0 ) {
          if ( m_LogLevel >= logthresh )
            cout << myname << "ERROR: Invalid fit parameter name: " << spar << endl;
          continue;
        }
        val = pf->GetParameter(ipar);
        if ( doRat ) {
          double mean = pf->GetParameter("Mean");
          val *= (mean == 0.0 ? 0.0 : 1.0/mean);
        }
      } else {
        cout << myname << "Invalid variable type " << vartype << " for " << hnam << endl;
        break;
      }
      float dval = 0.0;
      bool haveErr = true;
      if ( errtype == "none" ) {
        haveErr = false;
      } else if ( errtype == "zero" ) {
        dval = 0.0;
      } else if ( errtype.substr(0,3) == "rms" ) {
        dval = phvar->GetRMS();
        if ( dval < errmin ) dval = errmin;
      } else if ( errtype == "meanError" ) {
        dval = phvar->GetMeanError();
        if ( dval < errmin ) dval = errmin;
      } else if ( errtype == "rmsError" ) {
        dval = phvar->GetRMSError();
        if ( dval < errmin ) dval = errmin;
      } else if ( errtype.substr(0,3) == "fit" ) {
        Index nfun = phvar->GetListOfFunctions()->GetEntries();
        TF1* pf = nfun ? dynamic_cast<TF1*>(phvar->GetListOfFunctions()->At(0)) : nullptr;
        if ( phvar == nullptr ) {
          if ( m_LogLevel >= logthresh )
            cout << myname << "Unable to find find fit for sum hist " << hnam << endl;
          continue;
        }
        Name spar = errtype.substr(3);
        int ipar = pf->GetParNumber(spar.c_str());
        if ( ipar < 0 ) {
          cout << myname << "ERROR: Invalid fit parameter name: " << spar << endl;
        } else {
          dval = pf->GetParameter(ipar);
        }
      } else {
        cout << myname << "Invalid error type: " << errtype << endl;
        abort();
      }
      if ( ph->GetEntries() == 0 ) {
        TAxis* paxis = isDist ? ph->GetXaxis() : ph->GetYaxis();
        Name vlabOld = paxis->GetTitle();
        Name vlabNew = AdcChannelStringTool::build(m_adcStringBuilder, acd, vlabOld);
        if ( vlabNew != vlabOld ) {
          if ( m_LogLevel >= 3 ) cout << myname << "Setting variable label for " << hnam
                                      << " to \"" << vlabNew << "\"." << endl;
          paxis->SetTitle(vlabNew.c_str());
        }
        Name httlOld = ph->GetTitle();
        Name httlNew = AdcChannelStringTool::build(m_adcStringBuilder, acd, httlOld);
        if ( httlNew != httlOld ) {
          if ( m_LogLevel >= 3 ) cout << myname << "Setting title for " << hnam << " to \""
                                      << httlNew << "\"." << endl;
          ph->SetTitle(httlNew.c_str());
        }
      }
      if ( isDist ) {
        ph->Fill(val);
      } else {
        ph->SetBinContent(ichaBin, val);
        if ( haveErr ) ph->SetBinError(ichaBin, dval);
      }
      ++nchaGood;
    }
    if ( nchaGood < ncha ) {
      cout << myname << "WARNING: Only filled " << nchaGood << " of " << ncha
           << " channels for channel summary histogram " << hnam << endl;
    }
  }
}

//**********************************************************************

void AdcRoiViewer::writeChanSumHists() const {
  const string myname = "AdcRoiViewer::writeChanSumHists: ";
  if ( m_ChanSumRootFileName.size() == 0 ) return;
  if ( getState().chanSumHists.size() == 0 ) {
    cout << myname << "No channel summary histograms found." << endl;
    return;
  }
  Name ofrname = m_ChanSumRootFileName;
  TDirectory* savdir = gDirectory;
  TFile* pfile = TFile::Open(ofrname.c_str(), "UPDATE");
  if ( m_LogLevel >= 1 ) cout << myname << "Writing channel summary histograms. Count is "
                              << getState().chanSumHists.size() << "." << endl;
  for ( HistMap::value_type ihst : getState().chanSumHists ) {
    TH1* ph = ihst.second;
    TH1* phnew = dynamic_cast<TH1*>(ph->Clone());
    phnew->Write();
    if ( m_LogLevel >= 2 ) cout << myname << "  Wrote " << phnew->GetName() << endl;
  }
  if ( pfile != nullptr ) pfile->Close();
  delete pfile;
  if ( m_LogLevel >= 1 ) cout << myname << "Closed summary histogram file " << ofrname << endl;
  savdir->cd();
}

//**********************************************************************

void AdcRoiViewer::writeChanSumPlots() const {
  const string myname = "AdcRoiViewer::writeChanSumPlots: ";
  for ( HistMap::value_type ihst : getState().chanSumHists ) {
    TH1* ph = ihst.second;
    Name hnam = ph->GetName();
    Name pnam = getState().getChanSumPlotName(hnam);
    if ( pnam.size() == 0 ) continue;
    if ( m_LogLevel >= 2 ) cout << myname << "Hist:Plot name: " << hnam << ":" << pnam << endl;
    bool isDist = getState().chanSumHistTypes[hnam];
    Index wpadx = isDist ? 700 : 1400;
    Index wpady = 500;
    TPadManipulator* pman = new TPadManipulator(wpadx, wpady);
    Name plotOpt = isDist ? "hist" : "p";
    pman->add(ph, plotOpt, false);
    pman->addAxis();
    pman->showUnderflow();
    pman->showOverflow();
    if ( ! isDist ) {
      bool doRange = false;
      float ymin = ph->GetMinimum();
      float ymax = ph->GetMaximum();
      Name yopt;
      if ( getState().chanSumPlotYMins.find(hnam) != getState().chanSumPlotYMins.end() ) {
        ymin = getState().chanSumPlotYMins[hnam];
        doRange = true;
      }
      if ( getState().chanSumPlotYMaxs.find(hnam) != getState().chanSumPlotYMaxs.end() ) {
        ymax = getState().chanSumPlotYMaxs[hnam];
        doRange = true;
      }
      if ( doRange ) {
        Name yopt = getState().chanSumPlotYOpts[hnam];
        double yfac = 0.0;
        if ( yopt == "pamp") {
          const RunData& rdat = getState().runData;
          if ( rdat.havePulserAmplitude() ) {
            yfac = rdat.pulserAmplitude();
          } else {
            yfac = 10.0;
            cout << myname << "WARNING: Scaling option pamp requested without run data." << endl;
          }
        } else if ( yopt == "pampg14") {
          const RunData& rdat = getState().runData;
          if ( rdat.havePulserAmplitude() && rdat.haveGain() ) {
            yfac = rdat.pulserAmplitude()*rdat.gain()/14.0;
          } else {
            yfac = 10.0;
            cout << myname << "WARNING: Scaling option pampg14 requested without run data." << endl;
          }
        } else if ( yopt == "nevt") {
          yfac = getState().eventCallCount.size();
        }
        if ( yfac != 0.0 ) {
          ymin *= yfac;
          ymax *= yfac;
        }
        // If ymin > ymax and ymin > 0, then we center histogram on median and use width = ymin.
        if ( ymin > ymax && ymin > 0.0 ) {
          int nval = ph->GetNbinsX();
          FloatVector tmpvals(nval);
          for ( int ival=0; ival<nval; ++ival ) tmpvals[ival] = ph->GetBinContent(ival+1);
          std::sort(tmpvals.begin(), tmpvals.end());
          float ymed = 0.5*(tmpvals[(nval-1)/2] + tmpvals[nval/2]);
          float width = ymin;
          ymin = ymed - 0.5*width;
          ymax = ymin + width;
        }
        if ( ymax > ymin ) {
          if ( m_LogLevel >= 2 ) cout << myname << "Setting plot range to (" << ymin << ", " << ymax << ")" << endl;
          pman->setRangeY(ymin, ymax);
          TH1* php = pman->hist();
          double del = 1.e-4*(ymax -ymin);
          // Put points on the page.
          for ( int ibin=1; ibin<php->GetNbinsX(); ++ibin ) {
            if ( php->GetBinContent(ibin) > ymax ) php->SetBinContent(ibin, ymax-del);
            if ( php->GetBinContent(ibin) < ymin ) php->SetBinContent(ibin, ymin+del);
          }
        }
      }
    }
    bool highlightBadChannels = ! isDist;
    if ( highlightBadChannels ) {
      TH1* php = pman->hist();
      LineColors cols;
      TGraph* pgb = new TGraph;
      TGraph* pgn = new TGraph;
      pgb->SetMarkerStyle(4);
      pgb->SetMarkerColor(cols.red());
      pgn->SetMarkerStyle(4);
      pgn->SetMarkerColor(cols.brown());
      Index icha0 = php->GetBinLowEdge(1);
      Index nbad = 0;
      Index nnoi = 0;
      for ( int ibin=1; ibin<php->GetNbinsX(); ++ibin ) {
        Index icha = icha0 + Index(ibin) - 1;
        Index chanStat = getState().getChannelStatus(icha);
        if ( chanStat == AdcChannelStatusBad ) {
          pgb->SetPoint(nbad++, icha+0.5, php->GetBinContent(ibin));
        }
        if ( chanStat == AdcChannelStatusNoisy ) {
          pgn->SetPoint(nnoi++, icha+0.5, php->GetBinContent(ibin));
        }
      }
      if ( nnoi ) pman->add(pgn, "P");
      if ( nbad ) pman->add(pgb, "P");
      delete pgb;
      delete pgn;
    }
    if ( m_ChannelLineModulus ) {
      for ( Index icha : m_ChannelLinePattern ) {
        pman->addVerticalModLines(m_ChannelLineModulus, icha);
      }
    } else {
      for ( Index icha : m_ChannelLinePattern ) {
        pman->addVerticalLine(icha, 1.0, 3);
      }
    }
    if ( m_LogLevel >= 1 ) cout << myname << "Plotting channel summary " << pnam << endl;
    pman->print(pnam);
    delete pman;
  }
}

//**********************************************************************

void AdcRoiViewer::setPlotLabels(Name& sttl) const {
  StringManipulator sman(sttl, false);
  for ( NameMap::value_type isub : m_plotLabelSubs ) {
    sman.replace(isub.first, isub.second);
  }
  sttl = sman.str();
}

//**********************************************************************

DEFINE_ART_CLASS_TOOL(AdcRoiViewer)
