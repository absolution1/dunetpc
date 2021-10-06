// AdcChannelDftPlotter_tool.cc

#include "AdcChannelDftPlotter.h"
#include "dune/DuneInterface/Tool/AdcChannelStringTool.h"
#include "dune/DuneCommon/Utility/TPadManipulator.h"
#include "dune/DuneCommon/Utility/StringManipulator.h"
#include "dune/DuneCommon/Utility/LineColors.h"
#include "dune/ArtSupport/DuneToolManager.h"
#include "canvas/Utilities/Exception.h"
#include <iostream>
#include <sstream>
#include <vector>
#include <set>
#include <iomanip>
#include "TH1F.h"
#include "TGraph.h"

using std::string;
using std::cout;
using std::cin;
using std::endl;
using std::vector;
using std::ostringstream;
using std::setw;
using std::fixed;
using Index = AdcChannelDftPlotter::Index;
using IndexSet = std::set<Index>;
using IntVector = DataMap::IntVector;
using AcdVector = AdcChannelDftPlotter::AcdVector;

//**********************************************************************
// Local functions.
//**********************************************************************

namespace {

Index channelCount(const IntVector& chans) {
  Index ncha = 0;
  for ( IntVector::const_iterator iven=chans.begin(); iven!=chans.end(); ++iven ) {
    if ( find(chans.begin(), iven, *iven) == iven ) ++ncha;
  }
  return ncha;
}

Index channelCount(const AcdVector& acds) {
  IndexSet ichas;
  for ( const AdcChannelData* pacd : acds ) {
    ichas.insert(pacd->channel());
  }
  return ichas.size();
}

}

//**********************************************************************
// Class methods.
//**********************************************************************

AdcChannelDftPlotter::AdcChannelDftPlotter(fhicl::ParameterSet const& ps)
: AdcMultiChannelPlotter(ps, "Plot"),
  m_Variable(ps.get<Name>("Variable")),
  m_ChannelStatusFlag(ps.get<Index>("ChannelStatusFlag")),
  m_ChannelSelection(ps.get<Name>("ChannelSelection")),
  m_SampleFreq(ps.get<float>("SampleFreq")),
  m_XMin(0.0),
  m_XMax(0.0),
  m_YMinLog(ps.get<float>("YMinLog")),
  m_NBinX(0),
  m_HistName(ps.get<Name>("HistName")),
  m_HistTitle(ps.get<Name>("HistTitle")),
  m_HistSummaryTitles(ps.get<NameVector>("HistSummaryTitles")),
  m_pstate(new State) {
  const string myname = "AdcChannelDftPlotter::ctor: ";
  bool doMag = m_Variable == "magnitude";
  bool doPha = m_Variable == "phase";
  bool doPwr = m_Variable == "power";
  bool doPwt = m_Variable == "power/tick";
  // Check variable and get optional fields.
  if ( doPwr || doPwt ) {
    m_XMin = ps.get<double>("XMin");
    m_XMax = ps.get<double>("XMax");
    m_YMax = ps.get<float>("YMax");
    m_NBinX = ps.get<Index>("NBinX");
  } else if ( doMag ) {
    m_YMax = ps.get<float>("YMax");
  } else if ( doPha ) {
  } else {
    cout << myname << "Invalid Variable: " << m_Variable << endl;
    throw art::Exception(art::errors::Configuration, "InvalidFclValue");
  }
  // Fetch renaming tools.
  string snameBuilder = "adcStringBuilder";
  DuneToolManager* ptm = DuneToolManager::instance();
  m_adcStringBuilder = ptm->getShared<AdcChannelStringTool>(snameBuilder);
  if ( m_adcStringBuilder == nullptr ) {
    cout << myname << "WARNING: AdcChannelStringTool not found: " << snameBuilder << endl;
  }
  // Derived config.
  m_skipBad = m_ChannelStatusFlag==1 || m_ChannelStatusFlag==3;
  m_skipNoisy = m_ChannelStatusFlag==2 || m_ChannelStatusFlag==3;
  m_shiftFreq0 = (doPwr || doPwt) && (m_XMin >= m_XMax);
  if ( m_ChannelSelection.size() ) {
    m_ptfsel.reset(new TFormula("AdcChannelDftPlotter", m_ChannelSelection.c_str()));
  }
  // Display the configuration.
  if ( getLogLevel() >= 1 ) {
    cout << myname << "Configuration: " << endl;
    cout << myname << "           Variable: " << m_Variable << endl;
    cout << myname << "  ChannelStatusFlag: " << m_ChannelStatusFlag;
    if ( m_skipBad ) {
      if ( m_skipNoisy ) cout << " (skip bad and noisy)";
      else cout << " (skip bad)";
    } else if ( m_skipNoisy ) cout << " (skip noisy)";
    cout << endl;
    cout << myname << "   ChannelSelection: " << m_ChannelSelection << endl;
    cout << myname << "         SampleFreq: " << m_SampleFreq << endl;
    if ( doMag || doPwr || doPwt ) cout << myname << "              NBinX: " << m_NBinX << endl;
    if ( doPwr || doPwt ) {
      cout << myname << "               XMin: " << m_XMin << endl;
      cout << myname << "               XMax: " << m_XMax << endl;
      cout << myname << "               YMax: " << m_YMax << endl;
    }
    cout << myname << "            YMinLog: " << m_YMinLog << endl;
    cout << myname << "           HistName: " << m_HistName << endl;
    cout << myname << "          HistTitle: " << m_HistTitle << endl;
    cout << myname << "  HistSummaryTitles: [";
    bool first = true;
    for ( Name name :  m_HistSummaryTitles ) {
      if ( first ) first = false;
      else cout << ", ";
      cout << name;
    }
    cout << "]" << endl;
  }
}

//**********************************************************************

AdcChannelDftPlotter::~AdcChannelDftPlotter() {
  const string myname = "AdcChannelDftPlotter::dtor: ";
  if ( getLogLevel() >= 2 ) {
    cout << myname << "Closing." << endl;
    if ( getChannelRangeNames().size() ) {
      cout << myname << "Channel ranges" << endl;
      cout << myname << "     CR name    count nch/evt" << endl;
      for ( Name crn : getChannelRangeNames() ) {
        Index count = getJobState().count(crn);
        double nchan = getJobState().nchan(crn);
        cout << myname << setw(15) << crn << ":"
             << setw(8) << count << setw(8) << nchan/count << endl;
      }
    } else {
      cout << myname << "No channel ranges specified." << endl;
    }
    if ( getChannelRangeNames().size() ) {
      for ( Name cgn : getChannelGroupNames() ) {
        cout << myname << "Channel group " << cgn << endl;
        cout << myname << "     CR name    count nch/evt" << endl;
        const IndexRangeGroup& crg = getChannelGroup(cgn);
        for ( const IndexRange& ran : crg.ranges ) {
          Name crn = ran.name;
          Index count = getJobState().count(crn);
          double nchan = getJobState().nchan(crn);
          cout << myname << setw(15) << crn << ":"
               << setw(8) << count << setw(8) << nchan/count << endl;
        }
      }
    } else {
      cout << myname << "No channel groups specified." << endl;
    }
  }
  viewSummary(0);
}

//**********************************************************************

int AdcChannelDftPlotter::
viewMapChannels(Name crn, const AcdVector& acds, TPadManipulator& man, Index ncr, Index icr) const {
  const string myname = "AdcChannelDftPlotter::viewMapChannels: ";
  DataMap chret = viewLocal(crn, acds);
  bool doState = true;
  bool doFill = acds.size();
  if ( doState ) {
    ++getJobState().count(crn);
    ++getEventState().count(crn);
    Index jobCount = getJobState().count(crn);
    Index evtCount = getEventState().count(crn);
    bool doPwr = m_Variable == "power";
    bool doPwt = m_Variable == "power/tick";
    if ( doPwr || doPwt ) {
      TH1* ph = chret.getHist("dftHist");
      if ( ph != nullptr ) {
        TH1*& phjob = getJobState().hist(crn);
        if ( phjob == nullptr ) {
          if ( jobCount == 1 ) {
            phjob = dynamic_cast<TH1*>(ph->Clone());
            phjob->SetDirectory(nullptr);
            phjob->SetStats(0);
          } else {
            cout << myname << "ERROR: Hist missing for job count " << jobCount << endl;
          }
        } else {
          if ( jobCount > 1 ) phjob->Add(ph);
          else cout << myname << "ERROR: Hist existing for job count " << jobCount << endl;
        }
        TH1*& phevt = getEventState().hist(crn);
        if ( phevt == nullptr ) {
          if ( evtCount == 1 ) {
            phevt = dynamic_cast<TH1*>(ph->Clone());
            phevt->SetDirectory(nullptr);
            phevt->SetStats(0);
          } else {
            cout << myname << "ERROR: Hist missing for event count " << evtCount << endl;
          }
        } else {
          if ( evtCount > 1 ) phevt->Add(ph);
          else cout << myname << "ERROR: Hist existing for event count " << evtCount << endl;
        }
        const IntVector& allchans = chret.getIntVector("dftChannels");
        Index nven = allchans.size();
        Index ncha = channelCount(allchans);
        getJobState().nchan(crn) += ncha;
        getJobState().nviewentry(crn) += nven;
        getEventState().nchan(crn) += ncha;
        getEventState().nviewentry(crn) += nven;
      } else {
        doFill = false;      // No data in CR
      }
    }
  }
  chret.setString("dftCRLabel", crn);
  chret.setInt("dftCRCount", 1);
  chret.setInt("dftCRIndex", icr);
  if ( doFill ) fillPad(chret, man);
  return 0;
}

//**********************************************************************

int AdcChannelDftPlotter::
viewMapSummary(Index ilev, Name cgn, Name crn, TPadManipulator& man, Index ncr, Index icr) const {
  const string myname = "AdcChannelDftPlotter::viewMapSummary: ";
  if ( getLogLevel() >= 2 ) {
    cout << myname << "Processing " << cgn << "/" << crn << " (" << icr << "/" << ncr << ")" << endl;
  }
  if ( ilev > 1 ) {
    cout << myname << "ERROR: Invalid level: " << ilev << endl;
    return 12;
  }
  if ( icr >= ncr ) {
    cout << myname << "ERROR: Too many plots: " << icr << " >= " << ncr << endl;
    return 11;
  }
  SubState& levState = getSubState(ilev);
  Index count = levState.count(crn);
  Index nchanTot = levState.nchan(crn);
  float nchanEvt = count > 0 ? double(nchanTot)/count : 0.0;
  Index nvenTot = levState.nviewentry(crn);
  float nvenEvt = count > 0 ? double(nvenTot)/count : 0.0;
  TH1* ph = nullptr;
  TH1* phin = levState.hist(crn);
  bool doFill = false;
  if ( phin != nullptr ) {
    if ( count == 0 ) return 12;
    ph = (phin == nullptr) ? nullptr : dynamic_cast<TH1*>(phin->Clone());
    if ( ph == nullptr ) return 13;
    ph->SetDirectory(nullptr);
    if ( m_HistSummaryTitles.size() ) {
      Name htitl = ilev < m_HistSummaryTitles.size()
                   ? m_HistSummaryTitles[ilev]
                   : m_HistSummaryTitles.back();
      StringManipulator smanTitl(htitl, false);
      Name cglab = getChannelGroup(cgn).label();
      if ( cglab.size() == 0 ) cglab = cgn;
      smanTitl.replace("%CGNAME%", cgn);
      smanTitl.replace("%CGLABEL%", cglab);
      smanTitl.replace("%CRNAME%", crn);
      smanTitl.replace("%RUN%", getBaseState().run());
      smanTitl.replace("%EVENT%", getBaseState().event);
      smanTitl.replace("%VIEW%", getDataView());
      ph->SetTitle(htitl.c_str());
    }
    double fac = 1.0/count;
    ph->Scale(fac);
    doFill = true;
  }
  DataMap dm;
  dm.setHist("dftHist", ph, true);
  dm.setInt("dftEventCount", count);
  dm.setFloat("dftChanPerEventCount", nchanEvt);
  dm.setFloat("dftViewEntryPerEventCount", nvenEvt);
  dm.setString("dftDopt", "hist");
  dm.setString("dftCRLabel", crn);
  dm.setInt("dftCRCount", ncr);
  dm.setInt("dftCRIndex", icr);
  if ( doFill ) fillPad(dm, man);
  //man.add(ph, "hist");
  //delete ph;
  return 0;
}

//**********************************************************************

DataMap AdcChannelDftPlotter::view(const AdcChannelData& acd) const {
  const string myname = "AdcChannelDftPlotter::view: ";
  AcdVector acds(1, &acd);
  DataMap chret = viewLocal("", acds);
  if ( getPlotName().size() ) {
    string pname = AdcChannelStringTool::build(m_adcStringBuilder, acd, getPlotName());
    TPadManipulator man;
    fillPad(chret, man);
    if ( getLogLevel() >= 3 ) cout << myname << "Printing " << pname << endl;
    man.print(pname);
  }
  return chret;
}


//**********************************************************************

DataMap AdcChannelDftPlotter::beginEvent(const DuneEventInfo&) const {
  DataMap ret;
  getEventState().clear();
  return ret;
}

//**********************************************************************

DataMap AdcChannelDftPlotter::endEvent(const DuneEventInfo&) const {
  DataMap ret;
  viewSummary(1);
  getEventState().clear();
  return ret;
}

//**********************************************************************

DataMap AdcChannelDftPlotter::viewLocal(Name crn, const AcdVector& acds) const {
  const string myname = "AdcChannelDftPlotter::viewLocal: ";
  if ( getLogLevel() >= 4 ) {
    cout << myname << "Filling CR " << crn << " with ADC channel data size " << acds.size() << endl;
  }
  DataMap ret;
  if ( acds.size() == 0 ) return ret;
  const AdcChannelData* pacd = acds.front();
  Index evt = pacd->event();
  // Check if there have been any other views of this channel range for this event.
  // For now, we implicity expect configurations that do not repeat ranges.
  // Later we might cache results from an earlier attempt.
  // For now, user can expect some plots to overcount events and maybe worse.
  // Avoid this problem by ensuring configuration does not include any channel
  // range more than once. This includes channel ranges in channel groups.
  if ( getState().setEventChannelRange(evt, crn) ) {
    cout << myname << "WARNING: Duplicate view of channel range " << crn
         << " in event " << evt << endl;
  }
  if ( pacd == nullptr ) {
    cout << myname << "ERROR: First channel has no data." << endl;
    return ret;
  }
  const AdcChannelData& acd = *pacd;
  bool doMag = m_Variable == "magnitude";
  bool doPha = m_Variable == "phase";
  bool doPwr = m_Variable == "power";
  bool doPwt = m_Variable == "power/tick";
  bool haveFreq = m_SampleFreq > 0.0;
  if ( ! doMag && !doPha && !doPwr && !doPwt ) {
    cout << myname << "ERROR: Invalid plot variable: " << m_Variable << endl;
    return ret.setStatus(1);
  }
  Index nmag = acd.dftmags.size();
  Index npha = acd.dftphases.size();
  Index nsam = nmag + npha - 1;
  if ( nmag == 0 ) {
    cout << myname << "ERROR: DFT is not present." << endl;
    return ret.setStatus(2);
  }
  if ( npha > nmag || nmag - npha > 1 ) {
    cout << myname << "ERROR: DFT is not valid." << endl;
    return ret.setStatus(3);
  }
  // Check consistency of input data.
  Index nDataMissing = 0;
  Index nBadMagCount = 0;
  Index nBadPhaCount = 0;
  for ( const AdcChannelData* pacde : acds ) {
    if ( pacde == nullptr ) {
      ++nDataMissing;
    } else { 
      if ( pacde->dftmags.size() != nmag ) ++nBadMagCount;
      if ( pacde->dftphases.size() != npha ) ++nBadPhaCount;
    }
  }
  if ( nDataMissing ) cout << myname << "ERROR: Missing data channel count is " << nDataMissing << endl;
  if ( nBadMagCount ) cout << myname << "ERROR: Inconsistent mag size channel count is " << nBadMagCount << endl;
  if ( nBadPhaCount ) cout << myname << "ERROR: Inconsistent pha size channel count is " << nBadPhaCount << endl;
  if ( nDataMissing || nBadMagCount || nBadPhaCount ) return ret.setStatus(4);
  string hname = AdcChannelStringTool::build(m_adcStringBuilder, acd, m_HistName);
  string htitl = AdcChannelStringTool::build(m_adcStringBuilder, acd, m_HistTitle);
  StringManipulator smanName(hname, false);
  smanName.replace("%CRNAME%", crn);
  smanName.replace("%VIEW%", getDataView());
  StringManipulator smanTitl(htitl, false);
  smanTitl.replace("%CRNAME%", crn);
  smanTitl.replace("%VIEW%", getDataView());
  //xx
  //sman.replace("%CRNAME%", ran.name);
  //sman.replace("%CRLABEL%", ran.label());
  //sman.replace("%CRLABEL1%", ran.label(1));
  //sman.replace("%CRLABEL2%", ran.label(2));
  float pi = acos(-1.0);
  double xFac = haveFreq ? m_SampleFreq/nsam : 1.0;
  float yValMax = 0.0;
  string dopt;
  string xtitl = haveFreq ? "Frequency [kHz]" : "Frequency index";
  if ( doMag || doPha ) {  
    if ( acds.size() != 1 ) {
      cout << myname << "ERROR: " << (doMag ? "Magnitude" : "Phase")
           << " may only be filled for a single channel." << endl;
        return ret.setStatus(5);
    }
    string ytitl = "Phase";
    if ( doMag ) {
      ytitl = AdcChannelStringTool::build(m_adcStringBuilder, acd, "Amplitude% [SUNIT]%");
    }
    TGraph* pg = new TGraph;
    pg->SetName(hname.c_str());
    pg->SetTitle(htitl.c_str());
    pg->SetMarkerStyle(2);
    pg->SetMarkerColor(602);
    //ph = new TH1F(hname.c_str(), htitl.c_str(), nbin, 0, nbin);
    //ph->SetDirectory(nullptr);
    //ph->SetLineWidth(2);
    string xtitl = haveFreq ? "Frequency [kHz]" : "Frequency index";
    pg->GetXaxis()->SetTitle(xtitl.c_str());
    pg->GetYaxis()->SetTitle(ytitl.c_str());
    for ( Index ipha=0; ipha<nmag; ++ipha ) {
      float mag = acd.dftmags[ipha];
      float pha = ipha < npha ? acd.dftphases[ipha] : 0.0;
      if ( mag < 0 ) {
        pha += ( pha > 0 ? -pi : pi);
        mag = -mag;
      }
      float x = ipha*xFac;
      float y = doPha ? pha : mag;
      if ( y > yValMax ) yValMax = y;
      pg->SetPoint(ipha, x, y);
    }
    ret.setGraph("dftGraph", pg);
    ret.setString("dftDopt", "P");
  } else if ( doPwr || doPwt ) {
    // Build list of retained channels.
    DataMap::IntVector dftChannels;
    AcdVector keepAcds;
    for ( const AdcChannelData* pacd : acds ) {
      if ( m_ChannelStatusFlag ) {
        if ( m_skipBad && pacd->channelStatus()==1 ) continue;
        if ( m_skipNoisy && pacd->channelStatus()==2 ) continue;
      }
      if ( skipChannel(*pacd) ) continue;
      dftChannels.push_back(pacd->channel());
      keepAcds.push_back(pacd);
    }
    ret.setIntVector("dftChannels", dftChannels);
    if ( getLogLevel() >= 4 ) {
      Index ncha = channelCount(acds);
      Index nchaKeep = channelCount(keepAcds);
      cout << myname << "  Retaining " << keepAcds.size() << " of " << acds.size() << " entries in "
           << nchaKeep << " of " << ncha << " channels." << endl;
    }
    string ytitl = doPwr ? "Power" : "Power/tick";
    if ( acd.sampleUnit.size() ) {
      ytitl = AdcChannelStringTool::build(m_adcStringBuilder, acd, ytitl + " [(%SUNIT%)^{2}]");
    }
    double xmin = m_XMin;
    double xmax = m_XMax;
    Index nbin = m_NBinX;
    // If min >= max, then show the full range.
    // And for power plots shift so the zero frequency component is
    // in the underflow bin and the highest is included in the last bin.
    if ( xmin >= xmax ) {
      xmin = 0.0;
      xmax = (nmag-1)*xFac;
      if ( m_shiftFreq0 ) {
        double delx = 0.01*xFac;
        xmin += delx;
        xmax += delx;
      }
    }
    if ( nbin == 0 ) {
      nbin = (xmax - xmin)/xFac + 0.1;
    }
    TH1* ph = new TH1F(hname.c_str(), htitl.c_str(), nbin, xmin, xmax);
    ph->SetDirectory(nullptr);
    ph->SetLineWidth(2);
    ph->GetXaxis()->SetTitle(xtitl.c_str());
    ph->GetYaxis()->SetTitle(ytitl.c_str());
    float pwrFac = 1.0/keepAcds.size();
    if ( ! doPwr ) pwrFac /= nsam;
    for ( Index imag=0; imag<nmag; ++imag ) {
      float x = imag*xFac;
      float y = 0.0;
      for ( const AdcChannelData* pacde : keepAcds ) {
        float mag = pacde->dftmags[imag];
        if ( getLogLevel() >= 5 ) {
          cout << myname << "    Channel " << pacde->channel() << " sqrt(pwr[" << imag << "]): " << mag << endl;
        }
        y += pwrFac*mag*mag;
      }
      ph->Fill(x, y);
    }
    if ( ph->GetBinContent(nbin+1) && xmax > (nmag-1)*xFac ) {
      cout << myname << "ERROR: Full range histogram has overflow." << endl;
    }
    for ( Index ibin=0; ibin<nbin; ++ibin ) {
      double y = ph->GetBinContent(ibin);
      if ( y > yValMax ) yValMax = y;
    }
    dopt = "hist";
    ret.setHist("dftHist", ph, true);
    ret.setString("dftDopt", "hist");
  }
  ret.setFloat("dftYValMax", yValMax);
  return ret;
}

//**********************************************************************

int AdcChannelDftPlotter::fillPad(DataMap& dm, TPadManipulator& man) const {
  const string myname = "AdcChannelDftPlotter::fillPad: ";
  const bool dbg = false;
  TGraph* pg = dm.getGraph("dftGraph");
  TH1* ph = dm.getHist("dftHist");
  float yValMax = dm.getFloat("dftYValMax");
  bool doPha = m_Variable == "phase";
  //bool doPwr = m_Variable == "power";
  bool doPwt = m_Variable == "power/tick";
  string dopt = dm.getString("dftDopt");
  bool logy = false;
  // Assign y limits.
  double ymin = 0.0;
  double ymax = 0.0;
  if ( doPha ) {
    ymax =  3.2;
    ymin = -ymax;
  } else {
    ymin = 0.0;
    if ( m_YMax > 0 ) ymax = m_YMax;
    else if ( m_YMax < 0 && -m_YMax > yValMax ) ymax = -m_YMax;
    else ymax = yValMax*1.02;
  }
  if ( m_YMinLog ) {
    ymin = m_YMinLog;
    logy = true;
  }
  double xmin = 0.0;
  double xmax = 0.0;
  Index ncr = dm.getInt("dftCRCount");
  Index icr = 0;
  bool manyCR = ncr > 1;  // Does this plot have multiple CRs?
  bool lastCR = true;     // Is this the last CR on this plot?
  if ( manyCR ) {
    icr = dm.getInt("dftCRIndex");
    lastCR = icr + 1 == ncr;
  }
  if ( pg != nullptr ) {
    xmax = pg->GetXaxis()->GetXmax();
    xmin = -0.02*xmax;
    xmax *= 1.02;
    if ( manyCR ) {
      int icol = LineColors::color(icr, ncr);
      pg->SetMarkerColor(icol);
    }
    man.add(pg, dopt);
  } else if ( ph != nullptr ) {
    if ( manyCR ) {
      if ( icr > 0 ) dopt += " same";
      int icol = LineColors::color(icr, ncr);
      ph->SetLineColor(icol);
      if ( dbg ) cout << myname << "DEBUG: Color[" << icr << "] = " << icol << ", dopt = " << dopt << endl;
    }
    man.add(ph, dopt);
  } else {
    cout << myname << "ERROR: Neither hist or graph is defined." << endl;
    return 1;
  }
  if ( dbg ) cout << myname << "DEBUG: CR " << icr << "/" << ncr << endl;
  // Build the descriptor strings.
  string snevt;
  if ( dm.haveInt("dftEventCount") ) {
    ostringstream ssout;
    ssout << "N_{ev} = " << dm.getInt("dftEventCount");
    snevt = ssout.str();
  }
  string sncha;  // # unique channels
  string snven;  // # view entries
  {
    ostringstream ssoutch;
    ostringstream ssoutve;
    ssoutch.precision(1);
    ssoutve.precision(1);
    if ( dm.haveFloat("dftChanPerEventCount") ) {
      ssoutch << "N_{ch} = " << fixed << dm.getFloat("dftChanPerEventCount");
      ssoutve << "N_{ve} = " << fixed << dm.getFloat("dftViewEntryPerEventCount");
    } else {
      using IntVector = DataMap::IntVector;
      const IntVector& allchans = dm.getIntVector("dftChannels");
      Index nven = allchans.size();
      Index ncha = 0;
      for ( IntVector::const_iterator iven=allchans.begin(); iven!=allchans.end(); ++iven ) {
        if ( find(allchans.begin(), iven, *iven) == iven ) ++ncha;
      }
      ssoutch << "N_{ch} = " << ncha;
      ssoutve << "N_{ve} = " << nven;
      // Display the view entry count if a view is defined
      // or if the view entry and channl counts differ.
      if ( getDataView().size() == 0 ) {
        if ( nven != ncha ) {
          cout << "ERROR: View entry count differs from channel count: "
               << nven << " != " << ncha << "." << endl;
          ssoutve << " !!!";
        } else {
          ssoutve.str("");
        }
      }
    }
    sncha = ssoutch.str();
    snven = ssoutve.str();
  }
  string spow;
  if ( doPwt ) {
    ostringstream ssout;
    double sum = ph->Integral(0, ph->GetNbinsX()+1);
    ssout.precision(3);
    ssout << "#sqrt{#Sigma} = " << fixed << setw(2) << sqrt(sum);
    spow = ssout.str();
  }
  // If this is the last object added to the plot.
  if ( lastCR ) {
    if ( dbg ) cout << myname << "DEBUG: Closing plot." << endl;
    man.addAxis();
    if ( xmax > xmin ) man.setRangeX(xmin, xmax);
    if ( ymax > ymin ) man.setRangeY(ymin, ymax);
    if ( m_shiftFreq0 ) man.showUnderflow();
    if ( logy ) man.setLogY();
    if ( logy ) man.setGridY();
    double textSize = 0.04;
    int textFont = 42;
    if ( ! manyCR ) {
      double xlab = 0.70;
      double ylab = 0.85;
      double dylab = 1.2*textSize;
      if ( spow.size() ) {
        TLatex* ptxt = new TLatex(xlab, ylab, spow.c_str());
        ptxt->SetNDC();
        ptxt->SetTextFont(textFont);
        ptxt->SetTextSize(textSize);
        man.add(ptxt);
        ylab -= dylab;
      }
      if ( sncha.size() ) {
        TLatex* ptxt = new TLatex(xlab, ylab, sncha.c_str());
        ptxt->SetNDC();
        ptxt->SetTextFont(textFont);
        ptxt->SetTextSize(textSize);
        man.add(ptxt);
        ylab -= dylab;
      }
      if ( snven.size() ) {
        TLatex* ptxt = new TLatex(xlab, ylab, snven.c_str());
        ptxt->SetNDC();
        ptxt->SetTextFont(textFont);
        ptxt->SetTextSize(textSize);
        man.add(ptxt);
        ylab -= dylab;
      }
      if ( snevt.size() ) {
        TLatex* ptxt = new TLatex(xlab, ylab, snevt.c_str());
        ptxt->SetNDC();
        ptxt->SetTextFont(textFont);
        ptxt->SetTextSize(textSize);
        man.add(ptxt);
        ylab -= dylab;
      }
    } else {
      double xlab = 0.35;
      double ylab = 0.85;
      if ( snevt.size() ) {
        TLatex* ptxt = new TLatex(xlab, ylab, snevt.c_str());
        ptxt->SetNDC();
        ptxt->SetTextFont(textFont);
        ptxt->SetTextSize(textSize);
        man.add(ptxt);
      }
    }
  }
  // Update legend.
  if ( manyCR ) {
    TObject* pobj = man.object();
    Name lopt = pg == nullptr ? "l" : "p";
    if ( icr == 0 ) {
      double xlmin = 0.55;
      double xlmax = 0.93;
      double ylmax = 0.90;
      double ylmin = ylmax - 0.05*(ncr+0.5);
      if ( ylmin < 0.40 ) ylmin = 0.40;
      man.addLegend(xlmin, ylmin, xlmax, ylmax);
    } else {
      pobj = man.object(icr);
    }
    TLegend* pleg = man.getLegend();
    Name slab = dm.getString("dftCRLabel");
    if ( pleg != nullptr ) {
      pleg->SetMargin(0.1);   // Fraction of box used for symbols
      if ( spow.size() ) slab += " " + spow;
      if ( sncha.size() ) slab += " " + sncha;
      pleg->AddEntry(pobj, slab.c_str(), lopt.c_str());
    }
  }
  return 0;
}

//**********************************************************************

bool AdcChannelDftPlotter::skipChannel(const AdcChannelData& acd) const {
  const string myname = "AdcChannelDftPlotter::skipChannel: ";
  if ( ! m_ptfsel ) return false;
  Index npar = m_ptfsel->GetNpar();
  vector<double> pars(npar);
  for ( Index ipar=0; ipar<npar; ++ipar ) {
    string spar = m_ptfsel->GetParName(ipar);
    if ( ! acd.hasAttribute(spar) ) {
      cout << myname << "WARNING: Skipping run/event/channel " << acd.run() << "/" << acd.event()
           << "/" << acd.channel() << " because of missing attribute " << spar << endl;
      return true;
    }
    pars[ipar] = acd.getAttribute(spar);
  }
  double val = m_ptfsel->EvalPar(nullptr, &pars[0]);
  return val == 0.0;
}

//**********************************************************************

DEFINE_ART_CLASS_TOOL(AdcChannelDftPlotter)
