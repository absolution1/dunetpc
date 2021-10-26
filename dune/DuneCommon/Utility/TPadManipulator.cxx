// TPadManipulator.cxx

#include "TPadManipulator.h"
#include "StringManipulator.h"
#include <iostream>
#include <sstream>
#include "TPad.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TList.h"
#include "TStyle.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TFrame.h"
#include "TLine.h"
#include "TF1.h"
#include "TPaletteAxis.h"
#include "TError.h"
#include "TSystem.h"
#include "TBuffer.h"
#include "TClass.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TExec.h"
#include "TROOT.h" // gROOT

using std::string;
using std::cout;
using std::endl;
using std::ostringstream;
using Index = unsigned int;

bool dbg = 0;

//**********************************************************************

namespace {

//void EmptyErrorHandler(Int_t, Bool_t, const char*, const char*) { }

}  // end unnamed namespace

//**********************************************************************
// Static members.
//**********************************************************************

TPadManipulator* TPadManipulator::get(Name onam, TDirectory* tdir) {
  Name myname = "TPadManipulator::get: ";
  TPadManipulator* ppad = nullptr;
  TDirectory* pdir = tdir == nullptr ? gDirectory : tdir;
  if ( pdir == nullptr ) {
    cout << myname << "ERROR: Root directory not found." << endl;
    return ppad;
  }
  pdir->GetObject(onam.c_str(), ppad);
  return ppad;
}

//**********************************************************************

TPadManipulator* TPadManipulator::read(Name fnam, Name onam) {
  Name myname = "TPadManipulator::read: ";
  TPadManipulator* ppad = nullptr;
  TFile* pfil = TFile::Open(fnam.c_str(), "READ");
  if ( pfil == nullptr || ! pfil->IsOpen() ) {
    cout << myname << "Unable to open file " << fnam << endl;
  } else {
    ppad = get(onam, pfil);
  }
  delete pfil;
  return ppad;
}

//**********************************************************************
// Non-static members.
//**********************************************************************

TPadManipulator::TPadManipulator()
: m_parent(nullptr), m_ppad(nullptr),
  m_canvasWidth(0), m_canvasHeight(0),
  m_exec("myexec", ""),
  m_marginLeft(-999), m_marginRight(-999),
  m_marginBottom(-999), m_marginTop(-999),
  m_ph(nullptr),
  m_pg(nullptr),
  m_fillColor(0), m_frameFillColor(0),
  m_gridX(false), m_gridY(false),
  m_logX(false), m_logY(false), m_logZ(false),
  m_tickLengthX(0.03), m_tickLengthY(0.0),
  m_ndivX(0), m_ndivY(0),
  m_labSizeX(0.0), m_labSizeY(0.0), m_ttlSize(0.0),
  m_flowHist(nullptr),
  m_showUnderflow(false), m_showOverflow(false),
  m_flowGraph(nullptr),
  m_gflowMrk(0), m_gflowCol(0),
  m_top(false), m_right(false), m_iobjLegend(0),
  m_axisTitleOpt(0) {
  const string myname = "TPadManipulator::ctor: ";
  if ( dbg ) cout << myname << this << endl;
  m_label.SetNDC();
  m_label.SetX(0.01);
  m_label.SetY(0.02);
  m_label.SetTextFont(42);
  m_label.SetTextSize(0.035);
}

//**********************************************************************

TPadManipulator::TPadManipulator(TVirtualPad* ppad)
: TPadManipulator() {
  const string myname = "TPadManipulator::ctor(pad): ";
  cout << myname << "Not implemented--will copy objects from pad." << endl;
  if ( dbg ) cout << myname << this << endl;
}

//**********************************************************************

TPadManipulator::TPadManipulator(Index wx, Index wy, Index nPadX, Index nPadY)
: TPadManipulator() {
  m_canvasWidth = wx;
  m_canvasHeight = wy;
  Index nx = nPadX;
  Index ny = nPadY > 0 ? nPadY : nx;
  if ( nx ) split(nx, ny);
}

//**********************************************************************

TPadManipulator::TPadManipulator(const TPadManipulator& rhs)
: TPadManipulator() {
  *this = rhs;
}

//**********************************************************************

TPadManipulator& TPadManipulator::operator=(const TPadManipulator& rhs) {
  if ( this == &rhs ) return *this;
  bool saveAddDirectory = TH1::AddDirectoryStatus();
  TH1::AddDirectory(false);
  if ( m_ppad != nullptr ) {
    m_ppad->Clear();
  }
  m_canvasWidth = rhs.m_canvasWidth;
  m_canvasHeight = rhs.m_canvasHeight;
  m_exec = rhs.m_exec,
  m_marginLeft = rhs.m_marginLeft;
  m_marginRight = rhs.m_marginRight;
  m_marginBottom = rhs.m_marginBottom;
  m_marginTop = rhs.m_marginTop;
  // Clone drawn objects.
  delete m_ph;
  m_ph = nullptr;
  delete m_pg;
  m_pg = nullptr;
  if ( rhs.haveHist() ) {
    m_ph = dynamic_cast<TH1*>(rhs.hist()->Clone());
  } else if ( rhs.m_pg ) {
    m_pg = dynamic_cast<TGraph*>(rhs.m_pg->Clone());
  }
  m_dopt = rhs.m_dopt;
  m_objs.clear();
  for ( TObject* pobj : rhs.m_objs ) {
    m_objs.push_back(pobj != nullptr ? pobj->Clone() : nullptr);
  }
  m_opts = rhs.m_opts;
  TH1::AddDirectory(saveAddDirectory);
  m_setBounds = rhs.m_setBounds;
  m_setBoundsLog = rhs.m_setBoundsLog;
  // Copy decoration instructions.
  m_fillColor = rhs.m_fillColor;
  m_frameFillColor = rhs.m_frameFillColor;
  m_gridX = rhs.m_gridX;
  m_gridY = rhs.m_gridY;
  m_logX = rhs.m_logX;
  m_logY = rhs.m_logY;
  m_logZ = rhs.m_logZ;
  m_tickLengthX = rhs.m_tickLengthX;
  m_tickLengthY = rhs.m_tickLengthY;
  m_ndivX = rhs.m_ndivX;
  m_ndivY = rhs.m_ndivY;
  m_labSizeX = rhs.m_labSizeX;
  m_labSizeY = rhs.m_labSizeY;
  m_ttlSize = rhs.m_ttlSize;
  m_showUnderflow = rhs.m_showUnderflow;
  m_showOverflow = rhs.m_showOverflow;
  m_gflowOpt = rhs.m_gflowOpt;
  m_gflowMrk = rhs.m_gflowMrk;
  m_gflowCol = rhs.m_gflowCol;
  m_top = rhs.m_top;
  m_right = rhs.m_right;
  m_iobjLegend = rhs.m_iobjLegend;
  rhs.m_title.Copy(m_title); m_title.SetNDC();
  rhs.m_label.Copy(m_label); m_label.SetNDC();
  m_histFuns = rhs.m_histFuns;
  m_hmlXmod = rhs.m_hmlXmod;
  m_hmlXoff = rhs.m_hmlXoff;
  m_hmlXStyle = rhs.m_hmlXStyle;
  m_hmlXLength = rhs.m_hmlXLength;
  m_vmlXmod = rhs.m_vmlXmod;
  m_vmlXoff = rhs.m_vmlXoff;
  m_vmlXStyle = rhs.m_vmlXStyle;
  m_vmlXLength = rhs.m_vmlXLength;
  m_slSlop = rhs.m_slSlop;
  m_slYoff = rhs.m_slYoff;
  m_slStyl = rhs.m_slStyl;
  clearLineObjects();
  m_binLabelsX = rhs.m_binLabelsX;
  m_binLabelsY = rhs.m_binLabelsY;
  m_timeOffset = rhs.m_timeOffset;
  m_timeFormatX = rhs.m_timeFormatX;
  m_timeFormatY = rhs.m_timeFormatY;
  m_subBounds = rhs.m_subBounds;
  m_subMans.clear();
  for ( const TPadManipulator& man : rhs.m_subMans ) {
    m_subMans.emplace_back(man);
  }
  m_iobjLegend = rhs.m_iobjLegend;
  m_axisTitleOpt = rhs.m_axisTitleOpt;
  setParents(true);
  update();
  return *this;
}

//**********************************************************************

TPadManipulator::~TPadManipulator() {
  //cout << "Destroying manipulator " << this << " with pad " << m_ppad
  //     << " and parent " << m_parent << endl;
  if ( m_ppad != nullptr && m_parent == nullptr ) {
    // May 2020: If the pad is not in the Root list of canvases, then it may
    // be that Root has already deleted it.
    if ( gROOT->GetListOfCanvases()->IndexOf(m_ppad) >= 0 ) {
      m_ppad->Close();
      gSystem->ProcessEvents();
      delete m_ppad;
      m_ppad = nullptr;
    }
  }
  for ( TObject* pobj : m_objs ) delete pobj;
  m_objs.clear();
  clearLineObjects();
}

//**********************************************************************

double TPadManipulator::xmin() const {
  if ( m_logX &&
       m_setBoundsLog.x1 > 0.0 && m_setBoundsLog.x2 > 0.0 &&
       m_setBoundsLog.x1 <= m_setBoundsLog.x2 ) return m_setBoundsLog.x1;
  if ( m_setBounds.x2 > m_setBounds.x1 ) return m_setBounds.x1;
  if ( ! havePad() ) return 0.0;
  double val = pad()->GetUxmin();
  return pad()->GetLogx() ? pow(10,val) : val;
}

//**********************************************************************

double TPadManipulator::xmax() const {
  if ( m_logX &&
       m_setBoundsLog.x1 > 0.0 && m_setBoundsLog.x2 > 0.0 &&
       m_setBoundsLog.x1 <= m_setBoundsLog.x2 ) return m_setBoundsLog.x2;
  if ( m_setBounds.x2 > m_setBounds.x1 ) return m_setBounds.x2;
  if ( ! havePad() ) return 0.0;
  double val = pad()->GetUxmax();
  return pad()->GetLogx() ? pow(10,val) : val;
}

//**********************************************************************

double TPadManipulator::ymin() const {
  if ( m_logY &&
       m_setBoundsLog.y1 > 0.0 && m_setBoundsLog.y2 > 0.0 &&
       m_setBoundsLog.y1 <= m_setBoundsLog.y2 ) return m_setBoundsLog.y1;
  if ( m_setBounds.y2 > m_setBounds.y1 ) return m_setBounds.y1;
  if ( ! havePad() ) return 0.0;
  double val = pad()->GetUymin();
  return pad()->GetLogy() ? pow(10,val) : val;
}

//**********************************************************************

double TPadManipulator::ymax() const {
  if ( m_logY &&
       m_setBoundsLog.y1 > 0.0 && m_setBoundsLog.y2 > 0.0 &&
       m_setBoundsLog.y1 <= m_setBoundsLog.y2 ) return m_setBoundsLog.y2;
  if ( m_setBounds.y2 > m_setBounds.y1 ) return m_setBounds.y2;
  if ( ! havePad() ) return 0.0;
  double val = pad()->GetUymax();
  return pad()->GetLogy() ? pow(10,val) : val;
}

//**********************************************************************

double TPadManipulator::zmin() const {
  if ( m_logZ &&
       m_setBoundsLog.z1 > 0.0 && m_setBoundsLog.z2 > 0.0 &&
       m_setBoundsLog.z1 <= m_setBoundsLog.z2 ) return m_setBoundsLog.z1;
  if ( m_setBounds.z2 > m_setBounds.z1 ) return m_setBounds.z1;
  return 0.0;
}

//**********************************************************************

double TPadManipulator::zmax() const {
  if ( m_logZ &&
       m_setBoundsLog.z1 > 0.0 && m_setBoundsLog.z2 > 0.0 &&
       m_setBoundsLog.z1 <= m_setBoundsLog.z2 ) return m_setBoundsLog.z2;
  if ( m_setBounds.z2 > m_setBounds.z1 ) return m_setBounds.z2;
  return 0.0;
}

//**********************************************************************

TCanvas* TPadManipulator::canvas() const {
  if ( m_parent ) return m_parent->canvas();
  return dynamic_cast<TCanvas*>(m_ppad);
}

//**********************************************************************

TCanvas* TPadManipulator::canvas(bool doDraw) {
  if ( m_parent ) return m_parent->canvas(doDraw);
  if ( m_ppad == nullptr && doDraw ) draw();
  return dynamic_cast<TCanvas*>(m_ppad);
}

//**********************************************************************

int TPadManipulator::put(Name onam, TDirectory* tdir) const {
  Name myname = "TPadManipulator::put: ";
  TDirectory* pdir = tdir == nullptr ? gDirectory : tdir;
  if ( pdir == nullptr ) {
    cout << myname << "ERROR: Root directory not found." << endl;
    return 2;
  }
  if ( ! pdir->IsWritable() ) {
    cout << myname << "Root directory " << pdir->GetName() << " is not writable." << endl;
    return 3;
  }
  pdir->WriteObject(this, onam.c_str());
  return 0;
}

//**********************************************************************

int TPadManipulator::write(Name fnam, Name onam) const {
  Name myname = "TPadManipulator::write: ";
  TFile* pfil = TFile::Open(fnam.c_str(), "UPDATE");
  if ( pfil == nullptr || ! pfil->IsOpen() ) {
    cout << myname << "Unable to open file " << fnam << endl;
    return 1;
  }
  int wstat = put(onam, pfil);
  delete pfil;
  return wstat;
}

//**********************************************************************

int TPadManipulator::printOnce(string fname) {
  // If fname has suffix .root or .tpad, we save this object.
  Name::size_type idot = fname.rfind(".");
  if ( idot != Name::npos ) {
    Name suf = fname.substr(idot + 1);
    if ( suf == "root" || suf == "tpad" ) {
      write(fname);
      return 0;
    }
  }
  TCanvas* pcan = canvas(false);
  // If canvas does not yet exist, draw in batch mode to avoid
  // display on screen.
  // We should (but don't yet) also remove the canvas after draw so
  // there are no problems if a draw to screen is attempted later.
  bool setBackToNonBatch = false;
  if ( pcan == nullptr ) {
    bool isBatch = gROOT->IsBatch();
    if ( ! isBatch ) {
      gROOT->SetBatch(true);
      setBackToNonBatch = true;
    }
    pcan = canvas(true);
  }
  if ( pcan == nullptr ) {
    if ( setBackToNonBatch ) gROOT->SetBatch(false);
    return 1;
  }
  // Suppress printing message from Root.
  int levelSave = gErrorIgnoreLevel;
  gErrorIgnoreLevel = 1001;
  // Block non-default (e.g. art) from handling the Root "error".
  // We switch to the Root default handler while making the call to Print.
  ErrorHandlerFunc_t pehSave = nullptr;
  ErrorHandlerFunc_t pehDefault = DefaultErrorHandler;
  if ( GetErrorHandler() != pehDefault ) {
    pehSave = SetErrorHandler(pehDefault);
  }
  pcan->Print(fname.c_str());
  if ( pehSave != nullptr ) SetErrorHandler(pehSave);
  gErrorIgnoreLevel = levelSave;
  if ( setBackToNonBatch ) gROOT->SetBatch(false);
  return 0;
}

//**********************************************************************

int TPadManipulator::print(string fnamepat, string spat) {
  StringManipulator sman(fnamepat, true);
  int stats = 0;
  for ( string fnam : sman.patternSplit(spat) ) stats += printOnce(fnam);
  return stats;
}

//**********************************************************************

TPadManipulator* TPadManipulator::progenitor() {
  if ( haveParent() ) return parent()->progenitor();
  return this;
}

//**********************************************************************

unsigned int TPadManipulator::npad() const {
  return m_subMans.size();
}

//**********************************************************************

TPadManipulator* TPadManipulator::man(unsigned int ipad) {
  if ( npad() == 0 && ipad == 0 ) return this;
  return ipad<npad() ? &m_subMans[ipad] : nullptr;
}

//**********************************************************************

TObject* TPadManipulator::object() const {
  if ( haveHist() ) return hist();
  return graph();
}

//**********************************************************************

TH1* TPadManipulator::getHist(unsigned int iobj) {
  if ( iobj >= objects().size() ) return nullptr;
  return dynamic_cast<TH1*>(objects()[iobj]);
}

//**********************************************************************

TObject* TPadManipulator::lastObject() const {
  Index nobj = objects().size();
  if ( nobj ) return objects().back();
  return object();
}

//**********************************************************************

TH1* TPadManipulator::getHist(string hnam) {
  if ( hist() != nullptr && hist()->GetName() == hnam ) return hist();
  for ( TObject* pobj : objects() ) {
    TH1* ph = dynamic_cast<TH1*>(pobj);
    if ( ph == nullptr ) continue;
    if ( ph->GetName() == hnam ) return ph;
  }
  return nullptr;
}

//**********************************************************************

int TPadManipulator::canvasPixelsX() const {
  if ( ! haveCanvas() ) return 0;
  return canvas()->GetWw();
}

//**********************************************************************

int TPadManipulator::canvasPixelsY() const {
  if ( ! haveCanvas() ) return 0;
  return canvas()->GetWh();
}

//**********************************************************************

double TPadManipulator::padPixelsX() const {
  if ( ! havePad() ) return 0.0;
  double frac = pad()->GetAbsWNDC();
  return frac*canvasPixelsX();
}

//**********************************************************************

double TPadManipulator::padPixelsY() const {
  if ( ! havePad() ) return 0.0;
  double frac = pad()->GetAbsHNDC();
  return frac*canvasPixelsY();
}

//**********************************************************************

TFrame* TPadManipulator::frame() const {
  return havePad() ? pad()->GetFrame() : nullptr;
}

//**********************************************************************

int TPadManipulator::framePixelsX() const {
  return haveFrame() ? frame()->GetBBox().fWidth : 0;
}

//**********************************************************************

int TPadManipulator::framePixelsY() const {
  return haveFrame() ? frame()->GetBBox().fHeight : 0;
}

//**********************************************************************

TH1* TPadManipulator::frameHist() const {
  if ( ! havePad() ) return nullptr;
  return dynamic_cast<TH1*>(pad()->FindObject("hframe"));
}

//**********************************************************************

void TPadManipulator::clearLineObjects() {
  for ( TLine* pline : m_lines ) delete pline;
  m_lines.clear();
}

//**********************************************************************

int TPadManipulator::addPad(double x1, double y1, double x2, double y2, int icol) {
  if ( x2 <= x1 ) return 1;
  if ( y2 <= y1 ) return 2;
  if ( x1 < 0.0 ) return 3;
  if ( x2 > 1.0 ) return 4;
  if ( y1 < 0.0 ) return 5;
  if ( y2 > 1.0 ) return 6;
  m_subBounds.emplace_back(x1, y1, x2, y2);
  m_subMans.emplace_back();
  m_subMans.back().m_parent = this;
  return 0;
}

//**********************************************************************

int TPadManipulator::split(Index nx, Index nyin) {
  Index ny = nyin;
  if ( nx < 1 ) return 1;
  if ( ny < 1 ) ny = nx;
  double dx = 1.0/nx;
  double dy = 1.0/ny;
  double y2 = 1.0;
  for ( Index iy=0; iy<ny; ++iy ) {
    double y1 = y2 - dy;
    double x1 = 0.0;
    for ( Index ix=0; ix<nx; ++ix ) {
      double x2 = x1 + dx;
      int rstat = addPad(x1, y1, x2, y2);
      if ( rstat ) return 100 + rstat;
      x1 = x2;
    }
    y2 = y1;
  }
  return 0;
}

//**********************************************************************

int TPadManipulator::split(Index nx) {
  return split(nx, nx);
}

//**********************************************************************

void TPadManipulator::setPalette(int pal) {
  string scom;
  if ( pal >= 0 ) {
    ostringstream sscom;
    sscom << "RootPalette::set(" << pal << ")";
    scom = sscom.str();
  }
  setTExec(scom);
}

//**********************************************************************

int TPadManipulator::add(Index ipad, TObject* pobj, string sopt, bool replace) {
  const string myname = "TPadManipulator::add: ";
  if ( dbg ) cout << myname << this << endl;
  if ( pobj == nullptr ) return 101;
  TPadManipulator* pman = man(ipad);
  if ( pman == nullptr ) return 102;
  TH1* ph = dynamic_cast<TH1*>(pobj);
  TGraph* pg = dynamic_cast<TGraph*>(pobj);
  bool isNotHistOrGraph = ph==nullptr && pg==nullptr;
  // If we already have the primary histogram or graph, or this is not one of those,
  // this is an overlaid object.
  if ( (!replace && haveHistOrGraph()) || isNotHistOrGraph ) {
    TObject* pobjc = pobj->Clone();
    TH1* phc = dynamic_cast<TH1*>(pobjc);
    if ( phc != nullptr ) phc->SetDirectory(nullptr);
    m_objs.push_back(pobjc);
    m_opts.push_back(sopt);
  // Otherwise, the passed object becomes the primary object and must be
  // a histogram or graph.
  } else {
    if ( ph == nullptr && pg == nullptr ) return 103;
    if ( pman->hist() != nullptr ) {
      if ( replace ) delete pman->m_ph;
      else return 104;
    }
    if ( pman->graph() != nullptr ) {
      if ( replace ) delete pman->m_pg;
      else return 105;
    }
    //if ( m_ppad == nullptr ) return 106;
    // Transfer the hist/graph title to the pad title.
    pman->setTitle(pobj->GetTitle());
    // Clone the primary object.
    if ( ph != nullptr ) {
      TH1* phc = (TH1*) ph->Clone();
      phc->SetDirectory(nullptr);
      phc->SetTitle("");
      pman->m_ph = phc;
      pman->m_dopt = sopt;
    } else {
      TGraph* pgc = (TGraph*) pg->Clone();
      pgc->SetTitle("");
      pman->m_pg = pgc;
      string soptOut;
      for ( Index ipos=0; ipos<sopt.size(); ++ipos ) {
        char ch = sopt[ipos];
        if ( ch == 'a' || ch == 'A' ) {
          cout << myname << "WARNING: Dropping \"" << ch << "\" from drawing option string for graph "
               << pgc->GetName() << "." << endl;
        } else {
          soptOut += ch;
        }
      }
      pman->m_dopt = soptOut;
    }
  }
  return update();
} 

//**********************************************************************

int TPadManipulator::add(TObject* pobj, string sopt, bool replace) {
  return add(0, pobj, sopt, replace);
}

//**********************************************************************

TLegend* TPadManipulator::addLegend(double x1, double y1, double x2, double y2) {
  TLegend leg(x1, y1, x2, y2);
  add(0, &leg, "");
  TLegend* pleg = dynamic_cast<TLegend*>(objects().back());
  m_iobjLegend = objects().size() - 1;
  pleg->SetBorderSize(0);
  pleg->SetFillStyle(0);
  return pleg;
}
  
//**********************************************************************

int TPadManipulator::setTitle(string sttl, float height) {
  m_title.SetTitle(sttl.c_str());
  m_title.SetNDC();
  m_title.SetTextAlign(22);
  m_title.SetTextFont(42);
  float tsiz = height;
  if ( tsiz <= 0.0 ) {
    tsiz = getTitleSize();
    if ( tsiz <= 0.0 ) tsiz = 0.035;
  }
  m_title.SetTextSize(tsiz);
  double xttl = 0.5;
  double yttl = 1.0 - 0.70*tsiz;
  m_title.SetText(xttl, yttl, sttl.c_str());
  return 0;
}

//**********************************************************************

int TPadManipulator::setLabel(string slab) {
  m_label.SetTitle(slab.c_str());
  return 0;
}

//**********************************************************************

int TPadManipulator::clear() {
  delete m_ph;
  m_ph = nullptr;
  delete m_pg;
  m_pg = nullptr;
  if ( pad() != nullptr && npad() == 0 ) pad()->Clear();
  for ( TPadManipulator& man : m_subMans ) man.clear();
  return 0;
}

//**********************************************************************

int TPadManipulator::update() {
  const string myname = "TPadManipulator::update: ";
  if ( dbg ) cout << myname << this << endl;
  int rstat = 0;
  // No action if canvas/pad is not yet created.
  // Call draw() to create the pad.
  if ( ! havePad() ) return 0;
  // Update the subpads.
  // If needed, create the TPad for each subpad.
  TVirtualPad* pPadSave = gPad;
  m_ppad->cd();
  if ( m_subMans.size() ) {
    for ( Index isub=0; isub<m_subMans.size(); ++isub ) {
      TPadManipulator& man = m_subMans[isub];
      if ( man.m_ppad == nullptr ) {
        if ( m_subBounds.size() <= isub ) {
          cout << myname << "ERROR: Too few bounds." << endl;
          break;
        }
        const Bounds& bnd = m_subBounds[isub];
        ostringstream ssnam;
        ssnam << m_ppad->GetName() << "_sub" << isub;
        string spnam = ssnam.str();
        string spttl = spnam;
        TPad* ppad = new TPad(spnam.c_str(), spttl.c_str(), bnd.x1, bnd.y1, bnd.x2, bnd.y2);
        ppad->Draw();
        man.m_ppad = ppad;
        man.m_parent = this;
      }
      int srstat = man.update();
      if ( srstat ) {
        cout << myname << "WARNING: Error " << srstat << " updating pad " << isub << endl;
        rstat += srstat;
      }
    }
  }
  // If frame is not yet drawn, use the primary object to draw it.
  // Note that we will later redraw the frame.
  if ( ! haveFrameHist() ) {
    // Fetch the set bounds.
    if ( haveHist() ) {
      hist()->Draw(m_dopt.c_str());
    } else if ( m_pg ) {
      // If the graph has no points, we add one because root (6.12/06) raises an
      // exception if we draw an empty graph.
      if ( m_pg->GetN() == 0 ) {
        double xmin = m_pg->GetXaxis()->GetXmin();
        double ymin = m_pg->GetYaxis()->GetXmin();
        m_pg->SetPoint(0, xmin, ymin);
      }
      m_pg->Draw("AP");
    }
  }
/*
  if ( ! haveHistOrGraph() ) {
    const TList* prims = m_ppad->GetListOfPrimitives();
    for ( TObjLink* plnk=prims->FirstLink(); plnk; plnk=plnk->Next() ) {
      TObject* pobj = plnk->GetObject();
      string sopt = plnk->GetOption();
      TH1* ph = dynamic_cast<TH1*>(pobj);
      if ( ph != nullptr )  {
        m_ph.reset(dynamic_cast<TH1*>(ph->Clone("hmanip0")));
        m_ph->SetDirectory(nullptr);
        m_dopt = sopt;
        break;
      }
      TGraph* pg = dynamic_cast<TGraph*>(pobj);
      if ( pg != nullptr )  {
        m_pg.reset(dynamic_cast<TGraph*>(pg->Clone("gmanip0")));
        m_dopt = sopt;
        break;
      }
    }
  }
*/
  if ( ! haveHistOrGraph() ) {
    if ( npad() == 0 && !haveParent() ) {
      cout << myname << "Top-level pad does not have a histogram or graph or subpads!" << endl;
    }
    // Add the title and labels.
    m_title.Draw();
    if ( getLabel().size() ) m_label.Draw();
    gPad = pPadSave;
    return 0;
  }
  // Run TExec.
  if ( string(m_exec.GetTitle()).size() ) {
    m_exec.Draw();
  }
  // Set margins the first time the histogram is found.
  // After this, user can override with pad()->SetRightMargin(...), ...
  m_ppad->Update();    // This is needed to get color palette for 2D hists.
  bool isTH = haveHist();
  bool isTH2 = dynamic_cast<TH2*>(hist()) != nullptr;
  bool isTH1 = isTH && !isTH2;
  double xm0 = 0.015;
  double wx = padPixelsX();
  double wy = padPixelsY();
  double asp = wx > 0 ? wy/wx : 1.0;
  double aspx = asp < 1.0 ? asp : 1.0;        // Font is proportional to this for asp < 1.0
  double aspy = asp > 1.0 ? 1.0/asp : 1.0;    // Font is proportional to this for asp > 1.0
  if ( false ) {
    aspx = 1.0;
    aspy = 1.0/asp;
  }
  double xml = xm0 + 0.120*aspx;
  double xmr = 0.05*aspx;
  double xmb =       0.100*aspy;
  double xmt =       0.070*aspy;
  double xlb = -0.028 + 0.038*aspy;
  double xlz = 0.005*aspx;
  double xttl = 1.2*aspy;
  double zttl = 1.5*aspx;
  double yttl = 0.17 + 1.8*aspx;
  //double httl = 1.0 - 0.5*xmt;
  if ( m_marginLeft >= 0.0 ) {
    // When left margin is changed, we leave the y-axis title at the edge of the pad.
    xml = m_marginLeft;
    //double scalefac = m_marginLeft/xml;
    //yttl *= scalefac;
  }
  if ( m_marginRight >= 0.0 ) xmr = m_marginRight;
  if ( m_marginBottom >= 0.0 ) {
    xmb = m_marginBottom;
    //double scalefac = m_marginBottom/xmb;
    //xttl *= scalefac;
  }
  if ( m_marginTop >= 0.0 ) xmt = m_marginTop;
  if ( isTH2 ) {
    TPaletteAxis* pax = dynamic_cast<TPaletteAxis*>(hist()->GetListOfFunctions()->FindObject("palette"));
    if ( pax != nullptr ) {
      if ( m_marginRight < 0 ) xmr = xm0 + 0.150*aspx;
      double xp1 = 1.0 - 1.00*xmr;
      double xp2 = 1.0 - 0.65*xmr;
      double yp1 = xmb;
      double yp2 = 1.0 - xmt;
      double tlen = 0.25*(xp2 - xp1);
      //cout << myname << "Palette axis found." << endl;
      pax->SetX1NDC(xp1);
      pax->SetX2NDC(xp2);
      pax->SetY1NDC(yp1);
      pax->SetY2NDC(yp2);
      hist()->GetZaxis()->SetTickLength(tlen);
    //} else {
    //  cout << myname << "Palette axis not found." << endl;
    //  hist()->GetListOfFunctions()->Print();
    }
  }
  m_ppad->SetRightMargin(xmr);
  m_ppad->SetLeftMargin(xml);
  m_ppad->SetTopMargin(xmt);
  m_ppad->SetBottomMargin(xmb);
  // Set the axis tick lengths.
  // If the Y-size is zero, then they are drawn to have the same pixel length as the X-axis.
  double ticklenx = m_tickLengthX;
  double tickleny = m_tickLengthY;
  if ( ticklenx < 0.0 ) ticklenx = 0.0;
  if ( tickleny <= 0.0 ) {
    if ( framePixelsY() ) {
      double wx = framePixelsX();
      double wy = framePixelsY();
      if ( wy > wx ) {
        tickleny = ticklenx;
        ticklenx = (wx/wy)*tickleny;
      } else {
        tickleny = (wy/wx)*ticklenx;
      }
    }
  }
  int nbin = isTH1 ? hist()->GetNbinsX() : 0;
  int flowcol = kAzure - 9;
  // Build over/underflow histogram.
  delete m_flowHist;
  m_flowHist = nullptr;
  if ( (m_showUnderflow || m_showOverflow) && nbin > 0 ) {
    if ( m_flowHist == nullptr ) {
      m_flowHist = dynamic_cast<TH1*>(hist()->Clone("hmaniptmp"));
      m_flowHist->SetDirectory(nullptr);
      m_flowHist->SetStats(0);
      m_flowHist->SetLineColor(flowcol);
      m_flowHist->SetFillColor(flowcol);
      m_flowHist->SetLineWidth(1);
    }
    m_flowHist->Reset();
    if ( m_showUnderflow ) {
      if ( haveHist() ) {
        double binWidth = hist()->GetBinWidth(1);
        int binDisp1 = hist()->FindBin(xmin()+0.499*binWidth);  // First displayed bin.
        int binUnder2 = binDisp1 ? binDisp1 - 1 : 0;    // Last bin below display.
        double yunder = hist()->Integral(0, binUnder2);
        m_flowHist->SetBinContent(binDisp1, yunder);
      }
    }
    if ( m_showOverflow ) {
      if ( haveHist() ) {
        int binOver2 = hist()->GetNbinsX() + 1;  // Last bin above display
        double binWidth = hist()->GetBinWidth(1);
        int binOver1 = hist()->FindBin(xmax()+0.501*binWidth);  // First bin above display.
        int binDisp2 = binOver1 - 1;                        // Last displayed bin.
        double yover = hist()->Integral(binOver1, binOver2);
        m_flowHist->SetBinContent(binDisp2, yover);
      }
    }
  }
  // Make frame so we have clean axis to overlay.
  double xa1 = xmin();
  double xa2 = xmax();
  bool doLogx = m_logX;
  if ( doLogx && xa1 <= 0.0 ) {
    cout << myname << "WARNING: Must have range above 0.0 for logarithmic x-axis." << endl;
    cout << myname << "WARNING: Use setLogRangeX(x1, x2) to set range." << endl;
    cout << myname << "WARNING: Plot is displayed with linear axis." << endl;
    doLogx = false;
  }
  double ya1 = ymin();
  double ya2 = ymax();
  bool doLogy = m_logY;
  if ( doLogy && ya1 <= 0.0 ) {
    cout << myname << "WARNING: Must have range above 0.0 for logarithmic y-axis." << endl;
    cout << myname << "WARNING: Use setLogRangeY(y1, y2) to set range." << endl;
    cout << myname << "WARNING: Plot is displayed with linear axis." << endl;
    doLogy = false;
  }
  double za1 = zmin();
  double za2 = zmax();
  bool doLogz = m_logZ;
  if ( doLogz && za1 <= 0.0 ) {
    cout << myname << "WARNING: Must have range above 0.0 for logarithmic z-axis." << endl;
    cout << myname << "WARNING: Use setLogRangeY(z1, z2) to set range." << endl;
    cout << myname << "WARNING: Plot is displayed with linear axis." << endl;
    doLogy = false;
  }
  string sattl = ";";
  sattl += getXaxis()->GetTitle();
  sattl += ";";
  sattl += getYaxis()->GetTitle();
  sattl += ";";
  // Redraw everything.
  m_ppad->Clear();
  // Set the TPad attributes.
  m_ppad->SetFillColor(m_fillColor);
  m_ppad->SetFrameFillColor(m_frameFillColor);
  m_ppad->SetGridx(m_gridX);
  m_ppad->SetGridy(m_gridY);
  m_ppad->SetLogx(doLogx);
  m_ppad->SetLogy(doLogy);
  m_ppad->SetLogz(doLogz);
  // Draw frame and set axis parameters.
  m_ppad->DrawFrame(xa1, ya1, xa2, ya2, sattl.c_str());
  double labSizeX = getLabelSizeX();
  if ( labSizeX > 0.0 ) {
    getXaxis()->SetLabelSize(labSizeX);
    getXaxis()->SetTitleSize(labSizeX);
    //xttl *= labSizeX/0.035;
  }
  double labSizeY = getLabelSizeY();
  if ( labSizeY > 0.0 ) {
    getYaxis()->SetLabelSize(labSizeY);
    getYaxis()->SetTitleSize(labSizeY);
    //yttl *= labSizeY/0.035;
    if ( getZaxis() != nullptr ) {
      getZaxis()->SetLabelSize(labSizeY);
      getZaxis()->SetTitleSize(labSizeY);
    }
  }
  getXaxis()->SetLabelOffset(xlb);
  getXaxis()->SetTitleOffset(xttl);
  if ( getZaxis() != nullptr ) getZaxis()->SetTitleOffset(zttl);
  getYaxis()->SetTitleOffset(yttl);
  getXaxis()->SetTickLength(ticklenx);
  getYaxis()->SetTickLength(tickleny);
  if ( m_axisTitleOpt ) {
    for ( TAxis* pax : { getXaxis(), getYaxis(), getZaxis() } ) {
      if ( pax == nullptr ) continue;
      pax->CenterTitle();
    }
  }
  if ( m_ndivX ) getXaxis()->SetNdivisions(m_ndivX);
  if ( m_ndivY ) getYaxis()->SetNdivisions(m_ndivY);
  if ( m_timeFormatX.size() ) {
    getXaxis()->SetTimeDisplay(1);
    getXaxis()->SetTimeFormat(m_timeFormatX.c_str());
    getXaxis()->SetTimeOffset(m_timeOffset, "gmt");
  } else {
    getXaxis()->SetTimeDisplay(0);
  }
  if ( m_timeFormatY.size() ) {
    getYaxis()->SetTimeDisplay(1);
    getYaxis()->SetTimeFormat(m_timeFormatY.c_str());
    getYaxis()->SetTimeOffset(m_timeOffset, "gmt");
  } else {
    getYaxis()->SetTimeDisplay(0);
  }
  // May 2019. Ensure frame axis has same binning as histogram.
  // And set bin labels.
  // July 2019: The calls to TAxis::Set cause problems if the drawing bounds are
  // set here. Make those calls iff bin labels are set.
  if ( haveHist() ) {
    if ( m_binLabelsX.size() ) {
      TAxis* pah = hist()->GetXaxis();
      getXaxis()->Set(pah->GetNbins(), pah->GetXmin(), pah->GetXmax());
      for ( Index ilab=0; ilab<m_binLabelsX.size(); ++ilab ) {
        getXaxis()->SetBinLabel(ilab+1, m_binLabelsX[ilab].c_str());
      }
    }
    if ( dynamic_cast<TH2*>(hist()) != nullptr && m_binLabelsY.size() ) {
      TAxis* pah = hist()->GetYaxis();
      getYaxis()->Set(pah->GetNbins(), pah->GetXmin(), pah->GetXmax());
      for ( Index ilab=0; ilab<m_binLabelsY.size(); ++ilab ) {
        getYaxis()->SetBinLabel(ilab+1, m_binLabelsY[ilab].c_str());
      }
    }
  }
  // Primary object.
  if ( haveHist() ) {
    if ( m_flowHist != nullptr ) {
      m_flowHist->Draw("same");
    }
    TAxis* paz = hist()->GetZaxis();
    if ( paz != nullptr && za2 > za1 ) {
      paz->SetRangeUser(za1, za2);
      paz->SetLabelOffset(xlz);
    }
    string dopt = m_dopt.size() ? m_dopt + " same" : "same";
    hist()->Draw(dopt.c_str());
    drawHistFuns();
  } else {
    m_pg->Draw(m_dopt.c_str());
  }
  // Secondary objects.
  for ( Index iobj=0; iobj<m_objs.size(); ++iobj ) {
    TObject* pobj = m_objs[iobj];
    string sopt = m_opts[iobj];
    if ( pobj != nullptr ) pobj->Draw(sopt.c_str());
  }
  // no longer needed?  m_ppad->Update();   // Need an update here to get correct results for xmin, ymin, xmax, ymax
  // Make overflow graph.
  if ( m_gflowOpt.size() ) {
    bool doBot = m_gflowOpt.find("B") != string::npos;
    bool doTop = m_gflowOpt.find("T") != string::npos;
    bool doLef = m_gflowOpt.find("L") != string::npos;
    bool doRig = m_gflowOpt.find("R") != string::npos;
    TGraph* pgout = new TGraph();
    pgout->SetMarkerStyle(m_gflowMrk);
    pgout->SetMarkerColor(m_gflowCol);
    std::vector<TGraph*> gras;
    if ( haveGraph() ) gras.push_back(graph());
    for ( TObject* pobj : m_objs ) {
      TGraph* pgra = dynamic_cast<TGraph*>(pobj);
      if ( pgra != nullptr ) gras.push_back(pgra);
    }
    for ( TGraph* pgin : gras ) {
      for ( int ipt=0; ipt<pgin->GetN(); ++ipt ) {
        double x, y;
        pgin->GetPoint(ipt, x, y);
        double xout = x;
        double yout = y;
        bool doAdd = false;
        if ( doRig && x > xmax() ) { doAdd = true; xout = xmax(); }
        if ( doLef && x < xmin() ) { doAdd = true; xout = xmin(); }
        if ( doBot && y < ymin() ) { doAdd = true; yout = ymin(); }
        if ( doTop && y > ymax() ) { doAdd = true; yout = ymax(); }
        if ( doAdd ) pgout->SetPoint(pgout->GetN(), xout, yout);
      }
    }
    if ( pgout->GetN() == 0 ) {
      delete pgout;
      pgout = nullptr;
    } else {
      delete m_flowGraph;
      m_flowGraph = pgout;
      m_flowGraph->Draw("P");
    }
  }
  drawLines();
  if ( m_top ) drawAxisTop();
  if ( m_right ) drawAxisRight();
  pad()->RedrawAxis();  // In case they are covered
  pad()->RedrawAxis("G");  // In case they are covered
  // Add the title and labels.
  // First fix the title position and size (April 2020).
  // Title is positioned halfway up the top margin with text height no more
  // than 90% of the margin.
  float marg = m_ppad->GetTopMargin();
  float ttlY = 1.0 - 0.5*marg;
  m_title.SetY(ttlY);
  float ttlSizeMax = 0.90*marg;
  float ttlSize = m_ttlSize > 0 ? m_ttlSize : 0.035;
  if ( ttlSize > ttlSizeMax ) ttlSize = ttlSizeMax;
  m_title.SetTextSize(ttlSize);
  m_title.Draw();
  if ( getLabel().size() ) m_label.Draw();
  gPad = pPadSave;
  return rstat;
}

//**********************************************************************

TAxis* TPadManipulator::getXaxis() const {
  if ( haveFrameHist() ) return frameHist()->GetXaxis();
  if ( hist() != nullptr ) return hist()->GetXaxis();
  if ( m_pg != nullptr ) return m_pg->GetXaxis();
  return nullptr;
}

//**********************************************************************

TAxis* TPadManipulator::getYaxis() const {
  if ( haveFrameHist() ) return frameHist()->GetYaxis();
  if ( hist() != nullptr ) return hist()->GetYaxis();
  if ( m_pg != nullptr ) return m_pg->GetYaxis();
  return nullptr;
}

//**********************************************************************

TAxis* TPadManipulator::getZaxis() const {
  //if ( haveFrameHist() ) return frameHist()->GetZaxis();
  if ( hist() != nullptr ) return hist()->GetZaxis();
  return nullptr;
}

//**********************************************************************

int TPadManipulator::setCanvasSize(int wx, int wy) {
  if ( haveParent() ) return parent()->setCanvasSize(wx, wy);
  m_canvasWidth = wx;
  m_canvasHeight = wy;
  if ( haveCanvas() ) {
    canvas()->SetWindowSize(wx, wy);
    return update();
  }
  return 0;
}

//**********************************************************************

double TPadManipulator::getLabelSizeX() const {
  if ( m_labSizeX ) return m_labSizeX;
  if ( haveParent() ) {
    const TPadManipulator* pman = parent();
    double h1 = padPixelsY();
    double h2 = pman->padPixelsY();
    double h = h1 ? h2/h1*pman->getLabelSizeX() : 0.0;
    return h;
  }
  return 0.0;
}

//**********************************************************************

double TPadManipulator::getLabelSizeY() const {
  if ( m_labSizeY ) return m_labSizeY;
  if ( haveParent() ) {
    const TPadManipulator* pman = parent();
    double h1 = padPixelsY();
    double h2 = pman->padPixelsY();
    double h = h1 ? h2/h1*pman->getLabelSizeY() : 0.0;
    return h;
  }
  return 0.0;
}

//**********************************************************************

double TPadManipulator::getTitleSize() const {
  if ( m_ttlSize ) return m_ttlSize;
  if ( haveParent() ) {
    const TPadManipulator* pman = parent();
    double h1 = padPixelsY();
    double h2 = pman->padPixelsY();
    double h = h1 ? h2/h1*pman->getTitleSize() : 0.0;
    return h;
  }
  return 0.0;
}

//**********************************************************************

int TPadManipulator::setRangeX(double x1, double x2) {
  m_setBounds.x1 = x1;
  m_setBounds.x2 = x2;
  return update();
}

//**********************************************************************

int TPadManipulator::setRangeY(double y1, double y2) {
  m_setBounds.y1 = y1;
  m_setBounds.y2 = y2;
  return update();
}

//**********************************************************************

int TPadManipulator::setRangeZ(double z1, double z2) {
  m_setBounds.z1 = z1;
  m_setBounds.z2 = z2;
  return update();
}

//**********************************************************************

int TPadManipulator::setRangeXY(double x1, double x2, double y1, double y2) {
  m_setBounds.x1 = x1;
  m_setBounds.x2 = x2;
  m_setBounds.y1 = y1;
  m_setBounds.y2 = y2;
  return update();
}

//**********************************************************************

int TPadManipulator::setLogRangeX(double x1, double x2) {
  m_setBoundsLog.x1 = x1;
  m_setBoundsLog.x2 = x2;
  return update();
}

//**********************************************************************

int TPadManipulator::setLogRangeY(double y1, double y2) {
  m_setBoundsLog.y1 = y1;
  m_setBoundsLog.y2 = y2;
  return update();
}

//**********************************************************************

int TPadManipulator::setLogRangeZ(double z1, double z2) {
  m_setBoundsLog.z1 = z1;
  m_setBoundsLog.z2 = z2;
  return update();
}

//**********************************************************************

int TPadManipulator::setTimeOffset(double toff) {
  m_timeOffset = toff;
  return 0;
}

//**********************************************************************

int TPadManipulator::setTimeFormatX(string sfmt) {
  m_timeFormatX = sfmt;
  return 0;
}

//**********************************************************************

int TPadManipulator::setTimeFormatY(string sfmt) {
  m_timeFormatY = sfmt;
  return 0;
}

//**********************************************************************

int TPadManipulator::addAxis(bool flag) {
  return addAxisTop(flag) + addAxisRight(flag);
}

//**********************************************************************

int TPadManipulator::addAxisTop(bool flag) {
  m_top = flag;
  return update();
}

//**********************************************************************

int TPadManipulator::addAxisRight(bool flag) {
  m_right = flag;
  return update();
}

//**********************************************************************

int TPadManipulator::addHistFun(unsigned int ifun) {
  m_histFuns.push_back(ifun);
  return update();
}

//**********************************************************************

int TPadManipulator::showUnderflow(bool show) {
  if ( show == m_showUnderflow ) return 0;
  m_showUnderflow = show;
  return update();
}

//**********************************************************************

int TPadManipulator::showOverflow(bool show) {
  if ( show == m_showOverflow ) return 0;
  m_showOverflow = show;
  return update();
}

//**********************************************************************

int TPadManipulator::showGraphOverflow(std::string sopt, int imrk, int icol) {
  delete m_flowGraph;
  m_flowGraph = nullptr;
  m_gflowOpt = sopt;
  m_gflowMrk = imrk;
  m_gflowCol = icol;
  return 0;
}

//**********************************************************************

int TPadManipulator::clearLines() {
  m_vmlXmod.clear();
  m_vmlXoff.clear();
  m_vmlXStyle.clear();
  m_vmlXLength.clear();
  m_hmlXmod.clear();
  m_hmlXoff.clear();
  m_hmlXStyle.clear();
  m_hmlXLength.clear();
  return update();
}

//**********************************************************************

int TPadManipulator::addVerticalLine(double xoff, double lenfrac, int isty) {
  m_vmlXmod.push_back(0.0);
  m_vmlXoff.push_back(xoff);
  m_vmlXStyle.push_back(isty);
  m_vmlXLength.push_back(lenfrac);
  return update();
}
//**********************************************************************

int TPadManipulator::addHorizontalLine(double xoff, double lenfrac, int isty) {
  m_hmlXmod.push_back(0.0);
  m_hmlXoff.push_back(xoff);
  m_hmlXStyle.push_back(isty);
  m_hmlXLength.push_back(lenfrac);
  drawLines();
  return 0;
}
//**********************************************************************

int TPadManipulator::addSlopedLine(double slop, double yoff, int isty) {
  m_slSlop.push_back(slop);
  m_slYoff.push_back(yoff);
  m_slStyl.push_back(isty);
  drawLines();
  return 0;
}
//**********************************************************************

int TPadManipulator::addVerticalModLines(double xmod, double xoff, double lenfrac, int isty) {
  m_vmlXmod.push_back(xmod);
  m_vmlXoff.push_back(xoff);
  m_vmlXStyle.push_back(isty);
  m_vmlXLength.push_back(lenfrac);
  drawLines();
  return 0;
}

//**********************************************************************

int TPadManipulator::addHorizontalModLines(double xmod, double xoff, double lenfrac, int isty) {
  m_hmlXmod.push_back(xmod);
  m_hmlXoff.push_back(xoff);
  m_hmlXStyle.push_back(isty);
  m_hmlXLength.push_back(lenfrac);
  drawLines();
  return 0;
}

//**********************************************************************

int TPadManipulator::setBinLabelsX(const NameVector& labs) {
  m_binLabelsX = labs;
  return 0;
}

//**********************************************************************

int TPadManipulator::setBinLabelsY(const NameVector& labs) {
  m_binLabelsY = labs;
  return 0;
}

//**********************************************************************

int TPadManipulator::fixFrameFillColor() {
  if ( gStyle == nullptr ) return 2;
  m_ppad->SetFrameFillColor(gStyle->GetColorPalette(0));
  return 0;
}

//**********************************************************************

int TPadManipulator::draw() {
  if ( haveParent() ) return parent()->draw();
  if ( m_ppad == nullptr ) {
    if ( ! haveHistOrGraph() && npad() == 0 ) return 1;
    TCanvas* pcan = new TCanvas;
    if ( m_canvasWidth > 0 && m_canvasHeight > 0 ) {
      string snam = pcan->GetName();
      string sttl = pcan->GetTitle();
      delete pcan;
      pcan = new TCanvas(snam.c_str(), sttl.c_str(), m_canvasWidth, m_canvasHeight);
      //pcan->SetCanvasSize(m_canvasWidth, m_canvasHeight);
      //pcan->SetWindowSize(m_canvasWidth, m_canvasHeight);
    }
    m_ppad = pcan;
  }
  return update();
}

//**********************************************************************

int TPadManipulator::erase() {
  if ( haveParent() ) return parent()->erase();
  if ( m_ppad == nullptr ) return 0;
  delete m_ppad;
  m_ppad = nullptr;
  return 0;
}

//**********************************************************************

int TPadManipulator::redraw() {
  if ( haveParent() ) return parent()->redraw();
  erase();
  return draw();
}

//**********************************************************************

int TPadManipulator::drawAxisTop() {
  if ( ! m_top ) return 0;
  if ( m_ppad == nullptr ) return 1;
  TVirtualPad* pPadSave = gPad;
  m_ppad->cd();
  double xminPad, yminPad, xmaxPad, ymaxPad;
  m_ppad->GetRangeAxis(xminPad, yminPad, xmaxPad, ymaxPad);
  string sopt = "-US";
  double x1 = xmin();
  double x2 = xmax();
  if ( m_ppad->GetLogx() ) {
    xminPad = pow(10.0, xminPad);
    xmaxPad = pow(10.0, xmaxPad);
    x1 = pow(10.0, x1);
    x2 = pow(10.0, x2);
    sopt += "G";
  }
  if ( m_ppad->GetLogy() ) {
    yminPad = pow(10.0, yminPad);
    ymaxPad = pow(10.0, ymaxPad);
  }
  TAxis* paxold = getXaxis();
  if ( paxold == 0 ) return 2;
  double ticksize = paxold->GetTickLength();
  int ndiv = paxold->GetNdivisions();
  TGaxis* paxnew = new TGaxis(xminPad, ymaxPad, xmaxPad, ymaxPad,
                              x1, x2, ndiv, sopt.c_str());
  if ( ticksize > 0 ) paxnew->SetTickLength(ticksize);
  string name = "TopAxis";
  paxnew->SetName(name.c_str());
  TList* pobjs = m_ppad->GetListOfPrimitives();
  for ( int iobj=pobjs->GetEntries()-1; iobj>=0; --iobj ) {
    TGaxis* paxold = dynamic_cast<TGaxis*>(pobjs->At(iobj));
    if ( paxold == nullptr ) continue;
    if ( paxold->GetName() == name ) {
      pobjs->RemoveAt(iobj);
      break;
    }
  }
  if ( m_timeFormatX.size() ) {
    //paxnew->SetTimeDisplay(1);
    paxnew->SetTimeFormat(m_timeFormatX.c_str());
    paxnew->SetTimeOffset(m_timeOffset, "gmt");
  } else {
    //paxnew->SetTimeDisplay(0);
  }
  paxnew->Draw();
  gPad = pPadSave;
  return 0;
}

//**********************************************************************

int TPadManipulator::drawAxisRight() {
  if ( ! m_right ) return 0;
  if ( m_ppad == nullptr ) return 1;
  TVirtualPad* pPadSave = gPad;
  m_ppad->cd();
  double xminPad, yminPad, xmaxPad, ymaxPad;
  m_ppad->GetRangeAxis(xminPad, yminPad, xmaxPad, ymaxPad);
  double y1 = ymin();
  double y2 = ymax();
  string sopt = "+US";
  if ( m_ppad->GetLogx() ) {
    xminPad = pow(10.0, xminPad);
    xmaxPad = pow(10.0, xmaxPad);
  }
  if ( m_ppad->GetLogy() ) {
    yminPad = pow(10.0, yminPad);
    ymaxPad = pow(10.0, ymaxPad);
    sopt += "G";
  }
  TAxis* paxold = getYaxis();
  if ( paxold == 0 ) return 2;
  double ticksize = paxold->GetTickLength();
  int ndiv = paxold->GetNdivisions();
  TGaxis* paxnew = new TGaxis(xmaxPad, yminPad, xmaxPad, ymaxPad,
                              y1, y2, ndiv, sopt.c_str());
  if ( ticksize > 0 ) paxnew->SetTickLength(ticksize);
  string name = "RightAxis";
  paxnew->SetName(name.c_str());
  TList* pobjs = m_ppad->GetListOfPrimitives();
  for ( int iobj=pobjs->GetEntries()-1; iobj>=0; --iobj ) {
    TGaxis* paxold = dynamic_cast<TGaxis*>(pobjs->At(iobj));
    if ( paxold == nullptr ) continue;
    if ( paxold->GetName() == name ) {
      pobjs->RemoveAt(iobj);
      break;
    }
  }
  paxnew->Draw("");
  gPad = pPadSave;
  return 0;
}

//**********************************************************************

int TPadManipulator::drawHistFuns() {
  if ( m_ppad == nullptr ) return 1;
  const TList* pfuns = nullptr;
  if      ( haveHist() ) pfuns = hist()->GetListOfFunctions();
  else if ( m_pg != nullptr ) pfuns = m_pg->GetListOfFunctions();
  else    return 1;
  const TList& funs = *pfuns;
  unsigned int nfun = funs.GetEntries();
  TVirtualPad* pPadSave = gPad;
  m_ppad->cd();
  for ( unsigned int ifun : m_histFuns ) {
    if ( ifun >= nfun ) continue;
    TF1* pfun = dynamic_cast<TF1*>(funs.At(ifun));
    if ( pfun != nullptr ) pfun->Draw("same");
  }
  gPad = pPadSave;
  return 0;
}

//**********************************************************************

int TPadManipulator::drawLines() {
  if ( m_ppad == nullptr ) return 1;
  TVirtualPad* pPadSave = gPad;
  m_ppad->cd();
  clearLineObjects();
  for ( Index iset=0; iset<m_vmlXmod.size(); ++iset ) {
    double xmod = m_vmlXmod[iset];
    double xoff = m_vmlXoff[iset];
    double xlen = m_vmlXLength[iset];
    int isty = m_vmlXStyle[iset];
    if ( xmod < 0.0 ) continue;
    double x = xoff;
    if ( xmod > 0.0 ) {
      while ( x >= xmin() ) x -= xmod;
      while ( x < xmin() ) x += xmod;
    }
    double ytop = ymax() + (xlen - 1.0)*(ymax() - ymin());
    while ( x <= xmax() ) {
      TLine* pline = new TLine(x, ymin(), x, ytop);
      pline->SetLineStyle(isty);
      m_lines.push_back(pline);
      pline->Draw();
      if ( xmod == 0.0 ) break;
      x += xmod;
    }
  }
  for ( Index iset=0; iset<m_hmlXmod.size(); ++iset ) {
    double ymod = m_hmlXmod[iset];
    double yoff = m_hmlXoff[iset];
    double ylen = m_hmlXLength[iset];
    int isty = m_hmlXStyle[iset];
    if ( ymod < 0.0 ) continue;
    double y = yoff;
    if ( ymod > 0.0 ) {
      while ( y >= ymin() ) y -= ymod;
      while ( y < ymin() ) y += ymod;
    }
    double xtop = xmax() + (ylen - 1.0)*(xmax() - xmin());
    while ( y <= ymax() ) {
      TLine* pline = new TLine(xmin(), y, xtop, y);
      pline->SetLineStyle(isty);
      m_lines.push_back(pline);
      pline->Draw();
      if ( ymod == 0.0 ) break;
      y += ymod;
    }
  }
  for ( Index iset=0; iset<m_slSlop.size(); ++iset ) {
    double slop = m_slSlop[iset];
    double yoff = m_slYoff[iset];
    double isty = m_slStyl[iset];
    double wx = xmax() - xmin();
    double wy = ymax() - ymin();
    if ( wx < 0.0 || wy < 0.0 ) continue;
    double xlef = xmin();
    double xrig = xmax();
    double ybot = ymin();
    double ytop = ymax();
    double ylef = slop*xlef + yoff;
    double yrig = slop*xrig + yoff;
    double xbot = slop == 0.0  ? xmin() - wx : (ybot - yoff)/slop;
    double xtop = slop == 0.0  ? xmax() + wx : (ytop - yoff)/slop;
    double x1 = 0.0;
    double y1 = 0.0;
    if ( ylef > ybot && ylef < ytop ) {
      x1 = xlef;
      y1 = ylef;
    } else if ( ylef <= ybot && slop > 0.0 && xbot < xrig ) {
      x1 = xbot;
      y1 = ybot;
    } else if ( ylef >= ytop && slop < 0.0 && xtop < xrig ) {
      x1 = xtop;
      y1 = ytop;
    } else {
      continue;
    }
    double x2 = 0.0;
    double y2 = 0.0;
    if ( slop > 0.0 && xtop < xrig ) {
      x2 = xtop;
      y2 = ytop;
    } else if ( slop < 0.0 && xbot < xrig ) {
      x2 = xbot;
      y2 = ybot;
    } else {
      x2 = xrig;
      y2 = yrig;
    }
    TLine* pline = new TLine(x1, y1, x2, y2);
    pline->SetLineStyle(isty);
    m_lines.push_back(pline);
    pline->Draw();
  }
  gPad = pPadSave;
  return 0;
}

//**********************************************************************

void TPadManipulator::Streamer(TBuffer& buf) {
  const string myname = "TPadManipulator::Streamer: ";
  int dbg = 0;
  TClass* pclass = TClass::GetClass("TPadManipulator");
  if ( buf.IsReading() ) {
    if ( pclass == nullptr ) {
      cout << myname << "Dictionary not found for read." << endl;
      return;
    }
    // Stream in.
    pclass->ReadBuffer(buf, this);
    // Set the parents. This only needs to be done for the top level object
    // but I don't know how to tell if this is the case.
    setParents(true);
    // Log message.
    string msg = haveHist() ? "histogram" : haveGraph() ? "graph" : "no primary";
    cout << myname << "Read pad with " << msg << " and " 
         << objects().size() << " extra objects and "
         << m_subMans.size() << " subpads."  << endl;
  } else {
    if ( pclass == nullptr ) {
      cout << myname << "Dictionary not found for write." << endl;
      return;
    }
    // Log message.
    if ( dbg ) {
      string msg = haveHist() ? "histogram" : haveGraph() ? "graph" : "no primary";
      cout << myname << "Writing pad with " << msg << " and "
           << objects().size() << " extra objects and "
           << m_subMans.size() << " subpads." << endl;
      if ( haveHistOrGraph() ) object()->Print();
    }
    // Stream out.
    pclass->WriteBuffer(buf, this);
  }
}

//**********************************************************************

void TPadManipulator::setParents(bool recurse) {
  for ( TPadManipulator& man : m_subMans ) {
    man.m_parent = this;
    if ( recurse ) man.setParents(true);
  }
 }

//**********************************************************************
