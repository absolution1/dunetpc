// RootPalette.cxx

#include "RootPalette.h"
#include "TStyle.h"

using std::vector;

//**********************************************************************

namespace {

std::vector<const RootPalette*>& fixedPals() {
  static std::vector<const RootPalette*> pals(2000, nullptr);
  return pals;
}

}  // end unnamed namespace

//**********************************************************************
// Static methods.
//**********************************************************************

const RootPalette* RootPalette::defaultPalette() {
  static bool first = true;
  if ( first ) {
    first = false;
    fixedPals()[0] = new RootPalette;
  }
  return fixedPals()[0];
}

//**********************************************************************

const RootPalette* RootPalette::find(unsigned int ipal) {
  if ( ipal >= fixedPals().size() ) return nullptr;
  if ( fixedPals()[ipal] != nullptr ) return fixedPals()[ipal];
  defaultPalette();
  if ( ipal == 0 ) {
    return defaultPalette();
  } else if ( ipal < 1000 ) {
    if ( fixedPals()[ipal] == nullptr ) {
      fixedPals()[ipal] = new RootPalette(ipal);
    }
  } else if ( ipal >= 1010 && ipal < 1020 ) {
    TStyle* pstyleSave = gStyle;
    gStyle = dynamic_cast<TStyle*>(gStyle->Clone("TmpStyle"));
    TStyle* pstyleTmp = gStyle;
    const double alpha = 1.0;
    const int nRGBs = 5;
    const int ncol = 200;
    Double_t stops[nRGBs] = { 0.00, 0.30, 0.44, 0.70, 1.00};
    Double_t red[nRGBs]   = { 1.00, 1.00, 1.00, 0.70, 0.00};
    Double_t green[nRGBs] = { 1.00, 0.75, 0.55, 0.20, 0.00};
    Double_t blue[nRGBs]  = { 1.00, 0.00, 0.00, 0.10, 0.00};
    TColor::CreateGradientColorTable(nRGBs, stops, red, green, blue, ncol, alpha);
    RootPalette* ppal = new RootPalette;
    vector<unsigned int> ncfxs = {0, 1, 2, 4, 6, 8, 10, 20, 30, 40};
    ppal->setFirstColors(ncfxs[ipal - 1010]);
    fixedPals()[ipal] = ppal;
    gStyle = pstyleSave;
    delete pstyleTmp;
  } else if ( ipal >= 1020 && ipal < 1030 ) {
    TStyle* pstyleSave = gStyle;
    gStyle = dynamic_cast<TStyle*>(gStyle->Clone("TmpStyle"));
    TStyle* pstyleTmp = gStyle;
    const double alpha = 1.0;
    const int nRGBs = 7;
    const int ncol = 400;
    Double_t stops[nRGBs] = { 0.00, 0.50, 0.50, 0.65, 0.72, 0.85, 1.00};
    Double_t red[nRGBs]   = { 0.09, 0.90, 1.00, 1.00, 1.00, 0.70, 0.00};
    Double_t green[nRGBs] = { 0.60, 0.95, 1.00, 0.75, 0.55, 0.20, 0.00};
    Double_t blue[nRGBs]  = { 0.48, 1.00, 1.00, 0.00, 0.00, 0.10, 0.00};
    TColor::CreateGradientColorTable(nRGBs, stops, red, green, blue, ncol, alpha);
    RootPalette* ppal = new RootPalette;
    vector<unsigned int> ncfxs = {0, 1, 2, 4, 6, 8, 10, 20, 30, 40};
    ppal->setCentralColors(ncfxs[ipal - 1020]);
    fixedPals()[ipal] = ppal;
    gStyle = pstyleSave;
    delete pstyleTmp;
  }
  return fixedPals()[ipal];
}

//**********************************************************************

bool RootPalette::set(int ipalin) {
  unsigned int ipal = ipalin >= 0 ? ipalin : -ipalin;
  defaultPalette();
  if ( ipal >= fixedPals().size() ) return false;
  // Use predfined map if it exists.
  const RootPalette* ppal = RootPalette::find(ipal);
  // Otherwise, treat this as a Root palette index.
  if ( ppal == nullptr && ipal > 0 && ipal < 1000 ) {
    ppal = new RootPalette(ipal);
    fixedPals()[ipal] = ppal;
  }
  if ( ppal != nullptr && ppal->setRootPalette() == 0 ) {
    if ( ipalin < 0 ) TColor::InvertPalette();
    return true;
  }
  return false;
}

//**********************************************************************
// Non-static methods.
//**********************************************************************

RootPalette::RootPalette() {
  defaultPalette();
  getRootPalette();
}

//**********************************************************************

RootPalette::RootPalette(TStyle* pstyle) {
  defaultPalette();
  getRootPalette(pstyle);
}

//**********************************************************************

RootPalette::RootPalette(int ipal) {
  defaultPalette();
  if ( ipal > 0 && ipal < 1000 ) m_rootPalette = ipal;
}

//**********************************************************************

int RootPalette::getRootPalette() {
  return getRootPalette(gStyle);
}

//**********************************************************************

int RootPalette::getRootPalette(TStyle* pstyle) {
  m_cols.clear();
  if ( pstyle == nullptr ) return 2;
  unsigned int nscol = pstyle->GetNumberOfColors();
  if ( nscol == 0 ) return 1;
  for ( unsigned int icol=0; icol<nscol; ++icol ) m_cols.push_back(pstyle->GetColorPalette(icol));
  return 0;
}

//**********************************************************************

int RootPalette::setRootPalette() const {
  return setRootPalette(gStyle);
}

//**********************************************************************

int RootPalette::setRootPalette(TStyle* pstyle) const {
  if ( pstyle == nullptr ) return 2;
  int ipal = rootPaletteIndex();
  if ( ipal >0 && ipal < 1000 ) {
    pstyle->SetPalette(ipal);
  } else if ( ncol() == 0 ) {
    return 1;
  } else {
    vector<int> cols(ncol());
    for ( unsigned int icol=0; icol<ncol(); ++icol ) cols[icol] = m_cols[icol];
    pstyle->SetPalette(ncol(), &cols[0]);
  }
  return 0;
}

//**********************************************************************

void RootPalette::setFirstColors(unsigned int ncfx, unsigned int icol) {
  for ( unsigned int icfx=0; icfx<ncfx; ++icfx ) mutableColorVector()[icfx] = icol;
}

//**********************************************************************

void RootPalette::setCentralColors(unsigned int ncfx, unsigned int icol) {
  unsigned int icz = ncol()/2;
  unsigned int icfx1 = icz > ncfx ? icz - ncfx : 0;
  unsigned int icfx2 = icz + ncfx;
  if ( 2*icz + ncol() ) icfx2 -= 1;  // Handle odd # colors.
  if ( icfx2 > ncol() ) icfx2 = ncol();
  for ( unsigned int icfx=icfx1; icfx<icz+ncfx; ++icfx ) mutableColorVector()[icfx] = icol;
}

//**********************************************************************

