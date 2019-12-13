// RootPalette.h
//
// Class to define color palettes for Root.
//
// Example of use:
//   RootPalette oldPalette;       // Save the old palette
//   RootPalette::set(1026);
//   TCanvas* pcan = new TCanvas;
//   ph2->Draw("colz");
//   pcan->Print();
//   oldPalette.setRootPalette();  // Restore the old palette
//
// Palette indices:
//  < 1000 uses the Root definitions (e.g. kRainBow or kBird).
//         See the Root doc page for TColor for some of the permitted values.
//  K+10 - yellow to red to black. For zmin = 0.
//  K+11 - Same as 1010 with first 0.5% of colors set white
//  K+12 - Same as 1010 with first 1% of colors set white
//  K+13 - Same as 1010 with first 2% of colors set white
//  K+14 - Same as 1010 with first 3% of colors set white
//  K+15 - Same as 1010 with first 4% of colors set white
//  K+16 - Same as 1010 with first 5% of colors set white
//  K+17 - Same as 1010 with first 10% of colors set white
//  K+18 - Same as 1010 with first 15% of colors set white
//  K+19 - Same as 1010 with first 20% of colors set white
//  K+20 - dark to light blue, yellow to red to black. For zmin = -zmax.
//  K+21 - Same as 1020 with central 0.5% of colors set white
//  K+22 - Same as 1020 with central 1% of colors set white
//  K+23 - Same as 1020 with central 2% of colors set white
//  K+24 - Same as 1020 with central 3% of colors set white
//  K+25 - Same as 1020 with central 4% of colors set white
//  K+26 - Same as 1020 with central 5% of colors set white
//  K+27 - Same as 1020 with central 10% of colors set white
//  K+28 - Same as 1020 with central 15% of colors set white
//  K+29 - Same as 1020 with central 20% of colors set white
//  K = 1000 - Original palette
//      2000 - Palette added Dec 2019

#ifndef RootPalette_H
#define RootPalette_H

#include <vector>

class TStyle;

class RootPalette {

public:

  // Default Root palette.
  // This is created when the first object of this type is created.
  static const RootPalette* defaultPalette();

  // Return a palette by index.
  //   0 - default Root palette.
  //   1-999 - Root definitions
  //   1001 - ADC blue-red-black
  static const RootPalette* find(unsigned int ipal);

  // Set palette index to |ipal| (see above).
  // If ipal <  0, the color map is inverted.
  // Returns true if the palette was set successfully.
  static bool set(int ipal);

  // Ctor from the current style.
  // Copies the current color map.
  RootPalette();

  // Ctor from a given style.
  RootPalette(TStyle* pstyle);

  // Ctor from a Root palette index.
  RootPalette(int ipal);

  // Getters.
  int rootPaletteIndex() const { return m_rootPalette; }
  unsigned int ncol() const { return m_cols.size(); }
  const int* colorArray() const { return ncol() ? &m_cols[0] : nullptr; }
  const std::vector<int>& colorVector() const { return m_cols; }

  // Fetch the Root palette from the current or specified style.
  int getRootPalette();
  int getRootPalette(TStyle* pstyle);

  // Set the Root palette for the current or specified style.
  // Returns 0 for success.
  int setRootPalette() const;
  int setRootPalette(TStyle* pstyle) const;

  // Return the color map mutable.
  // With this, user may change any of the colors.
  std::vector<int>& mutableColorVector() { return m_cols; }

  // Set the first ncol entries in the color map to color icol.
  void setFirstColors(unsigned int ncol, unsigned int icol =0);

  // Set the ncol entries above and below the center of the color map to color icol.
  void setCentralColors(unsigned int ncol, unsigned int icol =0);

private:

  // If this is nonzero, the palette is one of the Root defaults.
  int m_rootPalette =0;

  // Otherwise, the palette is this color map.
  std::vector<int> m_cols;

};

#endif
