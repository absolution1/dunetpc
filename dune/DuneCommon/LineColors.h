// LineColors.h
//
// David Adams
// February 2018
//
// Standard Root colors for lines and text, e.g. for drawing histograms.
//
// Developers may inherit from the class and redefine any color to
// their taste.

#include "Rtypes.h"

class LineColors {

public:

  using ColorType = int;
  using Index = unsigned int;

  // Return the colors.
  static ColorType     white() { return kWhite; };
  static ColorType     black() { return kBlack; };
  static ColorType      blue() { return kBlue + 2; }
  static ColorType       red() { return kRed + 1; }
  static ColorType     green() { return kGreen + 2; }
  static ColorType     brown() { return kOrange + 3; }
  static ColorType    violet() { return kMagenta + 1; }
  static ColorType    orange() { return kOrange - 3; }
  static ColorType      cyan() { return kCyan + 1; }
  static ColorType    yellow() { return kYellow - 7; }
  static ColorType      gray() { return kGray + 2; }
  static ColorType    purple() { return kViolet + 2; }
  static ColorType lightGray() { return kGray + 0; }

  // Return the # "good" colors in the table.
  static Index size() { return 8; }

  // Return the color for an index and color count.
  // The colors repeat every ncol values.
  static ColorType color(Index icolin, Index ncol =size()) {
    Index icol = icolin;
    if ( ncol > 0 ) icol = icolin%ncol;
    if ( icol ==  0 ) return   blue();
    if ( icol ==  1 ) return    red();
    if ( icol ==  2 ) return  green();
    if ( icol ==  3 ) return  brown();
    if ( icol ==  4 ) return violet();
    if ( icol ==  5 ) return orange();
    if ( icol ==  6 ) return   cyan();
    if ( icol ==  7 ) return   gray();
    if ( icol ==  8 ) return purple();
    if ( icol ==  9 ) return  black();
    if ( icol == 10 ) return yellow();
    if ( icol == 11 ) return lightGray();
    return white();
  }

  // Return the color for an index using the default color count.
  ColorType operator[](Index icolin) const { return color(icolin); }

};
