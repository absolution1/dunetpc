// test_AdcChannelData.cxx
//
// David Adams
// January 2020
//
// Test AdcChannelData.

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include "dune/DuneInterface/Data/AdcChannelData.h"

#undef NDEBUG
#include <cassert>

using std::string;
using std::cout;
using std::endl;
using std::vector;

using Index = AdcIndex;
using View = AdcChannelData::View;

//**********************************************************************

int test_AdcChannelData() {
  const string myname = "test_AdcChannelData: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = myname + "-----------------------------";

  AdcChannelData acdtop;
  cout << line << endl;
  cout << "Checking fill." << endl;
  acdtop.channel = 12345;
  acdtop.fembID = 101;
  for ( int i=0; i<20; ++i ) acdtop.raw.push_back(i);
  assert( acdtop.channel == 12345 );
  assert( acdtop.fembID == 101 );
  assert( acdtop.raw.size() == 20 );

  cout << "Checking copy." << endl;
  AdcChannelData acdtopCopy(acdtop);
  assert( acdtopCopy.channel == 12345 );
  assert( acdtopCopy.fembID == 101 );
  assert( acdtopCopy.raw.size() == 0 );

#if 0
  cout << "Checking assignment." << endl;
  AdcChannelData acdtopAssg;
  acdtopAssg = acdtop;
  assert( acdtopAssg.channel == 12345 );
  assert( acdtopAssg.fembID == 101 );
  assert( acdtopAssg.raw.size() == 20 );
#endif

  cout << line << endl;
  cout << myname << "Add view1." << endl;
  assert( acdtop.viewSize() == 0 );
  View& view1 = acdtop.updateView("view1");
  assert ( acdtop.hasView("view1") );
  assert ( ! acdtop.hasView("view1/view2") );
  assert( acdtop.viewSize() == 1 );
  Index ient1 = 0;
  for ( ient1=0; ient1<4; ++ient1 ) view1.push_back(acdtop);
  ient1 = 0;
  for ( AdcChannelData& acd1 : view1 ) {
    assert( acd1.channel == acdtop.channel);
    acd1.setMetadata("viewIndex1", ient1);
    ++ient1;
  }

  cout << line << endl;
  cout << myname << "Check view1." << endl;
  Index nent1 = acdtop.viewSize("view1");
  assert( nent1 == 4 );
  for ( ient1=0; ient1<nent1; ++ient1 ) {
    const AdcChannelData* pacd1 = acdtop.viewEntry("view1", ient1);
    assert( pacd1 != nullptr );
    const AdcChannelData& acd1 = *pacd1;
    assert( acd1.channel == acdtop.channel);
    assert( acd1.hasMetadata("viewIndex1") );
    assert( acd1.getMetadata("viewIndex1", 99) == ient1 );
    AdcChannelData* pacdm = acdtop.mutableViewEntry("view1", ient1);
    assert( pacdm != nullptr );
    const AdcChannelData& acdm = *pacdm;
    assert( acdm.channel == acdtop.channel);
    assert( acdm.hasMetadata("viewIndex1") );
    assert( acdm.getMetadata("viewIndex1", 99) == ient1 );
  }

  cout << line << endl;
  cout << myname << "Add view2." << endl;
  Index ient2 = 0;
  vector<Index> nent2s = { 1, 2, 2, 3 };
  vector<Index> ient1s;
  vector<Index> ient2s;
  Index nent12 = 0;
  Index ient12 = 0;
  for ( Index ient1=0; ient1<nent1; ++ient1 ) {
    AdcChannelData& acd1 = view1[ient1];
    View& view2 = acd1.updateView("view2");
    Index nent2 = nent2s[ient1];
    nent12 += nent2;
    for ( ient2=0; ient2<nent2; ++ient2 ) {
      view2.push_back(acdtop);
      AdcChannelData& acd2 = view2.back();
      acd2.setMetadata("viewIndex1", ient1);
      acd2.setMetadata("viewIndex2", ient2);
      acd2.setMetadata("viewIndex12", ient12);
      ient1s.push_back(ient1);
      ient2s.push_back(ient2);
      ++ient12;
    }
    assert ( acd1.hasView("view2") );
    assert( acd1.viewSize() == 1 );
    assert( ient12 = nent12 );
  }
  assert ( acdtop.hasView("view1") );
  assert ( acdtop.hasView("view1/view2") );
  assert( ient12 = nent12 );
    
  cout << line << endl;
  cout << myname << "Check view2." << endl;
  ient12 = 0;
  for ( ient1=0; ient1<nent1; ++ient1 ) {
    const AdcChannelData* pacd1 = acdtop.viewEntry("view1", ient1);
    assert( pacd1 != nullptr );
    const AdcChannelData& acd1 = *pacd1;
    Index nent2 = acd1.viewSize("view2");
    for ( ient2=0; ient2<nent2; ++ient2 ) {
      cout << myname << " " << ient1 << " " << ient2 << " " << ient12 << endl;
      const AdcChannelData* pacd2 = acd1.viewEntry("view2", ient2);
      assert( pacd2 != nullptr );
      const AdcChannelData& acd2 = *pacd2;
      assert( acd2.channel == 12345 );
      assert( acd2.fembID == 101 );
      assert( acd2.getMetadata("viewIndex1") == ient1 );
      assert( acd2.getMetadata("viewIndex2") == ient2 );
      assert( acd2.getMetadata("viewIndex12") == ient12 );
      ++ient12;
    }
  }

  cout << line << endl;
  cout << myname << "Check view12." << endl;
  assert( acdtop.viewSize("view1/view2") == nent12 );
  for ( ient12=0; ient12<nent12; ++ient12 ) {
    const AdcChannelData* pacd12 = acdtop.viewEntry("view1/view2", ient12);
    assert( pacd12 != nullptr );
    const AdcChannelData& acd12 = *pacd12;
    ient1 = ient1s[ient12];
    ient2 = ient2s[ient12];
    cout << myname << " " << ient1 << " " << ient2 << " " << ient12 << endl;
    assert( acd12.channel == 12345 );
    assert( acd12.fembID == 101 );
    assert( acd12.getMetadata("viewIndex1") == ient1 );
    assert( acd12.getMetadata("viewIndex2") == ient2 );
    assert( acd12.getMetadata("viewIndex12") == ient12 );
  }

  cout << line << endl;
  cout << myname << "Check mutable view12." << endl;
  assert( acdtop.viewSize("view1/view2") == nent12 );
  for ( ient12=0; ient12<nent12; ++ient12 ) {
    const AdcChannelData* pacd12 = acdtop.mutableViewEntry("view1/view2", ient12);
    assert( pacd12 != nullptr );
    const AdcChannelData& acd12 = *pacd12;
    ient1 = ient1s[ient12];
    ient2 = ient2s[ient12];
    cout << myname << " " << ient1 << " " << ient2 << " " << ient12;
    cout << " ?= "
         << acd12.getMetadata("viewIndex1", 999) << " "
         << acd12.getMetadata("viewIndex2", 999) << " "
         << acd12.getMetadata("viewIndex12", 999);
    cout << endl;
    assert( acd12.channel == 12345 );
    assert( acd12.fembID == 101 );
    assert( acd12.hasMetadata("viewIndex1") );
    assert( acd12.hasMetadata("viewIndex2") );
    assert( acd12.hasMetadata("viewIndex12") );
    assert( acd12.getMetadata("viewIndex1") == ient1 );
    assert( acd12.getMetadata("viewIndex2") == ient2 );
    assert( acd12.getMetadata("viewIndex12") == ient12 );
  }

  cout << line << endl;
  cout << myname << "All tests passed." << endl;
  return 0;
}

//**********************************************************************

int main(int argc, char* argv[]) {
  if ( argc > 1 ) {
    cout << "Usage: " << argv[0] << " [ARG]" << endl;
    return 0;
  }
  return test_AdcChannelData();
}

//**********************************************************************
