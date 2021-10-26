// Prototype program to try out the VD coldbox channel mapping.  Components to be made
// into a service
// compile with clang++ --std=c++17 -o maptest maptest.C

#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <vector>

struct VDCBChanInfo {
  int offlchan;        // in gdml and channel sorting convention
  int wib;             // slot number +1:    1, 2, 3, or 4
  int wibconnector;    // which board on the WIB  1, 2, 3, or 4
  int cebchan;         // cold electronics channel on board:  0 to 127
  int femb;            // FEMB ID:  1 to 14
  int asic;            // ASIC:   1 to 8
  int asicchan;        // ASIC channel:  0 to 15
  int connector;       // detector connector
  std::string stripid;  // strip name, with plane and number.  e.g. U79
  bool valid;          // true if valid, false if not
};

//  map so we can look up channel info by offline channel
std::unordered_map<int,VDCBChanInfo> chantoinfomap;
std::unordered_map<int,std::unordered_map<int,std::unordered_map<int,int> > > infotochanmap;

// function prototypes

VDCBChanInfo getChanInfoFromOfflChan(int offlchan);

// this uses conventions from Nitish's spreadsheet.  WIB: 1-3, wibconnector: 1-4, cechan: 0-127
int getOfflChanFromWIBConnectorInfo(int wib, int wibconnector, int cechan);

// this uses conventions from the DAQ WIB header, with two FEMBs per fiber
// on FELIX readout:  slot: 0-2, fiber=1 or 2, cehcan: 0-255
int getOfflChanFromSlotFiberChan(int slot, int fiber, int chan);

int main(int argc, char **argv)
{
  const bool debug=false;

  // when using this program in test mode, turn testmainchans on.  Looks for duplicates and
  // looks up all the channels we put in
  
  const bool testmainchans = true;

  // when running in add Disconnected chans mode, turn this switch on.  Will print out an addition
  // to the map corresponding to the disconnected channels
  
  const bool addDisconnectedChans = false;
  const int disconnectedChanOffset = 0; // 3456;
  std::vector<int> wibvec{0,1,1,1,1,2,2,2,2,3,3,3,3,4,4};
  std::vector<int> wibconnectorvec{0,1,2,3,4,1,2,3,4,1,2,3,4,1,2};

  //std::string fullname("vdcbce_chanmap_v1.txt");
  std::string fullname("vdcbce_chanmap_v1_dcchan0.txt");
  std::ifstream inFile(fullname, std::ios::in);
  std::string line;
  int numchans = 0;
  while (std::getline(inFile,line)) {

    VDCBChanInfo chinfo;
    std::stringstream linestream(line);
    linestream >>
      chinfo.offlchan >>
      chinfo.wib >>
      chinfo.wibconnector >>
      chinfo.cebchan >>
      chinfo.femb >>
      chinfo.asic >>
      chinfo.asicchan >>
      chinfo.connector >>
      chinfo.stripid;
    chinfo.valid = true;

    // see if we have already entered this channel
    
    int otest = getOfflChanFromWIBConnectorInfo(chinfo.wib,chinfo.wibconnector,chinfo.cebchan);
    if (otest >= 0)
      {
	std::cout << "Duplicate info found: " << chinfo.wib << " " << chinfo.wibconnector << " " << chinfo.cebchan << std::endl;
	std::cout << chinfo.offlchan << " " << otest << std::endl;
      }
    
    // std::cout << chinfo.offlchan << std::endl;
    chantoinfomap[chinfo.offlchan] = chinfo;
    infotochanmap[chinfo.wib][chinfo.wibconnector][chinfo.cebchan] = chinfo.offlchan;
    ++numchans;
  }
  inFile.close();
  //std::cout << "num chans: " << numchans << std::endl;
  
  if (testmainchans)
    {
      //VDCBChanInfo ciret = getChanInfoFromOfflChan(1289);
      //std::cout << "looked up offline channel 1289: " << ciret.offlchan << " " << ciret.wib << " " << ciret.wibconnector << " " << ciret.cebchan << " " << ciret.femb << " " << ciret.asic << " " << ciret.asicchan << " "  << ciret.connector << " " << ciret.stripid << " " << ciret.valid << std::endl;

      //int ctest = getOfflChanFromWIBConnectorInfo(1, 3, 25);
      //std::cout << "reverse lookup: " << ctest << std::endl;

      //int ctest2 = getOfflChanFromSlotFiberChan(0, 2, 25);
      //std::cout << "reverse lookup2: " << ctest2 << std::endl;

      // last channel in the actual detector is 3456.  Have 192 extras, possibly on the end
      
      for (int i=0;i<3648; ++i)
	{
	  VDCBChanInfo ciret2 = getChanInfoFromOfflChan(i);
	  if (ciret2.valid)
	    {
	      std::cout << "looked up offline channel: " << i << " " << ciret2.offlchan << " " << ciret2.wib << " " << ciret2.wibconnector << " " << ciret2.cebchan << " " << ciret2.femb << " " << ciret2.asic << " " << ciret2.asicchan << " "  << ciret2.connector << " " << ciret2.stripid << " " << ciret2.valid << std::endl;

	      int circ1 = getOfflChanFromWIBConnectorInfo(ciret2.wib, ciret2.wibconnector, ciret2.cebchan);
	      if (circ1 != i)
		{
		  std::cout << "circularity check 1 failed: " << circ1 << std::endl;
		}
	    }
       
	}
    }

  if (addDisconnectedChans)
    {
      int idc = 0;
      for (int ifemb=1;ifemb<15; ++ifemb)
	{
	  for (int ichan=0; ichan<128; ++ichan)
	    {
	      int wib=wibvec.at(ifemb);
	      int wibconnector=wibconnectorvec.at(ifemb);
	      int iofflchan = getOfflChanFromWIBConnectorInfo(wib,wibconnector,ichan);
	      int asic = 1 + ichan/16;
	      int asicchan = ichan % 16;
	      int connector = 0;
	      if (iofflchan < 0)
		{
		  std::string stripid="D";
		  stripid += std::to_string(idc);
                  std::cout << std::setw(8) << idc+disconnectedChanOffset << " " << std::setw(8) << wib << " " << std::setw(8) << wibconnector << " " << std::setw(8) << ichan << " " << std::setw(8) << ifemb << " " << std::setw(8) << asic << " " << std::setw(8) << asicchan << " " << std::setw(8) << connector << " " << std::setw(8) << stripid << std::endl;
      
		  ++idc;  // to get ready for the next one
		}
	    }
	}
    }
  
  return 0;
}

// if not found in the map, return a chan info struct filled with -1's and set the valid flag to false.

VDCBChanInfo getChanInfoFromOfflChan(int offlchan)
{
  VDCBChanInfo r;
  auto fm = chantoinfomap.find(offlchan);
  if (fm == chantoinfomap.end())
    {
      r.offlchan = -1;
      r.wib = -1;
      r.wibconnector = -1;
      r.cebchan = -1;
      r.femb = -1;
      r.asic = -1;
      r.asicchan = -1;
      r.connector = -1;
      r.stripid = "INVALID";
      r.valid = false;
    }
  else
    {
      r = fm->second;
    }
  return r;
}

// returns -1 if information not found in the map

int getOfflChanFromWIBConnectorInfo(int wib, int wibconnector, int cechan)
{
  int r = -1;
  auto fm1 = infotochanmap.find(wib);
  if (fm1 == infotochanmap.end()) return r;
  auto m1 = fm1->second;
  auto fm2 = m1.find(wibconnector);
  if (fm2 == m1.end()) return r;
  auto m2 = fm2->second;
  auto fm3 = m2.find(cechan);
  if (fm3 == m2.end()) return r;
  r = fm3->second;  
  return r;
}

int getOfflChanFromSlotFiberChan(int slot, int fiber, int chan)
{
  int wc = fiber*2 - 1;
  if (chan>127)
    {
      chan -= 128;
      wc++;
    }
  return getOfflChanFromWIBConnectorInfo(slot+1,wc,chan);
}
