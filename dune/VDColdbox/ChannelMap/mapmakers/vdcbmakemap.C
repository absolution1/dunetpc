#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>

int main(int argc, char **argv)
{
  const bool debug=false;

  std::string stripid;
  int connector = 0;
  int pin = 0;
  int board = 0;
  int cebchannel = 0;
  int asic = 0;
  int asicchan = 0;

  int wib=0;
  int wibconnector=0;

  std::vector<int> wibvec{0,1,1,1,1,2,2,2,2,3,3,3,3,4,4};
  std::vector<int> wibconnectorvec{0,1,2,3,4,1,2,3,4,1,2,3,4,1,2};
  
  while (true)
    {
      std::cin >> stripid >> connector >> pin >> board >> cebchannel >> asic >> asicchan;
      if (feof(stdin)) break;
      
      if (debug) std::cout << stripid << " " << connector << " " << pin << " " << board << " " << cebchannel << " " << asic << " " << asicchan << std::endl;

      if (board < 1 || board > 14)
	{
	  std::cout << "invalid CE board: " << board << std::endl;
	  return 1;
	}
      wib = wibvec.at(board);
      wibconnector = wibconnectorvec.at(board);

      // parse the strip ID
      std::string sp = stripid.substr(0,1);
      if (debug) std::cout << " strip ID prefix: " << sp << std::endl;
      std::string sn = stripid.substr(1);
      std::stringstream sns;
      sns << sn;
      int stripidnum;
      sns >> stripidnum;
      if (debug) std::cout << " strip ID number: " << stripidnum << std::endl;

      int chan = -1;
      if (sp == "U")
	{
	  chan = stripidnum + 864 - 1;
	  if (stripidnum > 256)
	    {
	      chan = stripidnum + 2720 - 257;
	    }
	}
      else if (sp == "Y")
	{
	  chan = stripidnum + 1120 - 1;
	  if (stripidnum > 320)
	    {
	      chan = stripidnum + 2848 - 321;
	    }
	}
      else if (sp == "Z")
	{
	  chan = stripidnum + 1440 - 1;
	  if (stripidnum > 288)
	    {
	      chan = stripidnum + 2289 - 289;
	    }
	}
      if (chan == -1)
	{
	  std::cout << "Could not assign channel.  Use the debug flag to get more info." << std::endl;
	}
      std::cout << std::setw(8) << chan << " " << std::setw(8) << wib << " " << std::setw(8) << wibconnector << " " << std::setw(8) << cebchannel << " " << std::setw(8) << board << " " << std::setw(8) << asic << " " << std::setw(8) << asicchan << " " << std::setw(8) << connector << " " << std::setw(8) << stripid << std::endl;
      
    }
  
  return 0;
}
