#include <iostream>
#include <string>
#include <iomanip>

// adjust WIB numbers from the Vertical Drift Coldbox BDE channel maps
// so they can be used for the vertical slice test data

int main(int argc, char **argv)
{
  const bool debug=false;

  std::string stripid;
  int connector = 0;
  int board = 0;
  int cebchannel = 0;
  int asic = 0;
  int asicchan = 0;

  int wib=0;
  int wibconnector=0;
  int chan = 0;

  while (true)
    {
      std::cin >> chan >> wib >> wibconnector >> cebchannel >> board >> asic >> asicchan >> connector >> stripid;

      if (feof(stdin)) break;

      // two WIBs were shifted over in the crate due to cooling issues,
      // just for the vertical slice integration test
      
      if (wib > 2) wib++;
      
       std::cout << std::setw(8) << chan << " " << std::setw(8) << wib << " " << std::setw(8) << wibconnector << " " << std::setw(8) << cebchannel << " " << std::setw(8) << board << " " << std::setw(8) << asic << " " << std::setw(8) << asicchan << " " << std::setw(8) << connector << " " << std::setw(8) << stripid << std::endl;
       
    }
}
