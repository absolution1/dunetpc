
#include "RawDecodingUtils.hh"

#include "artdaq-core/Data/ContainerFragment.hh"

#include "art/Framework/Principal/Handle.h"

#include <stdexcept>
//#include "cetlib/exception.h"

//#include <sstream>

artdaq::Fragments
dune::getHandleChecked(const art::Event& evt, const std::string& label, 
			const std::string& fragtype, const bool& expectFragsInContainer) {
  
  artdaq::Fragments fragments;

/*  
  art::Handle<artdaq::Fragments> handleFragments;

  if (expectFragsInContainer) {
    art::InputTag itag1(label, "Container");
    handleFragments = evt.getHandle<artdaq::Fragments>(itag1);
  } 
  else {
    art::inputTag itag2(label, fragtype);
    handleFragments = evt.getHandle<artdaq::Fragments>(itag2);
  }

  if(!handleFragments.isValid()) {
    std::cerr << "Run: " << evt.run()
	      << ", SubRun: " << evt.subRun()
	      << ", Event: " << evt.event()
	      << " is NOT VALID" << std::endl;

    if (expectFragsInContainer) {
      //throw cet::exception("getHandleChecked") << "getHandle call for fragments of type \"Container\" NOT VALID";
      throw std::runtime_error("(THIS SHOULD BE CHANGED TO A cet::exception) getHandle call for fragments NOT VALID");
    } else {
      //      throw cet::exception("getHandleChecked") << "getHandle call for fragments of type \"" << fragtype << "\" NOT VALID";
      throw std::runtime_error("(THIS SHOULD BE CHANGED TO A cet::exception) getHandle call for fragments NOT VALID");
    }
  }
  if (expectFragsInContainer) {
    
    for (auto cont : *handleFragments)
      {
	artdaq::ContainerFragment contf(cont);

	for (size_t ii = 0; ii < contf.block_count(); ++ii)
	  {
	    size_t fragSize = contf.fragSize(ii);
	    artdaq::Fragment thisfrag;
	    thisfrag.resizeBytes(fragSize);

	    std::cout << "Copying " << fragSize << " bytes from " << contf.at(ii) << " to " << thisfrag.headerAddress() << std::endl;
	    memcpy(thisfrag.headerAddress(), contf.at(ii), fragSize);

	    std::cout << "Putting new fragment into output vector" << std::endl;             
	    fragments.emplace_back(thisfrag);
	  }
      }
  } else {

    for(auto const& frag: *handleFragments){
      fragments.emplace_back( frag );
    }
  }
*/
  // C++11 will move, rather than copy, a vector on return
  return fragments;
}

