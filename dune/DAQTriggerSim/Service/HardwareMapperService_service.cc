// HardwareMapperService_service.cc
//
// Jonathan Davies j.p.davies@sussex.ac.uk
// August 2016
//
// Description: Dish out vectors of raw::ChannelID_t to the user for requested hardware element

#include "dune/DAQTriggerSim/Service/HardwareMapperService.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <sstream>
#include <set>
#include <memory>

//......................................................
HardwareMapperService::HardwareMapperService(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg)
  : fLogLevel(1), fNASICsPerBoard(1), fNBoardsPerAPA(1), fServiceName("HardwareMapperService") {
  pset.get_if_present<unsigned int>("LogLevel", fLogLevel);
  pset.get_if_present<unsigned int>("NASICsPerBoard", fNASICsPerBoard);
  pset.get_if_present<unsigned int>("NBoardsPerAPA", fNBoardsPerAPA);
  const std::string func_name = "HardwareMapperService";
  {
    mf::LogInfo loginfo(fServiceName);
    loginfo << func_name << "\n"
            << "fLogLevel: " << fLogLevel << "\n"
            << "fNASICsPerBoard: " << fNASICsPerBoard << "\n"
            << "fNBoardsPerAPA: " << fNBoardsPerAPA;
  }
  //jpd -- This is how to register a callback function to the activity registry
  //    -- This means that upon the preBeginRun (which I assume is immediately before
  //    -- begin run) this function (HardwareMapperService::preBeginRun) will be called

  reg.sPreBeginRun.watch(this, &HardwareMapperService::checkGeomVsFileDetectorName);
  fillHardwareMaps();
}

//......................................................
void HardwareMapperService::checkGeomVsFileDetectorName(art::Run const& run){
  const std::string func_name = "checkGeomVsFileDetectorName";
  if(fLogLevel>1) mf::LogInfo(fServiceName) << "In Function: " << func_name;

  //jpd -- Check that the Detector Name matches between the Geometry and the File
  //Inspired by Geometry_service.cc
  mf::LogInfo loginfo(fServiceName);

  std::vector< art::Handle<sumdata::RunData> > rdcol;
  run.getManyByType(rdcol);
  if (rdcol.empty()) {
    mf::LogError(fServiceName) << "Can't find sumdata::RunData for geometry vs. file detector name checks";
    return;
  }

  fDetectorNameFromFile = rdcol.front()->DetName();

  loginfo << "This file has the following detector name:   " <<  fDetectorNameFromFile     << "\n";
  loginfo << "The geometry service has this detector name: " <<  fDetectorNameFromGeometry << "\n";
  
  if(fDetectorNameFromGeometry != fDetectorNameFromFile){
    throw cet::exception("HardwareMapperService") << "Detector Name missmatch!\n"
                                                  << "Detector Name from Geometry Service: " << fDetectorNameFromGeometry << "\n"
                                                  << "Detector Name from File: " << fDetectorNameFromFile << "\n"
                                                  << "This means that the HardwareMapperService is filled with information for the wrong Geometry\n"
                                                  << "Throwing an exception as the user should make sure the Geometry in their fcl and that used to produce their file match\n";
  }
}

//......................................................
void HardwareMapperService::fillASICMap(){
  const std::string func_name = "fillASICMap";

  mf::LogInfo loginfo(fServiceName);

  if(fLogLevel>1||true) loginfo << "In Function: " << func_name << "\n";//FIXME

  loginfo << "Filling ASIC Map - fNASICsPerBoard: " << fNASICsPerBoard << "\n";

  if(fLogLevel>0) loginfo << "Looping over Boards" << "\n";

  Hardware::ID current_asic_id = 0;
  for(auto id_board_pair: fBoardMap){
    std::shared_ptr<Hardware::Board> board_ptr = id_board_pair.second;
    
    if(fLogLevel>0) loginfo << *board_ptr << "\n";
    
    unsigned int channels_in_this_asic = 0;
    unsigned int max_channels_in_asic = board_ptr->getNChannels() / fNASICsPerBoard;
    
    if(fLogLevel>1) loginfo << "current_asic_id: " << current_asic_id << "\n";
    
    board_ptr->addHardwareID(Hardware::ASIC(current_asic_id));
    std::shared_ptr<Hardware::ASIC> this_asic_ptr = std::shared_ptr<Hardware::ASIC>(new Hardware::ASIC(current_asic_id));
    
    if(fLogLevel>1) loginfo << "Creating new asic" 
                            << " registered to Board: " << board_ptr->getID()
                            << "- This Board has: " << board_ptr->getNHardwareIDs()
                            << " pieces of hardware" 
                            << "\n\n";
    
    if(fLogLevel>1) loginfo << "Looping over channels" << "\n";
    
    for(auto channel : board_ptr->getChannels()){
      if(channels_in_this_asic >= max_channels_in_asic){
        if(fLogLevel>1) loginfo << "Need new asic - current_asic_id: " << current_asic_id
                                << " channels_in_this_asic: " << channels_in_this_asic
                                << "/" << max_channels_in_asic
                                << "\n";
        
        
        fASICMap[this_asic_ptr->getID()] = this_asic_ptr;
        current_asic_id++;
        channels_in_this_asic=0;
        
        if(fLogLevel>1) loginfo << "Last ASIC: " << *this_asic_ptr << "\n\n";
        
        this_asic_ptr=std::shared_ptr<Hardware::ASIC>(new Hardware::ASIC(current_asic_id));
        
        board_ptr->addHardwareID(Hardware::ASIC(current_asic_id));
        if(fLogLevel>1) loginfo << "Creating new asic" 
                                << " registered to Board: " << board_ptr->getID()
                                << "- This Board has: " << board_ptr->getNHardwareIDs()
                                << " pieces of hardware" 
                                << "\n\n";
      }
      this_asic_ptr->addChannel(channel);
      channels_in_this_asic++;
      
    }//channel in this board
    if(fLogLevel>1) loginfo << "Finished channel loop\n"
                            << "current_asic_id: " << current_asic_id
                            << " channels_in_this_asic: " << channels_in_this_asic
                            << "/" << max_channels_in_asic
                            << "\n";
    fASICMap[this_asic_ptr->getID()] = this_asic_ptr;
    current_asic_id++;//start each Board with a new ASIC 
    if(fLogLevel>0) loginfo << "Finished this Board\n"
                            << *board_ptr << "\n"
                            << "Total ASICs in fASICMap: " << getNASICs() << "\n"
                            << "-----------------------\n\n";
  }//id_board_pair;
  loginfo << "Finished filling ASIC Map\n"
          << "Filled: " << getNASICs() << " ASICs";
  
  
}

//......................................................
void HardwareMapperService::fillBoardMap(){
 const std::string func_name = "fillBoardMap";
  if(fLogLevel>1||true) mf::LogInfo(fServiceName) << "In Function: " << func_name;//FIXME

  mf::LogInfo loginfo(fServiceName);
  loginfo << "Filling Board Map - fNBoardsPerAPA: " << fNBoardsPerAPA << "\n";

  Hardware::ID current_board_id = 0;
  if(fLogLevel>0) loginfo << "Looping over apas" << "\n";

  // --- Lets loop over the APA map
  for(auto id_apa_pair: fAPAMap){
    // --- Get a handle to this APA object.
    std::shared_ptr<Hardware::APA> apa_ptr = id_apa_pair.second;
    // --- How many channels should there be per board? 
    unsigned int max_channels_in_board = apa_ptr->getNChannels() / fNBoardsPerAPA;

    // Karl::Added to only run on the first APA.
    //if (apa_ptr->getID() > 1) continue;

    // --- Record the id of the current board that we are on.
    Hardware::ID FirstBoard  = current_board_id;
    // --- Make an array of boards for this APA, and then create them.
    std::vector< std::shared_ptr<Hardware::Board> > BoardsFromAPA;
    std::vector< TVector3 >      ZPosOfColWires;
    for (unsigned int cb=0; cb<fNBoardsPerAPA; ++cb) {
      // --- Add a board object to this APA
      apa_ptr->addHardwareID(Hardware::Board( current_board_id ));
      // --- Push back objects to my vectors.
      BoardsFromAPA .push_back( std::shared_ptr<Hardware::Board>(new Hardware::Board( current_board_id )) );
      ZPosOfColWires.push_back( {DBL_MAX, -DBL_MAX, -1} ); 
      std::cout << "Creating a new board with id " << FirstBoard+cb << " - " << BoardsFromAPA[ cb ]->getID() << " on APA " << apa_ptr->getID() << std::endl;
      // --- Incremment the current_board_id number.
      ++current_board_id;
    }
    
    // --- Now that I have a vector of boards for this APA, lets loop over the channels on this APA and assign them to each board! 
    if(fLogLevel>1) loginfo << "Looping over channels" << "\n";
    // I need to find out the range of z values which this APA covers, and the number of collection plane wires, in order to correctly assign the channels.
    // --- !!! I am going to explicitly assuming wrapped wires, and that both sides of the TPC are instrumented !!! --- 
    // --- !!! This is valid in all of the geometries that I am aware of 2nd June 2017                          !!! ---
    double ThisAPA_MinZ = DBL_MAX;
    double ThisAPA_MaxZ = -DBL_MAX;
    int    TotColWires  = 0;
    for(auto channel : apa_ptr->getChannels()) {
      if (fGeometryService->View(channel) != 2) continue; // Collection planes wires are view 2.
      ++TotColWires; // Got a collection plane wire.
      // Work out the Z position of this wire.
      double WireSt[3], WireEn[3];
      fGeometryService->WireEndPoints( fGeometryService->ChannelToWire(channel)[0], WireSt, WireEn );
      if (WireSt[2] < ThisAPA_MinZ ) ThisAPA_MinZ = WireSt[2];
      if (WireSt[2] > ThisAPA_MaxZ ) ThisAPA_MaxZ = WireSt[2];
    }
    // --- Work out the number of collection plane wires per board.
    const int NumColWiresPerBoard = TotColWires/fNBoardsPerAPA;
    if (TotColWires%fNBoardsPerAPA != 0) {
      throw cet::exception("HardwareMapperService") << "The number of collection plane wires doesn't divide evenly. There are " << TotColWires << " col wires,"
						    << " and you want " << fNBoardsPerAPA << ", leaving a remainder of " << TotColWires%fNBoardsPerAPA << " wires.";
    }
    std::cout << " There are " << TotColWires << " collection plane wires on this APA."
	      << " Each board will have " << NumColWiresPerBoard << " collection plane wires, and " << max_channels_in_board - TotColWires/fNBoardsPerAPA << " ind wires."
	      << " The wires go from " << ThisAPA_MinZ << " to " << ThisAPA_MaxZ << ".\n"
	      << std::endl;

    // --- Set the collection plane wires on our boards.
    int NumColPush = 0, PushToBoard = 0;
    for(auto channel : apa_ptr->getChannels()) {
      if (fGeometryService->View(channel) != 2) continue; // Collection planes wires are view 2.
      // Get the start and end points of the first wire segment for this channel.
      double WireSt[3], WireEn[3];
      fGeometryService->WireEndPoints( fGeometryService->ChannelToWire(channel)[0], WireSt, WireEn );
      BoardsFromAPA[ PushToBoard ]->addChannel(channel);
      ++NumColPush;
      if( WireSt[2] < ZPosOfColWires[ PushToBoard ][0] ) ZPosOfColWires[ PushToBoard ][0] = WireSt[2];
      if( WireEn[2] > ZPosOfColWires[ PushToBoard ][1] ) ZPosOfColWires[ PushToBoard ][1] = WireEn[2];
      ZPosOfColWires[ PushToBoard ][2] = fGeometryService->FindTPCAtPosition( WireSt ).TPC;
      /*
      std::cout << "   Looking at channel " << channel << " on plane " << fGeometryService->View(channel) << " TPC [[" << fGeometryService->FindTPCAtPosition( WireSt ) << " ]]."
		<< " Start pos was (" << WireSt[0] << ", " << WireSt[1] << ", " << WireSt[2] << "). "
		<< " Endin pos was (" << WireEn[0] << ", " << WireEn[1] << ", " << WireEn[2] << ").\n"
		<< "        Pushing it back to board " << PushToBoard << ", now pushed back " << NumColPush << " wires."
		<< " Extremes are " << ZPosOfColWires[ PushToBoard ].first << " - " << ZPosOfColWires[ PushToBoard ].second
		<< std::endl;
      */
      // --- Check that if I want to move to the next board.
      if ( NumColPush%NumColWiresPerBoard == 0 ) {
	++PushToBoard;
      }
    }

    for (unsigned int nn=0; nn<fNBoardsPerAPA; ++nn) {
      std::cout << " Looking at board " << BoardsFromAPA[ nn ]->getID() << ", in TPC " << ZPosOfColWires[ nn ][2]
		<< " the extremes are " << ZPosOfColWires[ nn ][0] << " - " << ZPosOfColWires[ nn ][1] << std::endl;
    }

    // --- Now set the induction plane wires on our boards.
    int NumIndPush = 0;
    
    for(auto channel : apa_ptr->getChannels()) {
      if (fGeometryService->View(channel) == 2) continue; // Collection planes wires are view 2.
      // --- When looking at induction plane wires whether I look at start or end wire pos depends on APA and TPC number.
      // --- When looking at even APA number - WireEnd   for even TPCs, WireStart for odd TPCs.
      // --- When looking at odd  APA number - WireStart for odd TPCs, WireEnd   for even TPCs.
      // Get the start and end points of the first wire segment for this channel.
      double WireSt[3], WireEn[3];
      fGeometryService->WireEndPoints( fGeometryService->ChannelToWire(channel)[0], WireSt, WireEn );
      /*
      std::cout << "\nLooking at channel " << channel << " on plane " << fGeometryService->View(channel) << " TPC [[" << fGeometryService->FindTPCAtPosition( WireSt ) << " ]]."
		<< " Start pos was (" << WireSt[0] << ", " << WireSt[1] << ", " << WireSt[2] << "). "
		<< " Endin pos was (" << WireEn[0] << ", " << WireEn[1] << ", " << WireEn[2] << ")."
		<< std::endl;
      */
      for (unsigned int nb=0; nb<fNBoardsPerAPA; ++nb) {
	if ( ZPosOfColWires[nb][2] != fGeometryService->FindTPCAtPosition( WireSt ).TPC ) continue;
	if ( (abs(WireSt[1]) > 600 && WireSt[2] > ZPosOfColWires[nb][0] && WireSt[2] < ZPosOfColWires[nb][1] ) ||
	     (abs(WireEn[1]) > 600 && WireEn[2] > ZPosOfColWires[nb][0] && WireEn[2] < ZPosOfColWires[nb][1] ) ) {
	  BoardsFromAPA[ nb ]->addChannel(channel);
	  /*
	  std::cout << "   Pushing this channel onto board: " << BoardsFromAPA[ nb ]->getID()
		    << ", Extremeities " << ZPosOfColWires[nb][0] << " - " << ZPosOfColWires[nb][1] << " in TPC " << ZPosOfColWires[nb][2]
		    << ". This board now has " << BoardsFromAPA[ nb ]->getNChannels() << std::endl; 
	  */
	  break;
	}
      }
      ++NumIndPush;
    }//channel in this apa
    
    // --- Now that I have looped over all of the channels, lets update my map!
    for (unsigned int cb=0; cb<fNBoardsPerAPA; ++cb) {
      unsigned int ThisID = BoardsFromAPA[ cb ]->getID();
      fBoardMap[ ThisID ] = BoardsFromAPA[ cb ];
      std::cout << "Adding board " << BoardsFromAPA[ cb ]->getID() << " to APA " << apa_ptr->getID() 
		<< ", it has " << BoardsFromAPA[ cb ]->getNChannels() << " - " << fBoardMap[ ThisID ]->getNChannels()
		<< std::endl;
    }
    if(fLogLevel>0) loginfo << "Finished this APA\n" << *apa_ptr << "\n" << "Total Boards in fBoardMap: " << getNBoards() << "\n-----------------------\n\n";
    // --- Start each APA with a new board
    //current_board_id++;
  }//id_apa_pair;
  
  loginfo<< "Finished filling Board Map\n Filled: " << getNBoards() << " Boards";
}

//......................................................
void HardwareMapperService::fillTPCMap(){
  const std::string func_name = "fillTPCMap";
  mf::LogInfo loginfo(fServiceName);
  if(fLogLevel>1) loginfo << "In Function: " << func_name << "\n";
  unsigned int Nchannels   = fGeometryService->Nchannels();
  loginfo  << "Filling TPC Map with " << Nchannels << " channels" << "\n";

  for(raw::ChannelID_t channel=0; channel<Nchannels ;channel++){
    std::vector<geo::WireID> const Wires = fGeometryService->ChannelToWire(channel);
    for(auto wire : Wires){
      auto tpc_id = wire.TPC;
      //jpd -- See if we have already created a TPC object for this tpc_id in our map
      auto find_result = fTPCMap.find(tpc_id);
      std::shared_ptr<Hardware::TPC> this_tpc;
      if(find_result != fTPCMap.end()){
        //jpd -- We already have one, add this channel
        this_tpc = (*find_result).second;
      }
      else{
        //jpd -- We don't have one. Create a new one and add this channel
        this_tpc = std::shared_ptr<Hardware::TPC>(new Hardware::TPC(tpc_id));
      }
      this_tpc->addChannel(channel);
      fTPCMap[tpc_id] = this_tpc;
    }//loop over wires
  }//loop over channels
  loginfo  << "Finished filling TPC Map\n"
           << "Filled: " << getNTPCs() << " TPCs";

}

//......................................................
void HardwareMapperService::fillAPAMap(){
  const std::string func_name = "fillAPAMap";

  mf::LogInfo loginfo(fServiceName);
  if(fLogLevel>1)  loginfo << "In Function: " << func_name << "\n";

  unsigned int Nchannels   = fGeometryService->Nchannels();
  loginfo  << "Filling APA Map with " << Nchannels << " channels" << "\n";

  for(raw::ChannelID_t channel=0; channel<Nchannels ;channel++){
    std::vector<geo::WireID> const Wires = fGeometryService->ChannelToWire(channel);
    for(auto wire : Wires){
      auto apa_id = wire.TPC / 2;
      //jpd -- See if we have already created a APA object for this apa_id in our map
      auto find_result = fAPAMap.find(apa_id);
      std::shared_ptr<Hardware::APA> this_apa;
      if(find_result != fAPAMap.end()){
        //jpd -- We already have one, add this channel
        this_apa = (*find_result).second;
      }
      else{
        //jpd -- We don't have one. Create a new one and add this channel
        this_apa = std::shared_ptr<Hardware::APA>(new Hardware::APA(apa_id));
      }
      this_apa->addChannel(channel);
      fAPAMap[apa_id] = this_apa;
    }//loop over wires
  }//loop over channels
  loginfo  << "Finished filling APA Map\n"
           << "Filled: " << getNAPAs() << " APAs";
}


//......................................................
void HardwareMapperService::fillHardwareMaps(){
  const std::string func_name = "fillHardwareMaps";
  if(fLogLevel>1) mf::LogInfo(fServiceName) << "In Function: " << func_name;
  fDetectorNameFromGeometry = fGeometryService->DetectorName();

  {
    mf::LogInfo loginfo(fServiceName);
    
    loginfo << "HardwareMapperService\n";
    loginfo << "Filling Hardware Maps using geometry service\n";
    loginfo << "DetectorName: " << fGeometryService->DetectorName() << "\n";
    loginfo << "Ncryostats:   " << fGeometryService->Ncryostats()   << "\n";
    loginfo << "TotalNTPC:    " << fGeometryService->TotalNTPC()    << "\n";
    loginfo << "Nchannels:    " << fGeometryService->Nchannels()    << "\n";
    loginfo << "NOpChannels:  " << fGeometryService->NOpChannels()  << "\n";
  }//using annonymous namespace to force mf::LogInfo destructor call to flush output

  fillTPCMap();
  fillAPAMap();
  fillBoardMap();
  fillASICMap();
  
  mf::LogInfo loginfo(fServiceName);
  loginfo << "Number of TPCs:   " << getNTPCs() << "\n"
          << "Number of APAs:   " << getNAPAs() << "\n"
          << "Number of Boards: " << getNBoards() << "\n"
          << "Number of ASICs:  " << getNASICs() << "\n";

}

//......................................................
void HardwareMapperService::printASICMap(unsigned int num_asics_to_print){
  const std::string func_name = "printASICMap";
  if(fLogLevel>1) mf::LogInfo(fServiceName) << "In Function: " << func_name;
  mf::LogInfo loginfo(fServiceName);
  loginfo << "Printing the first: " << num_asics_to_print << " ASICs\n";
  unsigned int total_channels = 0;
  unsigned int count = 0;
  for(auto this_pair : fASICMap ){
    if(count++ >= num_asics_to_print) break;
    std::shared_ptr<Hardware::ASIC> this_asic = this_pair.second;
    total_channels += this_asic->getNChannels();
    loginfo << *this_asic << "\n";
  }
  loginfo << "\n";
  loginfo << "Total Channels: " << total_channels << "\n";
  loginfo << "\n";
}

//......................................................
void HardwareMapperService::printBoardMap(unsigned int num_boards_to_print){
  const std::string func_name = "printBoardMap";
  if(fLogLevel>1) mf::LogInfo(fServiceName) << "In Function: " << func_name;
  mf::LogInfo loginfo(fServiceName);
  loginfo << "Printing the first: " << num_boards_to_print << " Boards\n";
  unsigned int total_channels = 0;
  unsigned int count = 0;
  for(auto this_pair : fBoardMap ){
    if(count++ >= num_boards_to_print) break;
    std::shared_ptr<Hardware::Board> this_board = this_pair.second;
    total_channels += this_board->getNChannels();
    loginfo << *this_board << "\n";
  }
  loginfo << "\n";
  loginfo << "Total Channels: " << total_channels << "\n";
  loginfo << "\n";
}

//......................................................
void HardwareMapperService::printTPCMap(unsigned int num_tpcs_to_print){
  const std::string func_name = "printTPCMap";
  if(fLogLevel>1) mf::LogInfo(fServiceName) << "In Function: " << func_name;
  mf::LogInfo loginfo(fServiceName);
  loginfo << "Printing the first: " << num_tpcs_to_print << " TPCs\n";
  unsigned int total_channels = 0;
  unsigned int count = 0;
  for(auto this_pair : fTPCMap ){
    if(count++ >= num_tpcs_to_print) break;
    std::shared_ptr<Hardware::TPC> this_tpc = this_pair.second;
    total_channels += this_tpc->getNChannels();
    loginfo << *this_tpc << "\n";
  }
  loginfo << "\n";
  loginfo << "Total Channels: " << total_channels << "\n";
  loginfo << "\n";
}

//......................................................
void HardwareMapperService::printAPAMap(unsigned int num_apas_to_print){
  const std::string func_name = "printAPAMap";
  if(fLogLevel>1) mf::LogInfo(fServiceName) << "In Function: " << func_name;
  mf::LogInfo loginfo(fServiceName);

  loginfo << "Printing the first: " << num_apas_to_print << " APAs\n";
  unsigned int total_channels = 0;
  unsigned int count = 0;
  for(auto this_pair : fAPAMap ){
    if(count++ >= num_apas_to_print) break;
    std::shared_ptr<Hardware::APA> this_apa = this_pair.second;
    total_channels += this_apa->getNChannels();
    loginfo << *this_apa << "\n";
  }
  loginfo << "\n";
  loginfo << "Total Channels: " << total_channels << "\n";
  loginfo << "\n";
}

//......................................................
void HardwareMapperService::printHardwareMaps(){
  const std::string func_name = "printHardwareMaps";
  if(fLogLevel>1) mf::LogInfo(fServiceName) << "In Function: " << func_name;
  mf::LogInfo(fServiceName)  
    << "fDetectorName: " << fDetectorNameFromGeometry << " taken from Geometry Service at time of fillHardwareMaps()" 
    << "\n"
    << "N.B. the geometry service can reload a new geometry on beginRun(), which would be bad"
    << "\n";

  printASICMap(getNASICs());
  printBoardMap(getNBoards());
  printTPCMap(getNTPCs());
  printAPAMap(getNAPAs());


  mf::LogInfo(fServiceName)
    << "Num ASICs: " << getNASICs() << "\n"
    << "Num Boards: " << getNBoards() << "\n"
    << "Num TPCs: " << getNTPCs() << "\n"
    << "Num APAs: " << getNAPAs();

}

//......................................................
void HardwareMapperService::printGeometryInfo(){
  const std::string func_name = "printGeometryInfo";
  if(fLogLevel>1) mf::LogInfo(fServiceName) << "In Function: " << func_name;
  {
    mf::LogInfo loginfo(fServiceName);
    loginfo << "Geometry Information with current Geometry Service\n";
    loginfo << "DetectorName: " << fGeometryService->DetectorName() << "\n";
    loginfo << "TotalMass: " << fGeometryService->TotalMass() << "\n";
    unsigned int Ncryostats  = fGeometryService->Ncryostats();
    unsigned int TotalNTPC   = fGeometryService->TotalNTPC();
    unsigned int Nchannels   = fGeometryService->Nchannels();
    unsigned int NOpChannels = fGeometryService->NOpChannels();
    loginfo << "Ncryostats:   " << Ncryostats   << "\n";
    loginfo << "TotalNTPC:    " << TotalNTPC    << "\n";
    loginfo << "Nchannels:    " << Nchannels    << "\n";
    loginfo << "NOpChannels:  " << NOpChannels  << "\n";
    loginfo << "\n";
  }//using annonymous namespace to force mf::LogInfo destructor (which makes flushes the output)

  printHardwareMaps();  

  return;
}

//......................................................
std::vector<raw::ChannelID_t> const& HardwareMapperService::getTPCChannels(Hardware::ID tpc_id){
  const std::string func_name = "getTPCChannels";
  mf::LogInfo loginfo(fServiceName);

  if(fLogLevel>1) loginfo << "In Function: " << func_name << "\n";
  if(fLogLevel>1) loginfo << "Finding channels for TPC: " << tpc_id << "\n";

  auto find_result = fTPCMap.find(tpc_id);
  std::shared_ptr<Hardware::TPC> this_tpc;
  
  if(find_result == fTPCMap.end()){
    mf::LogError(fServiceName) << "Failed to find this TPC: " << tpc_id << "\n";
    mf::LogError(fServiceName) << "Returning an empty vector of channels\n";
    static std::vector<raw::ChannelID_t> emptyVector;
    return emptyVector;
  }
  this_tpc = (*find_result).second;
  if(fLogLevel>1) loginfo << "Found " << *this_tpc << "\n";
  if(fLogLevel>1) loginfo << "\n";
  return this_tpc->getChannels();  
}

//......................................................
std::vector<raw::ChannelID_t> const& HardwareMapperService::getAPAChannels(Hardware::ID apa_id){
  const std::string func_name = "getAPAChannels";
  if(fLogLevel>1) mf::LogInfo(fServiceName) << "In Function: " << func_name;
  mf::LogInfo loginfo(fServiceName);
  if(fLogLevel>1) loginfo << "Finding channels for APA: " << apa_id << "\n";

  auto find_result = fAPAMap.find(apa_id);
  std::shared_ptr<Hardware::APA> this_apa;
  
  if(find_result == fAPAMap.end()){
    mf::LogError(fServiceName) << "Failed to find this APA: " << apa_id << "\n";
    mf::LogError(fServiceName) << "Returning an empty vector of channels\n";
    static std::vector<raw::ChannelID_t> emptyVector;
    return emptyVector;
  }
  this_apa = (*find_result).second;
  if(fLogLevel>1) loginfo << "Found " << *this_apa << "\n";
  if(fLogLevel>1) loginfo << "\n";
  return this_apa->getChannels();  
}

//......................................................
std::set<raw::ChannelID_t> const& HardwareMapperService::getTPCChannelsSet(Hardware::ID tpc_id){
  const std::string func_name = "getTPCChannelsSet";
  if(fLogLevel>1) mf::LogInfo(fServiceName) << "In Function: " << func_name;
  mf::LogInfo loginfo(fServiceName);
  if(fLogLevel>1) loginfo << "Finding channels for TPC: " << tpc_id << "\n";

  auto find_result = fTPCMap.find(tpc_id);
  std::shared_ptr<Hardware::TPC> this_tpc;
  
  if(find_result == fTPCMap.end()){
    mf::LogError(fServiceName) << "Failed to find this TPC: " << tpc_id << "\n";
    mf::LogError(fServiceName) << "Returning an empty set of channels\n";
    static std::set<raw::ChannelID_t> emptySet;
    return emptySet;
  }
  this_tpc = (*find_result).second;
  if(fLogLevel>1) loginfo << "Found " << *this_tpc << "\n";
  if(fLogLevel>1) loginfo << "\n";
  return this_tpc->getChannelsSet();  
}

//......................................................
std::set<raw::ChannelID_t> const& HardwareMapperService::getAPAChannelsSet(Hardware::ID apa_id){
  const std::string func_name = "getAPAChannelsSet";
  if(fLogLevel>1) mf::LogInfo(fServiceName) << "In Function: " << func_name;
  mf::LogInfo loginfo(fServiceName);
  if(fLogLevel>1) loginfo << "Finding channels for APA: " << apa_id << "\n";

  auto find_result = fAPAMap.find(apa_id);
  std::shared_ptr<Hardware::APA> this_apa;
  
  if(find_result == fAPAMap.end()){
    mf::LogError(fServiceName) << "Failed to find this APA: " << apa_id << "\n";
    mf::LogError(fServiceName) << "Returning an empty set of channels\n";
    static std::set<raw::ChannelID_t> emptySet;
    return emptySet;
  }
  this_apa = (*find_result).second;
  if(fLogLevel>1) loginfo << "Found " << *this_apa << "\n";
  if(fLogLevel>1) loginfo << "\n";
  return this_apa->getChannelsSet();  
}

//......................................................
DEFINE_ART_SERVICE(HardwareMapperService)

