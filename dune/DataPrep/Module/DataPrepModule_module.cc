// DataPrepModule_module.cc

// David Adams
// July 2016
//
// Module that reads RawData and writes Wire and their associations.
// It uses RawDigitPrepService to build the wires.

//#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "dune/DuneInterface/RawDigitPrepService.h"

using std::cout;
using std::endl;
using std::string;
using art::ServiceHandle;

//**********************************************************************

class DataPrepModule : public art::EDProducer {

public:
    
  // Ctor.
  explicit DataPrepModule(fhicl::ParameterSet const& pset); 

  // Dtor.
  ~DataPrepModule();
    
  // Producer methods.
  void produce(art::Event& evt);
  void beginJob();
  void endJob();
  void reconfigure(fhicl::ParameterSet const& p);
    
private:
    
  // Configuration parameters.
  int m_LogLevel;
  std::string m_DigitLabel;  ///< Full label for the input digit container, e.g. daq:
  std::string m_WireName;    ///< Second field in full label for the output wire container.

  // Split label into producer and name: PRODUCER or PRODUCER:NAME
  std::string m_DigitProducer;
  std::string m_DigitName;

  bool m_DoAssns = false;

};

DEFINE_ART_MODULE(DataPrepModule)
  
//**********************************************************************

DataPrepModule::DataPrepModule(fhicl::ParameterSet const& pset) {
  this->reconfigure(pset);
  produces< std::vector<recob::Wire> >(m_WireName);
  if ( m_DoAssns ) {
    produces<art::Assns<raw::RawDigit, recob::Wire>>(m_WireName);
  }
}
  
//**********************************************************************

DataPrepModule::~DataPrepModule() { }

//**********************************************************************

void DataPrepModule::reconfigure(fhicl::ParameterSet const& pset) {
  const string myname = "DataPrepModule::reconfigure: ";
  m_LogLevel   = pset.get<int>("LogLevel");
  m_DigitLabel = pset.get<std::string>("DigitLabel", "daq");
  m_WireName   = pset.get<std::string>("WireName", "");
  m_DoAssns    = pset.get<bool >("DoAssns", "");
               
  size_t ipos = m_DigitLabel.find(":");
  if ( ipos == std::string::npos ) {
    m_DigitProducer = m_DigitLabel;
  } else {
    m_DigitProducer = m_DigitLabel.substr(0, ipos);
    m_DigitName = m_DigitLabel.substr(ipos + 1);
  }
  if ( m_LogLevel >= 1 ) {
    cout << myname << "    LogLevel: " << m_LogLevel << endl;
    cout << myname << "  DigitLabel: " << m_DigitLabel << " (" << m_DigitProducer
                   << ", " << m_DigitName << ")" << endl;
    cout << myname << "    WireName: " << m_WireName << endl;
    cout << myname << "     DoAssns: " << m_DoAssns << endl;
  }
}

//**********************************************************************

void DataPrepModule::beginJob() { }

//**********************************************************************

void DataPrepModule::endJob() { }
  
//**********************************************************************

void DataPrepModule::produce(art::Event& evt) {      
  const string myname = "DataPrepModule::produce: ";

  // Read in the digits. 
  art::Handle< std::vector<raw::RawDigit> > hdigits;
  evt.getByLabel(m_DigitProducer, m_DigitName, hdigits);
  if ( hdigits->size() == 0 ) mf::LogWarning("DataPrepModule") << "Input digit container is empty";

  // Process the digits.
  ServiceHandle<RawDigitPrepService> hrdp;
  AdcChannelDataMap acds;
  std::unique_ptr<std::vector<recob::Wire>> pwires(new std::vector<recob::Wire>);
  pwires->reserve(hdigits->size());
  int rstat = hrdp->prepare(*hdigits, acds, pwires.get());
  if ( rstat != 0 ) mf::LogWarning("DataPrepModule") << "Data preparation srvice returned error " << rstat;
  if ( pwires->size() == 0 ) mf::LogWarning("DataPrepModule") << "No wires made for this event.";

  // Build associations.
  std::unique_ptr<art::Assns<raw::RawDigit,recob::Wire>> passns(new art::Assns<raw::RawDigit,recob::Wire>);
  if ( m_DoAssns ) {
    for ( const AdcChannelDataMap::value_type& iacd : acds ) {
      const AdcChannelData& acd = iacd.second;
      AdcIndex idig = acd.digitIndex;
      if ( idig == AdcChannelData::badIndex )
        throw art::Exception(art::errors::InsertFailure) << "Digit index is not set.";
      AdcIndex iwir = acd.wireIndex;
      if ( iwir == AdcChannelData::badIndex ) continue;
      art::Ptr<raw::RawDigit> pdig(hdigits, idig);
      bool success = util::CreateAssn(*this, evt, *pwires, pdig, *passns, m_WireName, iwir);
      if ( !success ) throw art::Exception(art::errors::InsertFailure)
                            << "Can't associate wire " << iwir << " with raw digit " << idig;
    }
  }

  if ( m_LogLevel >= 2 ) {
    cout << myname << "Created wire count: " << pwires->size() << endl;
  }

  // Record wires and associations in the event.
  evt.put(std::move(pwires), m_WireName);
  if ( m_DoAssns ) {
    evt.put(std::move(passns), m_WireName);
  }
  return;
}
