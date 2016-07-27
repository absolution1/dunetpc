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
  std::string m_DigitLabel;  ///< Full label for the input digit container, e.g. daq:

  // Split label into producer and name: PRODUCER or PRODUCER:NAME
  std::string m_DigitProducer;
  std::string m_DigitName;
    
};

DEFINE_ART_MODULE(DataPrepModule)
  
//**********************************************************************

DataPrepModule::DataPrepModule(fhicl::ParameterSet const& pset) {
  this->reconfigure(pset);
  produces< std::vector<recob::Wire> >(m_DigitName);
  produces<art::Assns<raw::RawDigit, recob::Wire>>(m_DigitName);
}
  
//**********************************************************************

DataPrepModule::~DataPrepModule() { }

//**********************************************************************

void DataPrepModule::reconfigure(fhicl::ParameterSet const& pset) {

  m_DigitLabel = pset.get< std::string >("DigitModuleLabel", "daq");
               
  size_t ipos = m_DigitLabel.find(":");
  if ( ipos == std::string::npos ) {
    m_DigitProducer = m_DigitLabel;
  } else {
    m_DigitProducer = m_DigitLabel.substr(0, ipos);
    m_DigitName = m_DigitLabel.substr(ipos + 1);
  }
}

//**********************************************************************

void DataPrepModule::beginJob() { }

//**********************************************************************

void DataPrepModule::endJob() { }
  
//**********************************************************************

void DataPrepModule::produce(art::Event& evt) {      

  // Read in the digits. 
  art::Handle< std::vector<raw::RawDigit> > hdigits;
  evt.getByLabel(m_DigitProducer, m_DigitName, hdigits);
  if ( hdigits->size() == 0 ) mf::LogWarning("DataPrepModule") << "Input digit container is empty";

  // Process the digits.
  ServiceHandle<RawDigitPrepService> hrdp;
  AdcChannelDataMap acds;
  std::unique_ptr<std::vector<recob::Wire>> pwires(new std::vector<recob::Wire>);
  int rstat = hrdp->prepare(*hdigits, acds, pwires.get());
  if ( rstat != 0 ) mf::LogWarning("DataPrepModule") << "Data preparation srvice returned error " << rstat;
  if ( pwires->size() == 0 ) mf::LogWarning("DataPrepModule") << "No wires made for this event.";

  // Build associations.
  std::unique_ptr<art::Assns<raw::RawDigit,recob::Wire>>
    passns(new art::Assns<raw::RawDigit,recob::Wire>);
  for ( const AdcChannelDataMap::value_type& iacd : acds ) {
    const AdcChannelData& acd = iacd.second;
    AdcIndex idig = acd.digitIndex;
    if ( idig == AdcChannelData::badIndex )
      throw art::Exception(art::errors::InsertFailure) << "Digit index is not set.";
    AdcIndex iwir = acd.wireIndex;
    if ( iwir == AdcChannelData::badIndex )
      throw art::Exception(art::errors::InsertFailure) << "Wire index is not set.";
    art::Ptr<raw::RawDigit> pdig(hdigits, idig);
    int astat = util::CreateAssn(*this, evt, *pwires, pdig, *passns, m_DigitName, iwir);
    if ( astat != 0 ) throw art::Exception(art::errors::InsertFailure)
                           << "Can't associate wire " << iwir << " with raw digit " << idig;
  }

  // Record wires and associations in the event.
  evt.put(std::move(pwires), m_DigitName);
  evt.put(std::move(passns), m_DigitName);

  return;
}
