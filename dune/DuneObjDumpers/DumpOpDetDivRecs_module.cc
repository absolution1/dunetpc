/**
 * @file   DumpOpDetDivRecs_module.cc
 * @brief  Module dumping OpDetDivRecs information on screen
 * @date   August 8, 2016
 * @author Jason Stock (jason.stock@mines.sdsmt.edu)
 * Module based closely on the DumpOpDetBacktrackerRecords by Lucas Mendes Santos
 *
 */


// nutools libraries
#include "dune/DuneObj/OpDetDivRec.h"

// framework libraries
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
#include "messagefacility/MessageLogger/MessageLogger.h"


namespace sim {
  class DumpOpDetDivRecs;
} // namespace sim

namespace {
  using namespace fhicl;

  /// Collection of configuration parameters for the module
  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;

    fhicl::Atom<art::InputTag> InputOpDetDivRecs {
      Name("InputOpDetDivRecs"),
      Comment("data product with the OpDetDivRecs to be dumped")
      };

    fhicl::Atom<std::string> OutputCategory {
      Name("OutputCategory"),
      Comment("name of the output stream (managed by the message facility)"),
      "DumpOpDetDivRecs" /* default value */
      };

  }; // struct Config


} // local namespace


class sim::DumpOpDetDivRecs: public art::EDAnalyzer {
    public:
  // type to enable module parameters description by art
  using Parameters = art::EDAnalyzer::Table<Config>;

  /// Configuration-checking constructor
  explicit DumpOpDetDivRecs(Parameters const& config);

  // Plugins should not be copied or assigned.
  DumpOpDetDivRecs(DumpOpDetDivRecs const&) = delete;
  DumpOpDetDivRecs(DumpOpDetDivRecs &&) = delete;
  DumpOpDetDivRecs& operator = (DumpOpDetDivRecs const&) = delete;
  DumpOpDetDivRecs& operator = (DumpOpDetDivRecs &&) = delete;


  // Operates on the event
  void analyze(art::Event const& event) override;


  /**
   * @brief Dumps the content of the specified OpDetDivRecs in the output stream
   * @tparam Stream the type of output stream
   * @param out the output stream
   * @param simchannel the OpDetDivRecs to be dumped
   * @param indent base indentation string (default: none)
   * @param bIndentFirst if first output line should be indented (default: yes)
   *
   * The indent string is prepended to every line of output, with the possible
   * exception of the first one, in case bIndentFirst is true.
   *
   * The output starts on the current line, and the last line is NOT broken.
   */
  template <typename Stream>
  void DumpOpDetDivRec(
    Stream&& out, sim::OpDetDivRec const& simchannel,
    std::string indent = "", bool bIndentFirst = true
    ) const;


    private:

  art::InputTag fInputChannels; ///< name of OpDetDivRecs's data product
  std::string fOutputCategory; ///< name of the stream for output

}; // class sim::DumpOpDetDivRecs


//------------------------------------------------------------------------------
//---  module implementation
//---
//------------------------------------------------------------------------------
sim::DumpOpDetDivRecs::DumpOpDetDivRecs(Parameters const& config)
  : EDAnalyzer(config)
  , fInputChannels(config().InputOpDetDivRecs())
  , fOutputCategory(config().OutputCategory())
{}


//------------------------------------------------------------------------------
template <typename Stream>
void sim::DumpOpDetDivRecs::DumpOpDetDivRec(
  Stream&& out, sim::OpDetDivRec const& channel,
  std::string indent /* = "" */, bool bIndentFirst /* = true */
) const {
  if (bIndentFirst) out << indent;
  channel.Dump(out, indent);
} // sim::DumpOpDetDivRecs::DumpOpDetDivRecs()


//------------------------------------------------------------------------------
void sim::DumpOpDetDivRecs::analyze(art::Event const& event) {

  // get the particles from the event
  auto const& OpDetDivRecs
    = *(event.getValidHandle<std::vector<sim::OpDetDivRec>>(fInputChannels));

  mf::LogVerbatim(fOutputCategory) << "Event " << event.id()
    << " : data product '" << fInputChannels.encode() << "' contains "
    << OpDetDivRecs.size() << " OpDetDivRecs";

  unsigned int iOpDetDivRec = 0;
  for (sim::OpDetDivRec const& div_rec: OpDetDivRecs) {

    // a bit of a header
    mf::LogVerbatim log(fOutputCategory);
    log << "[#" << (iOpDetDivRec++) << "] ";
    DumpOpDetDivRec(log, div_rec, "  ", false);

  } // for
  mf::LogVerbatim(fOutputCategory) << "\n";

} // sim::DumpOpDetDivRecs::analyze()


//------------------------------------------------------------------------------
DEFINE_ART_MODULE(sim::DumpOpDetDivRecs)

//------------------------------------------------------------------------------
