/**
 * @file   FilterGENIELeptonType.cc
 * @brief  A simple filter to select events with specified lepton types (+/-11, +/-12, +/-13, +/-14, etc.) from GENIE
 * @author Yun-Tse Tsai (yuntse@slac.stanford.edu)
 * @date   October 14, 2018
 *
 */


// LArSoft libraries
#include "nusimdata/SimulationBase/MCTruth.h"


// framework libraries
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h" // art::ValidHandle<>
#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"

// C++ libraries
#include <vector>
#include <set>
#include <algorithm> // std::copy()
#include <iterator> // std::inserter()


namespace {
  
  template <typename T>
  std::set<T> vectorToSet(std::vector<T> const& v) {
    std::set<T> s;
    std::copy(v.begin(), v.end(), std::inserter(s, s.begin()));
    return s;
  } // vectorToSet()
  
} // local namespace


namespace bdm {

  class FilterGENIELeptonType : public art::EDFilter
  {
  public:

    /// Algorithm configuration
    struct Config {

      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom< art::InputTag > MCTruthLabel{
        Name("MCTruthLabel"),
        Comment("the producer of the MCTruth used in the filter.")
      };

      fhicl::Sequence< int > LeptonPdgCode{
        Name("LeptonPdgCode"),
        Comment("the lepton pdgcodes for the events to be kept.")
      };

    }; // Config

    using Parameters = art::EDFilter::Table<Config>;
    
    explicit FilterGENIELeptonType( Parameters const& config );

    virtual bool filter(art::Event&) override;

  private:
    
    art::InputTag fMCTruthLabel;
    std::set< int > fKeptLeptons;
  };

} // namespace bdm

// -------------------------------------------------------------------------
// --- FilterGENIELeptonType
// ---
bdm::FilterGENIELeptonType::FilterGENIELeptonType( Parameters const& config )
  : fMCTruthLabel( config().MCTruthLabel() )
  , fKeptLeptons( ::vectorToSet(config().LeptonPdgCode()) )
{
  consumes< std::vector< simb::MCTruth > >( fMCTruthLabel );
}

bool bdm::FilterGENIELeptonType::filter( art::Event& evt )
{
  
  auto MCTruthHandle = evt.getValidHandle< std::vector< simb::MCTruth > >( fMCTruthLabel );
  auto const& MCTruthObjs = *MCTruthHandle;

  bool isDesired = false;

  for (auto const& MCTruthObj: MCTruthObjs) {
    auto const& lepton = MCTruthObj.GetNeutrino().Lepton();
    auto const leptonPdgCode = lepton.PdgCode();
    
    isDesired = fKeptLeptons.count( leptonPdgCode ) > 0;
    if ( isDesired ) break;

  }
  return isDesired;
}

//------------------------------------------------------------------------------
DEFINE_ART_MODULE(bdm::FilterGENIELeptonType)
