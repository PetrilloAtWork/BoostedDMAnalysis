/**
 * @file   FilterTauDecayType.cc
 * @brief  A simple filter to select events with specified tau decay types (leptonic or hadronic)
 * @author Yun-Tse Tsai (yuntse@slac.stanford.edu)
 * @date   August 1, 2019
 *
 */


// LArSoft libraries
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"


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
#include "cetlib_except/exception.h"

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

  class FilterTauDecayType : public art::EDFilter
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
      
      fhicl::Atom< bool > vetoMode{
        Name("VetoMode"),
        Comment("whether to veto the events containing the specified decay products.")
      };
      
      fhicl::Sequence< int > SpecifiedParticlePdgCode{
        Name("SpecifiedParticlePdgCode"),
        Comment("the pdgcodes of the specified decay products for the events.")
      };

    }; // Config

    using Parameters = art::EDFilter::Table<Config>;
    
    explicit FilterTauDecayType( Parameters const& config );

    virtual bool filter(art::Event&) override;

  private:
    
    art::InputTag fMCTruthLabel;
    bool fVetoMode;
    std::set< int > fSpecifiedParticles;
    int fNEvents;
    int fNDesiredEvents;
  };

} // namespace bdm

// -------------------------------------------------------------------------
// --- FilterTauDecayType
// ---
bdm::FilterTauDecayType::FilterTauDecayType( Parameters const& config )
  : fMCTruthLabel( config().MCTruthLabel() )
  , fVetoMode( config().vetoMode() )
  , fSpecifiedParticles( ::vectorToSet(config().SpecifiedParticlePdgCode()) )
  , fNEvents( 0 )
  , fNDesiredEvents( 0 )
{
  consumes< std::vector< simb::MCTruth > >( fMCTruthLabel );
}

bool bdm::FilterTauDecayType::filter( art::Event& evt )
{
  
  auto MCTruthHandle = evt.getValidHandle< std::vector< simb::MCTruth > >( fMCTruthLabel );
  auto const& MCTruthObjs = *MCTruthHandle;
  
  bool isDesired = fVetoMode;
  ++fNEvents;

  for ( size_t iMCTruth = 0; iMCTruth < MCTruthObjs.size(); ++iMCTruth ) {
    simb::MCTruth const& MCTruthObj = MCTruthObjs[iMCTruth];
    auto const& lepton = MCTruthObj.GetNeutrino().Lepton();
    auto const leptonPdgCode = lepton.PdgCode();
    
    if ( abs( leptonPdgCode ) != 15 ) continue;
    
    int tauTrackID = lepton.TrackId();
    
    for ( int iMCParticle = 0; iMCParticle < MCTruthObj.NParticles(); ++iMCParticle ) {
      const simb::MCParticle& thisMCParticle = MCTruthObj.GetParticle( iMCParticle );
      if ( thisMCParticle.Mother() != tauTrackID ) continue;
      int thisPdgCode = thisMCParticle.PdgCode();
      // std::cout << "MCParticle " << iMCParticle << " with PDGCode: " << thisPdgCode << std::endl;
      if ( fSpecifiedParticles.count( thisPdgCode ) > 0 ) {
        isDesired = !( fVetoMode );
        break;
      }
    }
    
    if ( isDesired != fVetoMode ) break;
    
  }
  if ( isDesired ) ++fNDesiredEvents;
  std::cout << "Number of events: " << fNEvents << ", number of filtered events: " << fNDesiredEvents
            << std::endl;
  return isDesired;
}

//------------------------------------------------------------------------------
DEFINE_ART_MODULE(bdm::FilterTauDecayType)
