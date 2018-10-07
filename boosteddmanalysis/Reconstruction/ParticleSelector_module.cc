/**
 * @file   ParticleSelector_module.cc
 * @brief  Module running `bdm::SelectFinalState` algorithm.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   October 6, 2018
 * 
 * Provides:
 * 
 * * `bdm::ParticleSelector` module
 * 
 */

// LArSoft libraries
#include "boosteddmanalysis/Reconstruction/SelectFinalState.h"
#include "nusimdata/SimulationBase/MCParticle.h"

// framework libraries
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h" // art::ValidHandle<>
#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"

// C/C++ standard libraries
#include <vector>
#include <memory> // std::make_unique()


namespace bdm {
  
  /**
   * @brief Module: applies a selection to generated particles.
   * 
   * This module creates a collection of pointers to simulated particles that
   * survive a selection.
   * 
   * Input
   * ------
   * 
   * A collection of `simb::MCParticle` is required.
   * 
   * 
   * Output
   * ------
   * 
   * A collection of _art_ pointers to some of the input particles is returned.
   * 
   * 
   * Configuration parameters
   * -------------------------
   * 
   * * *particles* (input tag, default: `largeant`): label of the data product
   *     with input simulated particles
   * * *selector* (parameter set): configuration passed to the selection
   *     algorithm
   * 
   */
  class ParticleSelector: public art::EDProducer {
    
      public:
    
    /// Module configuration data.
    struct Config {
      
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;
      
      fhicl::Atom<art::InputTag> particles {
        Name("particles"),
        Comment("the data product of simulated particles to be processed"),
        "largeant" // default
        };
      
      fhicl::Table<bdm::SelectFinalState::Config> selector{
        Name("selector"),
        Comment("configuration of the selector algorithm")
        };
      
    }; // Config
    
    /// Standard _art_ alias for module configuration table.
    using Parameters = art::EDProducer::Table<Config>;
    
    /// Constructor; see the class documentation for the configuration.
    explicit ParticleSelector(Parameters const& config);
    
    
    virtual void produce(art::Event& event) override;
    
    
      private:
    art::InputTag fParticleTag; ///< Label of the input data product.
    
    /// Selection algorithm.
    std::unique_ptr<bdm::SelectFinalState> fSelector;
    
    
  }; // class ParticleSelector
  
} // namespace bdm



//------------------------------------------------------------------------------
//--- ParticleSelector
//--- 
bdm::ParticleSelector::ParticleSelector(Parameters const& config)
  : fParticleTag(config().particles())
  , fSelector(std::make_unique<bdm::SelectFinalState>(config().selector()))
{
  
  consumes<std::vector<simb::MCParticle>>(fParticleTag);
  
  produces<std::vector<art::Ptr<simb::MCParticle>>>();
  
} // bdm::ParticleSelector::ParticleSelector()


//------------------------------------------------------------------------------
void bdm::ParticleSelector::produce(art::Event& event) {
  
  //
  // read the input
  //
  auto particleHandle
    = event.getValidHandle<std::vector<simb::MCParticle>>(fParticleTag);
  
  //
  // prepare the output structures
  //
  art::PtrMaker<simb::MCParticle> makePartPtr(event, particleHandle.id());
  
  //
  // set up the algorithm
  //
  fSelector->setup();
  
  //
  // run the algorithm
  //
  auto const selectedIndices = fSelector->select(*particleHandle);
  
  //
  // convert the result in the chosen output format: collection of art pointers
  //
  auto selectedParticles
    = std::make_unique<std::vector<art::Ptr<simb::MCParticle>>>();
  selectedParticles->reserve(selectedIndices.size());
  for (auto particleIndex: selectedIndices)
    selectedParticles->push_back(makePartPtr(particleIndex));
  
  //
  // store the data products into the event (and print a short summary)
  //
  mf::LogInfo("ParticleSelector")
    << "Selected " << selectedParticles->size() << "/"
    << particleHandle->size() << " particles from '"
    << fParticleTag.encode() << "'";
  
  event.put(std::move(selectedParticles));
  
} // bdm::ParticleSelector::produce()


//------------------------------------------------------------------------------
DEFINE_ART_MODULE(bdm::ParticleSelector)


//------------------------------------------------------------------------------
