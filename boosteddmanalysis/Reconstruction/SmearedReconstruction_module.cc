/**
 * @file   SmearedReconstruction_module.cc
 * @brief  Module running `bdm::SmearedReconstructionAlg` algorithm.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   October 6, 2018
 * 
 * Provides:
 * 
 * * `bdm::SmearedReconstruction` module
 * 
 */

// Boosted Dark Matter libraries
#include "boosteddmanalysis/Reconstruction/SmearedReconstructionAlg.h"
#include "boosteddmanalysis/DataObjects/SmearedMCParticle.h"

// LArSoft libraries
#include "nutools/RandomUtils/NuRandomService.h"
#include "nusimdata/SimulationBase/MCParticle.h"

// framework libraries
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h" // art::ValidHandle<>
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"

// CLHEP libraries
#include "CLHEP/Random/RandomEngine.h" // CLHEP::HepRandomEngine


// C/C++ standard libraries
#include <vector>
#include <algorithm> // std::count_if()
#include <memory> // std::make_unique()


namespace bdm {
  
  /**
   * @brief Module: reconstruct particles by means of parametrized smearing.
   * 
   * Input
   * ------
   * 
   * A collection of `simb::MCParticle` or `art::Ptr<simb::MCParticle>` is
   * required.
   * 
   * 
   * Output
   * ------
   * 
   * A collection of `bdm::SmearedMCParticle` is produced, one per input
   * particle.
   * 
   * 
   * Configuration parameters
   * -------------------------
   * 
   * * *particles* (input tag, default: `largeant`): label of the data product
   *     with input simulated particles
   * * *reconstruction* (parameter set): configuration passed to the
   *     reconstruction algorithm
   * 
   */
  class SmearedReconstruction: public art::EDProducer {
    
      public:
    
#ifdef USE_FHICLCPP_VALIDATION
    /// Module configuration data.
    struct Config {
      
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;
      
      fhicl::Atom<art::InputTag> particles {
        Name("particles"),
        Comment("the data product of simulated particles to be processed"),
        "largeant" // default
        };
      
      fhicl::Atom<fhicl::ParameterSet> reconstruction {
        Name("reconstruction"),
        Comment("configuration of the reconstruction algorithm")
        };
      
    }; // Config
    
    /// Standard _art_ alias for module configuration table.
    using Parameters = art::EDProducer::Table<Config>;
    
    /// Constructor; see the class documentation for the configuration.
    explicit SmearedReconstruction(Parameters const& config);
    
#else // USE_FHICLCPP_VALIDATION
    
    explicit SmearedReconstruction(fhicl::ParameterSet const& pset);
    
#endif // USE_FHICLCPP_VALIDATION

    
    virtual void produce(art::Event& event) override;
    
    
      private:
    art::InputTag fParticleTag; ///< Label of the input data product.
    
    CLHEP::HepRandomEngine* fRandomEngine = nullptr;
    
    /// Selection algorithm.
    std::unique_ptr<bdm::SmearedReconstructionAlg> fSmearAlg;
    
    /// Retrieves the input particles.
    std::vector<simb::MCParticle const*> fetchParticles
      (art::Event const& event) const;
    
    
  }; // class SmearedReconstruction
  
} // namespace bdm



//------------------------------------------------------------------------------
//--- SmearedReconstruction
//--- 
#ifdef USE_FHICLCPP_VALIDATION
bdm::SmearedReconstruction::SmearedReconstruction(Parameters const& config)
  : fParticleTag(config().particles())
  , fSmearAlg(std::make_unique<bdm::SmearedReconstructionAlg>
      (config().reconstruction())
      )
#else // USE_FHICLCPP_VALIDATION
bdm::SmearedReconstruction::SmearedReconstruction
  (fhicl::ParameterSet const& pset)
  : fParticleTag(pset.get<art::InputTag>("particles"))
  , fSmearAlg(std::make_unique<bdm::SmearedReconstructionAlg>
      (pset.get<fhicl::ParameterSet>("reconstruction"))
      )
#endif // USE_FHICLCPP_VALIDATION
{
  
  mayConsume<std::vector<simb::MCParticle>>(fParticleTag);
  mayConsume<std::vector<art::Ptr<simb::MCParticle>>>(fParticleTag);
  
  produces<std::vector<bdm::SmearedMCParticle>>();
  
  art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this);
  fRandomEngine
    = &(art::ServiceHandle<art::RandomNumberGenerator>()->getEngine());
  
} // bdm::SmearedReconstruction::SmearedReconstruction()


//------------------------------------------------------------------------------
std::vector<simb::MCParticle const*>
bdm::SmearedReconstruction::fetchParticles(art::Event const& event) const
{
  
  std::vector<simb::MCParticle const*> particles;
  
  art::Handle<std::vector<simb::MCParticle>> particleVec;
  if (event.getByLabel(fParticleTag, particleVec)) {
    return bdm::SmearedReconstructionAlg::toMCParticlePointers(*particleVec);
  }
  
  art::Handle<std::vector<art::Ptr<simb::MCParticle>>> particlePtrVec;
  if (event.getByLabel(fParticleTag, particlePtrVec)) {
    return bdm::SmearedReconstructionAlg::toMCParticlePointers(*particlePtrVec);
  }
  
  throw art::Exception(art::errors::ProductNotFound)
    << "Unable to find smeared particles data product '"
    << fParticleTag.encode() << "'\n"
    ;
  
} // bdm::SmearedReconstruction::fetchParticles()


//------------------------------------------------------------------------------
void bdm::SmearedReconstruction::produce(art::Event& event) {
  
  //
  // read the input
  //
  auto particles = fetchParticles(event);
  
  //
  // set up the algorithm
  //
  fSmearAlg->setup(*fRandomEngine);
  
  //
  // run the algorithm and store the result
  //
  auto reconstructedParticles
    = std::make_unique<std::vector<bdm::SmearedMCParticle>>
    (fSmearAlg->reconstruct(particles));
  
  unsigned int const nFailures = std::count_if(
    reconstructedParticles->cbegin(),
    reconstructedParticles->cend(),
    std::mem_fn(&bdm::SmearedMCParticle::isValid)
    );
  
  //
  // store the data products into the event (and print a short summary)
  //
  mf::LogInfo("SmearedReconstruction")
    << "Reconstructed " << nFailures << "/"
    << particles.size() << " particles from '"
    << fParticleTag.encode() << "'";
  
  event.put(std::move(reconstructedParticles));
  
} // bdm::SmearedReconstruction::produce()


//------------------------------------------------------------------------------
DEFINE_ART_MODULE(bdm::SmearedReconstruction)


//------------------------------------------------------------------------------
