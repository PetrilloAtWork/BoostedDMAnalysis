/**
 * @file   DumpSmearedMCParticles_module.cc
 * @brief  Dumps the content of `bdm::SmearedMCParticle` data product.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   October 6, 2018
 * 
 * Provides:
 * 
 * * `bdm::DumpSmearedMCParticles` module
 * 
 */

// Boosted Dark Matter libraries
#include "boosteddmanalysis/DataObjects/SmearedMCParticle.h"

// framework libraries
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h" // art::ValidHandle<>
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/Atom.h"

// C/C++ standard libraries
#include <vector>


namespace bdm {
  
  /**
   * @brief Module: dumps the content of `bdm::SmearedMCParticle` data product.
   * 
   * Input
   * ------
   * 
   * A collection of `bdm::SmearedMCParticle` is required.
   * 
   * 
   * Output
   * ------
   * 
   * The dump is sent to the messagefacility service.
   * 
   * 
   * Configuration parameters
   * -------------------------
   * 
   * See `lar --print-description DumpSmearedMCParticles`.
   * 
   */
  class DumpSmearedMCParticles: public art::EDAnalyzer {
    
      public:
    
    /// Collection of configuration parameters for the module
    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      
      fhicl::Atom<art::InputTag> InputParticles {
        Name("InputParticles"),
        Comment("data product with the particles to be dumped")
        };
      
      fhicl::Atom<std::string> OutputCategory {
        Name("OutputCategory"),
        Comment("name of the output stream (managed by the message facility)"),
        "DumpSmearedMCParticles" /* default value */
        };
      
    }; // struct Config
    
    
    /// Standard _art_ alias for module configuration table.
    using Parameters = art::EDAnalyzer::Table<Config>;
    
    /// Constructor.
    explicit DumpSmearedMCParticles(Parameters const& config);
    
    virtual void analyze(art::Event const& event) override;
    
      private:
    
    art::InputTag fInputParticles; ///< Name of SmearedMCParticle data product.
    std::string fOutputCategory; ///< Name of the stream for output.
    
    
  }; // class DumpSmearedMCParticles
  
} // namespace bdm



//------------------------------------------------------------------------------
//--- DumpSmearedMCParticles
//--- 
bdm::DumpSmearedMCParticles::DumpSmearedMCParticles(Parameters const& config)
  : EDAnalyzer(config)
  , fInputParticles(config().InputParticles())
  , fOutputCategory(config().OutputCategory())
{
  consumes<std::vector<bdm::SmearedMCParticle>>(fInputParticles);
}


//------------------------------------------------------------------------------
void bdm::DumpSmearedMCParticles::analyze(art::Event const& event) {
  
  // get the particles from the event
  auto const& Particles = *(
    event.getValidHandle<std::vector<bdm::SmearedMCParticle>>(fInputParticles)
    );
  
  mf::LogVerbatim(fOutputCategory) << "Event " << event.id()
    << ": data product '" << fInputParticles.encode() << "' contains "
    << Particles.size() << " SmearedMCParticle's";
  
  unsigned int iParticle = 0;
  for (auto const& particle: Particles) {
    // flush on every particle,
    // since the output buffer might grow too large otherwise
    mf::LogVerbatim log(fOutputCategory);
    
    // a bit of a header
    log << "[#" << (iParticle++) << "] ";
    particle.dump(log, bdm::SmearedMCParticle::MaxDumpVerbosity, "  ", "");
  } // for
  
  mf::LogVerbatim(fOutputCategory) << "\n";
  
} // bdm::DumpSmearedMCParticles::analyze()


//------------------------------------------------------------------------------
DEFINE_ART_MODULE(bdm::DumpSmearedMCParticles)


//------------------------------------------------------------------------------
