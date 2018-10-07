/**
 * @file   boosteddmanalysis/Reconstruction/SelectFinalState.h
 * @brief  Algorithm selecting the particles of hard scattering final state.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   October 6, 2018
 * @see    boosteddmanalysis/Reconstruction/SelectFinalState.cxx
 * 
 */

#ifndef BOOSTEDDMANALYSIS_RECONSTRUCTION_SELECTFINALSTATE_H
#define BOOSTEDDMANALYSIS_RECONSTRUCTION_SELECTFINALSTATE_H

// BDM libraries
// #include "boosteddmanalysis/DataObjects/SmearedMCParticle.h"

// nutools libraries
#include "nusimdata/SimulationBase/MCParticle.h"

// framework libraries
#include "canvas/Persistency/Common/PtrVector.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"



namespace bdm {
  
  
  /**
   * @brief Select the particles of hard scattering final state.
   * 
   * @todo This algorithm is incomplete!
   * 
   * Example of usage:
   * -------------------------------------------------------------------{.cpp}
   * bdm::SelectFinalState selectorAlg();
   * selectorAlg.setup();
   * 
   * std::vector<simb::MCParticle const*> selection
   *   = selectorAlg.select(mcParticles);
   * -------------------------------------------------------------------
   * (assuming `mcParticles` a collection of `simb::MCParticle`).
   * 
   * 
   * Configuration
   * ==============
   * 
   * * `keepNeutrons` (boolean, default: `true`): neutrons in the final state
   *     are kept
   * 
   * 
   * Dependencies and setup
   * =======================
   * 
   * Currently none.
   * 
   */
  class SelectFinalState {
    
      public:
    
    /// Type of collection of selected particles returned by the algorithm.
    using Selection_t = std::vector<std::size_t>;
    
    struct Config {
      
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      
      fhicl::Atom<bool> keepNeutrons{
        Name("keepNeutrons"),
        Comment("select neutrons"),
        true
        };
      
    }; // Config
    
    
    using Parameters = fhicl::Table<Config>;
    
    /// Constructor: reads the algorithm configuration.
    SelectFinalState(Config const& config)
      : fKeepNeutrons(config.keepNeutrons())
      {}
    
    /// Constructor: reads the algorithm configuration from a FHiCL table.
    SelectFinalState(Parameters const& config): SelectFinalState(config()) {}
    
    /// Destructor (does nothing, but it does it virtually).
    virtual ~SelectFinalState() = default;
    
    
    /// Set up the algorithm (currently no operation).
    virtual void setup() {}
    
    
    /**
     * @brief Returns the set of selected particles.
     * @param particles the simulated particles to apply a selection on
     * @return a collection of pointers to the selected ones among the particles
     * 
     * The returned list contains the indices in the original list of the
     * selected particles.
     */
    virtual Selection_t select
      (std::vector<simb::MCParticle> const& particles) const;
    
    
    // --- BEGIN Conversion functions ------------------------------------------
    /// @name Conversion functions
    /// @{
    /**
     * @brief Converts a selection into a list of C pointers to particles.
     * @param particles the _original_, complete particle list
     * @param selection indices of the selected particles
     * @return a vector of pointers to the selected particles
     */
    static std::vector<simb::MCParticle const*> toMCParticlePointers(
      std::vector<simb::MCParticle> const& particles,
      Selection_t const& selection
      );
    
    /// @}
    // --- END Conversion functions --------------------------------------------
    
      private:
    
    bool fKeepNeutrons = false; ///< Whether to keep neutrons in the selection.
    
  }; // class SelectFinalState
  
  
} // namespace bdm



#endif // BOOSTEDDMANALYSIS_RECONSTRUCTION_SELECTFINALSTATE_H
