/**
 * @file   boosteddmanalysis/Reconstruction/SmearedReconstructionAlg.h
 * @brief  Algorithm reconstructing particles by parametric smearing.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   October 6, 2018
 * @see    boosteddmanalysis/Reconstruction/SmearedReconstructionAlg.cxx
 * 
 */

#ifndef BOOSTEDDMANALYSIS_RECONSTRUCTION_SMEAREDRECONSTRUCTIONALG_H
#define BOOSTEDDMANALYSIS_RECONSTRUCTION_SMEAREDRECONSTRUCTIONALG_H

// BDM libraries
#include "boosteddmanalysis/DataObjects/SmearedMCParticle.h"

// nutools libraries
#include "nusimdata/SimulationBase/MCParticle.h"

// framework libraries
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/ParameterSet.h"

// C/C++ standard libraries
#include <string>
#include <vector>


// forward declarations
namespace CLHEP {
  class HepRandomEngine;
}


namespace bdm {
  
  
  /**
   * @brief Reconstructing particles by parametric smearing.
   * 
   * @todo This algorithm is incomplete!
   * 
   * 
   * If the reconstruction of a input particle fails, the smeared corresponding
   * particle is still returned, but flagged as invalid.
   * 
   * 
   * Example of usage:
   * -------------------------------------------------------------------{.cpp}
   * bdm::SmearedReconstructionAlg recoAlg();
   * recoAlg.setup();
   * 
   * std::vector<bdm::SmearedMCParticle> selection
   *   = recoAlg.reconstruct(mcParticles);
   * -------------------------------------------------------------------
   * (assuming `mcParticles` a collection of `simb::MCParticle`).
   * 
   * 
   * Configuration
   * ==============
   * 
   * * `particleParameters` (list of reconstruction tables):
   *   a list of reconstruction parameters, each one for a specific type of
   *   particles; each table may contain the following keys:
   *     * `id` (list of integers, _mandatory_): list of PDG particle IDs this
   *       table applies to; to specify a particle and its antiparticle, both
   *       need to appear explicitly in the configuration; the PDG ID 0 is
   *       special in that it matches all the particle types that are not
   *       matched otherwise
   *     * `efficiency` (real, default: 1.0): probability for a particle to
   *       be reconstructed at all. If a particle is not reconstructed,
   *       it will be flagged as invalid
   *       (see `bdm::SmearedMCParticle::isValid()`)
   * 
   * Dependencies and setup
   * =======================
   * 
   * Currently none.
   * 
   */
  class SmearedReconstructionAlg {
    
    /// Class holding all the reconstruction parameters for all particle types.
    class ReconstructionParameters {
        public:
      using PDGID_t = bdm::SmearedMCParticle::PDGID_t;
      
        public:
      
      /// Reconstruction parameters for one type of particle.
      struct ParticleParameters_t {
        
        /// Probability of reconstructing the particle (uniform)
        double fRecoEff = 1.0;
        
        /// Transversal smearing of reconstructed direction
        /// [fraction of momentum].
        double fDirectionSmearingAngle = 0.0;
        
        /// Smearing of reconstructed energy [fraction].
        double fEnergySmearingFraction = 0.0;
        
        /// Detection threshold by kinetic energy [GeV].
        double fKEThreshold = 0.0;
        
        /// Prints into the specified stream the content of these parameters.
        template <typename Stream>
        void dump(
          Stream&& out,
          std::string indent = "", std::string firstIndent = ""
          ) const;
        
      }; // ParticleParameters_t
      
      
      /// Initializes the parameter registry to host `maxParameters` tables.
      ReconstructionParameters(std::size_t maxParameters);
      
      //@{
      /// Associates the parameters in `params` to particles types in `ids`.
      void registerParams
        (ParticleParameters_t const& params, std::vector<PDGID_t> const& ids);
      void registerParams
        (ParticleParameters_t&& params, std::vector<PDGID_t> const& ids);
      //@}
      
      
      /**
       * @brief Returns the parameters pertaining the specified particle ID.
       * @param id the PDG id to be queried
       * @return the parameters pertaining the specified particle ID
       * @throw std::out_of_range if no parameter is present for that ID
       * 
       * Note that if parameters for a PDG ID `0` is registered, those are
       * returned if no parameter for `id` is present.
       */
      ParticleParameters_t const& operator[] (PDGID_t id) const;
      
      
      /// Prints into the specified stream the content of all the parameters.
      template <typename Stream>
      void dump(
        Stream&& out,
        std::string indent = "", std::string firstIndent = ""
        ) const;
      
        private:
      std::vector<ParticleParameters_t> fParameters;
      std::map<PDGID_t, ParticleParameters_t*> fAssociations;
      
      /// Associates a parameter record with all specified IDs
      void associateParams
        (ParticleParameters_t& params, std::vector<PDGID_t> const& ids);
      
    }; // class ReconstructionParameters
    
    
      public:
    
    /// Type of collection of reconstructed particles returned by the algorithm.
    using SmearedParticles_t = std::vector<bdm::SmearedMCParticle>;
    
    /// Type of PDG ID.
    using PDGID_t = bdm::SmearedMCParticle::PDGID_t;
    
    
#ifdef USE_FHICLCPP_VALIDATION
    struct ParticleParameterConfig {
      
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      
      fhicl::Sequence<PDGID_t> id {
        Name("id"),
        Comment("PDG IDs these parameters apply to")
        };
      
      fhicl::Atom<bool> efficiency {
        Name("efficiency"),
        Comment("probability of reconstructing the particle (uniform)"),
        1.0
        };
      
      fhicl::Atom<double> direction {
        Name("direction"),
        Comment("transversal displacement of reconstructed momentum [fraction]"),
        0.0
        };

      fhicl::Atom<double> energy {
        Name("energy"),
        Comment("energy resolution [fraction]"),
        0.0
      };
      fhicl::Atom<double> threshold {
        Name("threshold"),
        Comment("detection threshold in kinetic energy [GeV]"),
        0.0
      };
        
    }; // struct ParticleParameterConfig
    
    
    struct Config {
      
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      
      fhicl::Sequence<ParticleParameterConfig> particleParameters{
        Name("particleParameters"),
        Comment("select neutrons")
        };
      
    }; // Config
    
    
    using Parameters = fhicl::Table<Config>;
    
    /// Constructor: reads the algorithm configuration.
    SmearedReconstructionAlg(Config const& config);
    
    /// Constructor: reads the algorithm configuration from a FHiCL table.
    SmearedReconstructionAlg(Parameters const& config)
      : SmearedReconstructionAlg(config()) {}
    
#else // USE_FHICLCPP_VALIDATION
    
    
    /// Constructor: reads the algorithm configuration.
    SmearedReconstructionAlg(fhicl::ParameterSet const& pset);
    
    
    
#endif // USE_FHICLCPP_VALIDATION
    
    /// Destructor (does nothing, but it does it virtually).
    virtual ~SmearedReconstructionAlg() = default;
    
    
    /**
     * @brief Set up the algorithm.
     * @param rndEng random engine used for efficiency and smearing
     */
    virtual void setup(CLHEP::HepRandomEngine& rndEng);
    
    
    /**
     * @brief Returns the set of selected particles.
     * @param particles the simulated particles to apply a selection on
     * @return a collection of pointers to the selected ones among the particles
     * 
     * Each of the particles in the specified list is reconstructed by applying
     * parametrized smearing to the true physical characteristics.
     * 
     * Particles are specified as a vector of pointers.
     */
    virtual SmearedParticles_t reconstruct
      (std::vector<simb::MCParticle const*> const& particles) const;
    
    
    /**
     * @brief Returns the set of selected particles.
     * @param particles the simulated particles to apply a selection on
     * @return a collection of pointers to the selected ones among the particles
     * @see reconstruct(std::vector<simb::MCParticle const*> const&) const
     * 
     * Reconstruction is applied as described in
     * `reconstruct(std::vector<simb::MCParticle const*> const&) const`.
     * 
     * Particles are specified as a vector of particles.
     */
    SmearedParticles_t reconstruct
      (std::vector<simb::MCParticle> const& particles) const
      { return reconstruct(toMCParticlePointers(particles)); }
    
    
    // --- BEGIN Conversion functions ------------------------------------------
    /// @name Conversion functions
    /// @{
    /**
     * @brief Converts a selection into a list of C pointers to particles.
     * @param particles the _original_, complete particle list
     * @param selection indices of the selected particles
     * @return a vector of pointers to the selected particles
     */
    static std::vector<simb::MCParticle const*> toMCParticlePointers
      (std::vector<simb::MCParticle> const& particles);
    
    /**
     * @brief Converts a selection into a list of C pointers to particles.
     * @param particles the _original_, complete particle list
     * @param selection indices of the selected particles
     * @return a vector of pointers to the selected particles
     */
    static std::vector<simb::MCParticle const*> toMCParticlePointers
      (art::PtrVector<simb::MCParticle> const& particles);
    
    /**
     * @brief Converts a selection into a list of C pointers to particles.
     * @param particles the _original_, complete particle list
     * @param selection indices of the selected particles
     * @return a vector of pointers to the selected particles
     */
    static std::vector<simb::MCParticle const*> toMCParticlePointers
      (std::vector<art::Ptr<simb::MCParticle>> const& particles);
    
    /// @}
    // --- END Conversion functions --------------------------------------------
    
      private:
    
    /// Reconstruction parameters by particle ID.
    ReconstructionParameters fParameters;
    
    /// Random engine.
    CLHEP::HepRandomEngine* fRndEngine = nullptr;
    
    /// Parses configuration and fills the reconstruction parameter register.
#ifdef USE_FHICLCPP_VALIDATION
    static ReconstructionParameters fillParameters
      (std::vector<ParticleParameterConfig> const& params);
#else // USE_FHICLCPP_VALIDATION
    static ReconstructionParameters fillParameters
      (std::vector<fhicl::ParameterSet> const& params);
#endif // USE_FHICLCPP_VALIDATION
    
    
    
    //
    // It is important to remember to create the generator with a reference to
    // the engine rather than the pointer, or else CLHEP will grab ownership;
    // this function helps with that, e.g.:
    // 
    // auto flat = makeRandomGen<CLHEP::RandFlat>(-1.0, 1.0);
    //
    template <typename Generator, typename... Args>
    Generator makeRandomGen(Args&&... args) const
      { return { *fRndEngine, std::forward<Args>(args)... }; }
    
  }; // class SmearedReconstructionAlg
  
  
} // namespace bdm


//------------------------------------------------------------------------------
//--- template implementation
//------------------------------------------------------------------------------
template <typename Stream>
void bdm::SmearedReconstructionAlg::ReconstructionParameters::ParticleParameters_t::dump
(
  Stream&& out,
  std::string /* indent */ /* = "" */, std::string firstIndent /* = "" */
) const {
  out << firstIndent
    << "reconstruction efficiency: " << fRecoEff
    << ";  direction: x" << fDirectionSmearingAngle
    << ";  energy: +/- x" << fEnergySmearingFraction;
  
} // bdm::...::ReconstructionParameters::ParticleParameters_t::dump()


//------------------------------------------------------------------------------
template <typename Stream>
void bdm::SmearedReconstructionAlg::ReconstructionParameters::dump(
  Stream&& out,
  std::string indent /* = "" */, std::string firstIndent /* = "" */
) const {
  
  // find the default, if any
  auto iDefParams = fAssociations.find(0);
  ParticleParameters_t const* defParams
    = (iDefParams == fAssociations.end())? nullptr: iDefParams->second;
  
  out << firstIndent << fParameters.size() << " parameter sets covering "
    << (defParams? fAssociations.size() - 1: fAssociations.size())
    << " particles:";
  
  for (auto const& params: fParameters) {
    out << "\n" << indent << "particles:";
    unsigned int nFound = 0;
    for (auto const& p: fAssociations) { // p: { PDG ID, parameter pointer }
      if (p.second != &params) continue;
      
      auto const id = p.first;
      if (id == 0) continue;
      
      ++nFound;
      auto const* pPDGinfo = TDatabasePDG::Instance()->GetParticle(id);
      if (pPDGinfo) out << " " << pPDGinfo->GetName() << " (ID=" << id << ")";
      else          out << " <ID=" << id << ">";
    } // for p
    if (&params == defParams) out << " (default)";
    else if (nFound == 0) out << " none!";
    out << "\n";
    params.dump(std::forward<Stream>(out), indent, indent + "  ");
  } // for
  
} // bdm::SmearedReconstructionAlg::ReconstructionParameters::dump()


//------------------------------------------------------------------------------


#endif // BOOSTEDDMANALYSIS_RECONSTRUCTION_SMEAREDRECONSTRUCTIONALG_H
