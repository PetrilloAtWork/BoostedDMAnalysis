/**
 * @file   boosteddmanalysis/Reconstruction/SmearedReconstructionTableAlg.h
 * @brief  Algorithm reconstructing particles by smearing, where the parameters come from Shirley Li and Alex Friedland's table
 * @author Yun-Tse Tsai (yuntse@slac.stanford.edu)
 * @date   July 26, 2019
 * @see    boosteddmanalysis/Reconstruction/SmearedReconstructionTableAlg.cxx
 * 
 */

#ifndef BOOSTEDDMANALYSIS_RECONSTRUCTION_SMEAREDRECONSTRUCTIONTABLEALG_H
#define BOOSTEDDMANALYSIS_RECONSTRUCTION_SMEAREDRECONSTRUCTIONTABLEALG_H

// Base class
#include "boosteddmanalysis/Reconstruction/SmearedReconstructionAlg.h"

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
#include "cetlib_except/exception.h"

// C/C++ standard libraries
#include <string>
#include <vector>



namespace bdm {
  
  
  /**
   * @brief Reconstructing particles by smearing, where the parameters come from the table
   * 
   * @todo
   * 
   * 
   * If the reconstruction of a input particle fails, the smeared corresponding
   * particle is still returned, but flagged as invalid.
   * 
   * 
   * Example of usage:
   * -------------------------------------------------------------------{.cpp}
   * bdm::SmearedReconstructionTableAlg recoAlg();
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

  class SmearedReconstructionTableAlg: public SmearedReconstructionAlg {
  
  private:
    
    
    class ReconstructionTables {

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
        
        /// Energy smearing from the table
        std::vector< double > fTrueKEs;
        std::vector< double > fMeanDepEs;
        std::vector< double > fSTDDepEs;
        std::vector< double > fEResolutions;
        
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
      ReconstructionTables(std::size_t maxParameters);
      
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
      
    }; // class ReconstructionTables
    
  public:

    /// Constructor: reads the algorithm configuration from a FHiCL table.
    SmearedReconstructionTableAlg(fhicl::ParameterSet const& pset);

    
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
      (std::vector<simb::MCParticle const*> const& particles) const override;

  private:
    /// Reconstruction parameters by particle ID.
    ReconstructionTables fParameters;
    
    ReconstructionTables fillParameters (std::vector<fhicl::ParameterSet> const& params);

  }; // class SmearedReconstructionTableAlg
  
} // namespace bdm

//------------------------------------------------------------------------------
//--- template implementation
//------------------------------------------------------------------------------
template <typename Stream>
void bdm::SmearedReconstructionTableAlg::ReconstructionTables::ParticleParameters_t::dump
(
  Stream&& out,
  std::string indent /* = "" */, std::string firstIndent /* = "" */
) const {
  out << firstIndent
    << "reconstruction efficiency: " << fRecoEff
    << ";  direction: x" << fDirectionSmearingAngle
    << ";  KE threshold: " << fKEThreshold;
    
  for ( size_t iTrueKE = 0; iTrueKE < fTrueKEs.size(); ++iTrueKE ) {
    out << "\n" << indent
    << "true KE: " << fTrueKEs[iTrueKE] << " GeV"
    << "; average deposited energy: " << fMeanDepEs[iTrueKE] << " GeV"
    << "; std of deposited energy: " << fSTDDepEs[iTrueKE] << " GeV"
    << "; energy resolution: " << fEResolutions[iTrueKE]*100. << "%";
  }
  
} // bdm::...::ReconstructionTables::ParticleParameters_t::dump()


//------------------------------------------------------------------------------
template <typename Stream>
void bdm::SmearedReconstructionTableAlg::ReconstructionTables::dump(
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
  
} // bdm::SmearedReconstructionTableAlg::ReconstructionTables::dump()


//------------------------------------------------------------------------------

#endif // BOOSTEDDMANALYSIS_RECONSTRUCTION_SMEAREDRECONSTRUCTIONTABLEALG_H
