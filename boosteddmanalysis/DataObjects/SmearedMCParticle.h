/**
 * @file   boosteddmanalysis/DataObjects/SmearedMCParticle.h
 * @brief  Particle obtained from pseudo-reconstruction.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   October 5, 2018
 * @version 10
 * @see    boosteddmanalysis/DataObjects/SmearedMCParticle.h
 * @ingroup BoostedDarkMatterData
 * 
 */

#ifndef BOOSTEDDMANALYSIS_DATAOBJECTS_SMEAREDMCPARTICLE_H
#define BOOSTEDDMANALYSIS_DATAOBJECTS_SMEAREDMCPARTICLE_H


// LArSoft libraries
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h" // geo::Point_t...
#include "Math/GenVector/LorentzVector.h"


// ROOT libraries
#include "TDatabasePDG.h"
#include "TParticlePDG.h"

// C/C++ standard libraries
#include <ostream>
#include <string>
#include <utility> // std::move(), std::forward()

/**
 * @brief Namespace for Boosted Dark Matter analysis code.
 * @addtogroup BoostedDarkMatterData
 */
namespace bdm {
  
  // BEGIN BoostedDarkMatterData group -----------------------------------------
  /// @ingroup BoostedDarkMatterData
  /// @{
  
  /**
    * @brief Particle obtained from pseudo-reconstruction.
    * 
    * This class represents a simulated particle which undergo a parametrised
    * simulation and reconstruction.
    * 
    */
  class SmearedMCParticle {
    
      public:
    
    using Position_t = geo::Point_t; ///< Type to represent positions.
    using Momentum_t = geo::Vector_t; ///< Type to represent momenta.
    using PDGID_t = int; ///< Type of the particle ID.
    
    /// 4-momentum
    using FourMomentum_t
      = ROOT::Math::LorentzVector<ROOT::Math::Cartesian3D<double>>;
    
    /// Value of a particle ID that denotes it as invalid.
    static constexpr PDGID_t InvalidParticleID = 0;
    
    
    /// Default constructor, only for ROOT I/O (do not use it!).
    SmearedMCParticle() = default;
    
    /**
      * @brief Constructor from trajectory (stolen) and particle ID.
      * @param id type of particle, in PDG ID standard
      * @param momentum momentum vector [GeV/c]
      * @param energy energy [GeV]
      * @param time time of appearance
      * @param startPosition initial position [cm]
      * @param endPosition final position [cm]
      * 
      * The particle is initialized with the specified values.
      */
    SmearedMCParticle(
      PDGID_t id,
      Momentum_t momentum, double energy,
      double time, Position_t startPosition, Position_t endPosition
      );
    
    
    // --- BEGIN access to data ------------------------------------------------
    /// @{
    /// @name Access to data
    
    //@{
    /// Returns the initial momentum of the particle [GeV/c].
    Momentum_t const& momentum() const { return fMomentum; }
    Momentum_t const& cp() const { return momentum(); }
    //@}
    
    //@{
    /// Returns the initial energy of the particle [GeV].
    double energy() const { return fEnergy; }
    double E() const { return energy(); }
    //@}
    
    /// Returns the position where the particle was created, in world
    /// coordinates [cm].
    Position_t const& startPosition() const { return fStartingVertex; }
    
    /// Returns the position where the particle ends, in world coordinates [cm].
    Position_t const& endPosition() const { return fEndVertex; }
    
    //@{
    /// Returns reconstructed time of particle [ns] since hardware trigger time.
    double time() const { return fTime; }
    double t() const { return fTime; }
    //@}
    
    
    /// Returns whether the particle ID is valid.
    bool hasParticleId() const { return particleId() != InvalidParticleID; }
    
    /// Returns the particle ID, in PDG standard.
    PDGID_t particleId() const { return fPDGID; }
    
    /**
     * @brief Returns particle information from a database.
     * @return a pointer to the information (TParticlePDG) or `nullptr` if none
     * 
     * The information is provided by ROOT's particle database.
     * 
     * If the particle is not known to ROOT, or if there is no particle ID
     * (`!hasParticleId()`), a null pointer is returned instead.
     */
    TParticlePDG const* particleInfo() const;
    
    /// Returns a string with the name of the type of particle.
    std::string particleName() const;
    
    /// @}
    // --- END access to data --------------------------------------------------
    
    
    // --- BEGIN printing data -------------------------------------------------
    /// @{
    /// @name Printing data
    
    /// Default verbosity level.
    static constexpr unsigned int DefaultDumpVerbosity = 1U;
    
    /// Maximum verbosity level.
    static constexpr unsigned int MaxDumpVerbosity = 3U;
    
    
    //@{
    /**
      * @brief Prints the content of this object into an output stream.
      * @tparam Stream type of the output text stream
      * @param out the output text stream
      * @param verbosity the amount of information printed
      *                  (_default: `DefaultDumpVerbosity`_)
      * @param indent indentation string for all output except the first line
      *               (_default: none_)
      * @param firstIndent indentation string for the first line
      *                    (_default: as `indent`_)
      * 
      * Verbosity levels:
      * * 0: ID (compact), energy and momentum
      * * 1: ID in expanded form, starting position, energy and momentum
      * * 2: added end position (`endPosition()`)
      * * 3: added particle time (`time()`)
      * 
      */
    template <typename Stream>
    void dump(
      Stream&& out, unsigned int verbosity,
      std::string indent, std::string firstIndent
      ) const;
    template <typename Stream>
    void dump(
      Stream&& out, unsigned int verbosity = DefaultDumpVerbosity,
      std::string indent = ""
      ) const
      { dump(std::forward<Stream>(out), verbosity, indent, indent); }
    //@}
    
    ///@}
    
    // --- END printing data ---------------------------------------------------
    
      private:
    
    double     fEnergy = 0.0;   ///< Reconstructed energy of the particle [GeV].
    Momentum_t fMomentum; ///< Reconstructed momentum of the particle [GeV/c].
    double     fTime = 0.0;           ///< Reconstructed creation time [ns].
    Position_t fStartingVertex; ///< Reconstructed starting position [cm].
    Position_t fEndVertex;      ///< Reconstructed ending position [cm].
    
    PDGID_t fPDGID = InvalidParticleID; ///< Particle ID in PDG standard.
    
  }; // class SmearedMCParticle
  
  
  
  /// Prints the content of the track into a text stream.
  /// @related bdm::SmearedMCParticle
  /// @ingroup BoostedDarkMatterData
  inline std::ostream& operator<<
    (std::ostream& out, bdm::SmearedMCParticle const& particle)
    { particle.dump(out); return out; }
  
  /// @}
  // END BoostedDarkMatterData group -------------------------------------------
  
} // namespace bdm


//------------------------------------------------------------------------------
//---  template implementation
//------------------------------------------------------------------------------
template <typename Stream>
void bdm::SmearedMCParticle::dump(
  Stream&& out, unsigned int verbosity,
  std::string /* indent */, std::string firstIndent
  ) const
{
  // we could use ROOT's TDatabasePDG to get the name of the ID, but we'd rather
  // not depend on ROOT here...
  if (verbosity >= 1) {
    out << firstIndent
      << "particle: ";
    auto const* pPDGinfo = particleInfo();
    if (pPDGinfo) out << pPDGinfo->GetName() << " (ID=" << particleId() << ")";
    else          out << "ID " << particleId();
  }
  else out << firstIndent << "ID=" << particleId(); // <== verbosity: 0
  
  if (verbosity >= 1) out << " starting";
  if (verbosity >= 3) out << " on " << time() << " ns";
  if (verbosity >= 1) out << " at " << startPosition() << " cm";
  if (verbosity >= 2) out << " and ending at " << endPosition() << " cm";
  
  // verbosity: 0
  out << "; E=" << energy() << " GeV, cp=" << momentum() << " GeV/c";
  
} // lar::example::CheatTrack::dump()

//------------------------------------------------------------------------------


#endif // BOOSTEDDMANALYSIS_DATAOBJECTS_SMEAREDMCPARTICLE_H
