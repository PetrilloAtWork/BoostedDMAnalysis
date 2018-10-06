/**
 * @file   boosteddmanalysis/DataObjects/SmearedMCParticle.cxx
 * @brief  Particle obtained from pseudo-reconstruction (implementation).
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   October 5, 2018
 * @see    boosteddmanalysis/DataObjects/SmearedMCParticle.h
 * @ingroup BoostedDarkMatterData
 * 
 */

// library header
#include "boosteddmanalysis/DataObjects/SmearedMCParticle.h"


//------------------------------------------------------------------------------
bdm::SmearedMCParticle::SmearedMCParticle(
  PDGID_t id,
  Momentum_t momentum,
  double energy,
  double time,
  Position_t startPosition,
  Position_t endPosition
  )
  : fEnergy(energy)
  , fMomentum(momentum)
  , fTime(time)
  , fStartingVertex(startPosition)
  , fEndVertex(endPosition)
  , fPDGID(id)
  {}


//------------------------------------------------------------------------------
TParticlePDG const* bdm::SmearedMCParticle::particleInfo() const {
  return hasParticleId()
    ? TDatabasePDG::Instance()->GetParticle(particleId())
    : nullptr;
} // bdm::SmearedMCParticle::particleInfo()


//------------------------------------------------------------------------------
std::string bdm::SmearedMCParticle::particleName() const {
  if (!hasParticleId()) return "<none>";
  auto const* info = particleInfo();
  return info? info->GetName(): "<ID=" + std::to_string(particleId()) + ">";
} // bdm::SmearedMCParticle::particleName()


//------------------------------------------------------------------------------

