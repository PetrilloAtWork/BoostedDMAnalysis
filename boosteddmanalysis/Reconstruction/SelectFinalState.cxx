/**
 * @file   boosteddmanalysis/Reconstruction/SelectFinalState.cxx
 * @brief  Algorithm selecting the particles of hard scattering final state.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   October 6, 2018
 * @see    boosteddmanalysis/Reconstruction/SelectFinalState.h
 * 
 */

// algorithm header
#include "boosteddmanalysis/Reconstruction/SelectFinalState.h"

// C++ standard library
#include <iterator> // std::back_inserter()
#include <memory> // std::addressof()



//------------------------------------------------------------------------------
bdm::SelectFinalState::Selection_t bdm::SelectFinalState::select
  (std::vector<simb::MCParticle> const& particles) const
{
  static constexpr int kAbsNeutron = 2112; // neutron
  
  Selection_t selected; // indices of the selected particles
  
  auto const nParticles = particles.size();
  for (std::size_t index = 0U; index < nParticles; ++index) {
    auto const& particle = particles[index];
    
    // TODO
    
    // filter neutrons out
    if (!fKeepNeutrons && (std::abs(particle.PdgCode()) == kAbsNeutron))
      continue;
    
    // exclude the particles that are not "primary" (i.e. not generated)
    if (particle.Mother() > 0) continue;
    
    selected.push_back(index);
    
  } // for
  
  return selected;
} // bdm::SelectFinalState::select()


//------------------------------------------------------------------------------
std::vector<simb::MCParticle const*> bdm::SelectFinalState::toMCParticlePointers
(
  std::vector<simb::MCParticle> const& particles,
  Selection_t const& selection
) {
  
  std::vector<simb::MCParticle const*> pointers;
  pointers.reserve(selection.size());
  
  for (auto selectedIndex: selection)
    pointers.push_back(&(particles[selectedIndex]));
  
  return pointers;
} // bdm::SelectFinalState::toMCParticlePointers()


//------------------------------------------------------------------------------
