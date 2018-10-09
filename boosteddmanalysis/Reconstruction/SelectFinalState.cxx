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
#include <map>



//------------------------------------------------------------------------------
bdm::SelectFinalState::Selection_t bdm::SelectFinalState::select
  (std::vector<simb::MCParticle> const& particles) const
{
  static constexpr int kAbsNeutron = 2112; // neutron
  
  Selection_t selected; // indices of the selected particles
  
  std::map<int, std::size_t> ParticleIndices; // ( particle ID , particle index)
  for (std::size_t index = 0; index < particles.size(); ++index)
    ParticleIndices[particles[index].TrackId()] = index;
  
  std::set<std::size_t> InterestingParticles;
  
  auto const nParticles = particles.size();
  for (std::size_t index = 0U; index < nParticles; ++index) {
    auto const& particle = particles[index];
    
    bool isInteresting = (InterestingParticles.count(index) > 0);
    
    if (particle.Mother() == 0) isInteresting = true;
    
    if(!isInteresting) continue;
    
    switch ( std::abs( particle.PdgCode() ) ) {
      case 311:  // K_0
      case 130:  // K_L
      case 310:  // K_S
      case 111:  // pi0
        // I am interested in this particle's daughter!
        for ( int iDaughter = 0; iDaughter < particle.NumberDaughters(); ++iDaughter ) {
          // find the index of that daughter in the list
          auto daughterTrackId = particle.Daughter(iDaughter);
          InterestingParticles.insert(ParticleIndices.at(daughterTrackId));
        }
        continue; // not interested in this particle
      case kAbsNeutron: // neutrons
        if (!fKeepNeutrons) continue;
        break;
      case 2212:  // proton
      case 211:   // pi+/-
      case 13:    // muon+/-
      case 11:    // e+/-
      case 22:    // photon
        break; // go out and keep it
      default:
        continue;
    } // switch
        
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
