/**
 * @file   boosteddmanalysis/Reconstruction/SmearedReconstructionAlg.cxx
 * @brief  Algorithm reconstructing particles by parametric smearing.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   October 6, 2018
 * @see    boosteddmanalysis/Reconstruction/SmearedReconstructionAlg.h
 * 
 */

// algorithm header
#include "boosteddmanalysis/Reconstruction/SmearedReconstructionAlg.h"

// framework libraries
#include "messagefacility/MessageLogger/MessageLogger.h"

// CLHEP library
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandomEngine.h" // CLHEP::HepRandomEngine

// C++ standard library
#include <iterator> // std::back_inserter()
#include <utility> // std::move()
#include <stdexcept> // std::out_of_range



//------------------------------------------------------------------------------
//--- bdm::SmearedReconstructionAlg::ReconstructionParameters
//------------------------------------------------------------------------------
bdm::SmearedReconstructionAlg::ReconstructionParameters::ReconstructionParameters
  (std::size_t maxParameters)
{
  // since we are playing with pointers stored in a vector, we don't want the
  // vector to expand after we have taken its pointers; so we add restrictions
  fParameters.reserve(maxParameters);
} // ReconstructionParameters::ReconstructionParameters()


//------------------------------------------------------------------------------
void bdm::SmearedReconstructionAlg::ReconstructionParameters::registerParams
  (ParticleParameters_t const& params, std::vector<PDGID_t> const& ids)
{
  if (fParameters.size() == fParameters.capacity()) {
    throw std::logic_error(
      "registerParams(): Attempt to register more than the expected "
      + std::to_string(fParameters.capacity()) + " parameter tables"
      );
  }
  fParameters.push_back(params);
  associateParams(fParameters.back(), ids);
} // bdm::SmearedReconstructionAlg::ReconstructionParameters::registerParams()


void bdm::SmearedReconstructionAlg::ReconstructionParameters::registerParams
  (ParticleParameters_t&& params, std::vector<PDGID_t> const& ids)
{
  fParameters.push_back(std::move(params));
  associateParams(fParameters.back(), ids);
} // bdm::SmearedReconstructionAlg::ReconstructionParameters::registerParams()


//------------------------------------------------------------------------------
bdm::SmearedReconstructionAlg::ReconstructionParameters::ParticleParameters_t const& bdm::SmearedReconstructionAlg::ReconstructionParameters::operator[]
  (PDGID_t id) const
{
  auto iParams = fAssociations.find(id);
  if ((iParams == fAssociations.end()) && (id != 0))
    iParams = fAssociations.find(0);
  if (iParams == fAssociations.end()) {
    throw std::out_of_range(
      "No reconstruction parameters registered for particle ID="
      + std::to_string(id)
    );
  }
  return *(iParams->second);
} // bdm::SmearedReconstructionAlg::ReconstructionParameters::operator[]


//------------------------------------------------------------------------------
void bdm::SmearedReconstructionAlg::ReconstructionParameters::associateParams
  (ParticleParameters_t& params, std::vector<PDGID_t> const& ids)
{
  for (auto id: ids) fAssociations[id] = &params;
} // bdm::SmearedReconstructionAlg::ReconstructionParameters::associateParams


//------------------------------------------------------------------------------
//--- bdm::SmearedReconstructionAlg
//------------------------------------------------------------------------------
#ifdef USE_FHICLCPP_VALIDATION
bdm::SmearedReconstructionAlg::SmearedReconstructionAlg(Config const& config)
  : fParameters(fillParameters(config.particleParameters()))
{
  fParameters.dump(mf::LogInfo("SmearedReconstructionAlg"));
}

#else // USE_FHICLCPP_VALIDATION

bdm::SmearedReconstructionAlg::SmearedReconstructionAlg
  (fhicl::ParameterSet const& pset)
  : fParameters(fillParameters(pset.get<std::vector<fhicl::ParameterSet>>("particleParameters")))
{
  fParameters.dump(mf::LogInfo("SmearedReconstructionAlg"));
}



#endif // USE_FHICLCPP_VALIDATION

//------------------------------------------------------------------------------
void bdm::SmearedReconstructionAlg::setup
  (CLHEP::HepRandomEngine& rndEng)
{
  fRndEngine = &rndEng;
} // bdm::SmearedReconstructionAlg::setup()


//------------------------------------------------------------------------------
bdm::SmearedReconstructionAlg::SmearedParticles_t
bdm::SmearedReconstructionAlg::reconstruct
  (std::vector<simb::MCParticle const*> const& particles) const
{
  
  SmearedParticles_t reconstructed; // reconstructed particles
  reconstructed.reserve(particles.size());
  
  auto rndUniform = makeRandomGen<CLHEP::RandFlat>();

  for (simb::MCParticle const* particle: particles) {
    
    // retrieve the reconstruction parameters for this type of particle
    auto const& recoParams = fParameters[particle->PdgCode()];
    
    
    unsigned int nPoints = particle->NumberTrajectoryPoints();
    
    bdm::SmearedMCParticle::Position_t startPosition {
      particle->Vx(),
      particle->Vy(),
      particle->Vz()
      };
    
    bdm::SmearedMCParticle::Position_t endPosition{
      particle->Vx(nPoints-1),
      particle->Vy(nPoints-1),
      particle->Vz(nPoints-1)
      };
    
    double time = particle->T();
    
    bdm::SmearedMCParticle::Momentum_t momentum
      { particle->Px(), particle->Py(), particle->Pz() };
    
    double energy = particle->E();
    
    // TODO apply smearing and thresholds
    bool bCouldReconstruct = rndUniform.fire() < recoParams.fRecoEff;
    
    reconstructed.emplace_back(
      particle->PdgCode(), // id
      momentum,
      energy,
      time,
      startPosition,
      endPosition,
      bCouldReconstruct
        ? bdm::SmearedMCParticle::defaultFlags()
        : bdm::SmearedMCParticle::invalidParticleFlags()
      );
    
  } // for
  
  return reconstructed;
} // bdm::SmearedReconstructionAlg::reconstruct()


//------------------------------------------------------------------------------
std::vector<simb::MCParticle const*>
bdm::SmearedReconstructionAlg::toMCParticlePointers
  (std::vector<simb::MCParticle> const& particles)
{
  std::vector<simb::MCParticle const*> pointers;
  pointers.reserve(particles.size());
  
  for (auto const& part: particles) pointers.push_back(&part);
  
  return pointers;
} // bdm::SmearedReconstructionAlg::toMCParticlePointers(vector)


//------------------------------------------------------------------------------
std::vector<simb::MCParticle const*>
bdm::SmearedReconstructionAlg::toMCParticlePointers
  (art::PtrVector<simb::MCParticle> const& particles)
{
  std::vector<simb::MCParticle const*> pointers;
  pointers.reserve(particles.size());
  
  for (auto const& part: particles) pointers.push_back(part.get());
  
  return pointers;
} // bdm::SmearedReconstructionAlg::toMCParticlePointers(PtrVector)


//------------------------------------------------------------------------------
std::vector<simb::MCParticle const*>
bdm::SmearedReconstructionAlg::toMCParticlePointers
  (std::vector<art::Ptr<simb::MCParticle>> const& particles)
{
  std::vector<simb::MCParticle const*> pointers;
  pointers.reserve(particles.size());
  
  for (auto const& part: particles) pointers.push_back(part.get());
  
  return pointers;
} // bdm::SmearedReconstructionAlg::toMCParticlePointers(PtrVector)


//------------------------------------------------------------------------------
#ifdef USE_FHICLCPP_VALIDATION
bdm::SmearedReconstructionAlg::ReconstructionParameters
bdm::SmearedReconstructionAlg::fillParameters
  (std::vector<ParticleParameterConfig> const& params)
{
  ReconstructionParameters parameters(params.size());
  
  for (ParticleParameterConfig const& particleConfig: params) {
    
    ReconstructionParameters::ParticleParameters_t partParams;
    partParams.fRecoEff = particleConfig.efficiency();
    partParams.fDirectionSmearingAngle = particleConfig.direction();
    
    parameters.registerParams(std::move(partParams), particleConfig.id());
    
  } // for
  
  return parameters;
} // bdm::SmearedReconstructionAlg::fillParameters()

#else // USE_FHICLCPP_VALIDATION
bdm::SmearedReconstructionAlg::ReconstructionParameters
bdm::SmearedReconstructionAlg::fillParameters
  (std::vector<fhicl::ParameterSet> const& params)
{
  ReconstructionParameters parameters(params.size());
  
  for (fhicl::ParameterSet const& particleConfig: params) {
    
    ReconstructionParameters::ParticleParameters_t partParams;
    partParams.fRecoEff
      = particleConfig.get<double>("efficiency", 1.0);
    partParams.fDirectionSmearingAngle
      = particleConfig.get<double>("direction", 0.0);
    
    parameters.registerParams
      (std::move(partParams), particleConfig.get<std::vector<PDGID_t>>("id"));
    
  } // for
  
  return parameters;
} // bdm::SmearedReconstructionAlg::fillParameters()


#endif // USE_FHICLCPP_VALIDATION


//------------------------------------------------------------------------------
