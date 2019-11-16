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

// LArSoft libraries
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h" // util::pi()
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h" // geo::Zaxis()

// framework libraries
#include "messagefacility/MessageLogger/MessageLogger.h"

// CLHEP library
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandomEngine.h" // CLHEP::HepRandomEngine

// ROOT libraries
#include "Math/GenVector/RotationY.h"
#include "Math/GenVector/RotationZ.h"

// C++ standard library
#include <iterator> // std::back_inserter()
#include <utility> // std::move()
#include <stdexcept> // std::out_of_range
#include <cmath> // std::abs()


//------------------------------------------------------------------------------
namespace {
  
  /// Returns the square of the specified value.
  template <typename T>
  T sqr(T value) { return value * value; }
  
  
  /// Returns a random value distributed within a truncated Gaussian between
  /// `min` and `max`.
  double cappedGaus(CLHEP::RandGauss& gaus, double min, double max)
    {
      while (true) {
        double const value = gaus.fire();
        if ((value >= min) && (value <= max)) return value;
      } // while
    }
    
  /// Returns a random value distributed according to a Gaussian distribution
  /// with specified `mean` and `sigma`. The deviation from `mean` must be
  /// contained within `minCapSigma` and `maxCapSigma` times the standard
  /// deviation `sigma`, _and_ the value must be larger than `min`.
  /// The value is extracted over and over until these requirements are
  /// satisfied.
  double truncatedCappedGaus(
    CLHEP::RandGauss& gaus,
    double mean, double sigma, double minCapSigma, double maxCapSigma, double min
    )
    {
      while (true) {
        double const nSigmas = cappedGaus(gaus, minCapSigma, maxCapSigma);
        double const value = mean + sigma * nSigmas;
        if (value >= min) return value;        
      } // while     
    } // cappedEnergyGaus()
  
  
  geo::Vector_t smearDirection
    (geo::Vector_t const& v, double dtheta, double dphi)
    {
      //
      // The idea is to add a displacement to the original vector;
      // the displacement must be added orthogonally to the original vector,
      // and it is generated on a plane: its modulus is proportional to the
      // smearing parameter, while the direction on the plane is uniformly
      // random.
      // 
      // The displacement vector is rotated to be orthogonal to the original
      // one, and then added. Given that the azimuthal angle is random, the
      // details of the direction of the displacement rotation around the
      // original vector are not relevant.
      //
      // In the end, we restore the modulus of the vector.
      //
      if (dtheta == 0.0) return v; // actually, no smearing
      double const vR = v.R();
      double const dr = dtheta* vR;
      geo::Vector_t dv { dr * std::cos(dphi), dr * std::sin(dphi), 0.0 };
      auto const rotation
        = ROOT::Math::RotationZ(v.Phi()) * ROOT::Math::RotationY(v.Theta());
      double const rescale = std::sqrt(1.0 / (1.0 + sqr(dr/vR)));
      return rescale * (v + rotation * dv);
    }
  
} // local namespace


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
  auto rndGaus = makeRandomGen<CLHEP::RandGauss>();

  for (simb::MCParticle const* particle: particles) {
    
    // retrieve the reconstruction parameters for this type of particle
    auto const& recoParams = fParameters[particle->PdgCode()];
    // std::cout << "TrackID: " << particle->TrackId() << ", PDGCode: " << particle->PdgCode() << std::endl;
    //--------------------------------------------------------------------------
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
    
    double const energy = particle->E();
    double const mass = particle->Mass();
    double const p = particle->P();
    
    //--------------------------------------------------------------------------
    
    bool bCouldReconstruct = rndUniform.fire() < recoParams.fRecoEff;
    // std::cout << " => passed reconstruction efficiency (" << recoParams.fRecoEff << ")? " << bCouldReconstruct << std::endl;
    
    // Smear the energy
    double thisESmearing = 0.0;
    if ( p > recoParams.fEnergySmearingLowEDef ) {
      thisESmearing = std::sqrt( sqr( recoParams.fEnergySmearingConstant ) 
        + sqr( recoParams.fEnergySmearingCoeffSqrtE )/energy );
    } else {
      thisESmearing = recoParams.fLowEnergySmearingConstant;
    }
  //  double smearedEnergy = std::max( mass, energy * (1.0 + thisESmearing * cappedGaus(rndGaus, -3.0, 3.0) ));
    double smearedEnergy = truncatedCappedGaus(rndGaus, energy, thisESmearing, -3.0, 3.0, mass);
      
    double smearedKEnergy = smearedEnergy - mass;
    bCouldReconstruct &= smearedKEnergy > recoParams.fKEThreshold;
    // std::cout << " => passed detection threshold (" << smearedKEnergy << " vs. " << recoParams.fKEThreshold << ")? " << bCouldReconstruct << std::endl;
    
    
    bdm::SmearedMCParticle::Momentum_t smearedMomentum;
    if (recoParams.fDirectionSmearingAngle == 0.0)
      smearedMomentum = momentum;
    else {
      double const phi = rndUniform() * 2.0 * util::pi() - util::pi();
      double const dtheta
        = recoParams.fDirectionSmearingAngle
        * std::abs(cappedGaus(rndGaus, -3.0, +3.0))
        ;
      smearedMomentum = smearDirection(momentum, dtheta, phi);
    } // if momentum smearing
    
    double smearedP = std::sqrt( sqr( smearedEnergy ) - sqr( mass ) );
    smearedMomentum *= smearedP / p;
    // std::cout << "PDGCode: " << particle->PdgCode() << ", Energy: " << energy << ", SmearedEnergy: " << smearedEnergy 
    //           << ", Momentum: " << p << ", SmearedMomentum: " << smearedP << std::endl;
    
    //--------------------------------------------------------------------------
    reconstructed.emplace_back(
      particle->PdgCode(), // id
      smearedMomentum,
      smearedEnergy,
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
    partParams.fEnergySmearingFraction = particleConfig.energy();
    partParams.fKEThreshold = particleConfig.threshold();
    
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
    partParams.fEnergySmearingConstant
      = particleConfig.get<double>("energy_constant", 0.0);
    partParams.fEnergySmearingCoeffSqrtE
      = particleConfig.get<double>("energy_coeff_sqrtE", 0.0);
    partParams.fEnergySmearingCoeffSqrtE
      = particleConfig.get<double>("low_energy_def", 0.0);
    partParams.fLowEnergySmearingConstant
      = particleConfig.get<double>("low_energy_constant", 0.0);
    partParams.fKEThreshold
      = particleConfig.get<double>("threshold", 0.0);
     
    
    parameters.registerParams
      (std::move(partParams), particleConfig.get<std::vector<PDGID_t>>("id"));
    
  } // for
  
  return parameters;
} // bdm::SmearedReconstructionAlg::fillParameters()


#endif // USE_FHICLCPP_VALIDATION


//------------------------------------------------------------------------------
