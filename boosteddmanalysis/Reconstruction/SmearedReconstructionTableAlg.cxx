/**
 * @file   boosteddmanalysis/Reconstruction/SmearedReconstructionTableAlg.cxx
 * @brief  Algorithm reconstructing particles by smearing, where the parameters come from a table.
 * @author Yun-Tse Tsai (yuntse@slac.stanford.edu)
 * @date   July 26, 2019
 * @see    boosteddmanalysis/Reconstruction/SmearedReconstructionTableAlg.h
 * 
 */

// algorithm header
#include "boosteddmanalysis/Reconstruction/SmearedReconstructionTableAlg.h"

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
//--- bdm::SmearedReconstructionTableAlg::ReconstructionTables
//------------------------------------------------------------------------------
bdm::SmearedReconstructionTableAlg::ReconstructionTables::ReconstructionTables
  (std::size_t maxParameters)
{
  // since we are playing with pointers stored in a vector, we don't want the
  // vector to expand after we have taken its pointers; so we add restrictions
  fParameters.reserve(maxParameters);
} // ReconstructionTables::ReconstructionTables()


//------------------------------------------------------------------------------
void bdm::SmearedReconstructionTableAlg::ReconstructionTables::registerParams
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
} // bdm::SmearedReconstructionTableAlg::ReconstructionTables::registerParams()


void bdm::SmearedReconstructionTableAlg::ReconstructionTables::registerParams
  (ParticleParameters_t&& params, std::vector<PDGID_t> const& ids)
{
  fParameters.push_back(std::move(params));
  associateParams(fParameters.back(), ids);
} // bdm::SmearedReconstructionTableAlg::ReconstructionTables::registerParams()


//------------------------------------------------------------------------------
bdm::SmearedReconstructionTableAlg::ReconstructionTables::ParticleParameters_t const& bdm::SmearedReconstructionTableAlg::ReconstructionTables::operator[]
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
} // bdm::SmearedReconstructionTableAlg::ReconstructionTables::operator[]


//------------------------------------------------------------------------------
void bdm::SmearedReconstructionTableAlg::ReconstructionTables::associateParams
  (ParticleParameters_t& params, std::vector<PDGID_t> const& ids)
{
  for (auto id: ids) fAssociations[id] = &params;
} // bdm::SmearedReconstructionTableAlg::ReconstructionTables::associateParams

//------------------------------------------------------------------------------
//--- bdm::SmearedReconstructionTableAlg
//------------------------------------------------------------------------------
bdm::SmearedReconstructionTableAlg::SmearedReconstructionTableAlg
  (fhicl::ParameterSet const& pset)
  : SmearedReconstructionAlg(pset)
  , fParameters(fillParameters(pset.get<std::vector<fhicl::ParameterSet>>("particleParameters")))
{
  fParameters.dump(mf::LogInfo("SmearedReconstructionTableAlg"));
}

//------------------------------------------------------------------------------

bdm::SmearedReconstructionTableAlg::SmearedParticles_t
bdm::SmearedReconstructionTableAlg::reconstruct
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
    double KEnergy = energy - mass;
    size_t thisKEBin = 0;
    for ( size_t iKEBin = 0; iKEBin < recoParams.fTrueKEs.size() - 1; ++iKEBin ) {
      if ( KEnergy >= recoParams.fTrueKEs[iKEBin] && KEnergy < recoParams.fTrueKEs[iKEBin+1] ) {
        thisKEBin = iKEBin;
        break;
      }
    }
    if ( KEnergy >= recoParams.fTrueKEs.back() ) thisKEBin = recoParams.fTrueKEs.size() - 1;
    
    double thisKESmearing = recoParams.fEResolutions[thisKEBin];
    // double smearedKEnergy = KEnergy * std::max( 0., 1.0 + thisKESmearing * cappedGaus( rndGaus, -3.0, 3.0 ) );
    double smearedKEnergy = truncatedCappedGaus( rndGaus, KEnergy, thisKESmearing, -3.0, 3.0, 0. );
    double smearedEnergy = smearedKEnergy + mass;
      
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

    // std::cout << "PDGCode: " << particle->PdgCode() << " thisBin: " << thisKEBin
    //           << ", KE: " << KEnergy << ", SmearedKE: " << smearedKEnergy
    //           << ", Energy: " << energy << ", SmearedEnergy: " << smearedEnergy 
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
} // bdm::SmearedReconstructionTableAlg::reconstruct()

//------------------------------------------------------------------------------
bdm::SmearedReconstructionTableAlg::ReconstructionTables
bdm::SmearedReconstructionTableAlg::fillParameters
  (std::vector<fhicl::ParameterSet> const& params)
{
  ReconstructionTables parameters(params.size());
  
  for (fhicl::ParameterSet const& particleConfig: params) {
    
    ReconstructionTables::ParticleParameters_t partParams;
    partParams.fRecoEff
      = particleConfig.get<double>("efficiency", 1.0);
    partParams.fDirectionSmearingAngle
      = particleConfig.get<double>("direction", 0.0);
    partParams.fTrueKEs
      = particleConfig.get< std::vector<double> >("trueKE");
    partParams.fMeanDepEs
      = particleConfig.get< std::vector<double> >("aveDepE");
    partParams.fSTDDepEs
      = particleConfig.get< std::vector<double> >("stdDepE");
    partParams.fKEThreshold
      = particleConfig.get<double>("threshold", 0.0);
    
    if ( partParams.fTrueKEs.size() != partParams.fMeanDepEs.size() 
      || partParams.fTrueKEs.size() != partParams.fSTDDepEs.size() ) {
      throw cet::exception("SmearedReconstructionTableAlg") << "The lengths of TrueKE ("
        << partParams.fTrueKEs.size() << "), AveDepE (" << partParams.fMeanDepEs.size()
        << "), and StdDepE (" << partParams.fSTDDepEs.size() << ") are not identical!\n";
    }
    
    for ( size_t iTrueKE = 1; iTrueKE < partParams.fTrueKEs.size(); ++iTrueKE ) {
      if ( partParams.fTrueKEs[iTrueKE] < partParams.fTrueKEs[iTrueKE-1] )
        throw cet::exception("SmearedReconstructionTableAlg") 
          << "TrueKE is not monotonously increasing, which will not make the algorithm.\n";
    }
    
    partParams.fEResolutions.resize( partParams.fTrueKEs.size() );
    
    for ( size_t iTrueKE = 0; iTrueKE < partParams.fTrueKEs.size(); ++iTrueKE ) {
      if ( partParams.fMeanDepEs[iTrueKE] < std::numeric_limits< double >::min() )
        partParams.fEResolutions[iTrueKE] = 1.;
      else
        partParams.fEResolutions[iTrueKE] = partParams.fSTDDepEs[iTrueKE] / partParams.fMeanDepEs[iTrueKE];  
    }
    
    parameters.registerParams
      (std::move(partParams), particleConfig.get<std::vector<PDGID_t>>("id"));
    
  } // for
  
  return parameters;
} // bdm::SmearedReconstructionAlg::fillParameters()
