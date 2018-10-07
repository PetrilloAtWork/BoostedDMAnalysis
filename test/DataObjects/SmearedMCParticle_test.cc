/**
 * @file   SmearedMCParticle_test.cc
 * @brief  "Test" of SmearedMCParticle object.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   October 5, 2018
 */

// Boosted Dark Matter libraries
#include "boosteddmanalysis/DataObjects/SmearedMCParticle.h"

// Boost libraries
#define BOOST_TEST_MODULE ( SmearedMCParticle_test )
#include <cetlib/quiet_unit_test.hpp> // BOOST_AUTO_TEST_CASE()
#include <boost/test/test_tools.hpp> // BOOST_CHECK(), BOOST_CHECK_EQUAL()

// C/C++ standard libraries
#include <type_traits> // add_const<>


//------------------------------------------------------------------------------
void test_SmearedMCParticle_defaultConstructor() {
  
  BOOST_TEST_MESSAGE("Testing default-constructed SmearedMCParticle");
  
  using Momentum_t = bdm::SmearedMCParticle::Momentum_t;
  using Position_t = bdm::SmearedMCParticle::Position_t;
  using PDGID_t    = bdm::SmearedMCParticle::PDGID_t;
  using Flags_t    = bdm::SmearedMCParticle::Flags_t;
  
  Momentum_t expected_momentum; // default-constructed
  double expected_energy = 0.0;
  Position_t expected_startPosition;
  Position_t expected_endPosition;
  double expected_time = 0.0;
  PDGID_t expected_id = bdm::SmearedMCParticle::InvalidParticleID;
  Flags_t expected_flags = bdm::SmearedMCParticle::invalidParticleFlags();
  TParticlePDG const* expected_info = nullptr;
  
  bdm::SmearedMCParticle const particle;
  
  BOOST_CHECK_EQUAL(particle.momentum(),      expected_momentum);
  BOOST_CHECK_EQUAL(particle.cp(),            expected_momentum);
  BOOST_CHECK_EQUAL(particle.energy(),        expected_energy);
  BOOST_CHECK_EQUAL(particle.E(),             expected_energy);
  BOOST_CHECK_EQUAL(particle.startPosition(), expected_startPosition);
  BOOST_CHECK_EQUAL(particle.endPosition(),   expected_endPosition);
  BOOST_CHECK_EQUAL(particle.time(),          expected_time);
  BOOST_CHECK_EQUAL(particle.t(),             expected_time);
  BOOST_CHECK      (!particle.hasParticleId());
  BOOST_CHECK_EQUAL(particle.particleId(),    expected_id);
  BOOST_CHECK_EQUAL(particle.particleInfo(),  expected_info);
  BOOST_CHECK_EQUAL(particle.particleName(),  "<none>");
  BOOST_CHECK_EQUAL(particle.flags(),         expected_flags);
  BOOST_CHECK      (!particle.isValid());
  BOOST_CHECK      (particle.isInvalid());
  
  std::cout << "Default-constructed particle: " << particle << std::endl;
  for (unsigned int verb = 0U; verb <= bdm::SmearedMCParticle::MaxDumpVerbosity;
    ++verb)
  {
    std::cout << "Default-constructed particle (verbosity: " << verb << "): ";
    particle.dump(std::cout, verb, " => ");
    std::cout << std::endl;
  } // for
  
} // test_SmearedMCParticle_defaultConstructor()


void test_SmearedMCParticle_directConstructor() {
  
  BOOST_TEST_MESSAGE("Testing SmearedMCParticle constructor");
  
  using Momentum_t = bdm::SmearedMCParticle::Momentum_t;
  using Position_t = bdm::SmearedMCParticle::Position_t;
  using PDGID_t    = bdm::SmearedMCParticle::PDGID_t;
  using Flags_t    = bdm::SmearedMCParticle::Flags_t;
  
  Momentum_t expected_momentum = { 3.0, 4.0, 12.0 }; // default-constructed
  double expected_energy = 13.0;
  Position_t expected_startPosition = { -1.0, -3.0, -2.0 };
  Position_t expected_endPosition = expected_startPosition + expected_momentum;
  double expected_time = 1234.0;
  PDGID_t expected_id = 11; // electron
  Flags_t expected_flags = bdm::SmearedMCParticle::flValid;
  TParticlePDG const* expected_info
    = TDatabasePDG::Instance()->GetParticle(expected_id);
  
  bdm::SmearedMCParticle const particle(
    expected_id, expected_momentum, expected_energy,
    expected_time, expected_startPosition, expected_endPosition
    );
  
  BOOST_CHECK_EQUAL(particle.momentum(),      expected_momentum);
  BOOST_CHECK_EQUAL(particle.cp(),            expected_momentum);
  BOOST_CHECK_EQUAL(particle.energy(),        expected_energy);
  BOOST_CHECK_EQUAL(particle.E(),             expected_energy);
  BOOST_CHECK_EQUAL(particle.startPosition(), expected_startPosition);
  BOOST_CHECK_EQUAL(particle.endPosition(),   expected_endPosition);
  BOOST_CHECK_EQUAL(particle.time(),          expected_time);
  BOOST_CHECK_EQUAL(particle.t(),             expected_time);
  BOOST_CHECK      (particle.hasParticleId());
  BOOST_CHECK_EQUAL(particle.particleId(),    expected_id);
  BOOST_CHECK_EQUAL(particle.particleInfo(),  expected_info);
  BOOST_CHECK_EQUAL(particle.particleName(),  "e-");
  BOOST_CHECK_EQUAL(particle.flags(),         expected_flags);
  BOOST_CHECK      (particle.isValid());
  BOOST_CHECK      (!particle.isInvalid());
  
  std::cout << "Particle: " << particle << std::endl;
  for (unsigned int verb = 0U; verb <= bdm::SmearedMCParticle::MaxDumpVerbosity;
    ++verb)
  {
    std::cout << "Particle (verbosity: " << verb << "): ";
    particle.dump(std::cout, verb, " => ");
    std::cout << std::endl;
  } // for
  
} // test_SmearedMCParticle_directConstructor()



// -----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(SmearedMCParticleTest) {
  test_SmearedMCParticle_defaultConstructor();
  test_SmearedMCParticle_directConstructor();
}


// -----------------------------------------------------------------------------
