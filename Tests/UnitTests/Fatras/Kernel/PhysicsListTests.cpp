// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include <random>

#include "Acts/Material/MaterialProperties.hpp"
#include "Acts/Tests/CommonHelpers/PredefinedMaterials.hpp"
#include "Acts/Utilities/Units.hpp"
#include "ActsFatras/EventData/Particle.hpp"
#include "ActsFatras/Kernel/PhysicsList.hpp"

using namespace Acts::UnitLiterals;

namespace {
/// Physics process that does not trigger a break
struct SterileProcess {
  int some_parameter = 0;

  /// call operator
  template <typename generator_t, typename detector_t, typename particle_t>
  bool operator()(generator_t &, const detector_t &, particle_t &,
                  std::vector<particle_t> &) const {
    return false;
  }
};

/// Physics process that DOES trigger a break
struct FatalProcess {
  /// call operator
  template <typename generator_t, typename detector_t, typename particle_t>
  bool operator()(generator_t &, const detector_t &, particle_t &,
                  std::vector<particle_t> &) const {
    return true;
  }
};

struct Fixture {
  std::ranlux48 generator{23};
  Acts::MaterialProperties slab =
      Acts::MaterialProperties(Acts::Test::makeBeryllium(), 1_mm);
  ActsFatras::Particle inputParticle;
};
}  // namespace

BOOST_AUTO_TEST_SUITE(FatrasPhysicsList)

BOOST_AUTO_TEST_CASE(Empty) {
  Fixture fix;
  ActsFatras::PhysicsList<> emptyList;
  std::vector<ActsFatras::Particle> outgoing;

  // w/o processes the list should never abort
  BOOST_TEST(
      not emptyList(fix.generator, fix.slab, fix.inputParticle, outgoing));
}

BOOST_AUTO_TEST_CASE(SingleSterile) {
  Fixture fix;
  ActsFatras::PhysicsList<SterileProcess> sterileList;
  std::vector<ActsFatras::Particle> outgoing;

  // set some process parameters
  sterileList.get<SterileProcess>().some_parameter = 2;
  BOOST_TEST(sterileList.get<SterileProcess>().some_parameter == 2);

  // sterile process should never abort
  BOOST_TEST(
      not sterileList(fix.generator, fix.slab, fix.inputParticle, outgoing));
}

BOOST_AUTO_TEST_CASE(SingleFatal) {
  Fixture fix;
  ActsFatras::PhysicsList<FatalProcess> fatalList;
  std::vector<ActsFatras::Particle> outgoing;

  // fatal process must always abort
  BOOST_TEST(fatalList(fix.generator, fix.slab, fix.inputParticle, outgoing));
  // unless we disable it
  fatalList.disable<FatalProcess>();
  BOOST_TEST(
      not fatalList(fix.generator, fix.slab, fix.inputParticle, outgoing));
}

BOOST_AUTO_TEST_CASE(SterileFatal) {
  Fixture fix;
  ActsFatras::PhysicsList<SterileProcess, FatalProcess> physicsList;
  std::vector<ActsFatras::Particle> outgoing;

  // the contained fatal process must always abort
  BOOST_TEST(physicsList(fix.generator, fix.slab, fix.inputParticle, outgoing));
  // with the fatal process disabled, it should go through again
  physicsList.disable<FatalProcess>();
  BOOST_TEST(
      not physicsList(fix.generator, fix.slab, fix.inputParticle, outgoing));
}

BOOST_AUTO_TEST_SUITE_END()
