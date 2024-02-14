// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialInteraction.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Material/SurfaceMaterialMapper.hpp"
#include "Acts/Navigation/DetectorNavigator.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Tests/CommonHelpers/CubicDetector.hpp"
#include "Acts/Tests/CommonHelpers/CubicTrackingGeometry.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Logger.hpp"

auto tContext = Acts::GeometryContext();

using namespace Acts::UnitLiterals;

// volumes:
// -3               0               3
/// |               |               |
// layers:
//  |    P0    P1   |    S0     S1  |
// surfaces:
//  |    S0    S1   |  S3|S4  S4/S5 |
//  |    -2    -1   |  e-1+e  e-2+e |

std::vector<Acts::ActsScalar> materialSteps = {
    -2.8_m, -2.3_m, -2.1_m, -1.9_m, -1.1_m, -0.9_m, 0.9_m, 0.95_m,
    1.05_m, 1.1_m,  1.9_m,  1.95_m, 2.3_m,  2.8_m,  2.9_m};

Acts::RecordedMaterialTrack createRecordedMaterialTrack(
    const std::vector<Acts::ActsScalar>& positions) {
  Acts::RecordedMaterialTrack rmt;
  rmt.first = {Acts::Vector3(-2.95_m, 0, 0), Acts::Vector3(1, 0, 0)};

  for (const auto& p : positions) {
    Acts::MaterialInteraction mint;
    mint.position = Acts::Vector3(p, 0, 0);
    mint.direction = Acts::Vector3(1, 0, 0);

    mint.materialSlab = Acts::MaterialSlab(
        Acts::Material::fromMassDensity(0.1, 0.1, 0.1, 0.1, 0.1), 0.1_mm);
    rmt.second.materialInteractions.push_back(mint);
  }

  return rmt;
}

BOOST_AUTO_TEST_SUITE(SurfaceMaterialMapperTests)

BOOST_AUTO_TEST_CASE(SurfaceMaterialMapperTrackingGeometry) {
  auto cube = Acts::Test::CubicTrackingGeometry(tContext);
  auto tGeometry = cube();

  Acts::Navigator navigator({tGeometry});
  Acts::StraightLineStepper stepper;
  Acts::SurfaceMaterialMapper::StraightLineTGPropagator propagator(
      stepper, std::move(navigator),
      Acts::getDefaultLogger("StraightLinePropagator", Acts::Logging::VERBOSE));

  Acts::SurfaceMaterialMapper::Config smmConfig;

  Acts::SurfaceMaterialMapper smMapper(
      smmConfig, propagator,
      Acts::getDefaultLogger("SurfaceMaterialMapper", Acts::Logging::VERBOSE));

  auto state = smMapper.createState(tContext, tContext, *tGeometry);

  Acts::SurfaceMaterialMapper::State* surfaceState =
      static_cast<Acts::SurfaceMaterialMapper::State*>(state.get());

  BOOST_CHECK_NE(surfaceState, nullptr);
  BOOST_CHECK_EQUAL(surfaceState->accumulatedMaterial.size(), 6u);

  // Run the mapping
  auto mTrack = createRecordedMaterialTrack(materialSteps);
  auto [mapped, unmapped] = smMapper.mapMaterialTrack(*surfaceState, mTrack);
  // Check that all material is mapped
  BOOST_CHECK_EQUAL(mapped.second.materialInteractions.size(),
                    mTrack.second.materialInteractions.size());
  BOOST_CHECK_EQUAL(unmapped.second.materialInteractions.size(), 0u);
  // Check that the material is collected
  CHECK_CLOSE_ABS(mapped.second.materialInX0,
                  mTrack.second.materialInteractions.size(), 1e-10);
  CHECK_CLOSE_ABS(mapped.second.materialInL0,
                  mTrack.second.materialInteractions.size(), 1e-10);

  // Create a new cache and run the mapping again - this time intersections are
  // not necessary
  auto restate = smMapper.createState(tContext, tContext, *tGeometry);
  auto [remapped, reunmapped] = smMapper.mapMaterialTrack(*restate, mapped);

  BOOST_CHECK_EQUAL(mapped.second.materialInteractions.size(),
                    remapped.second.materialInteractions.size());
  CHECK_CLOSE_ABS(mapped.second.materialInX0, remapped.second.materialInX0,
                  1e-10);
  CHECK_CLOSE_ABS(mapped.second.materialInL0, remapped.second.materialInL0,
                  1e-10);
  BOOST_CHECK_EQUAL(unmapped.second.materialInteractions.size(),
                    reunmapped.second.materialInteractions.size());
}

BOOST_AUTO_TEST_CASE(SurfaceMaterialMapperTrackingGeometryWithVeto) {
  // Helepr class for positive X veto
  struct PositiveXVeto {
    bool operator()(const Acts::MaterialInteraction& mi) const {
      return mi.position.x() > 0;
    }
  };

  auto cube = Acts::Test::CubicTrackingGeometry(tContext);
  auto tGeometry = cube();

  Acts::Navigator navigator({tGeometry});
  Acts::StraightLineStepper stepper;
  Acts::SurfaceMaterialMapper::StraightLineTGPropagator propagator(
      stepper, std::move(navigator),
      Acts::getDefaultLogger("StraightLinePropagator", Acts::Logging::VERBOSE));

  Acts::SurfaceMaterialMapper::Config smmConfig;
  smmConfig.veto = PositiveXVeto{};

  Acts::SurfaceMaterialMapper smMapper(
      smmConfig, propagator,
      Acts::getDefaultLogger("VetoSurfaceMaterialMapper",
                             Acts::Logging::VERBOSE));

  auto state = smMapper.createState(tContext, tContext, *tGeometry);
  auto mTrack = createRecordedMaterialTrack(materialSteps);
  auto [mapped, unmapped] = smMapper.mapMaterialTrack(*state, mTrack);

  std::size_t negativeSide = 0;
  std::for_each(materialSteps.begin(), materialSteps.end(),
                [&](const Acts::ActsScalar& p) {
                  if (p < 0) {
                    negativeSide++;
                  }
                });

  // Check that all material is mapped
  BOOST_CHECK_EQUAL(mapped.second.materialInteractions.size(), negativeSide);
  BOOST_CHECK_EQUAL(unmapped.second.materialInteractions.size(),
                    mTrack.second.materialInteractions.size() - negativeSide);
}

BOOST_AUTO_TEST_CASE(SurfaceMaterialMapperTrackingGeometryRemappingWithVeto) {
  // Helepr class for positive X veto
  struct PositiveXVeto {
    bool operator()(const Acts::MaterialInteraction& mi) const {
      return mi.position.x() > 0;
    }
  };

  auto cube = Acts::Test::CubicTrackingGeometry(tContext);
  auto tGeometry = cube();

  // First the plain mapper
  Acts::Navigator pnavigator({tGeometry});
  Acts::StraightLineStepper pstepper;
  Acts::SurfaceMaterialMapper::StraightLineTGPropagator ppropagator(
      pstepper, std::move(pnavigator),
      Acts::getDefaultLogger("StraightLinePropagator", Acts::Logging::VERBOSE));

  Acts::SurfaceMaterialMapper::Config psmmConfig;

  Acts::SurfaceMaterialMapper psmMapper(
      psmmConfig, ppropagator,
      Acts::getDefaultLogger("PlainSurfaceMaterialMapper",
                             Acts::Logging::VERBOSE));

  auto pstate = psmMapper.createState(tContext, tContext, *tGeometry);
  auto mTrack = createRecordedMaterialTrack(materialSteps);
  auto [mapped, unmapped] = psmMapper.mapMaterialTrack(*pstate, mTrack);

  // Check that all material is mapped
  BOOST_CHECK_EQUAL(mapped.second.materialInteractions.size(),
                    mTrack.second.materialInteractions.size());

  // Now the veto mapper
  Acts::Navigator navigator({tGeometry});
  Acts::StraightLineStepper stepper;
  Acts::SurfaceMaterialMapper::StraightLineTGPropagator propagator(
      stepper, std::move(navigator),
      Acts::getDefaultLogger("StraightLinePropagator", Acts::Logging::VERBOSE));

  Acts::SurfaceMaterialMapper::Config smmConfig;
  smmConfig.veto = PositiveXVeto{};

  Acts::SurfaceMaterialMapper smMapper(
      smmConfig, propagator,
      Acts::getDefaultLogger("VetoSurfaceMaterialMapper",
                             Acts::Logging::VERBOSE));

  auto restate = smMapper.createState(tContext, tContext, *tGeometry);
  auto [remapped, reunmapped] = smMapper.mapMaterialTrack(*restate, mapped);

  std::size_t negativeSide = 0;
  std::for_each(materialSteps.begin(), materialSteps.end(),
                [&](const Acts::ActsScalar& p) {
                  if (p < 0)
                    negativeSide++;
                });
  // Check that all material is mapped
  BOOST_CHECK_EQUAL(remapped.second.materialInteractions.size(), negativeSide);
  BOOST_CHECK_EQUAL(reunmapped.second.materialInteractions.size(),
                    mTrack.second.materialInteractions.size() - negativeSide);
}

BOOST_AUTO_TEST_CASE(SurfaceMaterialMapperDetector) {
  auto cube = Acts::Test::CubicDetector(tContext);
  auto detector = cube();

  auto dnCfg = Acts::Experimental::DetectorNavigator::Config();
  dnCfg.detector = detector.get();

  Acts::Experimental::DetectorNavigator navigator(dnCfg);
  Acts::StraightLineStepper stepper;
  Acts::SurfaceMaterialMapper::StraightLineDetPropagator propagator(
      stepper, std::move(navigator),
      Acts::getDefaultLogger("StraightLinePropagator", Acts::Logging::VERBOSE));

  Acts::SurfaceMaterialMapper::Config smmConfig;

  Acts::SurfaceMaterialMapper smMapper(
      smmConfig, propagator,
      Acts::getDefaultLogger("SurfaceMaterialMapper", Acts::Logging::VERBOSE));

  auto state = smMapper.createState(tContext, tContext, *detector);

  Acts::SurfaceMaterialMapper::State* surfaceState =
      static_cast<Acts::SurfaceMaterialMapper::State*>(state.get());

  BOOST_CHECK_NE(surfaceState, nullptr);
  BOOST_CHECK_EQUAL(surfaceState->accumulatedMaterial.size(), 6u);

  // Run the mapping
  auto mTrack = createRecordedMaterialTrack(materialSteps);
  auto [mapped, unmapped] = smMapper.mapMaterialTrack(*surfaceState, mTrack);
  // Check that all material is mapped
  BOOST_CHECK_EQUAL(mapped.second.materialInteractions.size(),
                    mTrack.second.materialInteractions.size());
  BOOST_CHECK_EQUAL(unmapped.second.materialInteractions.size(), 0u);
  // Check that the material is collected
  CHECK_CLOSE_ABS(mapped.second.materialInX0,
                  mTrack.second.materialInteractions.size(), 1e-10);
  CHECK_CLOSE_ABS(mapped.second.materialInL0,
                  mTrack.second.materialInteractions.size(), 1e-10);

  // Create a new cache and run the mapping again - this time intersections are
  // not necessary
  auto restate = smMapper.createState(tContext, tContext, *detector);
  auto [remapped, reunmapped] = smMapper.mapMaterialTrack(*restate, mapped);

  BOOST_CHECK_EQUAL(mapped.second.materialInteractions.size(),
                    remapped.second.materialInteractions.size());
  CHECK_CLOSE_ABS(mapped.second.materialInX0, remapped.second.materialInX0,
                  1e-10);
  CHECK_CLOSE_ABS(mapped.second.materialInL0, remapped.second.materialInL0,
                  1e-10);
  BOOST_CHECK_EQUAL(unmapped.second.materialInteractions.size(),
                    reunmapped.second.materialInteractions.size());
}

BOOST_AUTO_TEST_CASE(SurfaceMaterialMapperDetectorWithVeto) {
  // Helepr class for positive X veto
  struct PositiveXVeto {
    bool operator()(const Acts::MaterialInteraction& mi) const {
      return mi.position.x() > 0;
    }
  };

  auto cube = Acts::Test::CubicDetector(tContext);
  auto detector = cube();

  auto dnCfg = Acts::Experimental::DetectorNavigator::Config();
  dnCfg.detector = detector.get();

  Acts::Experimental::DetectorNavigator navigator(dnCfg);
  Acts::StraightLineStepper stepper;
  Acts::SurfaceMaterialMapper::StraightLineDetPropagator propagator(
      stepper, std::move(navigator),
      Acts::getDefaultLogger("StraightLinePropagator", Acts::Logging::VERBOSE));

  Acts::SurfaceMaterialMapper::Config smmConfig;
  smmConfig.veto = PositiveXVeto{};

  Acts::SurfaceMaterialMapper smMapper(
      smmConfig, propagator,
      Acts::getDefaultLogger("VetoSurfaceMaterialMapper",
                             Acts::Logging::VERBOSE));

  auto state = smMapper.createState(tContext, tContext, *detector);
  auto mTrack = createRecordedMaterialTrack(materialSteps);
  auto [mapped, unmapped] = smMapper.mapMaterialTrack(*state, mTrack);

  std::size_t negativeSide = 0;
  std::for_each(materialSteps.begin(), materialSteps.end(),
                [&](const Acts::ActsScalar& p) {
                  if (p < 0)
                    negativeSide++;
                });

  // Check that all material is mapped
  BOOST_CHECK_EQUAL(mapped.second.materialInteractions.size(), negativeSide);
  BOOST_CHECK_EQUAL(unmapped.second.materialInteractions.size(),
                    mTrack.second.materialInteractions.size() - negativeSide);
}

BOOST_AUTO_TEST_CASE(SurfaceMaterialMapperDetectorRemappingWithVeto) {
  // Helepr class for positive X veto
  struct PositiveXVeto {
    bool operator()(const Acts::MaterialInteraction& mi) const {
      return mi.position.x() > 0;
    }
  };

  auto cube = Acts::Test::CubicDetector(tContext);
  auto detector = cube();

  auto dnCfg = Acts::Experimental::DetectorNavigator::Config();
  dnCfg.detector = detector.get();

  // First the plain mapper
  Acts::Experimental::DetectorNavigator pnavigator(dnCfg);
  Acts::StraightLineStepper pstepper;
  Acts::SurfaceMaterialMapper::StraightLineDetPropagator ppropagator(
      pstepper, std::move(pnavigator),
      Acts::getDefaultLogger("StraightLinePropagator", Acts::Logging::VERBOSE));

  Acts::SurfaceMaterialMapper::Config psmmConfig;

  Acts::SurfaceMaterialMapper psmMapper(
      psmmConfig, ppropagator,
      Acts::getDefaultLogger("PlainSurfaceMaterialMapper",
                             Acts::Logging::VERBOSE));

  auto pstate = psmMapper.createState(tContext, tContext, *detector);
  auto mTrack = createRecordedMaterialTrack(materialSteps);
  auto [mapped, unmapped] = psmMapper.mapMaterialTrack(*pstate, mTrack);

  // Check that all material is mapped
  BOOST_CHECK_EQUAL(mapped.second.materialInteractions.size(),
                    mTrack.second.materialInteractions.size());

  // Now the veto mapper
  Acts::Experimental::DetectorNavigator navigator(dnCfg);
  Acts::StraightLineStepper stepper;
  Acts::SurfaceMaterialMapper::StraightLineDetPropagator propagator(
      stepper, std::move(navigator),
      Acts::getDefaultLogger("StraightLinePropagator", Acts::Logging::VERBOSE));

  Acts::SurfaceMaterialMapper::Config smmConfig;
  smmConfig.veto = PositiveXVeto{};

  Acts::SurfaceMaterialMapper smMapper(
      smmConfig, propagator,
      Acts::getDefaultLogger("VetoSurfaceMaterialMapper",
                             Acts::Logging::VERBOSE));

  auto restate = smMapper.createState(tContext, tContext, *detector);
  auto [remapped, reunmapped] = smMapper.mapMaterialTrack(*restate, mapped);

  std::size_t negativeSide = 0;
  std::for_each(materialSteps.begin(), materialSteps.end(),
                [&](const Acts::ActsScalar& p) {
                  if (p < 0)
                    negativeSide++;
                });
  // Check that all material is mapped
  BOOST_CHECK_EQUAL(remapped.second.materialInteractions.size(), negativeSide);
  BOOST_CHECK_EQUAL(reunmapped.second.materialInteractions.size(),
                    mTrack.second.materialInteractions.size() - negativeSide);
}

BOOST_AUTO_TEST_SUITE_END()
