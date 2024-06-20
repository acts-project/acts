// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Propagator/SurfaceCollector.hpp"
#include "Acts/Propagator/TryAllNavigator.hpp"
#include "Acts/Tests/CommonHelpers/CylindricalTrackingGeometry.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"

namespace bdata = boost::unit_test::data;
using namespace Acts::UnitLiterals;
using Acts::VectorHelpers::perp;

namespace Acts::Test {

// Create a test context
GeometryContext tgContext = GeometryContext();
MagneticFieldContext mfContext = MagneticFieldContext();

CylindricalTrackingGeometry cGeometry(tgContext);
auto tGeometry = cGeometry();

const double Bz = 2_T;
auto bField = std::make_shared<ConstantBField>(Vector3{0, 0, Bz});

using SurfaceCollector = SurfaceCollector<SurfaceSelector>;

std::vector<GeometryIdentifier> collectRelevantGeoIds(
    const SurfaceCollector::result_type& surfaceHits) {
  std::vector<GeometryIdentifier> geoIds;
  for (const auto& surfaceHit : surfaceHits.collected) {
    auto geoId = surfaceHit.surface->geometryId();
    auto material = surfaceHit.surface->surfaceMaterial();
    if (geoId.sensitive() == 0 && material == nullptr) {
      continue;
    }
    geoIds.push_back(geoId);
  }
  return geoIds;
}

/// the actual test method that runs the test can be used with several
/// propagator types
///
/// @tparam propagator_t is the actual propagator type
///
/// @param prop is the propagator instance
/// @param start start parameters for propagation
/// @param debugMode toggle debug mode
template <typename propagator_t>
void runSelfConsistencyTest(const propagator_t& prop,
                            const CurvilinearTrackParameters& start,
                            bool debugMode) {
  // Action list and abort list
  using ActionListType = ActionList<SurfaceCollector>;
  using AbortListType = AbortList<>;
  using Options = PropagatorOptions<ActionListType, AbortListType>;

  // forward surface test
  Options fwdOptions(tgContext, mfContext);
  fwdOptions.pathLimit = 25_cm;

  // get the surface collector and configure it
  auto& fwdSurfaceCollector =
      fwdOptions.actionList.template get<SurfaceCollector>();
  fwdSurfaceCollector.selector.selectSensitive = true;
  fwdSurfaceCollector.selector.selectMaterial = true;
  fwdSurfaceCollector.selector.selectPassive = true;

  if (debugMode) {
    std::cout << ">>> Forward Propagation : start." << std::endl;
  }
  auto fwdResult = prop.propagate(start, fwdOptions).value();
  auto fwdSurfaceHits =
      fwdResult.template get<SurfaceCollector::result_type>().collected;
  auto fwdSurfaces = collectRelevantGeoIds(
      fwdResult.template get<SurfaceCollector::result_type>());

  // get the forward output to the screen
  if (debugMode) {
    // check if the surfaces are free
    std::cout << ">>> Surface hits found on ..." << std::endl;
    for (const auto& fwdSteps : fwdSurfaces) {
      std::cout << "--> Surface with " << fwdSteps << std::endl;
    }
    std::cout << ">>> Forward Propagation : end." << std::endl;
  }

  // backward surface test
  Options bwdOptions(tgContext, mfContext);
  bwdOptions.pathLimit = 25_cm;
  bwdOptions.direction = Direction::Backward;

  // get the surface collector and configure it
  auto& bwdMSurfaceCollector =
      bwdOptions.actionList.template get<SurfaceCollector>();
  bwdMSurfaceCollector.selector.selectSensitive = true;
  bwdMSurfaceCollector.selector.selectMaterial = true;
  bwdMSurfaceCollector.selector.selectPassive = true;

  const auto& startSurface = start.referenceSurface();

  if (debugMode) {
    std::cout << ">>> Backward Propagation : start." << std::endl;
  }
  auto bwdResult =
      prop.propagate(*fwdResult.endParameters, startSurface, bwdOptions)
          .value();
  auto bwdSurfaceHits =
      bwdResult.template get<SurfaceCollector::result_type>().collected;
  auto bwdSurfaces = collectRelevantGeoIds(
      bwdResult.template get<SurfaceCollector::result_type>());

  // get the backward output to the screen
  if (debugMode) {
    // check if the surfaces are free
    std::cout << ">>> Surface hits found on ..." << std::endl;
    for (auto& bwdSteps : bwdSurfaces) {
      std::cout << "--> Surface with " << bwdSteps << std::endl;
    }
    std::cout << ">>> Backward Propagation : end." << std::endl;
  }

  // forward-backward compatibility test
  {
    std::reverse(bwdSurfaces.begin(), bwdSurfaces.end());
    BOOST_CHECK_EQUAL_COLLECTIONS(bwdSurfaces.begin(), bwdSurfaces.end(),
                                  fwdSurfaces.begin(), fwdSurfaces.end());
  }

  // stepping from one surface to the next
  // now go from surface to surface and check
  Options fwdStepOptions(tgContext, mfContext);

  // get the surface collector and configure it
  auto& fwdStepSurfaceCollector =
      fwdOptions.actionList.template get<SurfaceCollector>();
  fwdStepSurfaceCollector.selector.selectSensitive = true;
  fwdStepSurfaceCollector.selector.selectMaterial = true;
  fwdStepSurfaceCollector.selector.selectPassive = true;

  std::vector<GeometryIdentifier> fwdStepSurfaces;

  // move forward step by step through the surfaces
  BoundTrackParameters sParameters = start;
  std::vector<BoundTrackParameters> stepParameters;
  for (auto& fwdSteps : fwdSurfaceHits) {
    if (debugMode) {
      std::cout << ">>> Forward step : "
                << sParameters.referenceSurface().geometryId() << " --> "
                << fwdSteps.surface->geometryId() << std::endl;
    }

    // make a forward step
    auto fwdStep =
        prop.propagate(sParameters, *fwdSteps.surface, fwdStepOptions).value();

    auto fwdStepSurfacesTmp = collectRelevantGeoIds(
        fwdStep.template get<SurfaceCollector::result_type>());
    fwdStepSurfaces.insert(fwdStepSurfaces.end(), fwdStepSurfacesTmp.begin(),
                           fwdStepSurfacesTmp.end());

    if (fwdStep.endParameters.has_value()) {
      // make sure the parameters do not run out of scope
      stepParameters.push_back(*fwdStep.endParameters);
      sParameters = stepParameters.back();
    }
  }
  // final destination surface
  const Surface& dSurface = fwdResult.endParameters->referenceSurface();
  if (debugMode) {
    std::cout << ">>> Forward step : "
              << sParameters.referenceSurface().geometryId() << " --> "
              << dSurface.geometryId() << std::endl;
  }
  auto fwdStepFinal =
      prop.propagate(sParameters, dSurface, fwdStepOptions).value();
  auto fwdStepSurfacesTmp = collectRelevantGeoIds(
      fwdStepFinal.template get<SurfaceCollector::result_type>());
  fwdStepSurfaces.insert(fwdStepSurfaces.end(), fwdStepSurfacesTmp.begin(),
                         fwdStepSurfacesTmp.end());

  // TODO forward-forward step compatibility test

  // stepping from one surface to the next : backwards
  // now go from surface to surface and check
  Options bwdStepOptions(tgContext, mfContext);
  bwdStepOptions.direction = Direction::Backward;

  // get the surface collector and configure it
  auto& bwdStepSurfaceCollector =
      bwdOptions.actionList.template get<SurfaceCollector>();
  bwdStepSurfaceCollector.selector.selectSensitive = true;
  bwdStepSurfaceCollector.selector.selectMaterial = true;
  bwdStepSurfaceCollector.selector.selectPassive = true;

  std::vector<GeometryIdentifier> bwdStepSurfaces;

  // move forward step by step through the surfaces
  sParameters = *fwdResult.endParameters;
  for (auto& bwdSteps : bwdSurfaceHits) {
    if (debugMode) {
      std::cout << ">>> Backward step : "
                << sParameters.referenceSurface().geometryId() << " --> "
                << bwdSteps.surface->geometryId() << std::endl;
    }

    // make a forward step
    auto bwdStep =
        prop.propagate(sParameters, *bwdSteps.surface, bwdStepOptions).value();

    auto bwdStepSurfacesTmp = collectRelevantGeoIds(
        bwdStep.template get<SurfaceCollector::result_type>());
    bwdStepSurfaces.insert(bwdStepSurfaces.end(), bwdStepSurfacesTmp.begin(),
                           bwdStepSurfacesTmp.end());

    if (bwdStep.endParameters.has_value()) {
      // make sure the parameters do not run out of scope
      stepParameters.push_back(*bwdStep.endParameters);
      sParameters = stepParameters.back();
    }
  }
  // final destination surface
  const Surface& dbSurface = start.referenceSurface();
  if (debugMode) {
    std::cout << ">>> Backward step : "
              << sParameters.referenceSurface().geometryId() << " --> "
              << dSurface.geometryId() << std::endl;
  }
  auto bwdStepFinal =
      prop.propagate(sParameters, dbSurface, bwdStepOptions).value();
  auto bwdStepSurfacesTmp = collectRelevantGeoIds(
      bwdStepFinal.template get<SurfaceCollector::result_type>());
  bwdStepSurfaces.insert(bwdStepSurfaces.end(), bwdStepSurfacesTmp.begin(),
                         bwdStepSurfacesTmp.end());

  // TODO backward-backward step compatibility test

  std::reverse(bwdStepSurfaces.begin(), bwdStepSurfaces.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(bwdStepSurfaces.begin(), bwdStepSurfaces.end(),
                                fwdStepSurfaces.begin(), fwdStepSurfaces.end());
}

/// the actual test method that runs the test can be used with several
/// propagator types
///
/// @tparam propagator_probe_t is the probe propagator type
/// @tparam propagator_ref_t is the reference propagator type
///
/// @param propProbe is the probe propagator instance
/// @param propRef is the reference propagator instance
/// @param start start parameters for propagation
/// @param debugMode toggle debug mode
template <typename propagator_probe_t, typename propagator_ref_t>
void runConsistencyTest(const propagator_probe_t& propProbe,
                        const propagator_ref_t& propRef,
                        const CurvilinearTrackParameters& start,
                        bool debugMode) {
  // Action list and abort list
  using ActionListType = ActionList<SurfaceCollector>;
  using AbortListType = AbortList<>;
  using Options = PropagatorOptions<ActionListType, AbortListType>;

  auto run = [&](const auto& prop) {
    // forward surface test
    Options fwdOptions(tgContext, mfContext);
    fwdOptions.pathLimit = 25_cm;
    fwdOptions.maxStepSize = 1_cm;

    // get the surface collector and configure it
    auto& fwdSurfaceCollector =
        fwdOptions.actionList.template get<SurfaceCollector>();
    fwdSurfaceCollector.selector.selectSensitive = true;
    fwdSurfaceCollector.selector.selectMaterial = true;
    fwdSurfaceCollector.selector.selectPassive = true;

    auto fwdResult = prop.propagate(start, fwdOptions).value();
    auto fwdSurfaces = collectRelevantGeoIds(
        fwdResult.template get<SurfaceCollector::result_type>());

    // get the forward output to the screen
    if (debugMode) {
      // check if the surfaces are free
      std::cout << ">>> Surface hits found on ..." << std::endl;
      for (const auto& fwdSteps : fwdSurfaces) {
        std::cout << "--> Surface with " << fwdSteps << std::endl;
      }
    }

    return fwdSurfaces;
  };

  if (debugMode) {
    std::cout << ">>> Probe Propagation : start." << std::endl;
  }
  const auto& probeSurfaces = run(propProbe);
  if (debugMode) {
    std::cout << ">>> Probe Propagation : end." << std::endl;
  }

  if (debugMode) {
    std::cout << ">>> Reference Propagation : start." << std::endl;
  }
  const auto& refSurfaces = run(propRef);
  if (debugMode) {
    std::cout << ">>> Reference Propagation : end." << std::endl;
  }

  // probe-ref compatibility test
  BOOST_CHECK_EQUAL_COLLECTIONS(probeSurfaces.begin(), probeSurfaces.end(),
                                refSurfaces.begin(), refSurfaces.end());
}

const int nTestsSelfConsistency = 500;
const int nTestsRefConsistency = 500;
int skip = 0;
bool debugMode = false;

using EigenStepper = Acts::EigenStepper<>;
using EigenPropagator = Propagator<EigenStepper, Navigator>;
using StraightLinePropagator = Propagator<StraightLineStepper, Navigator>;
using Reference1EigenPropagator = Propagator<EigenStepper, TryAllNavigator>;
using Reference1StraightLinePropagator =
    Propagator<StraightLineStepper, TryAllNavigator>;
using Reference2EigenPropagator =
    Propagator<EigenStepper, TryAllOverstepNavigator>;
using Reference2StraightLinePropagator =
    Propagator<StraightLineStepper, TryAllOverstepNavigator>;

EigenStepper estepper(bField);
StraightLineStepper slstepper;

EigenPropagator epropagator(estepper,
                            Navigator({tGeometry, true, true, false},
                                      getDefaultLogger("e_nav", Logging::INFO)),
                            getDefaultLogger("e_prop", Logging::INFO));
StraightLinePropagator slpropagator(slstepper,
                                    Navigator({tGeometry, true, true, false},
                                              getDefaultLogger("sl_nav",
                                                               Logging::INFO)),
                                    getDefaultLogger("sl_prop", Logging::INFO));

Reference1EigenPropagator refepropagator1(
    estepper,
    TryAllNavigator({tGeometry, true, true, false, BoundaryCheck(false)},
                    getDefaultLogger("ref1_e_nav", Logging::INFO)),
    getDefaultLogger("ref1_e_prop", Logging::INFO));
Reference1StraightLinePropagator refslpropagator1(
    slstepper,
    TryAllNavigator({tGeometry, true, true, false},
                    getDefaultLogger("ref1_sl_nav", Logging::INFO)),
    getDefaultLogger("ref1_sl_prop", Logging::INFO));

Reference2EigenPropagator refepropagator2(
    estepper,
    TryAllOverstepNavigator({tGeometry, true, true, false,
                             BoundaryCheck(false)},
                            getDefaultLogger("ref2_e_nav", Logging::INFO)),
    getDefaultLogger("ref2_e_prop", Logging::INFO));
Reference2StraightLinePropagator refslpropagator2(
    slstepper,
    TryAllOverstepNavigator({tGeometry, true, true, false},
                            getDefaultLogger("ref2_sl_nav", Logging::INFO)),
    getDefaultLogger("ref2_sl_prop", Logging::INFO));

auto eventGen =
    bdata::random((bdata::engine = std::mt19937(), bdata::seed = 20,
                   bdata::distribution = std::uniform_real_distribution<double>(
                       0.5_GeV, 10_GeV))) ^
    bdata::random((bdata::engine = std::mt19937(), bdata::seed = 21,
                   bdata::distribution =
                       std::uniform_real_distribution<double>(-M_PI, M_PI))) ^
    bdata::random((bdata::engine = std::mt19937(), bdata::seed = 22,
                   bdata::distribution = std::uniform_real_distribution<double>(
                       1.0, M_PI - 1.0))) ^
    bdata::random(
        (bdata::engine = std::mt19937(), bdata::seed = 23,
         bdata::distribution = std::uniform_int_distribution<int>(0, 1)));

BOOST_DATA_TEST_CASE(NavigatorSelfConsistency,
                     eventGen ^ bdata::xrange(nTestsSelfConsistency), pT, phi,
                     theta, charge, index) {
  if (index < skip) {
    return;
  }

  double p = pT / std::sin(theta);
  double q = -1 + 2 * charge;
  CurvilinearTrackParameters start(Vector4(0, 0, 0, 0), phi, theta, q / p,
                                   std::nullopt, ParticleHypothesis::pion());

  if (debugMode) {
    std::cout << ">>> Run navigation tests with pT = " << pT
              << "; phi = " << phi << "; theta = " << theta
              << "; charge = " << charge << "; index = " << index << ";"
              << std::endl;
  }

  if (debugMode) {
    std::cout << ">>> Test self consistency epropagator" << std::endl;
  }
  runSelfConsistencyTest(epropagator, start, debugMode);
  if (debugMode) {
    std::cout << ">>> Test self consistency slpropagator" << std::endl;
  }
  runSelfConsistencyTest(slpropagator, start, debugMode);
}

BOOST_DATA_TEST_CASE(NavigatorRef1Consistency,
                     eventGen ^ bdata::xrange(nTestsRefConsistency), pT, phi,
                     theta, charge, index) {
  if (index < skip) {
    return;
  }

  double p = pT / std::sin(theta);
  double q = -1 + 2 * charge;
  CurvilinearTrackParameters start(Vector4(0, 0, 0, 0), phi, theta, q / p,
                                   std::nullopt, ParticleHypothesis::pion());

  if (debugMode) {
    std::cout << ">>> Run navigation tests with pT = " << pT
              << "; phi = " << phi << "; theta = " << theta
              << "; charge = " << charge << "; index = " << index << ";"
              << std::endl;
  }

  if (debugMode) {
    std::cout << ">>> Test reference 1 consistency epropagator" << std::endl;
  }
  runConsistencyTest(epropagator, refepropagator1, start, debugMode);
  if (debugMode) {
    std::cout << ">>> Test reference 1 consistency slpropagator" << std::endl;
  }
  runConsistencyTest(slpropagator, refslpropagator1, start, debugMode);
}

BOOST_DATA_TEST_CASE(NavigatorRef2Consistency,
                     eventGen ^ bdata::xrange(nTestsRefConsistency), pT, phi,
                     theta, charge, index) {
  if (index < skip) {
    return;
  }

  double p = pT / std::sin(theta);
  double q = -1 + 2 * charge;
  CurvilinearTrackParameters start(Vector4(0, 0, 0, 0), phi, theta, q / p,
                                   std::nullopt, ParticleHypothesis::pion());

  if (debugMode) {
    std::cout << ">>> Run navigation tests with pT = " << pT
              << "; phi = " << phi << "; theta = " << theta
              << "; charge = " << charge << "; index = " << index << ";"
              << std::endl;
  }

  if (debugMode) {
    std::cout << ">>> Test reference 2 consistency epropagator" << std::endl;
  }
  runConsistencyTest(epropagator, refepropagator2, start, debugMode);
  if (debugMode) {
    std::cout << ">>> Test reference 2 consistency slpropagator" << std::endl;
  }
  runConsistencyTest(slpropagator, refslpropagator2, start, debugMode);
}

}  // namespace Acts::Test
