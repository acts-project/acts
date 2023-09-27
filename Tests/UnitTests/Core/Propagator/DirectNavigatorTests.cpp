// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Propagator/ActionList.hpp"
#include "Acts/Propagator/DirectNavigator.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/MaterialInteractor.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/SurfaceCollector.hpp"
#include "Acts/Tests/CommonHelpers/CylindricalTrackingGeometry.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"

#include <memory>

namespace bdata = boost::unit_test::data;
namespace tt = boost::test_tools;
using namespace Acts::UnitLiterals;

namespace Acts {
namespace Test {

// Create a test context
GeometryContext tgContext = GeometryContext();
MagneticFieldContext mfContext = MagneticFieldContext();

CylindricalTrackingGeometry cGeometry(tgContext);
auto tGeometry = cGeometry();

// Create a navigator for this tracking geometry
Navigator navigator({tGeometry});
DirectNavigator dnavigator;

using BField = ConstantBField;
using Stepper = EigenStepper<>;
using ReferencePropagator = Propagator<Stepper, Navigator>;
using DirectPropagator = Propagator<Stepper, DirectNavigator>;

const double Bz = 2_T;
auto bField = std::make_shared<BField>(Vector3{0, 0, Bz});
Stepper estepper(bField);
Stepper dstepper(bField);

ReferencePropagator rpropagator(std::move(estepper), std::move(navigator));
DirectPropagator dpropagator(std::move(dstepper), std::move(dnavigator));

const int ntests = 1000;
const int skip = 0;
bool debugMode = false;
bool referenceTiming = false;
bool oversteppingTest = false;
double oversteppingMaxStepSize = 1_mm;

/// The actual test nethod that runs the test
/// can be used with several propagator types
///
/// @tparam rpropagator_t is the reference propagator type
/// @tparam dpropagator_t is the direct propagator type
///
/// @param rprop is the reference propagator instance
/// @param dprop is the direct propagator instance
/// @param pT the transverse momentum
/// @param phi the azimuthal angle of the track at creation
/// @param theta the polar angle of the track at creation
/// @param charge is the charge of the particle
/// @param index is the run index from the test
template <typename rpropagator_t, typename dpropagator_t>
void runTest(const rpropagator_t& rprop, const dpropagator_t& dprop, double pT,
             double phi, double theta, int charge, double time, int index) {
  double dcharge = -1 + 2 * charge;

  if (index < skip) {
    return;
  }

  // Define start parameters from ranom input
  double p = pT / sin(theta);
  CurvilinearTrackParameters start(Vector4(0, 0, 0, time), phi, theta,
                                   dcharge / p);

  using EndOfWorld = EndOfWorldReached;

  // Action list and abort list
  using RefereceActionList = ActionList<MaterialInteractor, SurfaceCollector<>>;
  using ReferenceAbortList = AbortList<EndOfWorld>;

  // Options definition
  using Options = PropagatorOptions<RefereceActionList, ReferenceAbortList>;
  Options pOptions(tgContext, mfContext);
  if (oversteppingTest) {
    pOptions.maxStepSize = oversteppingMaxStepSize;
  }

  // Surface collector configuration
  auto& sCollector = pOptions.actionList.template get<SurfaceCollector<>>();
  sCollector.selector.selectSensitive = true;
  sCollector.selector.selectMaterial = true;

  // Result is immediately used, non-valid result would indicate failure
  const auto& pResult = rprop.propagate(start, pOptions).value();
  auto& cSurfaces = pResult.template get<SurfaceCollector<>::result_type>();
  auto& cMaterial = pResult.template get<MaterialInteractor::result_type>();
  const Surface& destination = pResult.endParameters->referenceSurface();

  std::cout << " - the standard navigator yielded "
            << cSurfaces.collected.size() << " collected surfaces" << std::endl;

  if (not referenceTiming) {
    // Create the surface sequence
    std::vector<const Surface*> surfaceSequence;
    surfaceSequence.reserve(cSurfaces.collected.size());
    for (auto& cs : cSurfaces.collected) {
      surfaceSequence.push_back(cs.surface);
    }

    // Action list for direct navigator with its initalizer
    using DirectActionList = ActionList<DirectNavigator::Initializer,
                                        MaterialInteractor, SurfaceCollector<>>;

    // Direct options definition
    using DirectOptions = PropagatorOptions<DirectActionList, AbortList<>>;
    DirectOptions dOptions(tgContext, mfContext);
    // Set the surface sequence
    auto& dInitializer =
        dOptions.actionList.get<DirectNavigator::Initializer>();
    dInitializer.navSurfaces = surfaceSequence;
    // Surface collector configuration
    auto& dCollector = dOptions.actionList.template get<SurfaceCollector<>>();
    dCollector.selector.selectSensitive = true;
    dCollector.selector.selectMaterial = true;

    // Now redo the propagation with the direct propagator
    const auto& ddResult =
        dprop.propagate(start, destination, dOptions).value();
    auto& ddSurfaces = ddResult.template get<SurfaceCollector<>::result_type>();
    auto& ddMaterial = ddResult.template get<MaterialInteractor::result_type>();

    // CHECK if you have as many surfaces collected as the default navigator
    BOOST_CHECK_EQUAL(cSurfaces.collected.size(), ddSurfaces.collected.size());
    CHECK_CLOSE_REL(cMaterial.materialInX0, ddMaterial.materialInX0, 1e-3);

    // Now redo the propagation with the direct propagator - without destination
    const auto& dwResult = dprop.propagate(start, dOptions).value();
    auto& dwSurfaces = dwResult.template get<SurfaceCollector<>::result_type>();

    // CHECK if you have as many surfaces collected as the default navigator
    BOOST_CHECK_EQUAL(cSurfaces.collected.size(), dwSurfaces.collected.size());
  }
}

// This test case checks that no segmentation fault appears
// - this tests the collection of surfaces
BOOST_DATA_TEST_CASE(
    test_direct_navigator,
    bdata::random((bdata::seed = 20,
                   bdata::distribution =
                       std::uniform_real_distribution<>(0.15_GeV, 10_GeV))) ^
        bdata::random((bdata::seed = 21,
                       bdata::distribution =
                           std::uniform_real_distribution<>(-M_PI, M_PI))) ^
        bdata::random((bdata::seed = 22,
                       bdata::distribution =
                           std::uniform_real_distribution<>(1.0, M_PI - 1.0))) ^
        bdata::random(
            (bdata::seed = 23,
             bdata::distribution = std::uniform_int_distribution<>(0, 1))) ^
        bdata::random(
            (bdata::seed = 24,
             bdata::distribution = std::uniform_int_distribution<>(0, 100))) ^
        bdata::xrange(ntests),
    pT, phi, theta, charge, time, index) {
  // Run the test
  runTest(rpropagator, dpropagator, pT, phi, theta, charge, time, index);
}

}  // namespace Test
}  // namespace Acts
