// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///  Boost include(s)
#define BOOST_TEST_MODULE MaterialCollection Tests

#include <boost/test/included/unit_test.hpp>
// leave blank line

#include <boost/test/data/test_case.hpp>
// leave blank line

#include <boost/test/output_test_stream.hpp>
// leave blank line

#include <memory>
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Extrapolator/MaterialInteractor.hpp"
#include "Acts/Extrapolator/Navigator.hpp"
#include "Acts/Extrapolator/SurfaceCollector.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialCollector.hpp"
#include "Acts/Propagator/ActionList.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Propagator/detail/DebugOutputActor.hpp"
#include "Acts/Propagator/detail/IntersectionCorrector.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"
#include "ExtrapolatorTestGeometry.hpp"

namespace bdata = boost::unit_test::data;
namespace tt    = boost::test_tools;

namespace Acts {

namespace Test {

  // Global definitions
  // The path limit abort
  typedef detail::PathLimitReached path_limit;

  std::vector<std::unique_ptr<const Surface>> stepState;
  auto tGeometry = testGeometry<ModuleSurface>(stepState);

  // create a navigator for this tracking geometry
  Navigator navigatorES(tGeometry);
  Navigator navigatorSL(tGeometry);

  typedef ConstantBField                BField;
  typedef detail::IntersectionCorrector StepCorrector;
  typedef EigenStepper<BField, StepCorrector>        EigenStepper;
  typedef Propagator<EigenStepper, Navigator>        EigenPropagator;
  typedef Propagator<StraightLineStepper, Navigator> StraightLinePropagator;

  const double    Bz = 2. * units::_T;
  BField          bField(0, 0, Bz);
  EigenStepper    estepper(bField);
  EigenPropagator epropagator(std::move(estepper), std::move(navigatorES));

  StraightLineStepper    slstepper;
  StraightLinePropagator slpropagator(std::move(slstepper),
                                      std::move(navigatorSL));

  const int ntests           = 10;  // test with 0.4 lower boundary - 952;
  const int skip             = 0;   // test with 0.4 lower boundary - 951;
  bool      debugModeFwd     = false;
  bool      debugModeBwd     = false;
  bool      debugModeFwdStep = false;
  bool      debugModeBwdStep = false;

  /// the actual test nethod that runs the test
  /// can be used with several propagator types
  /// @tparam Propagator_type is the actual propagator type
  ///
  /// @param prop is the propagator instance
  /// @param pT the transverse momentum
  /// @param phi the azimuthal angle of the track at creation
  /// @param theta the polar angle of the track at creation
  /// @parm charge is the charge of the particle
  /// @param index is the run index from the test
  template <typename Propagator_type>
  void
  runTest(const Propagator_type& prop,
          double                 pT,
          double                 phi,
          double                 theta,
          int                    charge,
          int                    index)
  {
    double dcharge = -1 + 2 * charge;

    if (index < skip) return;

    // define start parameters
    double                x  = 0;
    double                y  = 0;
    double                z  = 0;
    double                px = pT * cos(phi);
    double                py = pT * sin(phi);
    double                pz = pT / tan(theta);
    double                q  = dcharge;
    Vector3D              pos(x, y, z);
    Vector3D              mom(px, py, pz);
    CurvilinearParameters start(nullptr, pos, mom, q);

    typedef detail::DebugOutputActor DebugOutput;

    // Action list and abort list
    typedef ActionList<MaterialCollector, DebugOutput> ActionList_type;
    typedef AbortList<> AbortConditions_type;

    typename Propagator_type::template Options<ActionList_type,
                                               AbortConditions_type>
        fwdOptions;

    fwdOptions.maxStepSize = 25. * units::_cm;
    fwdOptions.pathLimit   = 25 * units::_cm;
    fwdOptions.debug       = debugModeFwd;

    // get the material collector and configure it
    auto& fwdMaterialCollector
        = fwdOptions.actionList.template get<MaterialCollector>();
    fwdMaterialCollector.detailedCollection = true;

    // forward material test
    const auto& fwdResult = prop.propagate(start, fwdOptions);
    auto&       fwdMaterial
        = fwdResult.template get<MaterialCollector::result_type>();

    double fwdStepMaterialInX0 = 0.;
    double fwdStepMaterialInL0 = 0.;
    // check that the collected material is not zero
    BOOST_TEST(fwdMaterial.materialInX0 != 0.);
    BOOST_TEST(fwdMaterial.materialInL0 != 0.);
    // check that the sum of all steps is the total material
    for (auto& materialHit : fwdMaterial.collected) {
      auto material = materialHit.material;
      fwdStepMaterialInX0 += materialHit.pathLength / material.X0();
      fwdStepMaterialInL0 += materialHit.pathLength / material.L0();
    }
    BOOST_CHECK_CLOSE(fwdMaterial.materialInX0, fwdStepMaterialInX0, 1e-3);
    BOOST_CHECK_CLOSE(fwdMaterial.materialInL0, fwdStepMaterialInL0, 1e-3);

    // get the forward output to the screen
    if (debugModeFwd) {
      const auto& fwdOutput
          = fwdResult.template get<DebugOutput::result_type>();
      std::cout << ">>> Forward Propgation & Navigation output " << std::endl;
      std::cout << fwdOutput.debugString << std::endl;
      // check if the surfaces are free
      std::cout << ">>> Material steps found on ..." << std::endl;
      for (auto& fwdStepsC : fwdMaterial.collected) {
        std::cout << "--> Surface with "
                  << fwdStepsC.surface->geoID().toString() << std::endl;
      }
    }

    // backward material test
    typename Propagator_type::template Options<ActionList_type,
                                               AbortConditions_type>
        bwdOptions;
    bwdOptions.maxStepSize = -25 * units::_cm;
    bwdOptions.pathLimit   = -25 * units::_cm;
    bwdOptions.direction   = backward;
    bwdOptions.debug       = debugModeBwd;

    // get the material collector and configure it
    auto& bwdMaterialCollector
        = bwdOptions.actionList.template get<MaterialCollector>();
    bwdMaterialCollector.detailedCollection = true;

    const auto& startSurface = start.referenceSurface();
    const auto& bwdResult    = prop.propagate(
        *fwdResult.endParameters.template get(), startSurface, bwdOptions);
    auto& bwdMaterial
        = bwdResult.template get<MaterialCollector::result_type>();

    double bwdStepMaterialInX0 = 0.;
    double bwdStepMaterialInL0 = 0.;

    // check that the collected material is not zero
    BOOST_TEST(bwdMaterial.materialInX0 != 0.);
    BOOST_TEST(bwdMaterial.materialInL0 != 0.);
    // check that the sum of all steps is the total material
    for (auto& materialHit : bwdMaterial.collected) {
      auto material = materialHit.material;
      bwdStepMaterialInX0 += materialHit.pathLength / material.X0();
      bwdStepMaterialInL0 += materialHit.pathLength / material.L0();
    }

    BOOST_CHECK_CLOSE(bwdMaterial.materialInX0, bwdStepMaterialInX0, 1e-3);
    BOOST_CHECK_CLOSE(bwdMaterial.materialInL0, bwdStepMaterialInL0, 1e-3);

    // get the backward output to the screen
    if (debugModeBwd) {
      const auto& bwd_output
          = bwdResult.template get<DebugOutput::result_type>();
      std::cout << ">>> Backward Propgation & Navigation output " << std::endl;
      std::cout << bwd_output.debugString << std::endl;
      // check if the surfaces are free
      std::cout << ">>> Material steps found on ..." << std::endl;
      for (auto& bwdStepsC : bwdMaterial.collected) {
        std::cout << "--> Surface with "
                  << bwdStepsC.surface->geoID().toString() << std::endl;
      }
    }

    // forward-backward compatibility test
    BOOST_TEST(bwdMaterial.collected.size() == fwdMaterial.collected.size());

    BOOST_CHECK_CLOSE(bwdMaterial.materialInX0, fwdMaterial.materialInX0, 1e-3);
    BOOST_CHECK_CLOSE(bwdMaterial.materialInL0, bwdMaterial.materialInL0, 1e-3);

    // stepping from one surface to the next
    // now go from surface to surface and check
    typename Propagator_type::template Options<ActionList_type,
                                               AbortConditions_type>
        fwdStepOptions;

    fwdStepOptions.maxStepSize = 25. * units::_cm;
    fwdStepOptions.pathLimit   = 25 * units::_cm;
    fwdStepOptions.debug       = debugModeFwdStep;

    // get the material collector and configure it
    auto& fwdStepMaterialCollector
        = fwdStepOptions.actionList.template get<MaterialCollector>();
    fwdStepMaterialCollector.detailedCollection = true;

    double fwdStepStepMaterialInX0 = 0.;
    double fwdStepStepMaterialInL0 = 0.;

    if (debugModeFwdStep) {
      // check if the surfaces are free
      std::cout << ">>> Forward steps to be processed sequentially ..."
                << std::endl;
      for (auto& fwdStepsC : fwdMaterial.collected) {
        std::cout << "--> Surface with "
                  << fwdStepsC.surface->geoID().toString() << std::endl;
      }
    }

    // move forward step by step through the surfaces
    const TrackParameters*              sParameters = &start;
    std::vector<const TrackParameters*> stepParameters;
    for (auto& fwdSteps : fwdMaterial.collected) {
      if (debugModeFwdStep)
        std::cout << ">>> Forward step : "
                  << sParameters->referenceSurface().geoID().toString()
                  << " --> " << fwdSteps.surface->geoID().toString()
                  << std::endl;

      // make a forward step
      const auto& fwdStep
          = prop.propagate(*sParameters, (*fwdSteps.surface), fwdStepOptions);
      // get the backward output to the screen
      if (debugModeFwdStep) {
        const auto& fwdStepOutput
            = fwdStep.template get<DebugOutput::result_type>();
        std::cout << fwdStepOutput.debugString << std::endl;
      }

      auto& fwdStepMaterial
          = fwdStep.template get<MaterialCollector::result_type>();
      fwdStepStepMaterialInX0 += fwdStepMaterial.materialInX0;
      fwdStepStepMaterialInL0 += fwdStepMaterial.materialInL0;

      if (fwdStep.endParameters != nullptr) {
        sParameters = fwdStep.endParameters->clone();
        // make sure the parameters do not run out of scope
        stepParameters.push_back(sParameters);
      }
    }
    // final destination surface
    const Surface& dSurface = fwdResult.endParameters->referenceSurface();

    if (debugModeFwdStep)
      std::cout << ">>> Forward step : "
                << sParameters->referenceSurface().geoID().toString() << " --> "
                << dSurface.geoID().toString() << std::endl;

    const auto& fwdStepFinal
        = prop.propagate(*sParameters, dSurface, fwdStepOptions);

    auto& fwdStepMaterial
        = fwdStepFinal.template get<MaterialCollector::result_type>();
    fwdStepStepMaterialInX0 += fwdStepMaterial.materialInX0;
    fwdStepStepMaterialInL0 += fwdStepMaterial.materialInL0;

    // forward-forward step compatibility test
    BOOST_CHECK_CLOSE(fwdStepStepMaterialInX0, fwdStepMaterialInX0, 1e-3);
    BOOST_CHECK_CLOSE(fwdStepStepMaterialInL0, fwdStepMaterialInL0, 1e-3);

    // get the backward output to the screen
    if (debugModeFwdStep) {
      const auto& fwdStepOutput
          = fwdStepFinal.template get<DebugOutput::result_type>();
      std::cout << ">>> Forward final step propgation & navigation output "
                << std::endl;
      std::cout << fwdStepOutput.debugString << std::endl;
    }

    // stepping from one surface to the next : backwards
    // now go from surface to surface and check
    typename Propagator_type::template Options<ActionList_type,
                                               AbortConditions_type>
        bwdStepOptions;

    bwdStepOptions.maxStepSize = -25 * units::_cm;
    bwdStepOptions.pathLimit   = -25 * units::_cm;
    bwdStepOptions.direction   = backward;
    bwdStepOptions.debug       = debugModeBwdStep;

    // get the material collector and configure it
    auto& bwdStepMaterialCollector
        = bwdStepOptions.actionList.template get<MaterialCollector>();
    bwdStepMaterialCollector.detailedCollection = true;

    double bwdStepStepMaterialInX0 = 0.;
    double bwdStepStepMaterialInL0 = 0.;

    if (debugModeBwdStep) {
      // check if the surfaces are free
      std::cout << ">>> Backeard steps to be processed sequentially ..."
                << std::endl;
      for (auto& bwdStepsC : bwdMaterial.collected) {
        std::cout << "--> Surface with "
                  << bwdStepsC.surface->geoID().toString() << std::endl;
      }
    }

    // move forward step by step through the surfaces
    sParameters = fwdResult.endParameters.template get();
    for (auto& bwdSteps : bwdMaterial.collected) {
      if (debugModeBwdStep)
        std::cout << ">>> Backward step : "
                  << sParameters->referenceSurface().geoID().toString()
                  << " --> " << bwdSteps.surface->geoID().toString()
                  << std::endl;
      // make a forward step
      const auto& bwdStep
          = prop.propagate(*sParameters, (*bwdSteps.surface), bwdStepOptions);
      // get the backward output to the screen
      if (debugModeBwdStep) {
        const auto& bwdStepOutput
            = bwdStep.template get<DebugOutput::result_type>();
        std::cout << bwdStepOutput.debugString << std::endl;
      }

      auto& bwdStepMaterial
          = bwdStep.template get<MaterialCollector::result_type>();
      bwdStepStepMaterialInX0 += bwdStepMaterial.materialInX0;
      bwdStepStepMaterialInL0 += bwdStepMaterial.materialInL0;

      if (bwdStep.endParameters != nullptr) {
        sParameters = bwdStep.endParameters->clone();
        // make sure the parameters do not run out of scope
        stepParameters.push_back(sParameters);
      }
    }
    // final destination surface
    const Surface& dbSurface = start.referenceSurface();

    if (debugModeBwdStep)
      std::cout << ">>> Backward step : "
                << sParameters->referenceSurface().geoID().toString() << " --> "
                << dSurface.geoID().toString() << std::endl;

    const auto& bwdStepFinal
        = prop.propagate(*sParameters, dbSurface, bwdStepOptions);

    auto& bwdStepMaterial
        = bwdStepFinal.template get<MaterialCollector::result_type>();
    bwdStepStepMaterialInX0 += bwdStepMaterial.materialInX0;
    bwdStepStepMaterialInL0 += bwdStepMaterial.materialInL0;

    // forward-forward step compatibility test
    BOOST_CHECK_CLOSE(bwdStepStepMaterialInX0, bwdStepMaterialInX0, 1e-3);
    BOOST_CHECK_CLOSE(bwdStepStepMaterialInL0, bwdStepMaterialInL0, 1e-3);

    // get the backward output to the screen
    if (debugModeBwdStep) {
      const auto& bwdStepOutput
          = bwdStepFinal.template get<DebugOutput::result_type>();
      std::cout << ">>> Backward final step propgation & navigation output "
                << std::endl;
      std::cout << bwdStepOutput.debugString << std::endl;
    }
  }
  // This test case checks that no segmentation fault appears
  // - this tests the collection of surfaces
  BOOST_DATA_TEST_CASE(
      test_material_collector,
      bdata::random((bdata::seed = 20,
                     bdata::distribution
                     = std::uniform_real_distribution<>(0.5 * units::_GeV,
                                                        10. * units::_GeV)))
          ^ bdata::random((bdata::seed = 21,
                           bdata::distribution
                           = std::uniform_real_distribution<>(-M_PI, M_PI)))
          ^ bdata::random((bdata::seed = 22,
                           bdata::distribution
                           = std::uniform_real_distribution<>(1.0, M_PI - 1.0)))
          ^ bdata::random((bdata::seed = 23,
                           bdata::distribution
                           = std::uniform_int_distribution<>(0, 1)))
          ^ bdata::xrange(ntests),
      pT,
      phi,
      theta,
      charge,
      index)
  {
    runTest(epropagator, pT, phi, theta, charge, index);
    runTest(slpropagator, pT, phi, theta, charge, index);
  }

}  // namespace Test
}  // namespace Acts
