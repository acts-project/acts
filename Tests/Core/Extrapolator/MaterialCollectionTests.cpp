// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// clang-format off
#define BOOST_TEST_MODULE MaterialCollection Tests
#include <boost/test/included/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/output_test_stream.hpp>
// clang-format on

#include <memory>

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Extrapolator/MaterialInteractor.hpp"
#include "Acts/Extrapolator/Navigator.hpp"
#include "Acts/Extrapolator/SurfaceCollector.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Propagator/ActionList.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Propagator/detail/DebugOutputActor.hpp"
#include "Acts/Propagator/detail/RelativePathCorrector.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Tests/CommonHelpers/CylindricalTrackingGeometry.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"

namespace bdata = boost::unit_test::data;
namespace tt    = boost::test_tools;

namespace Acts {

namespace Test {

  // Global definitions
  // The path limit abort
  using path_limit = detail::PathLimitReached;

  CylindricalTrackingGeometry cGeometry;
  auto                        tGeometry = cGeometry();

  // create a navigator for this tracking geometry
  Navigator navigatorES(tGeometry);
  Navigator navigatorSL(tGeometry);

  using BField                 = ConstantBField;
  using StepCorrector          = detail::RelativePathCorrector;
  using EigenStepper           = EigenStepper<BField, StepCorrector>;
  using EigenPropagator        = Propagator<EigenStepper, Navigator>;
  using StraightLinePropagator = Propagator<StraightLineStepper, Navigator>;

  const double    Bz = 2. * units::_T;
  BField          bField(0, 0, Bz);
  EigenStepper    estepper(bField);
  EigenPropagator epropagator(std::move(estepper), std::move(navigatorES));

  StraightLineStepper    slstepper;
  StraightLinePropagator slpropagator(std::move(slstepper),
                                      std::move(navigatorSL));
  const int ntests           = 500;
  const int skip             = 0;
  bool      debugModeFwd     = false;
  bool      debugModeBwd     = false;
  bool      debugModeFwdStep = false;
  bool      debugModeBwdStep = false;

  /// the actual test nethod that runs the test
  /// can be used with several propagator types
  /// @tparam propagator_t is the actual propagator type
  ///
  /// @param prop is the propagator instance
  /// @param pT the transverse momentum
  /// @param phi the azimuthal angle of the track at creation
  /// @param theta the polar angle of the track at creation
  /// @param charge is the charge of the particle
  /// @param index is the run index from the test
  template <typename propagator_t>
  void
  runTest(const propagator_t& prop,
          double              pT,
          double              phi,
          double              theta,
          int                 charge,
          int                 index)
  {
    double dcharge = -1 + 2 * charge;

    if (index < skip) {
      return;
    }

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

    using DebugOutput = detail::DebugOutputActor;

    // Action list and abort list
    using ActionListType = ActionList<MaterialInteractor, DebugOutput>;
    using AbortListType  = AbortList<>;

    using Options = PropagatorOptions<ActionListType, AbortListType>;
    Options fwdOptions;

    fwdOptions.maxStepSize = 25. * units::_cm;
    fwdOptions.pathLimit   = 25 * units::_cm;
    fwdOptions.debug       = debugModeFwd;

    // get the material collector and configure it
    auto& fwdMaterialInteractor
        = fwdOptions.actionList.template get<MaterialInteractor>();
    fwdMaterialInteractor.recordInteractions = true;
    fwdMaterialInteractor.energyLoss         = false;
    fwdMaterialInteractor.multipleScattering = false;

    if (debugModeFwd) {
      std::cout << ">>> Forward Propagation : start." << std::endl;
    }
    // forward material test
    const auto& fwdResult = prop.propagate(start, fwdOptions);
    auto&       fwdMaterial
        = fwdResult.template get<MaterialInteractor::result_type>();

    double fwdStepMaterialInX0 = 0.;
    double fwdStepMaterialInL0 = 0.;
    // check that the collected material is not zero
    BOOST_CHECK_NE(fwdMaterial.materialInX0, 0.);
    BOOST_CHECK_NE(fwdMaterial.materialInL0, 0.);
    // check that the sum of all steps is the total material
    for (auto& mInteraction : fwdMaterial.materialInteractions) {
      fwdStepMaterialInX0 += mInteraction.materialProperties.thicknessInX0();
      fwdStepMaterialInL0 += mInteraction.materialProperties.thicknessInL0();
    }
    CHECK_CLOSE_REL(fwdMaterial.materialInX0, fwdStepMaterialInX0, 1e-3);
    CHECK_CLOSE_REL(fwdMaterial.materialInL0, fwdStepMaterialInL0, 1e-3);

    // get the forward output to the screen
    if (debugModeFwd) {
      const auto& fwdOutput
          = fwdResult.template get<DebugOutput::result_type>();
      std::cout << ">>> Forward Propagation & Navigation output " << std::endl;
      std::cout << fwdOutput.debugString << std::endl;
      // check if the surfaces are free
      std::cout << ">>> Material steps found on ..." << std::endl;
      for (auto& fwdStepsC : fwdMaterial.materialInteractions) {
        std::cout << "--> Surface with "
                  << fwdStepsC.surface->geoID().toString() << std::endl;
      }
    }

    // backward material test
    Options bwdOptions;
    bwdOptions.maxStepSize = -25 * units::_cm;
    bwdOptions.pathLimit   = -25 * units::_cm;
    bwdOptions.direction   = backward;
    bwdOptions.debug       = debugModeBwd;

    // get the material collector and configure it
    auto& bwdMaterialInteractor
        = bwdOptions.actionList.template get<MaterialInteractor>();
    bwdMaterialInteractor.recordInteractions = true;
    bwdMaterialInteractor.energyLoss         = false;
    bwdMaterialInteractor.multipleScattering = false;

    const auto& startSurface = start.referenceSurface();

    if (debugModeBwd) {
      std::cout << ">>> Backward Propagation : start." << std::endl;
    }
    const auto& bwdResult = prop.propagate(
        *fwdResult.endParameters.get(), startSurface, bwdOptions);

    if (debugModeBwd) {
      std::cout << ">>> Backward Propagation : end." << std::endl;
    }

    auto& bwdMaterial
        = bwdResult.template get<typename MaterialInteractor::result_type>();

    double bwdStepMaterialInX0 = 0.;
    double bwdStepMaterialInL0 = 0.;

    // check that the collected material is not zero
    BOOST_CHECK_NE(bwdMaterial.materialInX0, 0.);
    BOOST_CHECK_NE(bwdMaterial.materialInL0, 0.);
    // check that the sum of all steps is the total material
    for (auto& mInteraction : bwdMaterial.materialInteractions) {
      bwdStepMaterialInX0 += mInteraction.materialProperties.thicknessInX0();
      bwdStepMaterialInL0 += mInteraction.materialProperties.thicknessInL0();
    }

    CHECK_CLOSE_REL(bwdMaterial.materialInX0, bwdStepMaterialInX0, 1e-3);
    CHECK_CLOSE_REL(bwdMaterial.materialInL0, bwdStepMaterialInL0, 1e-3);

    // get the backward output to the screen
    if (debugModeBwd) {
      const auto& bwd_output
          = bwdResult.template get<DebugOutput::result_type>();
      std::cout << ">>> Backward Propagation & Navigation output " << std::endl;
      std::cout << bwd_output.debugString << std::endl;
      // check if the surfaces are free
      std::cout << ">>> Material steps found on ..." << std::endl;
      for (auto& bwdStepsC : bwdMaterial.materialInteractions) {
        std::cout << "--> Surface with "
                  << bwdStepsC.surface->geoID().toString() << std::endl;
      }
    }

    // forward-backward compatibility test
    BOOST_CHECK_EQUAL(bwdMaterial.materialInteractions.size(),
                      fwdMaterial.materialInteractions.size());

    CHECK_CLOSE_REL(bwdMaterial.materialInX0, fwdMaterial.materialInX0, 1e-3);
    CHECK_CLOSE_REL(bwdMaterial.materialInL0, bwdMaterial.materialInL0, 1e-3);

    // stepping from one surface to the next
    // now go from surface to surface and check
    Options fwdStepOptions;
    fwdStepOptions.maxStepSize = 25. * units::_cm;
    fwdStepOptions.pathLimit   = 25 * units::_cm;
    fwdStepOptions.debug       = debugModeFwdStep;

    // get the material collector and configure it
    auto& fwdStepMaterialInteractor
        = fwdStepOptions.actionList.template get<MaterialInteractor>();
    fwdStepMaterialInteractor.recordInteractions = true;
    fwdStepMaterialInteractor.energyLoss         = false;
    fwdStepMaterialInteractor.multipleScattering = false;

    double fwdStepStepMaterialInX0 = 0.;
    double fwdStepStepMaterialInL0 = 0.;

    if (debugModeFwdStep) {
      // check if the surfaces are free
      std::cout << ">>> Forward steps to be processed sequentially ..."
                << std::endl;
      for (auto& fwdStepsC : fwdMaterial.materialInteractions) {
        std::cout << "--> Surface with "
                  << fwdStepsC.surface->geoID().toString() << std::endl;
      }
    }

    // move forward step by step through the surfaces
    const TrackParameters*              sParameters = &start;
    std::vector<const TrackParameters*> stepParameters;
    for (auto& fwdSteps : fwdMaterial.materialInteractions) {
      if (debugModeFwdStep) {
        std::cout << ">>> Forward step : "
                  << sParameters->referenceSurface().geoID().toString()
                  << " --> " << fwdSteps.surface->geoID().toString()
                  << std::endl;
      }

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
          = fwdStep.template get<typename MaterialInteractor::result_type>();
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

    if (debugModeFwdStep) {
      std::cout << ">>> Forward step : "
                << sParameters->referenceSurface().geoID().toString() << " --> "
                << dSurface.geoID().toString() << std::endl;
    }

    const auto& fwdStepFinal
        = prop.propagate(*sParameters, dSurface, fwdStepOptions);

    auto& fwdStepMaterial
        = fwdStepFinal.template get<typename MaterialInteractor::result_type>();
    fwdStepStepMaterialInX0 += fwdStepMaterial.materialInX0;
    fwdStepStepMaterialInL0 += fwdStepMaterial.materialInL0;

    // forward-forward step compatibility test
    CHECK_CLOSE_REL(fwdStepStepMaterialInX0, fwdStepMaterialInX0, 1e-3);
    CHECK_CLOSE_REL(fwdStepStepMaterialInL0, fwdStepMaterialInL0, 1e-3);

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
    Options bwdStepOptions;

    bwdStepOptions.maxStepSize = -25 * units::_cm;
    bwdStepOptions.pathLimit   = -25 * units::_cm;
    bwdStepOptions.direction   = backward;
    bwdStepOptions.debug       = debugModeBwdStep;

    // get the material collector and configure it
    auto& bwdStepMaterialInteractor
        = bwdStepOptions.actionList.template get<MaterialInteractor>();
    bwdStepMaterialInteractor.recordInteractions = true;
    bwdStepMaterialInteractor.multipleScattering = false;
    bwdStepMaterialInteractor.energyLoss         = false;

    double bwdStepStepMaterialInX0 = 0.;
    double bwdStepStepMaterialInL0 = 0.;

    if (debugModeBwdStep) {
      // check if the surfaces are free
      std::cout << ">>> Backward steps to be processed sequentially ..."
                << std::endl;
      for (auto& bwdStepsC : bwdMaterial.materialInteractions) {
        std::cout << "--> Surface with "
                  << bwdStepsC.surface->geoID().toString() << std::endl;
      }
    }

    // move forward step by step through the surfaces
    sParameters = fwdResult.endParameters.get();
    for (auto& bwdSteps : bwdMaterial.materialInteractions) {
      if (debugModeBwdStep) {
        std::cout << ">>> Backward step : "
                  << sParameters->referenceSurface().geoID().toString()
                  << " --> " << bwdSteps.surface->geoID().toString()
                  << std::endl;
      }
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
          = bwdStep.template get<typename MaterialInteractor::result_type>();
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

    if (debugModeBwdStep) {
      std::cout << ">>> Backward step : "
                << sParameters->referenceSurface().geoID().toString() << " --> "
                << dSurface.geoID().toString() << std::endl;
    }

    const auto& bwdStepFinal
        = prop.propagate(*sParameters, dbSurface, bwdStepOptions);

    auto& bwdStepMaterial
        = bwdStepFinal.template get<typename MaterialInteractor::result_type>();
    bwdStepStepMaterialInX0 += bwdStepMaterial.materialInX0;
    bwdStepStepMaterialInL0 += bwdStepMaterial.materialInL0;

    // forward-forward step compatibility test
    CHECK_CLOSE_REL(bwdStepStepMaterialInX0, bwdStepMaterialInX0, 1e-3);
    CHECK_CLOSE_REL(bwdStepStepMaterialInL0, bwdStepMaterialInL0, 1e-3);

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
