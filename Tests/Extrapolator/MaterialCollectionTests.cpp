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
#include "ACTS/EventData/TrackParameters.hpp"
#include "ACTS/Extrapolator/MaterialInteractor.hpp"
#include "ACTS/Extrapolator/Navigator.hpp"
#include "ACTS/Extrapolator/SurfaceCollector.hpp"
#include "ACTS/MagneticField/ConstantBField.hpp"
#include "ACTS/Material/Material.hpp"
#include "ACTS/Material/MaterialCollector.hpp"
#include "ACTS/Propagator/ActionList.hpp"
#include "ACTS/Propagator/EigenStepper.hpp"
#include "ACTS/Propagator/Propagator.hpp"
#include "ACTS/Propagator/StraightLineStepper.hpp"
#include "ACTS/Propagator/detail/DebugOutputActor.hpp"
#include "ACTS/Surfaces/CylinderSurface.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/Units.hpp"
#include "ExtrapolatorTestGeometry.hpp"

namespace bdata = boost::unit_test::data;
namespace tt    = boost::test_tools;

namespace Acts {

namespace Test {

  // Global definitions
  // The path limit abort
  typedef detail::PathLimitReached path_limit;

  typedef ConstantBField                  BField;
  typedef EigenStepper<BField>            EigenStepper;
  typedef Propagator<EigenStepper>        EigenPropagator;
  typedef Propagator<StraightLineStepper> StraightLinePropagator;

  const double    Bz         = 0.;  // 2. * units::_T;
  const double    stepFactor = 1.;  // avoid overstepping
  BField          bField(0, 0, Bz);
  EigenStepper    estepper(bField);
  EigenPropagator epropagator(std::move(estepper));

  StraightLineStepper    slstepper;
  StraightLinePropagator slpropagator(std::move(slstepper));

  std::vector<std::unique_ptr<const Surface>> stepState;
  auto tGeometry = testGeometry<ModuleSurface>(stepState);

  const int ntests              = 1000;
  const int skip                = 0;
  bool      debug_mode_fwd      = false;
  bool      debug_mode_bwd      = false;
  bool      debug_mode_fwd_step = false;
  bool      debug_mode_bwd_step = false;

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
    typedef ActionList<Navigator, MaterialCollector, DebugOutput>
                        ActionList_type;
    typedef AbortList<> AbortConditions_type;

    typename Propagator_type::template Options<ActionList_type,
                                               AbortConditions_type>
        fwd_navigator_options;

    fwd_navigator_options.maxStepSize   = 25. * units::_cm;
    fwd_navigator_options.maxPathLength = 25 * units::_cm;
    fwd_navigator_options.debug         = debug_mode_fwd;

    // get the navigator and provide the TrackingGeometry
    auto& fwd_navigator
        = fwd_navigator_options.actionList.template get<Navigator>();
    fwd_navigator.trackingGeometry  = tGeometry;
    fwd_navigator.initialStepFactor = stepFactor;
    fwd_navigator.debug             = debug_mode_fwd;

    // get the material collector and configure it
    auto& fwd_materialCollector
        = fwd_navigator_options.actionList.template get<MaterialCollector>();
    fwd_materialCollector.detailedCollection = true;
    fwd_materialCollector.debug              = debug_mode_fwd;

    // forward material test
    const auto& fwd_result = prop.propagate(start, fwd_navigator_options);
    auto&       fwd_material
        = fwd_result.template get<MaterialCollector::result_type>();

    double fwd_step_materialInX0 = 0.;
    double fwd_step_materialInL0 = 0.;
    // check that the collected material is not zero
    BOOST_TEST(fwd_material.materialInX0 != 0.);
    BOOST_TEST(fwd_material.materialInL0 != 0.);
    // check that the sum of all steps is the total material
    for (auto& materialHit : fwd_material.collected) {
      auto material = materialHit.material;
      fwd_step_materialInX0 += materialHit.pathLength / material.X0();
      fwd_step_materialInL0 += materialHit.pathLength / material.L0();
    }
    BOOST_CHECK_CLOSE(fwd_material.materialInX0, fwd_step_materialInX0, 1e-5);
    BOOST_CHECK_CLOSE(fwd_material.materialInL0, fwd_step_materialInL0, 1e-5);

    // get the forward output to the screen
    if (debug_mode_fwd) {
      const auto& fwd_output
          = fwd_result.template get<DebugOutput::result_type>();
      std::cout << ">>> Forward Propgation & Navigation output " << std::endl;
      std::cout << fwd_output.debugString << std::endl;
      // check if the surfaces are free
      std::cout << ">>> Material steps found on ..." << std::endl;
      for (auto& fwd_steps_o : fwd_material.collected) {
        std::cout << "--> Surface with "
                  << fwd_steps_o.surface->geoID().toString() << std::endl;
      }
    }

    // backward material test
    typename Propagator_type::template Options<ActionList_type,
                                               AbortConditions_type>
        bwd_navigator_options;
    bwd_navigator_options.maxStepSize   = 25. * units::_cm;
    bwd_navigator_options.maxPathLength = 25 * units::_cm;
    bwd_navigator_options.direction     = backward;
    bwd_navigator_options.debug         = debug_mode_bwd;

    // get the backward navigator and provide the TrackingGeometry - for a
    // different logger
    auto& bwd_navigator
        = bwd_navigator_options.actionList.template get<Navigator>();
    bwd_navigator.trackingGeometry  = tGeometry;
    bwd_navigator.initialStepFactor = stepFactor;
    bwd_navigator.debug             = debug_mode_bwd;

    // get the material collector and configure it
    auto& bwd_materialCollector
        = bwd_navigator_options.actionList.template get<MaterialCollector>();
    bwd_materialCollector.detailedCollection = true;
    bwd_materialCollector.debug              = debug_mode_bwd;

    const auto& startSurface = start.referenceSurface();
    const auto& bwd_result
        = prop.propagate(*fwd_result.endParameters.template get(),
                         startSurface,
                         bwd_navigator_options);
    auto& bwd_material
        = bwd_result.template get<MaterialCollector::result_type>();

    double bwd_step_materialInX0 = 0.;
    double bwd_step_materialInL0 = 0.;

    // check that the collected material is not zero
    BOOST_TEST(bwd_material.materialInX0 != 0.);
    BOOST_TEST(bwd_material.materialInL0 != 0.);
    // check that the sum of all steps is the total material
    for (auto& materialHit : bwd_material.collected) {
      auto material = materialHit.material;
      bwd_step_materialInX0 += materialHit.pathLength / material.X0();
      bwd_step_materialInL0 += materialHit.pathLength / material.L0();
    }

    BOOST_CHECK_CLOSE(bwd_material.materialInX0, bwd_step_materialInX0, 1e-5);
    BOOST_CHECK_CLOSE(bwd_material.materialInL0, bwd_step_materialInL0, 1e-5);

    // get the backward output to the screen
    if (debug_mode_bwd) {
      const auto& bwd_output
          = bwd_result.template get<DebugOutput::result_type>();
      std::cout << ">>> Backward Propgation & Navigation output " << std::endl;
      std::cout << bwd_output.debugString << std::endl;
      // check if the surfaces are free
      std::cout << ">>> Material steps found on ..." << std::endl;
      for (auto& bwd_steps_o : bwd_material.collected) {
        std::cout << "--> Surface with "
                  << bwd_steps_o.surface->geoID().toString() << std::endl;
      }
    }

    // forward-backward compatibility test
    BOOST_TEST(bwd_material.collected.size() == fwd_material.collected.size());

    BOOST_CHECK_CLOSE(
        bwd_material.materialInX0, fwd_material.materialInX0, 1e-5);
    BOOST_CHECK_CLOSE(
        bwd_material.materialInL0, bwd_material.materialInL0, 1e-5);

    // stepping from one surface to the next
    // now go from surface to surface and check
    typename Propagator_type::template Options<ActionList_type,
                                               AbortConditions_type>
        fwdstep_navigator_options;

    fwdstep_navigator_options.maxStepSize   = 25. * units::_cm;
    fwdstep_navigator_options.maxPathLength = 25 * units::_cm;
    fwdstep_navigator_options.debug         = debug_mode_fwd_step;

    // get the navigator and provide the TrackingGeometry
    auto& fwdstep_navigator
        = fwdstep_navigator_options.actionList.template get<Navigator>();
    fwdstep_navigator.trackingGeometry  = tGeometry;
    fwdstep_navigator.initialStepFactor = stepFactor;
    fwdstep_navigator.debug             = debug_mode_fwd_step;

    // get the material collector and configure it
    auto& fwdstep_materialCollector = fwdstep_navigator_options.actionList
                                          .template get<MaterialCollector>();
    fwdstep_materialCollector.detailedCollection = true;
    fwdstep_materialCollector.debug              = debug_mode_fwd_step;

    double fwdstep_step_materialInX0 = 0.;
    double fwdstep_step_materialInL0 = 0.;

    if (debug_mode_fwd_step) {
      // check if the surfaces are free
      std::cout << ">>> Steps to be processed sequentially ..." << std::endl;
      for (auto& fwd_steps_o : fwd_material.collected) {
        std::cout << "--> Surface with "
                  << fwd_steps_o.surface->geoID().toString() << std::endl;
      }
    }

    // move forward step by step through the surfaces
    const TrackParameters*              sParameters = &start;
    std::vector<const TrackParameters*> stepParameters;
    for (auto& fwd_steps : fwd_material.collected) {
      if (debug_mode_bwd_step)
        std::cout << ">>> Step : "
                  << sParameters->referenceSurface().geoID().toString()
                  << " --> " << fwd_steps.surface->geoID().toString()
                  << std::endl;

      // make a forward step
      const auto& fwd_step = prop.propagate(
          *sParameters, (*fwd_steps.surface), fwdstep_navigator_options);
      // get the backward output to the screen
      if (debug_mode_fwd_step) {
        const auto& fwdstep_output
            = fwd_step.template get<DebugOutput::result_type>();
        std::cout << fwdstep_output.debugString << std::endl;
      }

      auto& fwdstep_material
          = fwd_step.template get<MaterialCollector::result_type>();
      fwdstep_step_materialInX0 += fwdstep_material.materialInX0;
      fwdstep_step_materialInL0 += fwdstep_material.materialInL0;

      if (fwd_step.endParameters != nullptr) {
        sParameters = fwd_step.endParameters->clone();
        // make sure the parameters do not run out of scope
        stepParameters.push_back(sParameters);
      }
    }
    // final destination surface
    const Surface& dSurface = fwd_result.endParameters->referenceSurface();

    if (debug_mode_fwd_step)
      std::cout << ">>> Step : "
                << sParameters->referenceSurface().geoID().toString() << " --> "
                << dSurface.geoID().toString() << std::endl;

    const auto& fwdstep_final
        = prop.propagate(*sParameters, dSurface, fwdstep_navigator_options);

    auto& fwdstep_material
        = fwdstep_final.template get<MaterialCollector::result_type>();
    fwdstep_step_materialInX0 += fwdstep_material.materialInX0;
    fwdstep_step_materialInL0 += fwdstep_material.materialInL0;

    // forward-forward step compatibility test
    BOOST_CHECK_CLOSE(fwdstep_step_materialInX0, fwd_step_materialInX0, 1e-5);
    BOOST_CHECK_CLOSE(fwdstep_step_materialInL0, fwd_step_materialInL0, 1e-5);

    // get the backward output to the screen
    if (debug_mode_fwd_step) {
      const auto& fwdstep_output
          = fwdstep_final.template get<DebugOutput::result_type>();
      std::cout << ">>> Forward Final Step Propgation & Navigation output "
                << std::endl;
      std::cout << fwdstep_output.debugString << std::endl;
    }

    // stepping from one surface to the next : backwards
    // now go from surface to surface and check
    typename Propagator_type::template Options<ActionList_type,
                                               AbortConditions_type>
        bwdstep_navigator_options;

    bwdstep_navigator_options.maxStepSize   = 25. * units::_cm;
    bwdstep_navigator_options.maxPathLength = 25 * units::_cm;
    bwdstep_navigator_options.direction     = backward;
    bwdstep_navigator_options.debug         = debug_mode_bwd_step;

    // get the navigator and provide the TrackingGeometry
    auto& bwdstep_navigator
        = bwdstep_navigator_options.actionList.template get<Navigator>();
    bwdstep_navigator.trackingGeometry  = tGeometry;
    bwdstep_navigator.initialStepFactor = stepFactor;
    bwdstep_navigator.debug             = debug_mode_bwd_step;

    // get the material collector and configure it
    auto& bwdstep_materialCollector = bwdstep_navigator_options.actionList
                                          .template get<MaterialCollector>();
    bwdstep_materialCollector.detailedCollection = true;
    bwdstep_materialCollector.debug              = debug_mode_bwd_step;

    double bwdstep_step_materialInX0 = 0.;
    double bwdstep_step_materialInL0 = 0.;

    if (debug_mode_bwd_step) {
      // check if the surfaces are free
      std::cout << ">>> Steps to be processed sequentially ..." << std::endl;
      for (auto& bwd_steps_o : bwd_material.collected) {
        std::cout << "--> Surface with "
                  << bwd_steps_o.surface->geoID().toString() << std::endl;
      }
    }

    // move forward step by step through the surfaces
    sParameters = fwd_result.endParameters.template get();
    for (auto& bwd_steps : bwd_material.collected) {
      if (debug_mode_bwd_step)
        std::cout << ">>> Step : "
                  << sParameters->referenceSurface().geoID().toString()
                  << " --> " << bwd_steps.surface->geoID().toString()
                  << std::endl;
      // make a forward step
      const auto& bwd_step = prop.propagate(
          *sParameters, (*bwd_steps.surface), bwdstep_navigator_options);
      // get the backward output to the screen
      if (debug_mode_bwd_step) {
        const auto& bwdstep_output
            = bwd_step.template get<DebugOutput::result_type>();
        std::cout << bwdstep_output.debugString << std::endl;
      }

      auto& bwdstep_material
          = bwd_step.template get<MaterialCollector::result_type>();
      bwdstep_step_materialInX0 += bwdstep_material.materialInX0;
      bwdstep_step_materialInL0 += bwdstep_material.materialInL0;

      if (bwd_step.endParameters != nullptr) {
        sParameters = bwd_step.endParameters->clone();
        // make sure the parameters do not run out of scope
        stepParameters.push_back(sParameters);
      }
    }
    // final destination surface
    const Surface& dbSurface = start.referenceSurface();

    if (debug_mode_bwd_step)
      std::cout << ">>> Step : "
                << sParameters->referenceSurface().geoID().toString() << " --> "
                << dSurface.geoID().toString() << std::endl;

    const auto& bwdstep_final
        = prop.propagate(*sParameters, dbSurface, bwdstep_navigator_options);

    auto& bwdstep_material
        = bwdstep_final.template get<MaterialCollector::result_type>();
    bwdstep_step_materialInX0 += bwdstep_material.materialInX0;
    bwdstep_step_materialInL0 += bwdstep_material.materialInL0;

    // forward-forward step compatibility test
    BOOST_CHECK_CLOSE(bwdstep_step_materialInX0, bwd_step_materialInX0, 1e-5);
    BOOST_CHECK_CLOSE(bwdstep_step_materialInL0, bwd_step_materialInL0, 1e-5);

    // get the backward output to the screen
    if (debug_mode_bwd_step) {
      const auto& bwdstep_output
          = bwdstep_final.template get<DebugOutput::result_type>();
      std::cout << ">>> Forward Final Step Propgation & Navigation output "
                << std::endl;
      std::cout << bwdstep_output.debugString << std::endl;
    }
  }

  // This test case checks that no segmentation fault appears
  // - this tests the collection of surfaces
  BOOST_DATA_TEST_CASE(
      test_material_collector,
      bdata::random((bdata::seed = 20,
                     bdata::distribution
                     = std::uniform_real_distribution<>(0.4 * units::_GeV,
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
