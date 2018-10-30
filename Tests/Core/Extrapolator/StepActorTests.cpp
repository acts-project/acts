// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///  Boost include(s)
#define BOOST_TEST_MODULE StepActor Tests

#include <boost/test/included/unit_test.hpp>
// leave blank line

#include <fstream>
#include "Acts/Detector/TrackingGeometry.hpp"
#include "Acts/Extrapolator/Navigator.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Propagator/DefaultExtension.hpp"
#include "Acts/Propagator/DenseEnvironmentExtension.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/detail/Auctioneer.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "DetectorBuilder.hpp"

// TMP
#include "../../Integration/PropagationTestHelper.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"

namespace Acts {
namespace Test {

  using cstep = detail::ConstrainedStep;

  ///
  /// @brief Aborter for the case that a particle leaves the detector
  ///
  struct EndOfWorld
  {
    /// Maximum value in x-direction of the detector
    double maxX = 1. * units::_m;

    /// @brief Constructor
    EndOfWorld() = default;

    /// @brief Main call operator for the abort operation
    ///
    /// @tparam propagator_state_t State of the propagator
    /// @param [in] state State of the propagation
    /// @return Boolean statement if the particle is still in the detector
    template <typename propagator_state_t>
    bool
    operator()(propagator_state_t& state) const
    {
      if (std::abs(state.stepping.position().x()) >= maxX
          || std::abs(state.stepping.position().y()) >= 0.5 * units::_m
          || std::abs(state.stepping.position().z()) >= 0.5 * units::_m)
        return true;
      return false;
    }
  };

  ///
  /// @brief Data collector while propagation
  ///
  struct StepCollector
  {

    ///
    /// @brief Data container for result analysis
    ///
    struct this_result
    {
      // Position of the propagator after each step
      std::vector<Vector3D> position;
      // Momentum of the propagator after each step
      std::vector<Vector3D> momentum;
      // Covariance matrix of the propagator after each step
      std::vector<ActsMatrixD<7, 7>> jac;
      // Step sizes of the propagator after each step
      std::vector<cstep> stepSize;
    };

    using result_type = this_result;

    /// @brief Main call operator for the action list. It stores the data for
    /// analysis afterwards
    ///
    /// @tparam propagator_state_t Type of the propagator state
    /// @param [in] state State of the propagator
    /// @param [out] result Struct which is filled with the data
    template <typename propagator_state_t>
    void
    operator()(propagator_state_t& state, result_type& result) const
    {
      result.position.push_back(state.stepping.position());
      result.momentum.push_back(state.stepping.momentum());
      result.jac.push_back(state.stepping.jacTransport);
      result.stepSize.push_back(state.stepping.stepSize);
    }
  };

  ActsMatrixD<7, 7> id77 = ActsMatrixD<7, 7>::Identity();
  bool output = false;

  /// @brief This function tests the EigenStepper with the DefaultExtension and
  /// the DenseEnvironmentExtension. The focus of this tests lies in the
  /// choosing of the right extension for the individual use case. This is
  /// performed with three different detectors:
  /// a) Pure vaccuum -> DefaultExtension needs to act
  /// b) Pure Be -> DenseEnvironmentExtension needs to act
  /// c) Vacuum - Be - Vacuum -> Both should act and switch during the
  /// propagation
  BOOST_AUTO_TEST_CASE(step_actor_test)
  {
    // Test case a). The DenseEnvironmentExtension should state that it is not
    // valid in this case.
    {
      // Build detector
      std::shared_ptr<TrackingGeometry> vacuum = buildVacDetector();

      // Build navigator
      Navigator naviVac(vacuum);
      naviVac.resolvePassive   = true;
      naviVac.resolveMaterial  = true;
      naviVac.resolveSensitive = true;

      // Set initial parameters for the particle track
      ActsSymMatrixD<5> cov    = ActsSymMatrixD<5>::Identity();
      auto              covPtr = std::make_unique<const ActsSymMatrixD<5>>(cov);
      Vector3D startParams(0., 0., 0.), startMom(1. * units::_GeV, 0., 0.);
      SingleCurvilinearTrackParameters<ChargedPolicy> sbtp(
          std::move(covPtr), startParams, startMom, 1.);

      // Create action list for surface collection
      ActionList<StepCollector,
                 DefaultExtensionActor,
                 DenseEnvironmentExtensionActor>
                            aList;
      AbortList<EndOfWorld> abortList;

      // Set options for propagator
      Propagator<EigenStepper<ConstantBField,
                              VoidCorrector,
                              ExtensionList<DefaultExtension,
                                            DenseEnvironmentExtension>>,
                 Navigator>::Options<ActionList<StepCollector,
                                                DefaultExtensionActor,
                                                DenseEnvironmentExtensionActor>,
                                     AbortList<EndOfWorld>>
          propOpts;
      propOpts.actionList     = aList;
      propOpts.stopConditions = abortList;
      propOpts.maxSteps       = 1e6;
      propOpts.maxStepSize    = 0.5 * units::_m;

      // Build stepper and propagator
      ConstantBField bField(Vector3D(0., 0., 0.));
      EigenStepper<ConstantBField,
                   VoidCorrector,
                   ExtensionList<DefaultExtension, DenseEnvironmentExtension>>
          es(bField);
      Propagator<EigenStepper<ConstantBField,
                              VoidCorrector,
                              ExtensionList<DefaultExtension,
                                            DenseEnvironmentExtension>>,
                 Navigator>
          prop(es, naviVac);

      // Launch and collect results
      const auto&                       result = prop.propagate(sbtp, propOpts);
      const StepCollector::this_result& stepResult
          = result.get<typename StepCollector::result_type>();

      // Check that the propagation happend without interactions
      for (const auto& pos : stepResult.position) {
        BOOST_TEST(pos.y() == 0.);
        BOOST_TEST(pos.z() == 0.);
        if (pos == stepResult.position.back())
          BOOST_TEST(pos.x() == 1. * units::_m);
      }
      for (const auto& mom : stepResult.momentum) {
        BOOST_TEST(mom.x() == 1. * units::_GeV);
        BOOST_TEST(mom.y() == 0.);
        BOOST_TEST(mom.z() == 0.);
      }

      // Rebuild and check the choice of extension
      ActionList<StepCollector, DefaultExtensionActor> aListDef;

      // Set options for propagator
      Propagator<EigenStepper<ConstantBField,
                              VoidCorrector,
                              ExtensionList<DefaultExtension>>,
                 Navigator>::Options<ActionList<StepCollector,
                                                DefaultExtensionActor>,
                                     AbortList<EndOfWorld>>
          propOptsDef;
      propOptsDef.actionList     = aListDef;
      propOptsDef.stopConditions = abortList;
      propOptsDef.maxSteps       = 1e6;
      propOptsDef.maxStepSize    = 0.5 * units::_m;

      EigenStepper<ConstantBField,
                   VoidCorrector,
                   ExtensionList<DefaultExtension>>
          esDef(bField);
      Propagator<EigenStepper<ConstantBField,
                              VoidCorrector,
                              ExtensionList<DefaultExtension>>,
                 Navigator>
          propDef(esDef, naviVac);

      // Launch and collect results
      const auto& resultDef = propDef.propagate(sbtp, propOptsDef);
      const StepCollector::this_result& stepResultDef
          = resultDef.get<typename StepCollector::result_type>();

      // Check that the right extension was chosen
      // If chosen correctly, the number of elements should be identical
      BOOST_TEST(stepResult.position.size() == stepResultDef.position.size());
      for (unsigned int i = 0; i < stepResult.position.size(); i++) {
        BOOST_TEST(stepResult.position[i] == stepResultDef.position[i]);
      }
      BOOST_TEST(stepResult.momentum.size() == stepResultDef.momentum.size());
      for (unsigned int i = 0; i < stepResult.momentum.size(); i++) {
        BOOST_TEST(stepResult.momentum[i] == stepResultDef.momentum[i]);
      }

      IntegrationTest::covariance_curvilinear(
          prop, startMom.x(), 0., M_PI * 0.5, 1., 0.5 * units::_m, 0);
    }
    // Test case b). The DefaultExtension should state that it is invalid here.
    {
      // Build detector
      std::shared_ptr<TrackingGeometry> material = buildMatDetector();

      // Build navigator
      Navigator naviMat(material);
      naviMat.resolvePassive   = true;
      naviMat.resolveMaterial  = true;
      naviMat.resolveSensitive = true;

      // Set initial parameters for the particle track
      ActsSymMatrixD<5> cov    = ActsSymMatrixD<5>::Identity();
      auto              covPtr = std::make_unique<const ActsSymMatrixD<5>>(cov);
      Vector3D          startParams(0., 0., 0.),
          startMom(5. * units::_GeV, 0., 0.);  // TODO: modified mom
      SingleCurvilinearTrackParameters<ChargedPolicy> sbtp(
          std::move(covPtr), startParams, startMom, 1.);

      // Create action list for surface collection
      ActionList<StepCollector,
                 DefaultExtensionActor,
                 DenseEnvironmentExtensionActor>
                            aList;
      AbortList<EndOfWorld> abortList;

      // Set options for propagator
      Propagator<EigenStepper<ConstantBField,
                              VoidCorrector,
                              ExtensionList<DefaultExtension,
                                            DenseEnvironmentExtension>>,
                 Navigator>::Options<ActionList<StepCollector,
                                                DefaultExtensionActor,
                                                DenseEnvironmentExtensionActor>,
                                     AbortList<EndOfWorld>>
          propOpts;
      propOpts.actionList     = aList;
      propOpts.stopConditions = abortList;
      propOpts.maxSteps       = 1e6;
      propOpts.maxStepSize    = 0.5 * units::_m;
      propOpts.debug          = true;

      // Build stepper and propagator
      ConstantBField bField(Vector3D(0., 0., 0.));
      EigenStepper<ConstantBField,
                   VoidCorrector,
                   ExtensionList<DefaultExtension, DenseEnvironmentExtension>>
          es(bField);
      Propagator<EigenStepper<ConstantBField,
                              VoidCorrector,
                              ExtensionList<DefaultExtension,
                                            DenseEnvironmentExtension>>,
                 Navigator>
          prop(es, naviMat);

      // Launch and collect results
      const auto&                       result = prop.propagate(sbtp, propOpts);
      const StepCollector::this_result& stepResult
          = result.get<typename StepCollector::result_type>();

      // Check that there occured interaction
      for (const auto& pos : stepResult.position) {
        BOOST_TEST(pos.y() == 0.);
        BOOST_TEST(pos.z() == 0.);
        if (pos == stepResult.position.front()) {
          BOOST_TEST(pos.x() == 0.);
        } else {
          BOOST_TEST(pos.x() != 0.);
        }
      }
      for (const auto& mom : stepResult.momentum) {
        BOOST_TEST(mom.y() == 0.);
        BOOST_TEST(mom.z() == 0.);
        if (mom == stepResult.momentum.front()) {
          BOOST_TEST(mom.x() == 5. * units::_GeV);
        } else {
          BOOST_TEST(mom.x() < 5. * units::_GeV);
        }
      }

      // Rebuild and check the choice of extension
      // Set options for propagator
      Propagator<EigenStepper<ConstantBField,
                              VoidCorrector,
                              ExtensionList<DenseEnvironmentExtension>>,
                 Navigator>::Options<ActionList<StepCollector,
                                                DefaultExtensionActor,
                                                DenseEnvironmentExtensionActor>,
                                     AbortList<EndOfWorld>>
          propOptsDense;
      propOptsDense.actionList     = aList;
      propOptsDense.stopConditions = abortList;
      propOptsDense.maxSteps       = 1e6;
      propOptsDense.maxStepSize    = 0.5 * units::_m;
      propOptsDense.debug          = true;

      // Build stepper and propagator
      EigenStepper<ConstantBField,
                   VoidCorrector,
                   ExtensionList<DenseEnvironmentExtension>>
          esDense(bField);
      Propagator<EigenStepper<ConstantBField,
                              VoidCorrector,
                              ExtensionList<DenseEnvironmentExtension>>,
                 Navigator>
          propDense(esDense, naviMat);

      // Launch and collect results
      const auto& resultDense = propDense.propagate(sbtp, propOptsDense);
      const StepCollector::this_result& stepResultDense
          = resultDense.get<typename StepCollector::result_type>();

      // Check that the right extension was chosen
      // If chosen correctly, the number of elements should be identical
      BOOST_TEST(stepResult.position.size() == stepResultDense.position.size());
      for (unsigned int i = 0; i < stepResult.position.size(); i++) {
        BOOST_TEST(stepResult.position[i] == stepResultDense.position[i]);
      }
      BOOST_TEST(stepResult.momentum.size() == stepResultDense.momentum.size());
      for (unsigned int i = 0; i < stepResult.momentum.size(); i++) {
        BOOST_TEST(stepResult.momentum[i] == stepResultDense.momentum[i]);
      }

      IntegrationTest::covariance_curvilinear(
          prop, startMom.x(), 0., M_PI * 0.5, 1., 0.5 * units::_m, 0);
    } /**
       ////////////////////////////////////////////////////////////////////

           // Re-launch the configuration with magnetic field
           bField.setField(0., 1. * units::_T, 0.);
           EigenStepper<ConstantBField,
                        VoidCorrector,
                        ExtensionList<DefaultExtension,
     DenseEnvironmentExtension>>
               esB(bField);
           Propagator<EigenStepper<ConstantBField,
                                   VoidCorrector,
                                   ExtensionList<DefaultExtension,
                                                 DenseEnvironmentExtension>>,
                      Navigator>
               propB(esB, naviMat);

           const auto& resultB = propB.propagate(sbtp, propOpts);
           const StepCollector::this_result& stepResultB
               = resultB.get<typename StepCollector::result_type>();

           // Check that there occured interaction
           for (const auto& pos : stepResultB.position) {
             if (pos == stepResultB.position.front()) {
               BOOST_TEST(pos.x() == 0.);
               BOOST_TEST(pos.y() == 0.);
               BOOST_TEST(pos.z() == 0.);
             } else {
               BOOST_TEST(pos.x() != 0.);
               BOOST_TEST(pos.y() == 0.);
               BOOST_TEST(pos.z() != 0.);
             }
           }
           for (const auto& mom : stepResultB.momentum) {
             if (mom == stepResultB.momentum.front()) {
               BOOST_TEST(mom.x() == 5. * units::_GeV);
               BOOST_TEST(mom.y() == 0.);
               BOOST_TEST(mom.z() == 0.);
             } else {
               BOOST_TEST(mom.x() != 5. * units::_GeV);
               BOOST_TEST(mom.y() == 0.);
               BOOST_TEST(mom.z() != 0.);
             }
           }
           IntegrationTest::covariance_curvilinear(
               propB, startMom.x(), 0., M_PI * 0.5, 1., 0.5 * units::_m, 0);
         }
         // Test case c). Both should be involved in their part of the detector
         {
           // Build detector
           std::shared_ptr<TrackingGeometry> det = buildVacMatVacDetector();

           // Build navigator
           Navigator naviDet(det);
           naviDet.resolvePassive   = true;
           naviDet.resolveMaterial  = true;
           naviDet.resolveSensitive = true;

           // Set initial parameters for the particle track
           ActsSymMatrixD<5> cov;
           cov << 1. * units::_mm, 0., 0., 0., 0., 0., 1. * units::_mm, 0., 0.,
     0.,
               0., 0., 1., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 1.;
           auto     covPtr = std::make_unique<const ActsSymMatrixD<5>>(cov);
           Vector3D startParams(0., 0., 0.), startMom(1. * units::_GeV, 0., 0.);
           SingleCurvilinearTrackParameters<ChargedPolicy> sbtp(
               std::move(covPtr), startParams, startMom, 1.);

           // Create action list for surface collection
           ActionList<StepCollector,
                      DefaultExtensionActor,
                      DenseEnvironmentExtensionActor>
                                 aList;
           AbortList<EndOfWorld> abortList;
           abortList.get<EndOfWorld>().maxX = 3. * units::_m;

           // Set options for propagator
           Propagator<EigenStepper<ConstantBField,
                                   VoidCorrector,
                                   ExtensionList<DefaultExtension,
                                                 DenseEnvironmentExtension>>,
                      Navigator>::Options<ActionList<StepCollector,
                                                     DefaultExtensionActor,
                                                     DenseEnvironmentExtensionActor>,
                                          AbortList<EndOfWorld>>
               propOpts;
           propOpts.actionList     = aList;
           propOpts.stopConditions = abortList;
           propOpts.maxSteps       = 1e6;
           propOpts.maxStepSize    = 0.5 * units::_m;

           // Re-configure propagation with B-field
           ConstantBField bField(Vector3D(0., 1. * units::_T, 0.));
           EigenStepper<ConstantBField,
                        VoidCorrector,
                        ExtensionList<DefaultExtension,
     DenseEnvironmentExtension>>
               es(bField);
           Propagator<EigenStepper<ConstantBField,
                                   VoidCorrector,
                                   ExtensionList<DefaultExtension,
                                                 DenseEnvironmentExtension>>,
                      Navigator>
               prop(es, naviDet);

           // Launch and collect results
           //~ const auto&                       result = prop.propagate(sbtp,
     propOpts);
           //~ const StepCollector::this_result& stepResult
               //~ = result.get<typename StepCollector::result_type>();
     // TODO: cov checks need to be in integration setup
     //~ for(size_t i = 1; i < 10; i++)
     //~ {
       //~ std::cout << i << std::endl;
         //~ IntegrationTest::covariance_curvilinear(
               //~ prop, startMom.x(), 0., M_PI * 0.5, 1., 0.1 * i * units::_m,
     0, 1e-3);
     //~ }
           // Check that the propagation step size is constrained and released
           // properly
           //~ for (unsigned int i = 0; i < stepResult.stepSize.size(); i++) {
           //~ if (stepResult.position[i].x() < 1. * units::_m
           //~ && std::abs(stepResult.position[i].z()) < 0.5)
           //~ BOOST_TEST(stepResult.stepSize[i].value(cstep::user)
           //~ == propOpts.maxStepSize);
           //~ if (stepResult.position[i].x() > 1. * units::_m
           //~ && stepResult.position[i].x() < 2. * units::_m
           //~ && std::abs(stepResult.position[i].z()) < 0.5)
           //~ BOOST_TEST(stepResult.stepSize[i].value(cstep::user)
           //~ == aList.get<StepActor>().maxStepSize);
           //~ if (stepResult.position[i].x() > 2. * units::_m
           //~ && std::abs(stepResult.position[i].z()) < 0.5)
           //~ BOOST_TEST(stepResult.stepSize[i].value(cstep::user)
           //~ == propOpts.maxStepSize);
           //~ std::cout << "pos: " << stepResult.position[i].x() << "\t"
           //~ << stepResult.position[i].y() << "\t"
           //~ << stepResult.position[i].z() << std::endl;
           //~ }
           //~
     //////////////////////////////////////////////////////////////////
           //~ std::ofstream ofs("out.txt");
           //~ for(unsigned int ss = 1; ss < 100; ss++)
           //~ {
           //~ for(unsigned int mom = 1; mom < 10; mom++)
           //~ {
           //~ for(double bfield = 0.5; bfield < 2.5; bfield += 0.1)
           //~ {
           //~ covPtr = std::make_unique<const ActsSymMatrixD<5>>(cov);
           //~ startMom = {mom * units::_GeV, 0., 0.};
           //~ SingleCurvilinearTrackParameters<ChargedPolicy> sbtp2(
           //~ std::move(covPtr), startParams, startMom, 1.);

           //~ aList.get<StepActor>().maxStepSize = ss;
           //~ propOpts.actionList = aList;

           //~ ConstantBField               bField2(Vector3D(0., bfield *
     units::_T,
           // 0.));
           //~ EigenStepper<ConstantBField> es2(bField2);
           //~ Propagator<EigenStepper<ConstantBField>, Navigator> prop2(es2,
           //~ naviDet);

           //~ const auto& result2 = prop2.propagate(sbtp2, propOpts);
           //~ const StepCollector::this_result& stepResult2 =
     result2.get<typename
           // StepCollector::result_type>();

           //~ for(unsigned int i = stepResult2.position.size() - 1; i > 0; i--)
           //~ if(stepResult2.position[i].x() <= 3000.)
           //~ {
           //~ std::cout << "posx: " << stepResult2.position[i].x() <<
     std::endl;
           //~ ofs << ss << " " << mom << " " << bfield << " " <<
           // stepResult2.position[i].x() << std::endl;
           //~ break;
           //~ }
           //~ }
           //~ }
           //~ }
           //~ ofs.close();
           //~
     //////////////////////////////////////////////////////////////////
         }
         * */
  }
}  // namespace Test
}  // namespace Acts
