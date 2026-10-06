// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Direction.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/MagneticField/MagneticFieldError.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Propagator/ButcherTableau.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Propagator/GenericRungeKuttaStepper.hpp"
#include "Acts/Propagator/HelixStepper.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <cmath>
#include <functional>
#include <limits>
#include <memory>
#include <optional>
#include <ostream>
#include <stdexcept>
#include <vector>

using namespace Acts;
using namespace Acts::UnitLiterals;
using Acts::VectorHelpers::makeVector4;

// The test contexts print the name of the tableau
BOOST_TEST_DONT_PRINT_LOG_VALUE(std::shared_ptr<const Acts::ButcherTableau>)
BOOST_TEST_DONT_PRINT_LOG_VALUE(Acts::GenericRungeKuttaStepper::ErrorEstimation)

namespace ActsTests {

namespace {

const auto tgContext = GeometryContext::dangerouslyDefaultConstruct();
const MagneticFieldContext mfContext;

constexpr double infinity = std::numeric_limits<double>::infinity();

/// A smooth, non-uniform field
///
///   B = (s a y z, -s a x z, B0 + s a (x^2 + y^2)), s = B0 / L^2, a = 1/2
class SmoothField final : public MagneticFieldProvider {
 public:
  struct Cache {
    explicit Cache(const MagneticFieldContext& /*mctx*/) {}
  };

  SmoothField(double b0, double length) : m_b0(b0), m_length(length) {}

  MagneticFieldProvider::Cache makeCache(
      const MagneticFieldContext& mctx) const override {
    return MagneticFieldProvider::Cache(std::in_place_type<Cache>, mctx);
  }

  Result<Vector3> getField(
      const Vector3& p,
      MagneticFieldProvider::Cache& /*cache*/) const override {
    constexpr double a = 0.5;
    const double s = m_b0 / (m_length * m_length);
    return Result<Vector3>::success(
        Vector3(s * a * p.y() * p.z(), -s * a * p.x() * p.z(),
                m_b0 + s * a * (p.x() * p.x() + p.y() * p.y())));
  }

 private:
  double m_b0;
  double m_length;
};

/// A constant field that fails outside of |x| < xMax
class BoundedField final : public MagneticFieldProvider {
 public:
  struct Cache {
    explicit Cache(const MagneticFieldContext& /*mctx*/) {}
  };

  BoundedField(const Vector3& field, double xMax)
      : m_field(field), m_xMax(xMax) {}

  MagneticFieldProvider::Cache makeCache(
      const MagneticFieldContext& mctx) const override {
    return MagneticFieldProvider::Cache(std::in_place_type<Cache>, mctx);
  }

  Result<Vector3> getField(
      const Vector3& p,
      MagneticFieldProvider::Cache& /*cache*/) const override {
    if (std::abs(p.x()) >= m_xMax) {
      return Result<Vector3>::failure(MagneticFieldError::OutOfBounds);
    }
    return Result<Vector3>::success(m_field);
  }

 private:
  Vector3 m_field;
  double m_xMax;
};

FreeVector makeStart(const Vector3& pos, const Vector3& dir, double qop) {
  FreeVector start = FreeVector::Zero();
  start.segment<3>(eFreePos0) = pos;
  start[eFreeTime] = 1.;
  start.segment<3>(eFreeDir0) = dir.normalized();
  start[eFreeQOverP] = qop;
  return start;
}

/// Propagate the free parameters over the path length and return the state.
template <typename stepper_t>
typename stepper_t::State propagateFree(
    const stepper_t& stepper, typename stepper_t::Options options,
    const FreeVector& start, double pathLength, bool covTransport,
    const ParticleHypothesis& particleHypothesis = ParticleHypothesis::pion()) {
  typename stepper_t::State state = stepper.makeState(options);
  // Initialise from the start parameters, so that the direction columns of
  // the bound-to-free jacobian are orthogonal to the start direction.
  stepper.initialize(
      state,
      BoundTrackParameters::createCurvilinear(
          makeVector4(start.segment<3>(eFreePos0), start[eFreeTime]),
          start.segment<3>(eFreeDir0), start[eFreeQOverP],
          covTransport ? std::optional<BoundMatrix>(BoundMatrix::Identity())
                       : std::nullopt,
          particleHypothesis));
  // Keep the exact start values, the round trip through the bound
  // parameters changes them at the level of the round-off.
  state.pars = start;
  while (std::abs(pathLength - state.pathAccumulated) > 1e-12) {
    state.stepSize.release(ConstrainedStep::Type::Navigator);
    state.stepSize.update(pathLength - state.pathAccumulated,
                          ConstrainedStep::Type::Navigator);
    auto res = stepper.step(state, Direction::fromScalar(pathLength), nullptr);
    BOOST_REQUIRE(res.ok());
  }
  return state;
}

GenericRungeKuttaStepper::Options fixedStepOptions(double stepSize) {
  GenericRungeKuttaStepper::Options options(tgContext, mfContext);
  options.errorEstimation = GenericRungeKuttaStepper::ErrorEstimation::None;
  options.initialStepSize = stepSize;
  options.maxStepSize = stepSize;
  return options;
}

}  // namespace

BOOST_AUTO_TEST_SUITE(PropagatorSuite)

BOOST_AUTO_TEST_CASE(butcher_tableau_validation) {
  BOOST_CHECK_THROW(ButcherTableau("sizes", 1, 0, {0., 1.}, {{}}, {1.}, {}),
                    std::invalid_argument);
  BOOST_CHECK_THROW(
      ButcherTableau("node", 1, 0, {0., 1.}, {{}, {0.5}}, {0.5, 0.5}, {}),
      std::invalid_argument);
  BOOST_CHECK_THROW(ButcherTableau("first", 1, 0, {1.}, {{}}, {1.}, {}),
                    std::invalid_argument);
  BOOST_CHECK_NO_THROW(
      ButcherTableau("heun", 2, 1, {0., 1.}, {{}, {1.}}, {0.5, 0.5}, {1., 0.}));
}

BOOST_AUTO_TEST_CASE(generic_runge_kutta_stepper_without_tableau) {
  auto field = std::make_shared<ConstantBField>(Vector3(0., 0., 2_T));
  BOOST_CHECK_THROW(GenericRungeKuttaStepper(
                        GenericRungeKuttaStepper::Config{field, nullptr}),
                    std::invalid_argument);
}

/// Check the order conditions of the built-in tableaus up to order 5. The
/// Verner 9(8) coefficients were checked against all conditions up to order 9
/// with 60-digit arithmetic when they were added.
BOOST_AUTO_TEST_CASE(butcher_tableau_order_conditions) {
  using Vec = Eigen::VectorXd;
  using Mat = Eigen::MatrixXd;

  // Value of each condition for the weights w, with its order and its
  // expected value 1 / gamma(tree)
  struct Condition {
    unsigned order;
    double expected;
    std::function<double(const Vec&, const Vec&, const Mat&)> value;
  };
  const std::vector<Condition> conditions = {
      {1, 1., [](auto& w, auto&, auto&) { return w.sum(); }},
      {2, 1. / 2, [](auto& w, auto& c, auto&) { return w.dot(c); }},
      {3, 1. / 3,
       [](auto& w, auto& c, auto&) { return w.dot(Vec(c.cwiseProduct(c))); }},
      {3, 1. / 6, [](auto& w, auto& c, auto& a) { return w.dot(Vec(a * c)); }},
      {4, 1. / 4,
       [](auto& w, auto& c, auto&) {
         return w.dot(Vec(c.array().pow(3).matrix()));
       }},
      {4, 1. / 8,
       [](auto& w, auto& c, auto& a) {
         return w.dot(Vec(c.cwiseProduct(a * c)));
       }},
      {4, 1. / 12,
       [](auto& w, auto& c, auto& a) {
         return w.dot(Vec(a * c.cwiseProduct(c)));
       }},
      {4, 1. / 24,
       [](auto& w, auto& c, auto& a) { return w.dot(Vec(a * a * c)); }},
      {5, 1. / 5,
       [](auto& w, auto& c, auto&) {
         return w.dot(Vec(c.array().pow(4).matrix()));
       }},
      {5, 1. / 10,
       [](auto& w, auto& c, auto& a) {
         return w.dot(Vec(c.cwiseProduct(c).cwiseProduct(a * c)));
       }},
      {5, 1. / 20,
       [](auto& w, auto& c, auto& a) {
         return w.dot(Vec((a * c).cwiseProduct(a * c)));
       }},
      {5, 1. / 15,
       [](auto& w, auto& c, auto& a) {
         return w.dot(Vec(c.cwiseProduct(a * c.cwiseProduct(c))));
       }},
      {5, 1. / 30,
       [](auto& w, auto& c, auto& a) {
         return w.dot(Vec(c.cwiseProduct(a * a * c)));
       }},
      {5, 1. / 20,
       [](auto& w, auto& c, auto& a) {
         return w.dot(Vec(a * c.array().pow(3).matrix()));
       }},
      {5, 1. / 40,
       [](auto& w, auto& c, auto& a) {
         return w.dot(Vec(a * c.cwiseProduct(a * c)));
       }},
      {5, 1. / 60,
       [](auto& w, auto& c, auto& a) {
         return w.dot(Vec(a * a * c.cwiseProduct(c)));
       }},
      {5, 1. / 120,
       [](auto& w, auto& c, auto& a) { return w.dot(Vec(a * a * a * c)); }},
  };

  for (const auto& tableau :
       {ButcherTableau::classicalRk4(), ButcherTableau::dormandPrince54(),
        ButcherTableau::verner98()}) {
    const std::size_t s = tableau->stages();
    Vec c(s);
    Vec b(s);
    Vec bEmbedded = Vec::Zero(s);
    Mat a = Mat::Zero(s, s);
    for (std::size_t i = 0; i < s; ++i) {
      c[i] = tableau->c(i);
      b[i] = tableau->b(i);
      if (tableau->hasEmbedded()) {
        bEmbedded[i] = tableau->bEmbedded(i);
      }
      for (std::size_t j = 0; j < i; ++j) {
        a(i, j) = tableau->a(i, j);
      }
    }

    BOOST_TEST_CONTEXT(tableau->name()) {
      for (const auto& condition : conditions) {
        if (condition.order <= tableau->order()) {
          CHECK_CLOSE_ABS(condition.value(b, c, a), condition.expected, 1e-14);
        }
        if (tableau->hasEmbedded() &&
            condition.order <= tableau->embeddedOrder()) {
          CHECK_CLOSE_ABS(condition.value(bEmbedded, c, a), condition.expected,
                          1e-14);
        }
      }
      // The solution does not have a higher order than stated
      if (tableau->order() == 4) {
        BOOST_CHECK_GT(std::abs(conditions[8].value(b, c, a) - 1. / 5), 1e-6);
      }
      // The embedded solution differs from the solution at the next order
      if (tableau->hasEmbedded() && tableau->embeddedOrder() < 5) {
        bool differs = false;
        for (const auto& condition : conditions) {
          if (condition.order == tableau->embeddedOrder() + 1) {
            differs |= std::abs(condition.value(bEmbedded, c, a) -
                                condition.expected) > 1e-6;
          }
        }
        BOOST_CHECK(differs);
      }
    }
  }
}

/// In a constant field the reference stepper reproduces the exact helix.
BOOST_AUTO_TEST_CASE(generic_runge_kutta_stepper_constant_field) {
  auto field = std::make_shared<ConstantBField>(Vector3(0.1_T, -0.2_T, 2_T));
  const FreeVector start =
      makeStart(Vector3(1., 2., 3.), Vector3(4., -5., 6.), -1. / 0.7_GeV);
  const double pathLength = 1_m;

  HelixStepper helix(field);
  HelixStepper::Options helixOptions(tgContext, mfContext);
  auto helixState = propagateFree(helix, helixOptions, start, pathLength, true);

  GenericRungeKuttaStepper reference(field);
  GenericRungeKuttaStepper::Options options(tgContext, mfContext);
  options.stepTolerance = 1e-12;
  options.initialStepSize = 1_cm;
  auto state = propagateFree(reference, options, start, pathLength, true);

  CHECK_CLOSE_ABS(state.pars, helixState.pars, 1e-9);
  CHECK_CLOSE_ABS(state.derivative, helixState.derivative, 1e-9);

  // The reference projects the jacobian of the direction onto the unit sphere
  // after each step, so compare the curvilinear jacobians.
  const BoundMatrix refJac = reference.transportToCurvilinear(state);
  const BoundMatrix helixJac = helix.transportToCurvilinear(helixState);
  CHECK_CLOSE_OR_SMALL(refJac, helixJac, 1e-8, 1e-8);
}

namespace {

struct ConvergenceCase {
  std::shared_ptr<const ButcherTableau> tableau;
  double qop{};
  std::vector<int> stepCounts;

  friend std::ostream& operator<<(std::ostream& os,
                                  const ConvergenceCase& testCase) {
    return os << testCase.tableau->name();
  }
};

// The ninth-order error needs a strongly bent track to stay above the
// round-off while the step size halves
const std::vector<ConvergenceCase> convergenceCases = {
    {ButcherTableau::classicalRk4(), 1. / 1_GeV, {8, 16, 32}},
    {ButcherTableau::dormandPrince54(), 1. / 1_GeV, {8, 16, 32}},
    {ButcherTableau::verner98(), 1. / 0.1_GeV, {8, 12, 18}},
};

}  // namespace

/// The fixed-step global error falls with the order of the tableau.
BOOST_DATA_TEST_CASE(generic_runge_kutta_stepper_convergence,
                     boost::unit_test::data::make(convergenceCases), testCase) {
  auto field = std::make_shared<SmoothField>(2_T, 1_m);
  const FreeVector start =
      makeStart(Vector3(100., -50., 20.), Vector3(1., 0.3, 0.2), testCase.qop);
  const double pathLength = 1_m;

  // Every consistent tableau converges to the same solution, so a wrong order
  // still shows in the ratios.
  GenericRungeKuttaStepper::Config config{field, ButcherTableau::verner98()};
  const GenericRungeKuttaStepper reference(config);
  const FreeVector referencePars =
      propagateFree(reference, fixedStepOptions(pathLength / 1024), start,
                    pathLength, false)
          .pars;

  config.tableau = testCase.tableau;
  const GenericRungeKuttaStepper stepper(config);
  std::vector<double> errors;
  for (int nSteps : testCase.stepCounts) {
    const FreeVector end =
        propagateFree(stepper, fixedStepOptions(pathLength / nSteps), start,
                      pathLength, false)
            .pars;
    errors.push_back(
        (end - referencePars).segment<3>(eFreePos0).cwiseAbs().maxCoeff());
  }

  BOOST_TEST_CONTEXT(testCase.tableau->name()
                     << " errors " << errors[0] << " " << errors[1] << " "
                     << errors[2]) {
    for (std::size_t i = 0; i + 1 < errors.size(); ++i) {
      const double ratio = errors[i] / errors[i + 1];
      const double expectedRatio =
          std::pow(static_cast<double>(testCase.stepCounts[i + 1]) /
                       testCase.stepCounts[i],
                   testCase.tableau->order());
      BOOST_CHECK_GT(ratio, 0.6 * expectedRatio);
      BOOST_CHECK_LT(ratio, 1.6 * expectedRatio);
    }
  }
}

/// The chain-rule jacobian matches finite differences of the full
/// propagation, and it needs the field gradient to do so.
BOOST_AUTO_TEST_CASE(generic_runge_kutta_stepper_jacobian) {
  const FreeVector start =
      makeStart(Vector3(100., -50., 20.), Vector3(1., 0.3, 0.2), -1. / 1_GeV);
  const double pathLength = 1_m;
  const auto options = fixedStepOptions(pathLength / 16);

  auto field = std::make_shared<SmoothField>(2_T, 1_m);
  const GenericRungeKuttaStepper stepper(field);

  FreeMatrix numeric;
  for (std::size_t j = 0; j < eFreeSize; ++j) {
    // A larger delta keeps the round-off of positions of order 1 m small
    const double delta =
        j == eFreeQOverP ? 1e-4 * std::abs(start[eFreeQOverP]) : 1e-4;
    FreeVector plus = start;
    FreeVector minus = start;
    plus[j] += delta;
    minus[j] -= delta;
    numeric.col(j) =
        (propagateFree(stepper, options, plus, pathLength, false).pars -
         propagateFree(stepper, options, minus, pathLength, false).pars) /
        (2. * delta);
  }

  const auto analytic =
      propagateFree(stepper, options, start, pathLength, true);
  // Both normalise the end direction, so the projected jacobian compares
  // directly.
  CHECK_CLOSE_OR_SMALL(analytic.jacTransport, numeric, 1e-5, 1e-5);

  // Without the gradient term the jacobian is wrong
  auto noGradientOptions = options;
  noGradientOptions.includeFieldGradient = false;
  const auto noGradient =
      propagateFree(stepper, noGradientOptions, start, pathLength, true);
  BOOST_CHECK_GT((noGradient.jacTransport - numeric).cwiseAbs().maxCoeff(),
                 1e-3);
}

namespace {

using ErrorEstimation = GenericRungeKuttaStepper::ErrorEstimation;

struct AdaptiveCase {
  std::shared_ptr<const ButcherTableau> tableau;
  ErrorEstimation errorEstimation{};

  friend std::ostream& operator<<(std::ostream& os,
                                  const AdaptiveCase& testCase) {
    return os << testCase.tableau->name() << " "
              << (testCase.errorEstimation == ErrorEstimation::Embedded
                      ? "embedded"
                      : "step doubling");
  }
};

const std::vector<AdaptiveCase> adaptiveCases = {
    {ButcherTableau::dormandPrince54(), ErrorEstimation::Embedded},
    {ButcherTableau::verner98(), ErrorEstimation::Embedded},
    {ButcherTableau::classicalRk4(), ErrorEstimation::StepDoubling},
    {ButcherTableau::dormandPrince54(), ErrorEstimation::StepDoubling},
};

}  // namespace

/// The adaptive step size meets the tolerance.
BOOST_DATA_TEST_CASE(generic_runge_kutta_stepper_adaptive,
                     boost::unit_test::data::make(adaptiveCases), testCase) {
  auto field = std::make_shared<SmoothField>(2_T, 1_m);
  const FreeVector start =
      makeStart(Vector3(100., -50., 20.), Vector3(1., 0.3, 0.2), 1. / 1_GeV);
  const double pathLength = 1_m;

  GenericRungeKuttaStepper::Config config{field, ButcherTableau::verner98()};
  // There is no closed-form solution in this field. The reference is the
  // fixed-step Verner 9(8) solution with 1024 steps. Its error is many orders
  // of magnitude below the tolerances, also for the adaptive Verner 9(8)
  // case, which takes far fewer steps.
  const FreeVector referencePars =
      propagateFree(GenericRungeKuttaStepper(config),
                    fixedStepOptions(pathLength / 1024), start, pathLength,
                    false)
          .pars;

  config.tableau = testCase.tableau;
  const GenericRungeKuttaStepper stepper(config);
  GenericRungeKuttaStepper::Options options(tgContext, mfContext);
  options.errorEstimation = testCase.errorEstimation;
  options.initialStepSize = 10_m;
  double previousError = infinity;
  for (double tolerance : {1e-4, 1e-7, 1e-10}) {
    options.stepTolerance = tolerance;
    const auto state =
        propagateFree(stepper, options, start, pathLength, false);
    const double error = (state.pars - referencePars)
                             .segment<3>(eFreePos0)
                             .cwiseAbs()
                             .maxCoeff();
    BOOST_TEST_CONTEXT(testCase << " tolerance " << tolerance << " error "
                                << error << " steps " << state.nSteps) {
      BOOST_CHECK_GT(state.statistics.nRejectedSteps, 0u);
      // The global error is the sum of the local errors
      BOOST_CHECK_LT(error, 10. * tolerance * state.nSteps);
      BOOST_CHECK_LT(error, previousError);
    }
    previousError = error;
  }
}

/// The embedded error estimate needs a tableau with embedded weights.
BOOST_AUTO_TEST_CASE(generic_runge_kutta_stepper_embedded_without_weights) {
  auto field = std::make_shared<SmoothField>(2_T, 1_m);
  const GenericRungeKuttaStepper stepper(
      GenericRungeKuttaStepper::Config{field, ButcherTableau::classicalRk4()});
  GenericRungeKuttaStepper::Options options(tgContext, mfContext);
  BOOST_REQUIRE(options.errorEstimation == ErrorEstimation::Embedded);
  BOOST_CHECK_THROW(stepper.makeState(options), std::invalid_argument);

  // The options in the state can change after makeState
  options.errorEstimation = ErrorEstimation::StepDoubling;
  auto state = stepper.makeState(options);
  stepper.initialize(state, BoundTrackParameters::createCurvilinear(
                                Vector4::Zero(), Vector3::UnitX(), 1. / 1_GeV,
                                std::nullopt, ParticleHypothesis::pion()));
  state.options.errorEstimation = ErrorEstimation::Embedded;
  BOOST_CHECK_THROW(
      static_cast<void>(stepper.step(state, Direction::Forward(), nullptr)),
      std::invalid_argument);
}

/// Without rejections, step doubling equals fixed steps of h/2.
BOOST_AUTO_TEST_CASE(generic_runge_kutta_stepper_step_doubling_halves) {
  auto field = std::make_shared<SmoothField>(2_T, 1_m);
  const GenericRungeKuttaStepper stepper(
      GenericRungeKuttaStepper::Config{field, ButcherTableau::classicalRk4()});
  const FreeVector start =
      makeStart(Vector3(100., -50., 20.), Vector3(1., 0.3, 0.2), 1. / 1_GeV);

  GenericRungeKuttaStepper::Options options(tgContext, mfContext);
  options.errorEstimation = ErrorEstimation::StepDoubling;
  options.stepTolerance = infinity;
  options.initialStepSize = 10_cm;
  options.maxStepSize = 10_cm;
  const auto doubling = propagateFree(stepper, options, start, 1_m, true);
  const auto halves =
      propagateFree(stepper, fixedStepOptions(5_cm), start, 1_m, true);

  BOOST_CHECK_EQUAL(doubling.nSteps, 10u);
  BOOST_CHECK_EQUAL(halves.nSteps, 20u);
  BOOST_CHECK_EQUAL(doubling.statistics.nRejectedSteps, 0u);
  CHECK_CLOSE_ABS(doubling.pars, halves.pars, 1e-12);
  CHECK_CLOSE_OR_SMALL(doubling.jacTransport, halves.jacTransport, 1e-12,
                       1e-12);
}

/// The jacobian in the error estimate makes the adaptive steps smaller.
BOOST_DATA_TEST_CASE(generic_runge_kutta_stepper_jacobian_in_error_estimate,
                     boost::unit_test::data::make(std::vector<ErrorEstimation>{
                         ErrorEstimation::Embedded,
                         ErrorEstimation::StepDoubling}),
                     errorEstimation) {
  auto field = std::make_shared<SmoothField>(2_T, 1_m);
  const GenericRungeKuttaStepper stepper(field);
  GenericRungeKuttaStepper::Options options(tgContext, mfContext);
  options.errorEstimation = errorEstimation;
  options.stepTolerance = 1e-7;
  options.initialStepSize = 10_m;
  const FreeVector start =
      makeStart(Vector3(100., -50., 20.), Vector3(1., 0.3, 0.2), 1. / 1_GeV);

  const auto withoutJacobian =
      propagateFree(stepper, options, start, 1_m, true);
  options.jacobianInErrorEstimate = true;
  const auto withJacobian = propagateFree(stepper, options, start, 1_m, true);

  BOOST_TEST_CONTEXT("steps " << withoutJacobian.nSteps << " "
                              << withJacobian.nSteps) {
    BOOST_CHECK_GT(withJacobian.nSteps, withoutJacobian.nSteps);
  }
  CHECK_CLOSE_ABS(withJacobian.pars, withoutJacobian.pars, 1e-6);
}

/// The time of a neutral particle does not depend on q/p, because its charge
/// and so its momentum hypothesis are fixed.
BOOST_AUTO_TEST_CASE(generic_runge_kutta_stepper_neutral_time) {
  auto field = std::make_shared<ConstantBField>(Vector3::Zero());
  const GenericRungeKuttaStepper stepper(field);
  const auto options = fixedStepOptions(10_cm);
  const FreeVector start =
      makeStart(Vector3(1., 2., 3.), Vector3(1., 0.3, 0.2), 1. / 0.5_GeV);

  const auto neutral = propagateFree(stepper, options, start, 1_m, true,
                                     ParticleHypothesis::pion0());
  BOOST_CHECK_EQUAL(neutral.jacTransport(eFreeTime, eFreeQOverP), 0.);

  const auto charged = propagateFree(stepper, options, start, 1_m, true,
                                     ParticleHypothesis::pion());
  BOOST_CHECK_NE(charged.jacTransport(eFreeTime, eFreeQOverP), 0.);
}

/// A failed field lookup returns the error and leaves the state unchanged,
/// both for a stage and for the end of the step.
///
/// The explicit Euler method has a single stage at the start, so only the
/// lookup at the end of the step can fail.
BOOST_DATA_TEST_CASE(generic_runge_kutta_stepper_field_failure,
                     boost::unit_test::data::make(
                         std::vector<std::shared_ptr<const ButcherTableau>>{
                             ButcherTableau::classicalRk4(),
                             std::make_shared<const ButcherTableau>(
                                 "Euler", 1, 0, std::vector<double>{0.},
                                 std::vector<std::vector<double>>{{}},
                                 std::vector<double>{1.},
                                 std::vector<double>{})}),
                     tableau) {
  auto field = std::make_shared<BoundedField>(Vector3(0., 0., 2_T), 50_cm);
  const GenericRungeKuttaStepper stepper(
      GenericRungeKuttaStepper::Config{field, tableau});
  const auto options = fixedStepOptions(1_m);
  const FreeVector start =
      makeStart(Vector3::Zero(), Vector3(1., 0., 0.), 1. / 1_GeV);

  auto state = stepper.makeState(options);
  stepper.initialize(
      state, BoundTrackParameters::createCurvilinear(
                 makeVector4(start.segment<3>(eFreePos0), start[eFreeTime]),
                 start.segment<3>(eFreeDir0), start[eFreeQOverP],
                 BoundMatrix::Identity(), ParticleHypothesis::pion()));
  const FreeVector parsBefore = state.pars;
  const FreeMatrix jacBefore = state.jacTransport;

  BOOST_TEST_CONTEXT(tableau->name()) {
    auto res = stepper.step(state, Direction::Forward(), nullptr);
    BOOST_REQUIRE(!res.ok());
    BOOST_CHECK(res.error() == MagneticFieldError::OutOfBounds);
    BOOST_CHECK_EQUAL(state.pars, parsBefore);
    BOOST_CHECK_EQUAL(state.jacTransport, jacBefore);
    BOOST_CHECK_EQUAL(state.pathAccumulated, 0.);
    BOOST_CHECK_EQUAL(state.nSteps, 0u);
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
