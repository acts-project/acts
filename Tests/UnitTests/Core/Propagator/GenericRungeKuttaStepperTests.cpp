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
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Propagator/GenericRungeKuttaStepper.hpp"
#include "Acts/Propagator/HelixStepper.hpp"
#include "Acts/Propagator/detail/ButcherTableau.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <cmath>
#include <functional>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <vector>

using namespace Acts;
using namespace Acts::UnitLiterals;
using Acts::detail::ButcherTableau;
using Acts::VectorHelpers::makeVector4;

namespace ActsTests {

namespace {

const auto tgContext = GeometryContext::dangerouslyDefaultConstruct();
const MagneticFieldContext mfContext;

constexpr double infinity = std::numeric_limits<double>::infinity();

/// A smooth, non-uniform field with an analytic gradient
///
///   B = (s a y z, -s a x z, B0 + s a (x^2 + y^2)), s = B0 / L^2, a = 1/2
class SmoothField final : public MagneticFieldProvider {
 public:
  struct Cache {
    explicit Cache(const MagneticFieldContext& /*mctx*/) {}
  };

  SmoothField(double b0, double length, bool withGradient)
      : m_b0(b0), m_length(length), m_withGradient(withGradient) {}

  MagneticFieldProvider::Cache makeCache(
      const MagneticFieldContext& mctx) const override {
    return MagneticFieldProvider::Cache(std::in_place_type<Cache>, mctx);
  }

  Result<Vector3> getField(
      const Vector3& position,
      MagneticFieldProvider::Cache& /*cache*/) const override {
    return Result<Vector3>::success(evaluate(position).field);
  }

  bool providesFieldGradient() const override { return m_withGradient; }

  Result<FieldAndGradient> getFieldAndGradient(
      const Vector3& position,
      MagneticFieldProvider::Cache& /*cache*/) const override {
    if (!m_withGradient) {
      return Result<FieldAndGradient>::failure(
          MagneticFieldError::NotImplemented);
    }
    return Result<FieldAndGradient>::success(evaluate(position));
  }

 private:
  FieldAndGradient evaluate(const Vector3& p) const {
    constexpr double a = 0.5;
    const double s = m_b0 / (m_length * m_length);
    const double x = p.x();
    const double y = p.y();
    const double z = p.z();
    FieldAndGradient result;
    result.field =
        Vector3(s * a * y * z, -s * a * x * z, m_b0 + s * a * (x * x + y * y));
    // clang-format off
    result.gradient <<          0.,  s * a * z,  s * a * y,
                        -s * a * z,         0., -s * a * x,
                     2 * s * a * x, 2 * s * a * y,        0.;
    // clang-format on
    return result;
  }

  double m_b0;
  double m_length;
  bool m_withGradient;
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
typename stepper_t::State propagateFree(const stepper_t& stepper,
                                        typename stepper_t::Options options,
                                        const FreeVector& start,
                                        double pathLength, bool covTransport) {
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
          ParticleHypothesis::pion()));
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
  options.adaptiveStepSize = false;
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

/// The fixed-step global error falls with the order of the tableau.
BOOST_DATA_TEST_CASE(generic_runge_kutta_stepper_convergence,
                     boost::unit_test::data::make({0, 1, 2}), tableauIndex) {
  struct Case {
    std::shared_ptr<const ButcherTableau> tableau;
    double qop{};
    std::vector<int> stepCounts;
  };
  // The ninth-order error needs a strongly bent track to stay above the
  // round-off while the step size halves
  const Case testCase =
      tableauIndex == 0
          ? Case{ButcherTableau::classicalRk4(), 1. / 1_GeV, {8, 16, 32}}
      : tableauIndex == 1
          ? Case{ButcherTableau::dormandPrince54(), 1. / 1_GeV, {8, 16, 32}}
          : Case{ButcherTableau::verner98(), 1. / 0.1_GeV, {8, 12, 18}};

  auto field = std::make_shared<SmoothField>(2_T, 1_m, true);
  const FreeVector start =
      makeStart(Vector3(100., -50., 20.), Vector3(1., 0.3, 0.2), testCase.qop);
  const double pathLength = 1_m;

  // Every consistent tableau converges to the same solution, so a wrong order
  // still shows in the ratios.
  GenericRungeKuttaStepper::Config config{field, ButcherTableau::verner98()};
  const GenericRungeKuttaStepper reference(config);
  const FreeVector exact =
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
    errors.push_back((end - exact).segment<3>(eFreePos0).cwiseAbs().maxCoeff());
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

  auto analyticField = std::make_shared<SmoothField>(2_T, 1_m, true);
  const GenericRungeKuttaStepper stepper(analyticField);

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

  // A field without the analytic gradient falls back to finite differences
  auto numericField = std::make_shared<SmoothField>(2_T, 1_m, false);
  const GenericRungeKuttaStepper numericStepper(numericField);
  const auto fallback =
      propagateFree(numericStepper, options, start, pathLength, true);
  CHECK_CLOSE_OR_SMALL(fallback.jacTransport, analytic.jacTransport, 1e-6,
                       1e-8);

  // Without the gradient term the jacobian is wrong
  auto noGradientOptions = options;
  noGradientOptions.includeFieldGradient = false;
  const auto noGradient =
      propagateFree(stepper, noGradientOptions, start, pathLength, true);
  BOOST_CHECK_GT((noGradient.jacTransport - numeric).cwiseAbs().maxCoeff(),
                 1e-3);
}

/// The adaptive step size meets the tolerance.
BOOST_DATA_TEST_CASE(generic_runge_kutta_stepper_adaptive,
                     boost::unit_test::data::make({0, 1}), tableauIndex) {
  const auto tableau = tableauIndex == 0 ? ButcherTableau::dormandPrince54()
                                         : ButcherTableau::verner98();
  auto field = std::make_shared<SmoothField>(2_T, 1_m, true);
  const FreeVector start =
      makeStart(Vector3(100., -50., 20.), Vector3(1., 0.3, 0.2), 1. / 1_GeV);
  const double pathLength = 1_m;

  GenericRungeKuttaStepper::Config config{field, ButcherTableau::verner98()};
  const FreeVector exact = propagateFree(GenericRungeKuttaStepper(config),
                                         fixedStepOptions(pathLength / 1024),
                                         start, pathLength, false)
                               .pars;

  config.tableau = tableau;
  const GenericRungeKuttaStepper stepper(config);
  GenericRungeKuttaStepper::Options options(tgContext, mfContext);
  options.initialStepSize = 10_m;
  double previousError = infinity;
  for (double tolerance : {1e-4, 1e-7, 1e-10}) {
    options.stepTolerance = tolerance;
    const auto state =
        propagateFree(stepper, options, start, pathLength, false);
    const double error =
        (state.pars - exact).segment<3>(eFreePos0).cwiseAbs().maxCoeff();
    BOOST_TEST_CONTEXT(tableau->name()
                       << " tolerance " << tolerance << " error " << error
                       << " steps " << state.nSteps) {
      BOOST_CHECK_GT(state.statistics.nRejectedSteps, 0u);
      // The global error is the sum of the local errors
      BOOST_CHECK_LT(error, 10. * tolerance * state.nSteps);
      BOOST_CHECK_LT(error, previousError);
    }
    previousError = error;
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
