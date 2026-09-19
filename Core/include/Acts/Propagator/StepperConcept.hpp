// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/TrackParametersConcept.hpp"
#include "Acts/EventData/detail/CorrectedTransformationFreeToBound.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Propagator/StepperStatistics.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Concepts.hpp"
#include "Acts/Utilities/Logger.hpp"

namespace Acts {
namespace Concepts {

/// @brief Concept that is satisfied by both single- and multi-steppers.
template <typename Stepper, typename State = typename Stepper::State>
concept CommonStepper = requires {
  typename Stepper::BoundParameters;
  typename Stepper::State;
  typename Stepper::Jacobian;
  typename Stepper::Covariance;

  requires requires(const Stepper& s, State& t) {
    { s.transportToCurvilinear(t) } -> std::same_as<typename Stepper::Jacobian>;
    { s.prepareCurvilinearState(t) } -> std::same_as<bool>;

    requires requires(const typename Stepper::BoundParameters& par) {
      { s.initialize(t, par) } -> std::same_as<void>;
    };

    requires requires(const Surface& sf, const FreeToBoundCorrection& corr) {
      {
        s.transportToBound(t, sf, corr)
      } -> std::same_as<Result<typename Stepper::Jacobian>>;
    };

    requires requires(const Surface& sf, std::uint8_t ui, Direction d,
                      const BoundaryTolerance& bt, double sc,
                      ConstrainedStep::Type st, const Logger& l) {
      { s.updateSurfaceStatus(t, sf, ui, d, bt, sc, st, l) };
    };

    requires requires(const ConstrainedStep::Type st) {
      { s.releaseStepSize(t, st) } -> std::same_as<void>;

      requires requires(double d) {
        { s.updateStepSize(t, d, st) } -> std::same_as<void>;
      };
    };
  };

  requires requires(const Stepper& s, const State& t) {
    { s.position(t) } -> std::same_as<Vector3>;
    { s.direction(t) } -> std::same_as<Vector3>;
    { s.qOverP(t) } -> std::same_as<double>;
    { s.absoluteMomentum(t) } -> std::same_as<double>;
    { s.momentum(t) } -> std::same_as<Vector3>;
    { s.charge(t) } -> std::same_as<double>;
    { s.time(t) } -> std::same_as<double>;
    { s.outputStepSize(t) } -> std::same_as<std::string>;
    { s.pathLength(t) } -> std::same_as<double>;
    { s.hasCovariance(t) } -> std::same_as<bool>;
    { s.statistics(t) } -> std::same_as<const StepperStatistics&>;
    {
      s.curvilinearParameters(t)
    } -> std::same_as<typename Stepper::BoundParameters>;

    requires requires(const Surface& sf) {
      {
        s.boundParameters(t, sf)
      } -> std::same_as<Result<typename Stepper::BoundParameters>>;
    };

    requires requires(const ConstrainedStep::Type st) {
      { s.getStepSize(t, st) } -> std::same_as<double>;
    };
  };
};

/// @brief Concept that is satisfied by single-steppers.
template <typename Stepper, typename State = typename Stepper::State>
concept SingleStepper =
    CommonStepper<Stepper, State> && requires(const Stepper& s, State& t) {
      requires BoundTrackParametersConcept<typename Stepper::BoundParameters>;

      requires requires(const FreeVector& fv, const BoundVector& bv,
                        const BoundMatrix& bm, const Surface& sf) {
        { s.update(t, fv, bv, bm, sf) } -> std::same_as<void>;
      };

      requires requires(const Vector3& v1, const Vector3& v2, double d1,
                        double d2) {
        { s.update(t, v1, v2, d1, d2) } -> std::same_as<void>;
        { s.getField(t, v1) } -> std::same_as<Result<Vector3>>;
      };

      requires requires(const typename Stepper::Covariance& cov) {
        { s.setCovariance(t, cov) } -> std::same_as<void>;
      };

      requires requires(const Stepper& cs, const State& ct) {
        {
          cs.covariance(ct)
        } -> std::same_as<const typename Stepper::Covariance&>;
        { cs.stepSize(ct) } -> std::same_as<const ConstrainedStep&>;
      };
    };

/// @brief Concept that is satisfied by multi-steppers.
template <typename Stepper, typename State = typename Stepper::State>
concept MultiStepper = CommonStepper<Stepper, State> && requires {
  // TODO for now we do not check if the ComponentProxy does fulfill a concept
  typename Stepper::ComponentProxy;

  // TODO for now we do not check if the ConstComponentProxy does fulfill a
  // concept
  typename Stepper::ConstComponentProxy;

  requires requires(const Stepper& s, State& t) {
    { s.numberComponents(t) } -> std::same_as<std::size_t>;
    { s.clearComponents(t) } -> std::same_as<void>;
    { s.removeMissedComponents(t) } -> std::same_as<void>;
  };
};
}  // namespace Concepts

/// @brief Concept that is satisfied by steppers.
template <typename _Stepper, typename State = typename _Stepper::State>
concept StepperConcept = Concepts::SingleStepper<_Stepper, State> ||
                         Concepts::MultiStepper<_Stepper, State>;

}  // namespace Acts
