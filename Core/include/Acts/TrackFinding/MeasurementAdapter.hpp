// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/TrackFinding/MeasurementSelector.hpp"

#include <cstdint>

namespace Acts {

class Surface;

// TODO not too happy with the name `MeasurementAdapter` as it is still quite
// generic even tho this concept is only useful during combinatorial track
// finding. Maybe `CombinatorialMeasurementAdapter` would be better?

template <typename MeasurementAdapter>
concept MeasurementAdapterConcept = requires(const MeasurementAdapter& a,
                                             const Surface& s) {
  typename MeasurementAdapter::State;

  typename MeasurementAdapter::TrackStateContainer;
  typename MeasurementAdapter::TrackStateProxy;

  typename MeasurementAdapter::PreBoundSelectionRange;
  typename MeasurementAdapter::BoundSelectionRange;

  { a.makeState() } -> std::same_as<typename MeasurementAdapter::State>;

  { a.clearState(s) } -> std::same_as<void>;

  { a.isActive(s) } -> std::same_as<bool>;

  {
    a.selectPreBound(s)
  } -> std::same_as<typename MeasurementAdapter::PreBoundSelectionRange>;

  requires requires(
      const MeasurementAdapter::PreBoundSelectionRange& preBoundSelectionRange,
      const MeasurementAdapter::TrackStateContainer& trackStateContainer,
      const Surface& surface, const BoundVector& parameters,
      const BoundMatrix& covariance) {
    {
      a.selectPreCalibrationImpl(preBoundSelectionRange, surface, parameters,
                                 covariance, trackStateContainer)
    } -> std::same_as<typename MeasurementAdapter::BoundSelectionRange>;
  };
};

template <typename Derived, typename traj_t>
class MeasurementAdapterBase {
 public:
  using TrackStateContainer = traj_t;
  using TrackStateProxy = typename traj_t::TrackStateProxy;

  using BoundSelectionRange = std::pair<std::size_t, std::size_t>;

  struct State {};

  Derived::State makeState() const { return {}; }

  void clearState(Derived::State& state, const Surface& surface) const {
    (void)state;
    (void)surface;
  }

  bool isActive(Derived::State& state, const Surface& surface) const {
    (void)state;
    (void)surface;

    return true;
  }

  bool shouldCreateHoleImpl(
      Derived::State& state,
      const Derived::BoundSelectionRange& postCalibrationSelectionRange,
      const TrackStateContainer& trackStateContainer, const Surface& surface,
      const BoundVector& parameters, const BoundMatrix& covariance) const {
    (void)state;
    (void)postCalibrationSelectionRange;
    (void)trackStateContainer;
    (void)surface;
    (void)parameters;
    (void)covariance;

    return postCalibrationSelectionRange.first ==
           postCalibrationSelectionRange.second;
  }

  Derived::BoundSelectionRange selectPostBound(
      Derived::State& state,
      const Derived::PreBoundSelectionRange& preBoundSelectionRange,
      const Derived::TrackStateContainer& trackStateContainer,
      const Surface& surface, const BoundVector& parameters,
      const BoundMatrix& covariance) const {
    typename Derived::BoundSelectionRange range;

    // pre calibration selection and calibration
    {
      range.first = trackStateContainer.size();

      for (auto measurementIt = preBoundSelectionRange.first;
           measurementIt != preBoundSelectionRange.second; ++measurementIt) {
        const auto& measurement = *measurementIt;

        if (!self().doSelectPreCalibrationImpl(state, measurement, surface,
                                               parameters, covariance)) {
          continue;
        }

        TrackStateProxy trackState = trackStateContainer.makeTrackState(
            TrackStatePropMask::None, MultiTrajectoryTraits::kInvalid);

        self().calibrateImpl(state, trackState, measurement);

        if (!self().doSelectPostCalibrationImpl(state, trackState, surface,
                                                parameters, covariance)) {
          // TODO pop the last track state
          continue;
        }
      }

      range.second = trackStateContainer.size();
    }

    if (self().shouldCreateHoleImpl(state, range, trackStateContainer, surface,
                                    parameters, covariance)) {
      // TODO
    }

    return range;
  }

 private:
  const Derived& self() const { return static_cast<const Derived&>(*this); }
};

}  // namespace Acts
