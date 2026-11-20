// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/EventData/Types.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/TrackFitting/GsfComponent.hpp"
#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

#include <functional>
#include <vector>

namespace Acts {

struct ScatteringTailPrior;

/// @addtogroup track_finding
/// @{
namespace CkfTypes {

/// expected max number of track states that are expected to be added by
/// stateCandidateCreator
/// @note if the number of states exceeds this number dynamic memory allocation will occur.
///       the number is chosen to yield a container size of 64 bytes.

static constexpr std::size_t s_maxBranchesPerSurface = 10;
template <typename T>
using BranchVector = boost::container::small_vector<T, s_maxBranchesPerSurface>;

/// The predicted state of the trajectory on a surface
struct BoundState {
  /// The predicted bound parameters, including their reference surface
  BoundTrackParameters parameters;
  /// The jacobian from the previous track state to this surface
  BoundMatrix jacobian;
  /// The path length accumulated up to this surface
  double pathLength{};
};

/// The arguments of the track state creation on one surface.
///
/// Only the CKF makes this object, so it can get more fields without a change
/// of the track state creator signature.
///
/// @note It borrows everything it refers to for the duration of one call of
///       the track state creator. Do not store it.
///
/// @tparam traj_t the track state container backend
template <typename traj_t>
// NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
struct TrackStateCreatorArguments {
  /// The current geometry context
  std::reference_wrapper<const GeometryContext> geoContext;
  /// The current calibration context
  std::reference_wrapper<const CalibrationContext> calibrationContext;
  /// The predicted state on the surface the trajectory is extended onto
  std::reference_wrapper<const BoundState> boundState;
  /// The tip of the trajectory which is extended
  TrackIndexType prevTip{};
  /// The trajectory the new track states are added to
  std::reference_wrapper<traj_t> trajectory;
  /// The scattering tail of the prediction, nullptr if there is none
  const ScatteringTailPrior* scatteringTail{nullptr};
};

}  // namespace CkfTypes

/// Return type of the `BranchStopper` delegate for the
/// CombinatorialKalmanFilter
enum class CombinatorialKalmanFilterBranchStopperResult {
  Continue,
  StopAndDrop,
  StopAndKeep,
};

/// Extension struct which holds the delegates to customize the CKF behavior
template <typename track_container_t>
struct CombinatorialKalmanFilterExtensions {
  /// Type alias for track state container backend
  using traj_t = typename track_container_t::TrackStateContainerBackend;
  /// Type alias for track proxy from the container
  using TrackProxy = typename track_container_t::TrackProxy;
  /// Type alias for track state proxy from the container
  using TrackStateProxy = typename track_container_t::TrackStateProxy;

  /// Type alias for branch stopper result enumeration
  using BranchStopperResult = CombinatorialKalmanFilterBranchStopperResult;

  /// Type alias for Kalman filter updater delegate
  using Updater = typename KalmanFitterExtensions<traj_t>::Updater;
  /// Type alias for branch stopper delegate function
  using BranchStopper =
      Delegate<BranchStopperResult(const TrackProxy&, const TrackStateProxy&)>;

  /// The updater incorporates measurement information into the track parameters
  Updater updater{DelegateFuncTag<detail::voidFitterUpdater<traj_t>>{}};

  /// The branch stopper is called during the filtering by the Actor.
  BranchStopper branchStopper{DelegateFuncTag<voidBranchStopper>{}};

  /// @brief Delegate the extension of the trajectory onto the surface of
  ///        `args.boundState` to an external unit.
  ///
  /// @note Expected to create track states for measurements on that surface
  ///       which match the predicted state. The "filtered" data is not
  ///       expected to be set, but the outlier flag should be set for states
  ///       that are considered to be outlier. A track state may have its own
  ///       predicted parameters instead of sharing them.
  ///
  /// @param args the predicted state, the trajectory to extend, and the
  ///        contexts
  /// @param logger a logger for messages
  /// @return indices of new track states which extend the trajectory given by
  ///         `args.prevTip`. Each one must be a measurement or an outlier,
  ///         otherwise the CKF fails with
  ///         @ref CombinatorialKalmanFilterError::UnsupportedTrackState. An
  ///         empty list means that the CKF creates a hole, the error
  ///         @ref CombinatorialKalmanFilterError::NoMeasurementExpected that it
  ///         creates a state without an expected measurement.
  using TrackStateCreator =
      Delegate<Result<CkfTypes::BranchVector<TrackIndexType>>(
          const CkfTypes::TrackStateCreatorArguments<traj_t>& args,
          const Logger& logger)>;

  /// The delegate to create new track states.
  /// @note a reference implementation can be found in @ref TrackStateCreator
  ///   which makes uses of @ref MeasurementSelector and SourceLinkAccessor
  TrackStateCreator createTrackStates;

  // The following options are only relevant if a multi stepper is used

  /// Type alias for component reducer delegate function
  using ComponentReducer =
      Delegate<void(std::vector<GsfComponent>&, std::size_t, const Surface&)>;

  /// Takes a vector of components and reduces its number
  ComponentReducer mixtureReducer;

 private:
  /// Default branch stopper which will never stop
  /// @return false
  static BranchStopperResult voidBranchStopper(
      const TrackProxy& /*track*/, const TrackStateProxy& /*trackState*/) {
    return BranchStopperResult::Continue;
  }
};

/// @}

}  // namespace Acts
