// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Direction.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/AnyTrackStateProxy.hpp"
#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/EventData/TrackContainerFrontendConcept.hpp"
#include "Acts/EventData/TrackProxyCommon.hpp"
#include "Acts/EventData/TrackProxyConcept.hpp"
#include "Acts/EventData/TrackStateProxyConcept.hpp"
#include "Acts/EventData/TrackStateType.hpp"
#include "Acts/EventData/TransformationHelpers.hpp"
#include "Acts/EventData/Types.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Propagator/StandardAborters.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/KalmanFitterError.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

#include <optional>
#include <utility>
#include <vector>

namespace Acts {

/// Strategy for track extrapolation to target surface
enum class TrackExtrapolationStrategy {
  /// Use the first track state to reach target surface
  first,
  /// Use the last track state to reach target surface
  last,
  /// Use the first or last track state to reach target surface depending on the
  /// distance
  firstOrLast,
};

/// Error codes for track extrapolation operations
/// @ingroup errors
enum class TrackExtrapolationError {
  /// Did not find a compatible track state
  CompatibleTrackStateNotFound = 1,
  /// Provided reference surface is unreachable
  ReferenceSurfaceUnreachable = 2,
};

/// Create error code from TrackExtrapolationError
/// @param e The error code enum value
/// @return Standard error code
std::error_code make_error_code(TrackExtrapolationError e);

/// Find the first measurement state in a track
/// @param track The track to search
/// @return Result containing the first measurement state proxy or error
template <TrackProxyConcept track_proxy_t>
Result<typename track_proxy_t::ConstTrackStateProxy> findFirstMeasurementState(
    const track_proxy_t &track) {
  using TrackStateProxy = typename track_proxy_t::ConstTrackStateProxy;

  // TODO specialize if track is forward linked

  std::optional<TrackStateProxy> firstMeasurementOpt;

  for (const auto &trackState : track.trackStatesReversed()) {
    if (trackState.typeFlags().isMeasurement()) {
      firstMeasurementOpt = trackState;
    }
  }

  if (firstMeasurementOpt) {
    return Result<TrackStateProxy>::success(*firstMeasurementOpt);
  }

  return Result<TrackStateProxy>::failure(
      TrackExtrapolationError::CompatibleTrackStateNotFound);
}

/// Find the last measurement state in a track
/// @param track The track to search
/// @return Result containing the last measurement state proxy or error
template <TrackProxyConcept track_proxy_t>
Result<typename track_proxy_t::ConstTrackStateProxy> findLastMeasurementState(
    const track_proxy_t &track) {
  using TrackStateProxy = typename track_proxy_t::ConstTrackStateProxy;

  for (const auto &trackState : track.trackStatesReversed()) {
    if (trackState.typeFlags().isMeasurement()) {
      return TrackStateProxy{trackState};
    }
  }

  return Result<TrackStateProxy>::failure(
      TrackExtrapolationError::CompatibleTrackStateNotFound);
}

/// @brief Smooth a track using the gain matrix smoother
///
/// @tparam track_proxy_t The track proxy type
/// @tparam smoother_t The smoother type
///
/// @param geoContext The geometry context
/// @param track The track to smooth
/// @param logger The logger
/// @param smoother The smoother
///
/// @return The result of the smoothing
template <TrackProxyConcept track_proxy_t,
          typename smoother_t = GainMatrixSmoother>
Result<void> smoothTrack(
    const GeometryContext &geoContext, track_proxy_t &track,
    const Logger &logger = *getDefaultLogger("TrackSmoother", Logging::INFO),
    smoother_t smoother = GainMatrixSmoother()) {
  auto &trackContainer = track.container();
  auto &trackStateContainer = trackContainer.trackStateContainer();

  auto last = findLastMeasurementState(track);
  if (!last.ok()) {
    ACTS_DEBUG("no last track state found");
    return last.error();
  }

  auto smoothingResult =
      smoother(geoContext, trackStateContainer, last->index(), logger);

  if (!smoothingResult.ok()) {
    ACTS_DEBUG("Smoothing track " << track.index() << " failed with error "
                                  << smoothingResult.error());
    return smoothingResult.error();
  }

  return Result<void>::success();
}

/// @brief Smooth tracks using the gain matrix smoother
///
/// @tparam track_container_t The track container type
///
/// @param geoContext The geometry context
/// @param trackContainer The track container
/// @param logger The logger
///
/// @return The result of the smoothing
template <TrackContainerFrontend track_container_t>
Result<void> smoothTracks(
    const GeometryContext &geoContext, const track_container_t &trackContainer,
    const Logger &logger = *getDefaultLogger("TrackSmoother", Logging::INFO)) {
  Result<void> result = Result<void>::success();

  for (const auto &track : trackContainer) {
    auto smoothingResult = smoothTrack(geoContext, track, logger);

    // Only keep the first error
    if (!smoothingResult.ok() && result.ok()) {
      result = smoothingResult.error();
    }
  }

  return result;
}

/// @brief Condition the smoothed states of a segment on a new estimate of
/// one of its states
///
/// The segment runs from @p sharedIndex to @p tipIndex in filter order. Its
/// smoothed states hold the segment's own measurements. A second segment that
/// was filtered from the shared state as prior adds information to the shared
/// state; this function transports that change to every state of the first
/// segment without a new propagation. Conditional on the shared state the two
/// segments are independent, so the exact update in the linear-Gaussian model
/// is the regression of each state on the shared state:
/// `L_k = Cov(x_k, x_s) P_s^-1`, `x_k += L_k (x_s' - x_s)`,
/// `P_k += L_k (P_s' - P_s) L_k^T`.
///
/// The regression is computed in the form of the modified Bryson-Frazier
/// smoother, so no state covariance is inverted. The smoothed states are
/// affine in the prior mean of the segment. A forward sweep carries the
/// sensitivity of the filtered states to the prior mean,
/// `Dp_k = F_k Df_(k-1)`, `Df_k = (I - K_k H_k) Dp_k`, and a backward sweep
/// carries the sensitivity of the smoother adjoint,
/// `Psi~_k = H_k^T S_k^-1 H_k Dp_k + (I - K_k H_k)^T Psi^_k`,
/// `Psi^_(k-1) = F_k^T Psi~_k`, so that `T_k = Df_k - Pf_k Psi^_k` is the
/// sensitivity of the smoothed state. With `Z = Pp_s Ps_s^-1`, the move of
/// the posterior of the shared state per unit of its prior mean, the
/// regression is `L_k = T_k Z`. Only `S_k` and the shared covariance are
/// inverted.
///
/// States before the shared state in filter order get a plain smoothing step
/// from the corrected shared state. States after the last smoothed state, such
/// as holes after the last measurement, stay as they are.
///
/// @tparam traj_t The MultiTrajectory type
///
/// @param trajectory The trajectory that holds the segment
/// @param sharedIndex The index of the shared state
/// @param tipIndex The index of the last state of the segment in filter order
/// @param newMean The new estimate of the shared state
/// @param newCov Its covariance
/// @param logger The logger
///
/// @return Success, or @ref KalmanFitterError::SmoothFailed if the shared
///         state is not on the chain, a state lacks the needed components, or
///         a result is not finite
template <typename traj_t>
Result<void> conditionSmoothedStatesOnSharedState(
    traj_t &trajectory, std::size_t sharedIndex, std::size_t tipIndex,
    const BoundVector &newMean, const BoundMatrix &newCov,
    const Logger &logger = getDummyLogger()) {
  using TrackStateProxy = typename traj_t::TrackStateProxy;

  // collect the segment from the tip down to the shared state
  std::vector<std::size_t> indices;
  bool found = false;
  trajectory.applyBackwards(tipIndex, [&](TrackStateProxy ts) {
    indices.push_back(ts.index());
    found = (ts.index() == sharedIndex);
    return !found;
  });
  if (!found) {
    ACTS_DEBUG("Shared state " << sharedIndex << " is not on the chain of "
                               << tipIndex);
    return KalmanFitterError::SmoothFailed;
  }

  // states beyond the last smoothed state, such as holes after the last
  // measurement, carry no information and stay as they are
  std::size_t skip = 0;
  while (skip + 1 < indices.size() &&
         !trajectory.getTrackState(indices[skip]).hasSmoothed()) {
    ++skip;
  }
  indices.erase(indices.begin(), indices.begin() + skip);

  auto shared = trajectory.getTrackState(sharedIndex);
  if (!shared.hasSmoothed()) {
    ACTS_DEBUG("Shared state " << sharedIndex << " has no smoothed parameters");
    return KalmanFitterError::SmoothFailed;
  }

  const BoundVector sharedMean = shared.smoothed();
  const BoundMatrix sharedCov = shared.smoothedCovariance();
  const BoundVector deltaMean = subtractBoundParameters(newMean, sharedMean);
  const BoundMatrix deltaCov = newCov - sharedCov;
  const Eigen::LDLT<BoundMatrix> sharedLdlt(sharedCov);
  if (sharedLdlt.info() != Eigen::Success || !sharedLdlt.isPositive()) {
    ACTS_DEBUG("Shared state covariance is not positive definite");
    return KalmanFitterError::SmoothFailed;
  }
  // Z = Pp_s Ps_s^-1 turns the sensitivity to the prior mean into the
  // regression on the smoothed shared state; it acts on the change only
  const BoundMatrix sharedPredictedCov = shared.predictedCovariance();
  const BoundMatrix z = sharedLdlt.solve(sharedPredictedCov).transpose();
  const BoundVector u = z * deltaMean;
  const BoundMatrix w = z * deltaCov * z.transpose();

  // measurement pieces of a state: the filtered sensitivity from the
  // predicted one, and the adjoint sensitivity from the predicted one
  auto filteredSensitivity = [](const AnyConstTrackStateProxy &ts,
                                const BoundMatrix &dp) -> BoundMatrix {
    if (!ts.typeFlags().isMeasurement()) {
      return dp;
    }
    BoundMatrix df = dp;
    visit_measurement(
        ts.calibratedSize(),
        [&]<std::size_t N>(std::integral_constant<std::size_t, N>) {
          const auto H = ts.projectorSubspaceHelper<N>().projector();
          const BoundMatrix pp = ts.predictedCovariance();
          const Eigen::Matrix<double, N, N> s =
              H * pp * H.transpose() + ts.calibratedCovariance<N>();
          const Eigen::Matrix<double, eBoundSize, N> k =
              pp * H.transpose() * s.inverse();
          df = dp - k * (H * dp);
        });
    return df;
  };
  auto adjointSensitivity = [](const AnyConstTrackStateProxy &ts,
                               const BoundMatrix &dp,
                               const BoundMatrix &psiHat) -> BoundMatrix {
    if (!ts.typeFlags().isMeasurement()) {
      return psiHat;
    }
    BoundMatrix psiTilde = psiHat;
    visit_measurement(
        ts.calibratedSize(),
        [&]<std::size_t N>(std::integral_constant<std::size_t, N>) {
          const auto H = ts.projectorSubspaceHelper<N>().projector();
          const BoundMatrix pp = ts.predictedCovariance();
          const Eigen::Matrix<double, N, N> s =
              H * pp * H.transpose() + ts.calibratedCovariance<N>();
          const Eigen::Matrix<double, N, N> sInv = s.inverse();
          const Eigen::Matrix<double, eBoundSize, N> k =
              pp * H.transpose() * sInv;
          // H^T S^-1 H Dp + (I - K H)^T Psi^
          psiTilde = H.transpose() * (sInv * (H * dp)) + psiHat -
                     H.transpose() * (k.transpose() * psiHat);
        });
    return psiTilde;
  };

  // forward sweep from the shared state to the tip: the sensitivity of the
  // predicted and filtered states to the prior mean of the segment
  std::vector<std::pair<BoundMatrix, BoundMatrix>> sensitivities;
  sensitivities.resize(indices.size());
  BoundMatrix df = BoundMatrix::Identity();
  for (std::size_t k = indices.size(); k-- > 0;) {
    auto ts = trajectory.getTrackState(indices[k]);
    if (!ts.hasSmoothed() || !ts.hasPredicted() || !ts.hasJacobian()) {
      ACTS_DEBUG("State " << ts.index()
                          << " lacks smoothed, predicted or jacobian");
      return KalmanFitterError::SmoothFailed;
    }
    const bool isShared = (k + 1 == indices.size());
    const BoundMatrix dp =
        isShared ? BoundMatrix::Identity() : BoundMatrix(ts.jacobian() * df);
    df = filteredSensitivity(AnyConstTrackStateProxy{ts}, dp);
    sensitivities[k] = {dp, df};
  }

  // backward sweep from the tip to the shared state: the sensitivity of the
  // adjoint, the sensitivity of the smoothed state, and the update
  BoundMatrix psiHat = BoundMatrix::Zero();
  for (std::size_t k = 0; k < indices.size(); ++k) {
    auto ts = trajectory.getTrackState(indices[k]);
    const bool isShared = (k + 1 == indices.size());
    const auto &[dp, dfk] = sensitivities[k];

    if (!isShared) {
      // non-measurement states carry the prediction as filtered estimate
      const BoundMatrix filteredCov =
          ts.hasFiltered() ? BoundMatrix(ts.filteredCovariance())
                           : BoundMatrix(ts.predictedCovariance());
      const BoundMatrix t = dfk - filteredCov * psiHat;

      const BoundVector mean = normalizeBoundParameters(ts.smoothed() + t * u);
      const BoundMatrix cov = ts.smoothedCovariance() + t * w * t.transpose();
      if (!mean.allFinite() || !cov.allFinite()) {
        ACTS_DEBUG("Conditioned state " << ts.index() << " is not finite");
        return KalmanFitterError::SmoothFailed;
      }
      ts.smoothed() = mean;
      // keep the covariance symmetric against rounding
      ts.smoothedCovariance() = 0.5 * (cov + cov.transpose());

      const BoundMatrix psiTilde =
          adjointSensitivity(AnyConstTrackStateProxy{ts}, dp, psiHat);
      psiHat = ts.jacobian().transpose() * psiTilde;
    }
  }

  shared.smoothed() = newMean;
  shared.smoothedCovariance() = newCov;

  // states before the shared state: plain smoothing from the corrected state
  if (shared.hasPrevious()) {
    GainMatrixSmoother smoother;
    std::error_code error;
    auto later = shared;
    trajectory.applyBackwards(shared.previous(), [&](TrackStateProxy ts) {
      ts.addComponents(TrackStatePropMask::Smoothed);
      auto res = smoother.calculate(AnyMutableTrackStateProxy{ts},
                                    AnyConstTrackStateProxy{later}, logger);
      if (!res.ok()) {
        error = res.error();
        return false;
      }
      later = ts;
      return true;
    });
    if (error) {
      return Result<void>::failure(error);
    }
  }

  return Result<void>::success();
}

/// @brief Link a second tracklet behind a state of a first tracklet
///
/// The second tracklet was filtered from the shared state as prior. Its
/// states are reversed so that the chain reads from the tip of the first
/// tracklet through the shared state into the second tracklet. With
/// @p invertJacobians the jacobians of the second tracklet are inverted and
/// shifted, and the transport from the second tracklet onto the shared state
/// is stored on the shared state.
///
/// @tparam track_state_proxy_t The track state proxy type
/// @tparam track_proxy_t The track proxy type
///
/// @param sharedState The shared state of the first tracklet
/// @param second The second tracklet
/// @param invertJacobians Whether to invert and shift the jacobians
///
/// @return The transport from the second tracklet onto the shared state, if
///         @p invertJacobians and the second tracklet has one
template <TrackStateProxyConcept track_state_proxy_t,
          TrackProxyConcept track_proxy_t>
std::optional<BoundMatrix> linkTrackletBehind(track_state_proxy_t sharedState,
                                              track_proxy_t &second,
                                              bool invertJacobians = true) {
  std::optional<BoundMatrix> stitch =
      second.reverseTrackStates(invertJacobians);
  if (second.tipIndex() == kTrackIndexInvalid) {
    return stitch;
  }

  auto secondTip = second.outermostTrackState();
  sharedState.previous() = secondTip.index();
  secondTip.template component<TrackIndexType>(detail_tp::kNextKey) =
      sharedState.index();

  if (stitch.has_value()) {
    sharedState.jacobian() = *stitch;
  }

  return stitch;
}

/// @brief Find a track state for extrapolation
///
/// @tparam track_proxy_t The track proxy type
///
/// @param geoContext The geometry context
/// @param track The track
/// @param referenceSurface The reference surface
/// @param strategy The extrapolation strategy
/// @param logger The logger
///
/// @return The result of the search containing the track state
///         and the distance to the reference surface
template <TrackProxyConcept track_proxy_t>
Result<std::pair<typename track_proxy_t::ConstTrackStateProxy, double>>
findTrackStateForExtrapolation(
    const GeometryContext &geoContext, const track_proxy_t &track,
    const Surface &referenceSurface, TrackExtrapolationStrategy strategy,
    const Logger &logger = *getDefaultLogger("TrackExtrapolation",
                                             Logging::INFO)) {
  using TrackStateProxy = typename track_proxy_t::ConstTrackStateProxy;

  // Intersect the reference surface with the trajectory at a track state.
  // Returns `std::nullopt` if the state carries no parameters at all and can
  // therefore not be started from; an invalid intersection means the state has
  // parameters but does not reach the reference surface.
  auto intersect =
      [&](const TrackStateProxy &state) -> std::optional<Intersection3D> {
    if (!state.hasSmoothed() && !state.hasFiltered() && !state.hasPredicted()) {
      return std::nullopt;
    }

    // `parameters` picks smoothed over filtered over predicted, which is what
    // `TrackProxy::createParametersFromState` starts the propagation from. The
    // distance has to be measured on the same parameters.
    const FreeVector freeVector = transformBoundToFreeParameters(
        state.referenceSurface(), geoContext, state.parameters());

    return referenceSurface
        .intersect(geoContext, freeVector.template segment<3>(eFreePos0),
                   freeVector.template segment<3>(eFreeDir0),
                   BoundaryTolerance::None(), s_onSurfaceTolerance)
        .closest();
  };

  switch (strategy) {
    case TrackExtrapolationStrategy::first: {
      ACTS_VERBOSE("looking for first track state");

      auto first = findFirstMeasurementState(track);
      if (!first.ok()) {
        ACTS_DEBUG("no first track state found");
        return first.error();
      }

      std::optional<Intersection3D> intersection = intersect(*first);
      if (!intersection.has_value()) {
        ACTS_DEBUG("first track state carries no parameters");
        return Result<std::pair<TrackStateProxy, double>>::failure(
            TrackExtrapolationError::CompatibleTrackStateNotFound);
      }
      if (!intersection->isValid()) {
        ACTS_DEBUG("no intersection found");
        return Result<std::pair<TrackStateProxy, double>>::failure(
            TrackExtrapolationError::ReferenceSurfaceUnreachable);
      }

      ACTS_VERBOSE("found intersection at " << intersection->pathLength());
      return std::pair(*first, intersection->pathLength());
    }

    case TrackExtrapolationStrategy::last: {
      ACTS_VERBOSE("looking for last track state");

      auto last = findLastMeasurementState(track);
      if (!last.ok()) {
        ACTS_DEBUG("no last track state found");
        return last.error();
      }

      std::optional<Intersection3D> intersection = intersect(*last);
      if (!intersection.has_value()) {
        ACTS_DEBUG("last track state carries no parameters");
        return Result<std::pair<TrackStateProxy, double>>::failure(
            TrackExtrapolationError::CompatibleTrackStateNotFound);
      }
      if (!intersection->isValid()) {
        ACTS_DEBUG("no intersection found");
        return Result<std::pair<TrackStateProxy, double>>::failure(
            TrackExtrapolationError::ReferenceSurfaceUnreachable);
      }

      ACTS_VERBOSE("found intersection at " << intersection->pathLength());
      return std::pair(*last, intersection->pathLength());
    }

    case TrackExtrapolationStrategy::firstOrLast: {
      ACTS_VERBOSE("looking for first or last track state");

      auto first = findFirstMeasurementState(track);
      if (!first.ok()) {
        ACTS_DEBUG("no first track state found");
        return first.error();
      }

      auto last = findLastMeasurementState(track);
      if (!last.ok()) {
        ACTS_DEBUG("no last track state found");
        return last.error();
      }

      std::optional<Intersection3D> intersectionFirstOpt = intersect(*first);
      std::optional<Intersection3D> intersectionLastOpt = intersect(*last);

      if (!intersectionFirstOpt.has_value() &&
          !intersectionLastOpt.has_value()) {
        ACTS_DEBUG("neither first nor last track state carries parameters");
        return Result<std::pair<TrackStateProxy, double>>::failure(
            TrackExtrapolationError::CompatibleTrackStateNotFound);
      }

      // an end without parameters cannot be started from, so it loses the
      // comparison below through the infinite path length of an invalid
      // intersection
      Intersection3D intersectionFirst =
          intersectionFirstOpt.value_or(Intersection3D::Invalid());
      Intersection3D intersectionLast =
          intersectionLastOpt.value_or(Intersection3D::Invalid());

      double absDistanceFirst = std::abs(intersectionFirst.pathLength());
      double absDistanceLast = std::abs(intersectionLast.pathLength());

      if (intersectionFirst.isValid() && absDistanceFirst <= absDistanceLast) {
        ACTS_VERBOSE("using first track state with intersection at "
                     << intersectionFirst.pathLength());
        return std::pair(*first, intersectionFirst.pathLength());
      }

      if (intersectionLast.isValid() && absDistanceLast <= absDistanceFirst) {
        ACTS_VERBOSE("using last track state with intersection at "
                     << intersectionLast.pathLength());
        return std::pair(*last, intersectionLast.pathLength());
      }

      ACTS_DEBUG("no intersection found");
      return Result<std::pair<TrackStateProxy, double>>::failure(
          TrackExtrapolationError::ReferenceSurfaceUnreachable);
    }
  }

  // unreachable
  return Result<std::pair<TrackStateProxy, double>>::failure(
      TrackExtrapolationError::CompatibleTrackStateNotFound);
}

/// @brief Extrapolate a track to a reference surface
///
/// @tparam track_proxy_t The track proxy type
/// @tparam propagator_t The propagator type
/// @tparam propagator_options_t The propagator options type
///
/// @param track The track which is modified in-place
/// @param referenceSurface The reference surface
/// @param propagator The propagator
/// @param options The propagator options
/// @param strategy The extrapolation strategy
/// @param logger The logger
///
/// @return The result of the extrapolation
template <TrackProxyConcept track_proxy_t, typename propagator_t,
          typename propagator_options_t>
Result<void> extrapolateTrackToReferenceSurface(
    track_proxy_t &track, const Surface &referenceSurface,
    const propagator_t &propagator, propagator_options_t options,
    TrackExtrapolationStrategy strategy,
    const Logger &logger = *getDefaultLogger("TrackExtrapolation",
                                             Logging::INFO)) {
  auto findResult = findTrackStateForExtrapolation(
      options.geoContext, track, referenceSurface, strategy, logger);

  if (!findResult.ok()) {
    ACTS_DEBUG("failed to find track state for extrapolation");
    return findResult.error();
  }

  auto &[trackState, distance] = *findResult;

  options.direction = Direction::fromScalarZeroAsPositive(distance);

  BoundTrackParameters parameters = track.createParametersFromState(trackState);
  ACTS_VERBOSE("extrapolating track to reference surface at distance "
               << distance << " with direction " << options.direction
               << " with starting parameters " << parameters);

  auto propagateResult =
      propagator.template propagate<propagator_options_t, ForcedSurfaceReached>(
          parameters, referenceSurface, options);

  if (!propagateResult.ok()) {
    ACTS_DEBUG("failed to extrapolate track: " << propagateResult.error());
    return propagateResult.error();
  }

  track.setReferenceSurface(referenceSurface.getSharedPtr());
  track.parameters() = propagateResult->endParameters.value().parameters();
  track.covariance() =
      propagateResult->endParameters.value().covariance().value();

  return Result<void>::success();
}

/// @brief Extrapolate tracks to a reference surface
///
/// @tparam track_container_t The track container type
/// @tparam propagator_t The propagator type
/// @tparam propagator_options_t The propagator options type
///
/// @param trackContainer The track container which is modified in-place
/// @param referenceSurface The reference surface
/// @param propagator The propagator
/// @param options The propagator options
/// @param strategy The extrapolation strategy
/// @param logger The logger
///
/// @return The result of the extrapolation
template <TrackContainerFrontend track_container_t, typename propagator_t,
          typename propagator_options_t>
Result<void> extrapolateTracksToReferenceSurface(
    const track_container_t &trackContainer, const Surface &referenceSurface,
    const propagator_t &propagator, propagator_options_t options,
    TrackExtrapolationStrategy strategy,
    const Logger &logger = *getDefaultLogger("TrackExtrapolation",
                                             Logging::INFO)) {
  Result<void> result = Result<void>::success();

  for (const auto &track : trackContainer) {
    auto extrapolateResult = extrapolateTrackToReferenceSurface(
        track, referenceSurface, propagator, options, strategy, logger);

    // Only keep the first error
    if (!extrapolateResult.ok() && result.ok()) {
      result = extrapolateResult.error();
    }
  }

  return result;
}

/// Helper function to calculate a number of track level quantities and store
/// them on the track itself
/// @tparam track_proxy_t The track proxy type
/// @param track A mutable track proxy to operate on
template <TrackProxyConcept track_proxy_t>
void calculateTrackQuantities(track_proxy_t track)
  requires(!track_proxy_t::ReadOnly)
{
  track.chi2() = 0;
  track.nDoF() = 0;

  track.nHoles() = 0;
  track.nMeasurements() = 0;
  track.nSharedHits() = 0;
  track.nOutliers() = 0;

  // The split hit column is optional, so only count if the track has it
  unsigned int nSplitHits = 0;

  for (const auto &trackState : track.trackStatesReversed()) {
    ConstTrackStateTypeMap typeFlags = trackState.typeFlags();

    if (typeFlags.isHole()) {
      track.nHoles()++;
    } else if (typeFlags.isOutlier()) {
      track.nOutliers()++;
    } else if (typeFlags.isMeasurement()) {
      if (typeFlags.isSharedHit()) {
        track.nSharedHits()++;
      }
      if (typeFlags.isSplitHit()) {
        nSplitHits++;
      }
      track.nMeasurements()++;
      track.chi2() += trackState.chi2();
      track.nDoF() += trackState.calibratedSize();
    }
  }

  if constexpr (requires { track.nSplitHits(); }) {
    if (track.hasColumn(detail_tp::kSplitHitsKey)) {
      track.nSplitHits() = nSplitHits;
    }
  }
}

/// Helper function to trim track states from the front of a track
/// @tparam track_proxy_t the track proxy type
/// @param track the track to trim
/// @param trimHoles whether to trim holes
/// @param trimOutliers whether to trim outliers
/// @param trimMaterial whether to trim pure material states
/// @param trimOtherNoneMeasurement whether to trim other, non measurement, states
template <TrackProxyConcept track_proxy_t>
void trimTrackFront(track_proxy_t track, bool trimHoles, bool trimOutliers,
                    bool trimMaterial, bool trimOtherNoneMeasurement)
  requires(!track_proxy_t::ReadOnly)
{
  using TrackStateProxy = typename track_proxy_t::TrackStateProxy;

  // TODO specialize if track is forward linked

  std::optional<TrackStateProxy> front;

  for (TrackStateProxy trackState : track.trackStatesReversed()) {
    TrackStateTypeMap typeFlags = trackState.typeFlags();
    bool isHole = typeFlags.isHole();
    bool isOutlier = typeFlags.isOutlier();
    bool isMaterial = typeFlags.isMaterial();
    bool isOtherNoneMeasurement =
        !typeFlags.hasMeasurement() && !isHole && !isOutlier && !isMaterial;
    if (trimHoles && isHole) {
      continue;
    }
    if (trimOutliers && isOutlier) {
      continue;
    }
    if (trimMaterial && isMaterial) {
      continue;
    }
    if (trimOtherNoneMeasurement && isOtherNoneMeasurement) {
      continue;
    }

    front = trackState;
  }

  if (front.has_value()) {
    front.value().previous() = TrackStateProxy::kInvalid;
  }
}

/// Helper function to trim track states from the back of a track
/// @tparam track_proxy_t the track proxy type
/// @param track the track to trim
/// @param trimHoles whether to trim holes
/// @param trimOutliers whether to trim outliers
/// @param trimMaterial whether to trim pure material states
/// @param trimOtherNoneMeasurement whether to trim other, non measurement, states
template <TrackProxyConcept track_proxy_t>
void trimTrackBack(track_proxy_t track, bool trimHoles, bool trimOutliers,
                   bool trimMaterial, bool trimOtherNoneMeasurement)
  requires(!track_proxy_t::ReadOnly)
{
  using TrackStateProxy = typename track_proxy_t::TrackStateProxy;

  std::optional<TrackStateProxy> back;

  for (TrackStateProxy trackState : track.trackStatesReversed()) {
    back = trackState;

    TrackStateTypeMap typeFlags = trackState.typeFlags();
    bool isHole = typeFlags.isHole();
    bool isOutlier = typeFlags.isOutlier();
    bool isMaterial = typeFlags.isMaterial();
    bool isOtherNoneMeasurement =
        !typeFlags.hasMeasurement() && !isHole && !isOutlier && !isMaterial;
    if (trimHoles && isHole) {
      continue;
    }
    if (trimOutliers && isOutlier) {
      continue;
    }
    if (trimMaterial && isMaterial) {
      continue;
    }
    if (trimOtherNoneMeasurement && isOtherNoneMeasurement) {
      continue;
    }

    break;
  }

  if (back.has_value()) {
    track.tipIndex() = back.value().index();
  }
}

/// Helper function to trim track states from the front and back of a track
/// @tparam track_proxy_t the track proxy type
/// @param track the track to trim
/// @param trimHoles whether to trim holes
/// @param trimOutliers whether to trim outliers
/// @param trimMaterial whether to trim pure material states
/// @param trimOtherNoneMeasurement whether to trim other, non measurement, states
template <TrackProxyConcept track_proxy_t>
void trimTrack(track_proxy_t track, bool trimHoles, bool trimOutliers,
               bool trimMaterial, bool trimOtherNoneMeasurement)
  requires(!track_proxy_t::ReadOnly)
{
  trimTrackFront(track, trimHoles, trimOutliers, trimMaterial,
                 trimOtherNoneMeasurement);
  trimTrackBack(track, trimHoles, trimOutliers, trimMaterial,
                trimOtherNoneMeasurement);
}

/// Helper function to calculate the predicted residual and its covariance
/// @tparam nMeasurementDim the dimension of the measurement
/// @tparam track_state_proxy_t the track state proxy type
/// @param trackState the track state to calculate the residual from
/// @return a pair of the residual and its covariance
template <std::size_t nMeasurementDim,
          TrackStateProxyConcept track_state_proxy_t>
std::pair<Vector<nMeasurementDim>, SquareMatrix<nMeasurementDim>>
calculatePredictedResidual(track_state_proxy_t trackState) {
  using MeasurementVector = Vector<nMeasurementDim>;
  using MeasurementMatrix = SquareMatrix<nMeasurementDim>;

  if (!trackState.hasPredicted()) {
    throw std::invalid_argument("track state has no predicted parameters");
  }
  if (!trackState.hasCalibrated()) {
    throw std::invalid_argument("track state has no calibrated parameters");
  }

  auto subspaceHelper =
      trackState.template projectorSubspaceHelper<nMeasurementDim>();

  auto measurement = trackState.template calibrated<nMeasurementDim>();
  auto measurementCovariance =
      trackState.template calibratedCovariance<nMeasurementDim>();
  MeasurementVector predicted =
      subspaceHelper.projectVector(trackState.predicted());
  MeasurementMatrix predictedCovariance =
      subspaceHelper.projectMatrix(trackState.predictedCovariance());

  MeasurementVector residual = measurement - predicted;
  MeasurementMatrix residualCovariance =
      measurementCovariance + predictedCovariance;

  return {residual, residualCovariance};
}

/// Helper function to calculate the filtered residual and its covariance
/// @tparam nMeasurementDim the dimension of the measurement
/// @tparam track_state_proxy_t the track state proxy type
/// @param trackState the track state to calculate the residual from
/// @return a pair of the residual and its covariance
template <std::size_t nMeasurementDim,
          TrackStateProxyConcept track_state_proxy_t>
std::pair<Vector<nMeasurementDim>, SquareMatrix<nMeasurementDim>>
calculateFilteredResidual(track_state_proxy_t trackState) {
  using MeasurementVector = Vector<nMeasurementDim>;
  using MeasurementMatrix = SquareMatrix<nMeasurementDim>;

  if (!trackState.hasFiltered()) {
    throw std::invalid_argument("track state has no filtered parameters");
  }
  if (!trackState.hasCalibrated()) {
    throw std::invalid_argument("track state has no calibrated parameters");
  }

  auto subspaceHelper =
      trackState.template projectorSubspaceHelper<nMeasurementDim>();

  auto measurement = trackState.template calibrated<nMeasurementDim>();
  auto measurementCovariance =
      trackState.template calibratedCovariance<nMeasurementDim>();
  MeasurementVector filtered =
      subspaceHelper.projectVector(trackState.filtered());
  MeasurementMatrix filteredCovariance =
      subspaceHelper.projectMatrix(trackState.filteredCovariance());

  MeasurementVector residual = measurement - filtered;
  MeasurementMatrix residualCovariance =
      measurementCovariance - filteredCovariance;

  return {residual, residualCovariance};
}

/// Helper function to calculate the smoothed residual and its covariance
/// @tparam nMeasurementDim the dimension of the measurement
/// @tparam track_state_proxy_t the track state proxy type
/// @param trackState the track state to calculate the residual from
/// @return a pair of the residual and its covariance
template <std::size_t nMeasurementDim,
          TrackStateProxyConcept track_state_proxy_t>
std::pair<Vector<nMeasurementDim>, SquareMatrix<nMeasurementDim>>
calculateSmoothedResidual(track_state_proxy_t trackState) {
  using MeasurementVector = Vector<nMeasurementDim>;
  using MeasurementMatrix = SquareMatrix<nMeasurementDim>;

  if (!trackState.hasSmoothed()) {
    throw std::invalid_argument("track state has no smoothed parameters");
  }
  if (!trackState.hasCalibrated()) {
    throw std::invalid_argument("track state has no calibrated parameters");
  }

  auto subspaceHelper =
      trackState.template projectorSubspaceHelper<nMeasurementDim>();

  auto measurement = trackState.template calibrated<nMeasurementDim>();
  auto measurementCovariance =
      trackState.template calibratedCovariance<nMeasurementDim>();
  MeasurementVector smoothed =
      subspaceHelper.projectVector(trackState.smoothed());
  MeasurementMatrix smoothedCovariance =
      subspaceHelper.projectMatrix(trackState.smoothedCovariance());

  MeasurementVector residual = measurement - smoothed;
  MeasurementMatrix residualCovariance =
      measurementCovariance - smoothedCovariance;

  return {residual, residualCovariance};
}

/// Helper function to calculate the predicted chi2
/// @tparam track_state_proxy_t the track state proxy type
/// @param trackState the track state to calculate the chi2 from
/// @return the chi2
template <TrackStateProxyConcept track_state_proxy_t>
double calculatePredictedChi2(track_state_proxy_t trackState) {
  if (!trackState.hasPredicted()) {
    throw std::invalid_argument("track state has no predicted parameters");
  }
  if (!trackState.hasCalibrated()) {
    throw std::invalid_argument("track state has no calibrated parameters");
  }

  return visit_measurement(
      trackState.calibratedSize(),
      [&]<std::size_t measdim>(
          std::integral_constant<std::size_t, measdim>) -> double {
        auto [residual, residualCovariance] =
            calculatePredictedResidual<measdim>(trackState);

        return (residual.transpose() * residualCovariance.inverse() * residual)
            .eval()(0, 0);
      });
}

/// Helper function to calculate the filtered chi2
/// @tparam track_state_proxy_t the track state proxy type
/// @param trackState the track state to calculate the chi2 from
/// @return the chi2
template <TrackStateProxyConcept track_state_proxy_t>
double calculateFilteredChi2(track_state_proxy_t trackState) {
  if (!trackState.hasFiltered()) {
    throw std::invalid_argument("track state has no filtered parameters");
  }
  if (!trackState.hasCalibrated()) {
    throw std::invalid_argument("track state has no calibrated parameters");
  }

  return visit_measurement(
      trackState.calibratedSize(),
      [&]<std::size_t measdim>(
          std::integral_constant<std::size_t, measdim>) -> double {
        auto [residual, residualCovariance] =
            calculateFilteredResidual<measdim>(trackState);

        return (residual.transpose() * residualCovariance.inverse() * residual)
            .eval()(0, 0);
      });
}

/// Helper function to calculate the smoothed chi2
/// @tparam track_state_proxy_t the track state proxy type
/// @param trackState the track state to calculate the chi2 from
/// @return the chi2
template <TrackStateProxyConcept track_state_proxy_t>
double calculateSmoothedChi2(track_state_proxy_t trackState) {
  if (!trackState.hasSmoothed()) {
    throw std::invalid_argument("track state has no smoothed parameters");
  }
  if (!trackState.hasCalibrated()) {
    throw std::invalid_argument("track state has no calibrated parameters");
  }

  return visit_measurement(
      trackState.calibratedSize(),
      [&]<std::size_t measdim>(
          std::integral_constant<std::size_t, measdim>) -> double {
        auto [residual, residualCovariance] =
            calculateSmoothedResidual<measdim>(trackState);

        return (residual.transpose() * residualCovariance.inverse() * residual)
            .eval()(0, 0);
      });
}

/// Helper function to calculate the unbiased track parameters and their
/// covariance (i.e. fitted track parameters with this measurement removed)
/// using Eq.(12a)-Eq.(12c) of NIMA 262, 444 (1987)
/// @tparam track_state_proxy_t the track state proxy type
/// @param trackState the track state to calculate the unbiased parameters from
/// @return a pair of the unbiased parameters and their covariance
/// @deprecated Instantiating this template is very expensive in compiler memory
///   (it expands the Eigen expression templates over all measurement
///   dimensions) and it does so in every calling translation unit. Prefer the
///   non-template overload taking a type-erased @c AnyConstTrackStateProxy,
///   which is compiled once in the Acts core library:
///   @code
///   calculateUnbiasedParametersCovariance(Acts::AnyConstTrackStateProxy{state});
///   @endcode
template <TrackStateProxyConcept track_state_proxy_t>
[[deprecated(
    "Use calculateUnbiasedParametersCovariance(const AnyConstTrackStateProxy&) "
    "instead; the templated form instantiates expensive Eigen code in every "
    "translation unit.")]]
std::pair<BoundVector, BoundMatrix> calculateUnbiasedParametersCovariance(
    track_state_proxy_t trackState) {
  // Explicitly select the non-template overload taking a type-erased proxy.
  // A plain call here would re-resolve to this very template (the wrapped
  // AnyConstTrackStateProxy satisfies TrackStateProxyConcept and is an exact
  // by-value match), causing infinite recursion and a self-deprecation error.
  std::pair<BoundVector, BoundMatrix> (&impl)(const AnyConstTrackStateProxy &) =
      calculateUnbiasedParametersCovariance;
  return impl(AnyConstTrackStateProxy{trackState});
}

/// Calculate the unbiased track parameters and their covariance for a
/// type-erased track state proxy. See the templated overload above for the
/// underlying formula.
///
/// This is the preferred entry point. It is not a template, so the (very
/// expensive) Eigen expression templates are instantiated exactly once, in the
/// Acts core library (TrackHelpersUnbiased.cpp), instead of in every calling
/// translation unit. Callers holding a concrete track state proxy wrap it:
/// @code
/// calculateUnbiasedParametersCovariance(Acts::AnyConstTrackStateProxy{state});
/// @endcode
/// @param trackState the (type-erased) track state to calculate from
/// @return a pair of the unbiased parameters and their covariance
std::pair<BoundVector, BoundMatrix> calculateUnbiasedParametersCovariance(
    const AnyConstTrackStateProxy &trackState);

}  // namespace Acts

namespace std {
// register with STL
template <>
struct is_error_code_enum<Acts::TrackExtrapolationError> : std::true_type {};
}  // namespace std
