// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// Plugin include(s)
#include "Acts/Plugins/Traccc/Detail/AlgebraConversion.hpp"

// Acts include(s)
#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/TrackProxy.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"

// Detray include(s)
#include "detray/core/detector.hpp"
#include "detray/tracks/bound_track_parameters.hpp"

// Traccc include(s)
#include "traccc/edm/track_state.hpp"

// Boost include(s)
#include <boost/range/combine.hpp>

// System include(s)
#include <memory>
#include <stdexcept>
#include <variant>

namespace Acts::TracccPlugin {

/// @brief Creates a new Acts bound track parameters from detray bound track parameters.
/// @param dparams the detray bound track parameters.
/// @param detector the detray detector.
/// @param trackingGeometry the Acts tracking geometry.
/// @return An Acts BoundTrackParameters with data copied from a detray bound_track_parameters.
template <typename algebra_t, typename metadata_t, typename container_t>
inline auto newParams(const detray::bound_track_parameters<algebra_t>& dparams,
                      const detray::detector<metadata_t, container_t>& detector,
                      const Acts::TrackingGeometry& trackingGeometry) {
  Acts::ActsVector<6U> parameterVector =
      Detail::newVector<6U>(dparams.vector()[0]);
  typename Acts::BoundTrackParameters::CovarianceMatrix cov =
      Detail::newSqaureMatrix<6U>(dparams.covariance());
  Acts::ParticleHypothesis particleHypothesis =
      Acts::ParticleHypothesis::pion();

  auto geoID =
      Acts::GeometryIdentifier(detector.surface(dparams.surface_link()).source);

  auto surface = trackingGeometry.findSurface(geoID);

  if (surface == nullptr) {
    throw std::runtime_error(
        "Mismatch between Acts geometry and detray detector: Acts tracking "
        "geometry does not contain geometry ID " +
        std::to_string(geoID.value()));
  }

  Acts::BoundTrackParameters params(surface->getSharedPtr(), parameterVector,
                                    std::make_optional(std::move(cov)),
                                    particleHypothesis);

  return params;
}

/// @brief Copies data from a traccc fitting result to an Acts track proxy.
/// @param source the traccc track fitting result to copy from.
/// @param destination the Acts track proxy to copy to.
/// @param detector the detray detector of the traccc track fitting result.
/// @param trackingGeometry the Acts tracking geometry.
template <typename algebra_t, typename track_container_t, typename trajectory_t,
          template <typename> class holder_t, typename metadata_t,
          typename container_t>
inline void copyFittingResult(
    const traccc::fitting_result<algebra_t>& source,
    Acts::TrackProxy<track_container_t, trajectory_t, holder_t, false>&
        destination,
    const detray::detector<metadata_t, container_t>& detector,
    const Acts::TrackingGeometry& trackingGeometry) {
  const auto params = newParams(source.fit_params, detector, trackingGeometry);
  // track.tipIndex() = kalmanResult.lastMeasurementIndex;
  destination.parameters() = params.parameters();
  destination.covariance() = params.covariance().value();
  destination.setReferenceSurface(params.referenceSurface().getSharedPtr());
}

/// @brief Copies data from a traccc track state to a Acts track state proxy.
/// @param source the traccc track state to copy from.
/// @param destination the Acts track state proxy to copy to.
/// @param detector the detray detector of the traccc track track state.
/// @param trackingGeometry the Acts tracking geometry.
/// @note Does not set the uncalibrated source link and calibrated measurements.
template <typename algebra_t, typename metadata_t, typename container_t,
          typename trajectory_t, std::size_t M>
inline void copyTrackState(
    const traccc::track_state<algebra_t>& source,
    Acts::TrackStateProxy<trajectory_t, M, false>& destination,
    const detray::detector<metadata_t, container_t>& detector,
    const Acts::TrackingGeometry& trackingGeometry) {
  auto geoID =
      Acts::GeometryIdentifier(detector.surface(source.surface_link()).source);
  auto surface = trackingGeometry.findSurface(geoID)->getSharedPtr();
  destination.setReferenceSurface(surface);

  using Parameters =
      typename Acts::TrackStateProxy<trajectory_t, M, false>::Parameters;
  using Covariance =
      typename Acts::TrackStateProxy<trajectory_t, M, false>::Covariance;

  destination.predicted() =
      Parameters(Detail::newVector<6U>(source.predicted().vector()[0]).data());
  destination.predictedCovariance() = Covariance(
      Detail::newSqaureMatrix<6U>(source.predicted().covariance()).data());

  destination.smoothed() =
      Parameters(Detail::newVector<6U>(source.smoothed().vector()[0]).data());
  destination.smoothedCovariance() = Covariance(
      Detail::newSqaureMatrix<6U>(source.smoothed().covariance()).data());

  destination.filtered() =
      Parameters(Detail::newVector<6U>(source.filtered().vector()[0]).data());
  destination.filteredCovariance() = Covariance(
      Detail::newSqaureMatrix<6U>(source.filtered().covariance()).data());

  destination.jacobian() =
      Covariance(Detail::newSqaureMatrix<6U>(source.jacobian()).data());

  destination.chi2() = source.smoothed_chi2();

  auto typeFlags = destination.typeFlags();
  typeFlags.set(TrackStateFlag::ParameterFlag);
  if (surface->surfaceMaterial() != nullptr) {
    typeFlags.set(TrackStateFlag::MaterialFlag);
  }
  if (source.is_hole) {
    typeFlags.set(TrackStateFlag::HoleFlag);
  }
  typeFlags.set(TrackStateFlag::MeasurementFlag);
}

/// @brief Gets the fitting result of a traccc container element.
/// @param e the traccc container element.
/// @return the fitting result contained in the container element (i.e., e.header).
template <typename fitting_result_t, typename track_state_vector_t>
inline auto& getFittingResult(
    const traccc::container_element<fitting_result_t, track_state_vector_t>&
        tracccTrack) {
  return tracccTrack.header;
}

/// @brief Gets the track states of a traccc container element.
/// @param e the traccc container element.
/// @return the track states contained in the container element (i.e., e.items).
template <typename fitting_result_t, typename track_state_vector_t>
inline auto& getTrackStates(
    const traccc::container_element<fitting_result_t, track_state_vector_t>&
        tracccTrack) {
  return tracccTrack.items;
}

/// @brief Creates a new track in the Acts track container.
/// This new track will contain data copied from the traccc track container
/// element (track and track state data).
/// @param tracccTrack The traccc container element to copy from.
/// @param actsTrackContainer The Acts track container. This is the new track will be made in this container.
/// @param detector The detray detector.
/// @param trackingGeometry The Acts tracking geometry.
/// @note Does not set the uncalibrated source link and calibrated measurements.
template <typename fitting_result_t, typename track_state_vector_t,
          typename track_container_t, typename trajectory_t,
          template <typename> class holder_t, typename metadata_t,
          typename container_t>
inline auto makeTrack(
    const traccc::container_element<fitting_result_t, track_state_vector_t>&
        tracccTrack,
    Acts::TrackContainer<track_container_t, trajectory_t, holder_t>&
        actsTrackContainer,
    const detray::detector<metadata_t, container_t>& detector,
    const Acts::TrackingGeometry& trackingGeometry) {
  auto fittingResult = getFittingResult(tracccTrack);
  auto trackStates = getTrackStates(tracccTrack);

  auto track = actsTrackContainer.makeTrack();
  copyFittingResult(fittingResult, track, detector, trackingGeometry);

  // Make the track states.
  for (const auto& tstate : trackStates) {
    auto astate = track.appendTrackState();
    copyTrackState(tstate, astate, detector, trackingGeometry);
  }

  track.linkForward();

  return track;
}

/// @brief Creates a vector of Acts and traccc track state pairs.
/// The states in the pairs are paired 1:1 by index.
/// @tparam traccc_track_t type of traccc container_element.
/// @tparam acts_track_t type of Acts track proxy.
/// @param tracccTrack the traccc track.
/// @param actsTrack the Acts track.
/// @returns a vector of tuples.
/// @note The Acts track and traccc track must contain the same number of track states.
template <typename traccc_track_t, typename acts_track_t>
auto trackStateZipView(traccc_track_t& tracccTrack, acts_track_t& actsTrack) {
  auto& tracccTrackStates = getTrackStates(tracccTrack);
  auto actsTrackStates = actsTrack.trackStates();

  using TracccTrackStateType = typename std::iterator_traits<
      decltype(tracccTrackStates.begin())>::value_type;
  using ActsTrackStateType = typename std::iterator_traits<
      decltype(actsTrackStates.begin())>::value_type;

  std::vector<std::tuple<TracccTrackStateType, ActsTrackStateType>>
      zippedTrackStates;
  zippedTrackStates.reserve(tracccTrackStates.size());

  auto tracccIter = tracccTrackStates.begin();
  auto actsIter = actsTrackStates.begin();
  while (tracccIter != tracccTrackStates.end() &&
         actsIter != actsTrackStates.end()) {
    zippedTrackStates.emplace_back(*tracccIter, *actsIter);
    ++tracccIter;
    ++actsIter;
  }

  return zippedTrackStates;
}

/// @brief Sets the uncalibrated source link and calibrated measurement of the Acts track state.
/// @param actsTrackState the Acts track state.
/// @param measurement the Acts bound variant measurement.
/// @note The uncalibrated source link and calibrated measurement are set to the same measurement.
template <typename trajectory_t, std::size_t M>
void setSourceAndMeasurement(
    Acts::TrackStateProxy<trajectory_t, M, false>& actsTrackState,
    const Acts::BoundVariantMeasurement& measurement) {
  std::visit(
      [&actsTrackState](auto& m) {
        actsTrackState.setUncalibratedSourceLink(m.sourceLink());
        actsTrackState.setCalibrated(m);
      },
      measurement);
}

/// @brief Sets the uncalibrated source link and calibrated measurement of the Acts track state
/// using the measurement data in their traccc pair.
/// @param trackStatePairs a list of tuples: (traccc track state, Acts track state)
/// @param map the map from traccc measurement to acts measurement.
template <typename track_state_pairs_t>
void setSourceAndMeasurements(
    track_state_pairs_t& trackStatePairs,
    const std::map<traccc::measurement, Acts::BoundVariantMeasurement>& map) {
  for (auto pair : trackStatePairs) {
    // First item in pair is the traccc track.
    // Second item in pair is the acts track.
    auto& measurement = map.at(std::get<0>(pair).get_measurement());
    setSourceAndMeasurement(std::get<1>(pair), measurement);
  }
}

}  // namespace Acts::TracccPlugin
