// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/TrackStatePropMask.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/TrackHelpers.hpp"

#include <cstddef>
#include <functional>
#include <ostream>

#include <detray/geometry/identifier.hpp>
#include <detray/geometry/tracking_surface.hpp>
#include <traccc/definitions/primitives.hpp>
#include <traccc/edm/track_container.hpp>
#include <traccc/edm/track_fit_outcome.hpp>
#include <traccc/edm/track_parameters.hpp>

namespace ActsPlugins::Traccc {

/// traccc bound track parameters with the default algebra
using BoundTrackParameters =
    traccc::bound_track_parameters<traccc::default_algebra>;

/// traccc track container with the default algebra
using TrackContainer = traccc::edm::track_container<traccc::default_algebra>;

/// Convert the parameter vector of traccc bound track parameters
Acts::BoundVector convertBoundVector(const BoundTrackParameters& params);

/// Convert the covariance of traccc bound track parameters
Acts::BoundMatrix convertBoundCovariance(const BoundTrackParameters& params);

/// Check for degenerate (all zero) traccc bound track parameters, which
/// traccc reports for parameters that were never set
bool isDegenerate(const BoundTrackParameters& params);

/// Create the source link of a measurement track state
///
/// @param measurementIndex index of the measurement in the traccc
///        measurement collection the tracks refer to
/// @param surface the ACTS surface of the measurement
using SourceLinkCreator = std::function<Acts::SourceLink(
    unsigned int measurementIndex, const Acts::Surface& surface)>;

/// Options for the traccc to ACTS track conversion
struct TrackConversionOptions {
  /// Particle hypothesis assigned to the converted tracks
  Acts::ParticleHypothesis particleHypothesis =
      Acts::ParticleHypothesis::pion();
};

/// Bookkeeping of the traccc to ACTS track conversion
struct TrackConversionStatistics {
  /// Number of input traccc tracks
  std::size_t nInput = 0;
  /// Number of converted tracks
  std::size_t nConverted = 0;
  /// Tracks dropped because the fit was not successful
  std::size_t nFailedFit = 0;
  /// Tracks dropped because a surface could not be found
  std::size_t nMissingSurface = 0;
  /// Tracks dropped because of degenerate track or track state parameters
  std::size_t nDegenerate = 0;
};

std::ostream& operator<<(std::ostream& os,
                         const TrackConversionStatistics& stats);

/// Find the ACTS surface of a detray surface
///
/// The detray surface source is the ACTS geometry identifier, as set by the
/// ACTS to detray geometry conversion.
///
/// @param detector the detray detector converted from @p trackingGeometry
/// @param trackingGeometry the ACTS tracking geometry
/// @param surfaceLink the detray surface identifier
///
/// @return the surface or nullptr if it can not be found
template <typename detector_t>
const Acts::Surface* findSurface(
    const detector_t& detector, const Acts::TrackingGeometry& trackingGeometry,
    const detray::geometry::identifier& surfaceLink) {
  if (surfaceLink.is_invalid()) {
    return nullptr;
  }
  const detray::tracking_surface surface{detector, surfaceLink};
  return trackingGeometry.findSurface(
      Acts::GeometryIdentifier{surface.source()});
}

/// Convert traccc tracks into ACTS tracks
///
/// Only tracks with a successful fit are converted. Measurement track states
/// get the filtered and smoothed parameters, holes only their surface. The
/// track parameters are the traccc fitted parameters, on the surface of the
/// first measurement.
///
/// @param tracks the traccc tracks and track states (on the host)
/// @param detector the detray detector the tracks were reconstructed with,
///        converted from @p trackingGeometry
/// @param trackingGeometry the ACTS tracking geometry
/// @param makeSourceLink creates the source link of a measurement state
/// @param output the ACTS track container to add the tracks to
/// @param options conversion options
///
/// @return conversion bookkeeping
template <typename detector_t, typename track_container_t>
TrackConversionStatistics convertTracks(
    const TrackContainer::const_view& tracks, const detector_t& detector,
    const Acts::TrackingGeometry& trackingGeometry,
    const SourceLinkCreator& makeSourceLink, track_container_t& output,
    const TrackConversionOptions& options = {}) {
  const TrackContainer::const_device input(tracks);

  constexpr auto kStateMask =
      Acts::TrackStatePropMask::Filtered | Acts::TrackStatePropMask::Smoothed;

  TrackConversionStatistics stats;
  stats.nInput = input.tracks.size();

  for (unsigned int iTrack = 0; iTrack < input.tracks.size(); ++iTrack) {
    const auto track = input.tracks.at(iTrack);
    if (track.fit_outcome() != traccc::track_fit_outcome::SUCCESS) {
      ++stats.nFailedFit;
      continue;
    }
    if (isDegenerate(track.params())) {
      ++stats.nDegenerate;
      continue;
    }

    const Acts::Surface* referenceSurface = findSurface(
        detector, trackingGeometry, track.params().surface_link());
    if (referenceSurface == nullptr) {
      ++stats.nMissingSurface;
      continue;
    }

    auto actsTrack = output.makeTrack();
    actsTrack.setReferenceSurface(referenceSurface->getSharedPtr());
    actsTrack.parameters() = convertBoundVector(track.params());
    actsTrack.covariance() = convertBoundCovariance(track.params());
    actsTrack.setParticleHypothesis(options.particleHypothesis);

    bool missingSurface = false;
    bool degenerate = false;
    for (const auto& link : track.constituent_links()) {
      if (link.type != traccc::edm::track_constituent_link::track_state) {
        continue;
      }
      const auto state = input.states.at(link.index);
      const auto& smoothed = state.smoothed_params();
      const Acts::Surface* surface = findSurface(
          detector, trackingGeometry, smoothed.surface_link());
      if (surface == nullptr) {
        missingSurface = true;
        break;
      }

      if (state.is_hole()) {
        auto actsState =
            actsTrack.appendTrackState(Acts::TrackStatePropMask::None);
        actsState.setReferenceSurface(surface->getSharedPtr());
        actsState.typeFlags().setIsHole();
        continue;
      }

      if (isDegenerate(smoothed)) {
        degenerate = true;
        break;
      }

      auto actsState = actsTrack.appendTrackState(kStateMask);
      actsState.setReferenceSurface(surface->getSharedPtr());
      actsState.filtered() = convertBoundVector(state.filtered_params());
      actsState.filteredCovariance() =
          convertBoundCovariance(state.filtered_params());
      actsState.smoothed() = convertBoundVector(smoothed);
      actsState.smoothedCovariance() = convertBoundCovariance(smoothed);
      actsState.chi2() = state.smoothed_chi2();
      actsState.typeFlags().setIsMeasurement();
      actsState.setUncalibratedSourceLink(
          makeSourceLink(state.measurement_index(), *surface));
    }

    if (missingSurface || degenerate) {
      output.removeTrack(actsTrack.index());
      ++(missingSurface ? stats.nMissingSurface : stats.nDegenerate);
      continue;
    }

    Acts::calculateTrackQuantities(actsTrack);
    actsTrack.chi2() = track.chi2();
    actsTrack.nDoF() = static_cast<unsigned int>(track.ndf());
    ++stats.nConverted;
  }

  return stats;
}

}  // namespace ActsPlugins::Traccc
