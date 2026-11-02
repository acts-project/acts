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
#include "Acts/Utilities/Delegate.hpp"
#include "Acts/Utilities/TrackHelpers.hpp"

#include <cstddef>
#include <ostream>

#include <detray/geometry/identifier.hpp>
#include <detray/geometry/tracking_surface.hpp>
#include <traccc/definitions/primitives.hpp>
#include <traccc/edm/track_container.hpp>
#include <traccc/edm/track_fit_outcome.hpp>
#include <traccc/edm/track_parameters.hpp>

namespace ActsPlugins {

/// traccc bound track parameters with the default algebra
using TracccBoundTrackParameters =
    traccc::bound_track_parameters<traccc::default_algebra>;

/// traccc track container with the default algebra
using TracccTrackContainer =
    traccc::edm::track_container<traccc::default_algebra>;

/// Convert the parameter vector of traccc bound track parameters
Acts::BoundVector convertBoundVectorFromTraccc(
    const TracccBoundTrackParameters& params);

/// Convert the covariance of traccc bound track parameters
Acts::BoundMatrix convertBoundCovarianceFromTraccc(
    const TracccBoundTrackParameters& params);

namespace ActsPlugins::detail {
/// Check for degenerate (all zero) traccc bound track parameters, which
/// traccc reports for parameters that were never set
bool hasZeroBoundParameters(const TracccBoundTrackParameters& params);
}  // namespace ActsPlugins::detail

/// Create the source link of a measurement track state
///
/// @param measurementIndex index of the measurement in the traccc
///        measurement collection the tracks refer to
/// @param surface the ACTS surface of the measurement
using TracccSourceLinkCreator = Acts::Delegate<Acts::SourceLink(
    unsigned int measurementIndex, const Acts::Surface& surface)>;

/// Options for the traccc to ACTS track conversion
struct TracccTrackConversionOptions {
  /// Particle hypothesis assigned to the converted tracks
  Acts::ParticleHypothesis particleHypothesis =
      Acts::ParticleHypothesis::pion();
};

/// Bookkeeping of the traccc to ACTS track conversion
struct TracccTrackConversionStatistics {
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
                         const TracccTrackConversionStatistics& stats);

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
const Acts::Surface* findSurfaceFromDetray(
    const detector_t& detector, const Acts::TrackingGeometry& trackingGeometry,
    const detray::geometry::identifier& surfaceLink) {
  if (surfaceLink.is_invalid()) {
    throw std::runtime_error("find surface: could not find reference surface");
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
TracccTrackConversionStatistics convertTracksFromTraccc(
    const TracccTrackContainer::const_view& tracks, const detector_t& detector,
    const Acts::TrackingGeometry& trackingGeometry,
    const TracccSourceLinkCreator& makeSourceLink, track_container_t& output,
    const TracccTrackConversionOptions& options = {}) {
  const TracccTrackContainer::const_device input(tracks);

  constexpr auto kStateMask =
      Acts::TrackStatePropMask::Filtered | Acts::TrackStatePropMask::Smoothed;

  TracccTrackConversionStatistics stats;
  stats.nInput = input.tracks.size();

  for (unsigned int iTrack = 0; iTrack < input.tracks.size(); ++iTrack) {
    const auto track = input.tracks.at(iTrack);
    if (track.fit_outcome() != traccc::track_fit_outcome::SUCCESS) {
      ++stats.nFailedFit;
      continue;
    }
    if (ActsPlugins::detail::hasZeroBoundParameters(track.params())) {
      ++stats.nDegenerate;
      continue;
    }

    const Acts::Surface* referenceSurface = findSurfaceFromDetray(
        detector, trackingGeometry, track.params().surface_link());
    if (referenceSurface == nullptr) {
      ++stats.nMissingSurface;
      throw std::runtime_error(
          "TracccTrackToActsAlg: could not find reference surface for track "
          "parameters");
    }

    auto actsTrack = output.makeTrack();
    actsTrack.setReferenceSurface(referenceSurface->getSharedPtr());
    actsTrack.parameters() = convertBoundVectorFromTraccc(track.params());
    actsTrack.covariance() = convertBoundCovarianceFromTraccc(track.params());
    actsTrack.setParticleHypothesis(options.particleHypothesis);

    bool missingSurface = false;
    bool degenerate = false;
    for (const auto& link : track.constituent_links()) {
      if (link.type != traccc::edm::track_constituent_link::track_state) {
        continue;
      }
      const auto state = input.states.at(link.index);
      const auto& smoothed = state.smoothed_params();
      const Acts::Surface* surface = findSurfaceFromDetray(
          detector, trackingGeometry, smoothed.surface_link());
      if (surface == nullptr) {
        missingSurface = true;
        throw std::runtime_error(
            "TracccTrackToActsAlg: could not find reference surface for "
            "smoothed track parameters");
      }

      if (state.is_hole()) {
        auto actsState =
            actsTrack.appendTrackState(Acts::TrackStatePropMask::None);
        actsState.setReferenceSurface(surface->getSharedPtr());
        actsState.typeFlags().setIsHole();
        continue;
      }

      if (ActsPlugins::detail::hasZeroBoundParameters(smoothed)) {
        degenerate = true;
        break;
      }

      auto actsState = actsTrack.appendTrackState(kStateMask);
      actsState.setReferenceSurface(surface->getSharedPtr());
      actsState.filtered() =
          convertBoundVectorFromTraccc(state.filtered_params());
      actsState.filteredCovariance() =
          convertBoundCovarianceFromTraccc(state.filtered_params());
      actsState.smoothed() = convertBoundVectorFromTraccc(smoothed);
      actsState.smoothedCovariance() =
          convertBoundCovarianceFromTraccc(smoothed);
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

}  // namespace ActsPlugins
