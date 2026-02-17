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
#include "Acts/EventData/GenericBoundTrackParameters.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/TrackProxyConcept.hpp"
#include "Acts/EventData/TrackStatePropMask.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <stdexcept>

#include <edm4hep/MCParticle.h>
#include <edm4hep/MutableSimTrackerHit.h>
#include <edm4hep/MutableTrack.h>
#include <edm4hep/SimTrackerHit.h>
#include <edm4hep/Track.h>
#include <edm4hep/TrackState.h>
#include <podio/podioVersion.h>

#if podio_VERSION_MAJOR == 0 || \
    (podio_VERSION_MAJOR == 1 && podio_VERSION_MINOR <= 2)

template <>
struct std::hash<podio::ObjectID> {
  std::size_t operator()(const podio::ObjectID& id) const noexcept {
    auto hash_collectionID = std::hash<std::uint32_t>{}(id.collectionID);
    auto hash_index = std::hash<int>{}(id.index);

    return hash_collectionID ^ hash_index;
  }
};

#endif

/// @namespace ActsPlugins::EDM4hepUtil
/// @ingroup edm4hep_plugin

namespace ActsPlugins::EDM4hepUtil {

static constexpr std::int32_t EDM4HEP_ACTS_POSITION_TYPE = 42;

namespace detail {
struct Parameters {
  Acts::Vector<6> values{};
  // Dummy default
  Acts::ParticleHypothesis particleHypothesis =
      Acts::ParticleHypothesis::pion();
  std::optional<Acts::BoundMatrix> covariance;
  std::shared_ptr<const Acts::Surface> surface;
};

Acts::SquareMatrix<6> jacobianToEdm4hep(double theta, double qOverP, double Bz);

Acts::SquareMatrix<6> jacobianFromEdm4hep(double tanLambda, double omega,
                                          double Bz);

void unpackCovariance(const float* from, Acts::SquareMatrix<6>& to);
void packCovariance(const Acts::SquareMatrix<6>& from, float* to);

Parameters convertTrackParametersToEdm4hep(
    const Acts::GeometryContext& gctx, double Bz,
    const Acts::BoundTrackParameters& params);

Acts::BoundTrackParameters convertTrackParametersFromEdm4hep(
    double Bz, const Parameters& params);

}  // namespace detail

/// @addtogroup edm4hep_plugin
/// @{

/// Get the particle from a SimTrackerHit (compatibility with EDM4hep < 0.99 and
/// >= 0.99)
/// @param hit The SimTrackerHit
/// @return The associated MCParticle
edm4hep::MCParticle getParticle(const edm4hep::SimTrackerHit& hit);

/// Set the particle for a MutableSimTrackerHit (compatibility with EDM4hep <
/// 0.99 and >= 0.99)
/// @param hit The MutableSimTrackerHit
/// @param particle The MCParticle to set
void setParticle(edm4hep::MutableSimTrackerHit& hit,
                 const edm4hep::MCParticle& particle);

/// Write an Acts track to EDM4hep format
/// @param gctx The geometry context
/// @param track The Acts track to convert
/// @param to The EDM4hep track to write to
/// @param Bz The magnetic field z-component
/// @param logger The logger instance
template <Acts::TrackProxyConcept track_proxy_t>
void writeTrack(const Acts::GeometryContext& gctx, track_proxy_t track,
                edm4hep::MutableTrack to, double Bz,
                const Acts::Logger& logger = Acts::getDummyLogger()) {
  ACTS_VERBOSE("Converting track to EDM4hep");
  to.setChi2(track.chi2());
  to.setNdf(track.nDoF());

  std::vector<edm4hep::TrackState> outTrackStates;
  outTrackStates.reserve(track.nTrackStates());

  auto setParameters = [](edm4hep::TrackState& trackState,
                          const detail::Parameters& params) {
    trackState.D0 = params.values[0];
    trackState.Z0 = params.values[1];
    trackState.phi = params.values[2];
    trackState.tanLambda = params.values[3];
    trackState.omega = params.values[4];
    trackState.time = params.values[5];

    if (params.covariance) {
      detail::packCovariance(params.covariance.value(),
                             trackState.covMatrix.data());
    }
  };

  ACTS_VERBOSE("Converting " << track.nTrackStates() << " track states");

  for (const auto& state : track.trackStatesReversed()) {
    auto typeFlags = state.typeFlags();
    if (!typeFlags.isMeasurement()) {
      continue;
    }

    edm4hep::TrackState& trackState = outTrackStates.emplace_back();
    trackState.location = edm4hep::TrackState::AtOther;

    Acts::BoundTrackParameters params{state.referenceSurface().getSharedPtr(),
                                      state.parameters(), state.covariance(),
                                      track.particleHypothesis()};

    // Convert to LCIO track parametrization expected by EDM4hep
    detail::Parameters converted =
        detail::convertTrackParametersToEdm4hep(gctx, Bz, params);

    // Write the converted parameters to the EDM4hep track state
    setParameters(trackState, converted);
    ACTS_VERBOSE("- parameters: " << state.parameters().transpose() << " -> "
                                  << converted.values.transpose());
    ACTS_VERBOSE("- covariance: \n"
                 << state.covariance() << "\n->\n"
                 << converted.covariance.value());

    // Converted parameters are relative to an ad-hoc perigee surface created at
    // the hit location
    auto center = converted.surface->center(gctx);
    trackState.referencePoint.x = center.x();
    trackState.referencePoint.y = center.y();
    trackState.referencePoint.z = center.z();
    ACTS_VERBOSE("- ref surface ctr: " << center.transpose());
  }
  outTrackStates.front().location = edm4hep::TrackState::AtLastHit;
  outTrackStates.back().location = edm4hep::TrackState::AtFirstHit;

  // Add a track state that represents the IP parameters
  auto& ipState = outTrackStates.emplace_back();

  // Convert the track parameters at the IP
  Acts::BoundTrackParameters trackParams{
      track.referenceSurface().getSharedPtr(), track.parameters(),
      track.covariance(), track.particleHypothesis()};

  // Convert to LCIO track parametrization expected by EDM4hep
  auto converted =
      detail::convertTrackParametersToEdm4hep(gctx, Bz, trackParams);
  setParameters(ipState, converted);
  ipState.location = edm4hep::TrackState::AtIP;
  ACTS_VERBOSE("Writing track level quantities as IP track state");
  ACTS_VERBOSE("- parameters: " << track.parameters().transpose());
  ACTS_VERBOSE("           -> " << converted.values.transpose());
  ACTS_VERBOSE("- covariance: \n"
               << track.covariance() << "\n->\n"
               << converted.covariance.value());

  // Write the converted parameters to the EDM4hep track state
  // The reference point is at the location of the reference surface of the
  // track itself, but if that's not a perigee surface, another ad-hoc perigee
  // at the position will be created.
  auto center = converted.surface->center(gctx);
  ipState.referencePoint.x = center.x();
  ipState.referencePoint.y = center.y();
  ipState.referencePoint.z = center.z();

  ACTS_VERBOSE("- ref surface ctr: " << center.transpose());

  for (auto& trackState : outTrackStates) {
    to.addToTrackStates(trackState);
  }
}

/// Read an EDM4hep track into Acts format
/// @param from The EDM4hep track to read
/// @param track The Acts track proxy to fill
/// @param Bz The magnetic field z-component
/// @param logger The logger instance
template <Acts::TrackProxyConcept track_proxy_t>
void readTrack(const edm4hep::Track& from, track_proxy_t& track, double Bz,
               const Acts::Logger& logger = Acts::getDummyLogger()) {
  using namespace Acts;
  ACTS_VERBOSE("Reading track from EDM4hep");
  TrackStatePropMask mask = TrackStatePropMask::Smoothed;

  std::optional<edm4hep::TrackState> ipState;

  auto unpack =
      [](const edm4hep::TrackState& trackState) -> detail::Parameters {
    detail::Parameters params;
    params.covariance = BoundMatrix::Zero();
    params.values = BoundVector::Zero();
    detail::unpackCovariance(trackState.covMatrix.data(),
                             params.covariance.value());
    params.values[0] = trackState.D0;
    params.values[1] = trackState.Z0;
    params.values[2] = trackState.phi;
    params.values[3] = trackState.tanLambda;
    params.values[4] = trackState.omega;
    params.values[5] = trackState.time;

    Vector3 center = {
        trackState.referencePoint.x,
        trackState.referencePoint.y,
        trackState.referencePoint.z,
    };
    params.surface = Acts::Surface::makeShared<PerigeeSurface>(center);

    return params;
  };

  ACTS_VERBOSE("Reading " << from.trackStates_size()
                          << " track states (including IP state)");
  // We write the trackstates out outside in, need to reverse iterate to get the
  // same order
  for (std::size_t i = from.trackStates_size() - 1;
       i <= from.trackStates_size(); i--) {
    auto trackState = from.getTrackStates(i);
    if (trackState.location == edm4hep::TrackState::AtIP) {
      ipState = trackState;
      continue;
    }

    auto params = unpack(trackState);

    auto ts = track.appendTrackState(mask);
    ts.typeFlags().setIsMeasurement();

    auto converted = detail::convertTrackParametersFromEdm4hep(Bz, params);

    ts.smoothed() = converted.parameters();
    ts.smoothedCovariance() =
        converted.covariance().value_or(BoundMatrix::Zero());
    ts.setReferenceSurface(params.surface);
  }

  if (!ipState.has_value()) {
    ACTS_ERROR("Did not find IP state in edm4hep input");
    throw std::runtime_error{"Did not find IP state in edm4hep input"};
  }

  detail::Parameters params = unpack(ipState.value());

  auto converted = detail::convertTrackParametersFromEdm4hep(Bz, params);

  ACTS_VERBOSE("IP state parameters: " << converted.parameters().transpose());
  ACTS_VERBOSE("-> covariance:\n"
               << converted.covariance().value_or(BoundMatrix::Zero()));

  track.parameters() = converted.parameters();
  track.covariance() = converted.covariance().value_or(BoundMatrix::Zero());
  track.setReferenceSurface(params.surface);

  track.chi2() = from.getChi2();
  track.nDoF() = from.getNdf();
  track.nMeasurements() = track.nTrackStates();
}

/// @brief Helper class for associating simulation hits between EDM4hep and internal indices
class SimHitAssociation {
 public:
  /// Reserve space for associations
  /// @param size Number of associations to reserve
  void reserve(std::size_t size);

  /// Get number of associations
  /// @return Number of associations
  std::size_t size() const;

  /// Add association between internal index and EDM4hep hit
  /// @param internalIndex Internal hit index
  /// @param edm4hepHit EDM4hep hit object
  void add(std::size_t internalIndex, const edm4hep::SimTrackerHit& edm4hepHit);

  /// Look up EDM4hep hit by internal index
  /// @param internalIndex Internal hit index
  /// @return EDM4hep hit object
  [[nodiscard]]
  edm4hep::SimTrackerHit lookup(std::size_t internalIndex) const;

  /// Look up internal index by EDM4hep hit
  /// @param hit EDM4hep hit object
  /// @return Internal hit index
  std::size_t lookup(const edm4hep::SimTrackerHit& hit) const;

 private:
  std::vector<edm4hep::SimTrackerHit> m_internalToEdm4hep;
  std::unordered_map<podio::ObjectID, std::size_t> m_edm4hepToInternal;
};

/// @}

}  // namespace ActsPlugins::EDM4hepUtil
