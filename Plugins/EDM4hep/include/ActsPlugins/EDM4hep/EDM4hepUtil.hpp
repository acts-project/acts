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
#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/TrackProxyConcept.hpp"
#include "Acts/EventData/TrackStatePropMask.hpp"
#include "Acts/EventData/Types.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "ActsPodioEdm/MutableTrackerHitLocal.h"
#include "ActsPodioEdm/TrackerHitLocal.h"

#include <algorithm>
#include <functional>
#include <optional>
#include <span>
#include <stdexcept>
#include <vector>

#include <edm4hep/MCParticle.h>
#include <edm4hep/MutableSimTrackerHit.h>
#include <edm4hep/MutableTrack.h>
#include <edm4hep/SimTrackerHit.h>
#include <edm4hep/Track.h>
#include <edm4hep/TrackState.h>
#include <edm4hep/TrackerHit.h>
#include <edm4hep/Vector4f.h>
#include <edm4hep/Vertex.h>
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

#include "ActsPodioEdm/TrackerHitLocalCollection.h"
#include "ActsPodioEdm/TrackerHitLocalSimTrackerHitLinkCollection.h"

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

/// Default tracker-hit lookup used by @ref writeTrack: never associates a hit,
/// reproducing the historical behaviour where no tracker hits are written.
struct NoTrackerHitLookup {
  template <typename state_proxy_t>
  std::optional<edm4hep::TrackerHit> operator()(const state_proxy_t& /*state*/
  ) const {
    return std::nullopt;
  }
};

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

/// Callable returning the local magnetic field z-component (in Acts native
/// units) at a given global position.
///
/// This is used by @ref writeTrack to support spatially varying magnetic
/// fields: the LCIO/EDM4hep perigee parametrization (in particular the
/// curvature @c omega) depends on the local field, so for a non-uniform field
/// the conversion should use the field value at each track state's location
/// rather than a single global constant.
using LocalBzProvider = std::function<double(const Acts::Vector3& position)>;

/// Write an Acts track to EDM4hep format, using a spatially varying field.
///
/// This is the general form of @ref writeTrack. In addition to the track
/// summary quantities (chi2, ndf, number of holes) it writes one EDM4hep track
/// state per measurement plus a dedicated @c AtIP state. The perigee
/// conversion of each state evaluates the local field via @p bzAtPosition at
/// the global position of that state, which makes this a drop-in replacement
/// for bespoke ACTS->EDM4hep converters that rely on a
/// @c Acts::MagneticFieldProvider (e.g. k4ActsTracking's @c ACTS2edm4hep_track).
///
/// @note Resolving tracker hits requires application-specific
///       source-link/hit-container knowledge, so it is delegated to the
///       optional @p hitLookup callback. It is invoked once per measurement
///       track state (in the same order the states are written) with the Acts
///       track state proxy, and should return the associated
///       @c edm4hep::TrackerHit (any of the interfaced hit types) or
///       @c std::nullopt. Returned hits are attached via @c addToTrackerHits.
///
/// @param gctx The geometry context
/// @param track The Acts track to convert
/// @param to The EDM4hep track to write to
/// @param bzAtPosition Callable returning the local Bz at a global position
/// @param logger The logger instance
/// @param hitLookup Optional callback mapping a measurement track state to its
///                  EDM4hep tracker hit (defaults to writing no tracker hits)
template <Acts::TrackProxyConcept track_proxy_t,
          typename hit_lookup_t = detail::NoTrackerHitLookup>
void writeTrack(const Acts::GeometryContext& gctx, track_proxy_t track,
                edm4hep::MutableTrack to, const LocalBzProvider& bzAtPosition,
                const Acts::Logger& logger = Acts::getDummyLogger(),
                const hit_lookup_t& hitLookup = {}) {
  ACTS_VERBOSE("Converting track to EDM4hep");
  to.setChi2(track.chi2());
  to.setNdf(track.nDoF());
  to.setNholes(static_cast<std::int32_t>(track.nHoles()));

  std::vector<edm4hep::TrackState> outTrackStates;
  outTrackStates.reserve(track.nTrackStates());
  // Tracker hits resolved for each measurement state, kept in parallel with
  // outTrackStates so both can be emitted in the final order below.
  std::vector<std::optional<edm4hep::TrackerHit>> outHits;
  outHits.reserve(track.nTrackStates());

  auto setParameters = [](edm4hep::TrackState& trackState,
                          const detail::Parameters& params) {
    trackState.D0 = static_cast<float>(params.values[0]);
    trackState.Z0 = static_cast<float>(params.values[1]);
    trackState.phi = static_cast<float>(params.values[2]);
    trackState.tanLambda = static_cast<float>(params.values[3]);
    trackState.omega = static_cast<float>(params.values[4]);
    trackState.time = static_cast<float>(params.values[5]);

    if (params.covariance) {
      detail::packCovariance(params.covariance.value(),
                             trackState.covMatrix.data());
    }
  };

  ACTS_VERBOSE("Converting " << track.nTrackStates() << " track states");

  for (const auto& state : track.trackStatesReversed()) {
    if (!state.typeFlags().isMeasurement()) {
      continue;
    }

    // Resolve the associated tracker hit for this measurement (if any). Stored
    // in parallel with the track state and attached below in output order.
    outHits.push_back(hitLookup(state));

    edm4hep::TrackState& trackState = outTrackStates.emplace_back();
    trackState.location = edm4hep::TrackState::AtOther;

    Acts::BoundTrackParameters params{state.referenceSurface().getSharedPtr(),
                                      state.parameters(), state.covariance(),
                                      track.particleHypothesis()};

    // Evaluate the local field at the global position of this state
    double Bz = bzAtPosition(params.position(gctx));

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
    trackState.referencePoint.x = static_cast<float>(center.x());
    trackState.referencePoint.y = static_cast<float>(center.y());
    trackState.referencePoint.z = static_cast<float>(center.z());
    ACTS_VERBOSE("- ref surface ctr: " << center.transpose());
  }
  // At this point the measurement states are ordered outside-in (last hit
  // first, first hit last), following Acts' reverse track state iteration.
  if (!outTrackStates.empty()) {
    outTrackStates.front().location = edm4hep::TrackState::AtLastHit;
    outTrackStates.back().location = edm4hep::TrackState::AtFirstHit;
  }

  // Track state that represents the IP parameters
  edm4hep::TrackState ipState;

  // Convert the track parameters at the IP
  Acts::BoundTrackParameters trackParams{
      track.referenceSurface().getSharedPtr(), track.parameters(),
      track.covariance(), track.particleHypothesis()};

  // Evaluate the local field at the track reference (perigee) position
  double ipBz = bzAtPosition(trackParams.position(gctx));

  // Convert to LCIO track parametrization expected by EDM4hep
  auto converted =
      detail::convertTrackParametersToEdm4hep(gctx, ipBz, trackParams);
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
  ipState.referencePoint.x = static_cast<float>(center.x());
  ipState.referencePoint.y = static_cast<float>(center.y());
  ipState.referencePoint.z = static_cast<float>(center.z());

  ACTS_VERBOSE("- ref surface ctr: " << center.transpose());

  // Emit track states following the EDM4hep/LCIO convention used by
  // k4ActsTracking: the IP state first, then the measurement states ordered
  // inside-out (first hit ... last hit). Reverse the outside-in buffers so the
  // measurement states, and their tracker hits, come out in that order.
  std::reverse(outTrackStates.begin(), outTrackStates.end());
  std::reverse(outHits.begin(), outHits.end());

  to.addToTrackStates(ipState);
  for (std::size_t i = 0; i < outTrackStates.size(); ++i) {
    to.addToTrackStates(outTrackStates[i]);
    if (outHits[i].has_value()) {
      to.addToTrackerHits(outHits[i].value());
    }
  }
}

/// Write an Acts track to EDM4hep format, using a uniform magnetic field.
///
/// Convenience overload of @ref writeTrack for the common case of a constant
/// solenoidal field. Delegates to the @ref LocalBzProvider form with a constant
/// field lookup.
///
/// @param gctx The geometry context
/// @param track The Acts track to convert
/// @param to The EDM4hep track to write to
/// @param Bz The (uniform) magnetic field z-component in Acts native units
/// @param logger The logger instance
/// @param hitLookup Optional callback mapping a measurement track state to its
///                  EDM4hep tracker hit (defaults to writing no tracker hits)
template <Acts::TrackProxyConcept track_proxy_t,
          typename hit_lookup_t = detail::NoTrackerHitLookup>
void writeTrack(const Acts::GeometryContext& gctx, track_proxy_t track,
                edm4hep::MutableTrack to, double Bz,
                const Acts::Logger& logger = Acts::getDummyLogger(),
                const hit_lookup_t& hitLookup = {}) {
  writeTrack(
      gctx, std::move(track), to,
      LocalBzProvider{[Bz](const Acts::Vector3& /*position*/) { return Bz; }},
      logger, hitLookup);
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

  auto unpack = [](const edm4hep::TrackState& trackState) {
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
  // Track states are written IP-first, followed by the measurement states
  // ordered inside-out (first hit ... last hit). Iterate forward and skip the
  // IP state so the reconstructed Acts track state order matches the input that
  // was passed to writeTrack.
  for (std::size_t i = 0; i < from.trackStates_size(); ++i) {
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

namespace detail {
/// Support for both EDM4hep versions where the vertex position is a 3 or 4
/// vector
constexpr bool kEdm4hepVertexHasTime =
    std::is_same_v<edm4hep::Vector4f,
                   decltype(std::declval<edm4hep::Vertex>().getPosition())> &&
    std::is_same_v<edm4hep::CovMatrix4f,
                   decltype(std::declval<edm4hep::Vertex>().getCovMatrix())>;

}  // namespace detail

void writeVertex(const Acts::Vertex& vertex, edm4hep::MutableVertex to);

/// Write a measurement to an EDM4hep tracker hit
///
/// This function converts an ACTS measurement into the EDM4hep format. It
/// handles:
/// - Position conversion from local to global coordinates (in mm)
/// - Time storage (in ns)
/// - Measurement values and covariance matrix storage
/// - Encoding of measurement indices into a 32-bit integer:
///   - First 4 bits: number of indices (max
///     `ActsPodioEdm::detail::kMaxSubspaceSize`)
///   - Next 4 bits per index: measured bound parameter index (max
///     `ActsPodioEdm::detail::kMaxSubspaceIndex`)
///
/// The function will throw if:
/// - The number of indices exceeds `ActsPodioEdm::detail::kMaxSubspaceSize`
/// - Any index is larger than `ActsPodioEdm::detail::kMaxSubspaceIndex`
/// - There's a size mismatch between parameters and covariance matrix
///
/// @param gctx The geometry context
/// @param parameters The parameters of the measurement
/// @param covariance The covariance of the measurement
/// @param indices The indices of the measurement
/// @param cellId The cell ID of the measurement
/// @param surface The surface of the measurement
/// @param to The EDM4hep tracker hit to write to
void writeMeasurement(const Acts::GeometryContext& gctx,
                      const Eigen::Map<const Acts::DynamicVector>& parameters,
                      const Eigen::Map<const Acts::DynamicMatrix>& covariance,
                      std::span<const std::uint8_t> indices,
                      std::uint64_t cellId, const Acts::Surface& surface,
                      ActsPodioEdm::MutableTrackerHitLocal& to);

/// Data extracted when reading a measurement from EDM4hep
struct MeasurementData {
  /// Measurement parameters (local coordinates, full bound space)
  Acts::BoundVector parameters{Acts::BoundVector::Zero()};
  /// Covariance matrix of the measurement (full bound space)
  Acts::BoundMatrix covariance{Acts::BoundMatrix::Zero()};
  /// Indices of the measured parameters (subspace)
  std::vector<Acts::SubspaceIndex> indices;
  /// Cell ID of the measurement
  std::uint64_t cellId{0};
};

/// Read a measurement from an EDM4hep tracker hit
///
/// This function extracts measurement parameters, covariance, and indices from
/// an EDM4hep TrackerHitLocal. It is the inverse of writeMeasurement.
///
/// @param from The EDM4hep tracker hit to read from
/// @return The extracted measurement data (parameters, covariance, indices,
///         cellId)
MeasurementData readMeasurement(const ActsPodioEdm::TrackerHitLocal& from);

/// Callback type for sim hit lookup during TrackerHitLocal link writing.
/// Given a hit index (position in the TrackerHitLocalCollection), returns the
/// associated edm4hep SimTrackerHit, or std::nullopt if no association exists.
using SimHitForHitIndex =
    std::function<std::optional<edm4hep::SimTrackerHit>(std::size_t hitIndex)>;

/// Write sim-hit links for a TrackerHitLocal collection.
///
/// For each hit in @p hits, calls @p lookup with its position index. If the
/// callback returns a hit, a link entry is created in @p links. Hits with no
/// association are silently skipped, so the resulting link collection may be
/// sparse.
///
/// @param hits    The TrackerHitLocal collection to link from
/// @param links   The link collection to populate
/// @param lookup  Callback mapping hit index → optional SimTrackerHit
void writeTrackerHitSimHitLinks(
    const ActsPodioEdm::TrackerHitLocalCollection& hits,
    ActsPodioEdm::TrackerHitLocalSimTrackerHitLinkCollection& links,
    const SimHitForHitIndex& lookup);

}  // namespace ActsPlugins::EDM4hepUtil
