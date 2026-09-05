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
#include "Acts/EventData/AnyTrackStateProxy.hpp"
#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/TrackProxyConcept.hpp"
#include "Acts/EventData/TrackStatePropMask.hpp"
#include "Acts/EventData/Types.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
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

#include <edm4hep/Constants.h>
#include <edm4hep/CovMatrix6f.h>
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

/// Unpack an EDM4hep track state covariance into the parameter order used by
/// @c Parameters::values (d0, z0, phi, tanLambda, omega, time).
void unpackCovariance(const edm4hep::CovMatrix6f& from,
                      Acts::SquareMatrix<6>& to);
/// Pack a covariance in @c Parameters::values order into an EDM4hep track
/// state covariance. The element placement is delegated to
/// @c edm4hep::CovMatrix6f::setValue, so the on-disk layout always follows the
/// EDM4hep parameter order (@c edm4hep::TrackParams), which is a different
/// permutation from the internal one.
void packCovariance(const Acts::SquareMatrix<6>& from,
                    edm4hep::CovMatrix6f& to);

/// Look up the local magnetic field z-component at @p position.
/// @throws std::runtime_error if the field lookup fails
double localBz(const Acts::MagneticFieldProvider& magneticField,
               Acts::MagneticFieldProvider::Cache& fieldCache,
               const Acts::Vector3& position, const Acts::Logger& logger);

/// Global position an EDM4hep track state is expressed relative to. This is
/// also the position at which the local magnetic field is evaluated.
Acts::Vector3 referencePoint(const edm4hep::TrackState& trackState);

/// Unpack an EDM4hep track state into the internal parameter representation,
/// including an ad-hoc perigee surface placed at its reference point.
Parameters unpackTrackState(const edm4hep::TrackState& trackState);

/// Global position of the point of closest approach described by an unpacked
/// EDM4hep track state, i.e. where the particle actually is, as opposed to the
/// reference point the parameters are expressed relative to. The two differ by
/// d0/z0. This is the position @ref writeTrackState evaluates the local
/// magnetic field at, so the read side has to use it as well.
Acts::Vector3 pointOfClosestApproach(const Acts::GeometryContext& gctx,
                                     const Parameters& params);

Parameters convertTrackParametersToEdm4hep(
    const Acts::GeometryContext& gctx, double Bz,
    const Acts::BoundTrackParameters& params);

Acts::BoundTrackParameters convertTrackParametersFromEdm4hep(
    double Bz, const Parameters& params);

}  // namespace detail

/// @addtogroup edm4hep_plugin
/// @{

/// Get the particle from a SimTrackerHit
/// @param hit The SimTrackerHit
/// @return The associated MCParticle
edm4hep::MCParticle getParticle(const edm4hep::SimTrackerHit& hit);

/// Set the particle for a MutableSimTrackerHit
/// @param hit The MutableSimTrackerHit
/// @param particle The MCParticle to set
void setParticle(edm4hep::MutableSimTrackerHit& hit,
                 const edm4hep::MCParticle& particle);

/// Convert a single set of bound track parameters into an EDM4hep track state.
///
/// The EDM4hep track state uses a perigee parametrization
/// (D0, Z0, phi, tanLambda, omega, time) defined relative to a reference point.
/// If @p params are not already expressed on a perigee surface they are
/// re-expressed at an ad-hoc perigee created at their global position,
/// transporting both parameters and covariance, and the state's
/// @c referencePoint is set to that perigee centre. The curvature @c omega
/// depends on the local field, which is evaluated at the position of @p params.
///
/// This is the per-state building block used by @ref writeTrack. It is exposed
/// so that track states which do not come from a fitted track - a seed, or an
/// extrapolation to the calorimeter face - can be produced with exactly the
/// same conversion instead of a private reimplementation.
///
/// @param gctx The geometry context
/// @param location The EDM4hep track state location, e.g.
///                 @c edm4hep::TrackState::AtIP
/// @param params The bound track parameters to convert
/// @param magneticField The magnetic field provider
/// @param fieldCache Field lookup cache, reused across calls
/// @param logger The logger instance
/// @return The converted EDM4hep track state
edm4hep::TrackState writeTrackState(
    const Acts::GeometryContext& gctx, std::int32_t location,
    const Acts::BoundTrackParameters& params,
    const Acts::MagneticFieldProvider& magneticField,
    Acts::MagneticFieldProvider::Cache& fieldCache,
    const Acts::Logger& logger = Acts::getDummyLogger());

/// Convert a single set of bound track parameters into an EDM4hep track state.
///
/// Convenience overload of @ref writeTrackState that creates its own field
/// cache from @p mctx. Prefer the cache-taking overload when converting many
/// states in a row.
///
/// @param gctx The geometry context
/// @param mctx The magnetic field context
/// @param location The EDM4hep track state location, e.g.
///                 @c edm4hep::TrackState::AtIP
/// @param params The bound track parameters to convert
/// @param magneticField The magnetic field provider
/// @param logger The logger instance
/// @return The converted EDM4hep track state
edm4hep::TrackState writeTrackState(
    const Acts::GeometryContext& gctx, const Acts::MagneticFieldContext& mctx,
    std::int32_t location, const Acts::BoundTrackParameters& params,
    const Acts::MagneticFieldProvider& magneticField,
    const Acts::Logger& logger = Acts::getDummyLogger());

/// Callback resolving the EDM4hep tracker hit associated with a measurement
/// track state.
///
/// It is invoked by @ref writeTrack once per measurement track state and
/// should return the associated @c edm4hep::TrackerHit (any of the interfaced
/// hit types) or @c std::nullopt if there is none. The track state is passed
/// type-erased, so this is a plain callable independent of the trajectory
/// backend. An empty callback means no tracker hits are written.
using TrackerHitLookup = std::function<std::optional<edm4hep::TrackerHit>(
    const Acts::AnyConstTrackStateProxy& state)>;

/// Write an Acts track to EDM4hep format, using a magnetic field provider.
///
/// This is the general form of @ref writeTrack. In addition to the track
/// summary quantities (chi2, ndf, number of holes) it writes one EDM4hep track
/// state per measurement plus a dedicated @c AtIP state. The perigee
/// conversion of each state evaluates the local field via @p magneticField at
/// the global position of that state, which supports spatially varying fields
///
/// @note Resolving tracker hits requires application-specific
///       source-link/hit-container knowledge, so it is delegated to the
///       optional @p hitLookup callback. It is invoked exactly once per
///       measurement track state, in the order the states are emitted, i.e.
///       inside-out from @c AtFirstHit to @c AtLastHit, so a stateful callback
///       sees the same order as the output. Returned hits are attached via
///       @c addToTrackerHits.
///
/// @param gctx The geometry context
/// @param mctx The magnetic field context
/// @param track The Acts track to convert
/// @param to The EDM4hep track to write to
/// @param magneticField The magnetic field provider, evaluated at each state
/// @param hitLookup Optional callback mapping a measurement track state to its
///                  EDM4hep tracker hit (defaults to writing no tracker hits)
/// @param logger The logger instance
template <Acts::TrackProxyConcept track_proxy_t>
void writeTrack(const Acts::GeometryContext& gctx,
                const Acts::MagneticFieldContext& mctx, track_proxy_t track,
                edm4hep::MutableTrack to,
                const Acts::MagneticFieldProvider& magneticField,
                const TrackerHitLookup& hitLookup = {},
                const Acts::Logger& logger = Acts::getDummyLogger()) {
  ACTS_VERBOSE("Converting track to EDM4hep");
  to.setChi2(track.chi2());
  to.setNdf(track.nDoF());
  to.setNholes(static_cast<std::int32_t>(track.nHoles()));

  auto fieldCache = magneticField.makeCache(mctx);

  std::vector<edm4hep::TrackState> outTrackStates;
  outTrackStates.reserve(track.nTrackStates());
  // Source track states kept in parallel with outTrackStates, so that the hit
  // lookup below can be run in the same order the states are emitted.
  std::vector<Acts::AnyConstTrackStateProxy> outStates;
  outStates.reserve(track.nTrackStates());

  ACTS_VERBOSE("Converting " << track.nTrackStates() << " track states");

  for (auto state : track.trackStatesReversed()) {
    if (!state.typeFlags().isMeasurement()) {
      continue;
    }

    outStates.push_back(Acts::AnyConstTrackStateProxy{state});

    Acts::BoundTrackParameters params{state.referenceSurface().getSharedPtr(),
                                      state.parameters(), state.covariance(),
                                      track.particleHypothesis()};

    outTrackStates.push_back(writeTrackState(gctx, edm4hep::TrackState::AtOther,
                                             params, magneticField, fieldCache,
                                             logger));
  }
  // At this point the measurement states are ordered outside-in (last hit
  // first, first hit last), following Acts' reverse track state iteration.
  if (!outTrackStates.empty()) {
    outTrackStates.front().location = edm4hep::TrackState::AtLastHit;
    outTrackStates.back().location = edm4hep::TrackState::AtFirstHit;
  }

  // Track state that represents the IP parameters. The reference point ends up
  // at the track's own reference surface, or at an ad-hoc perigee created at
  // its position if that surface is not already a perigee.
  Acts::BoundTrackParameters trackParams{
      track.referenceSurface().getSharedPtr(), track.parameters(),
      track.covariance(), track.particleHypothesis()};

  ACTS_VERBOSE("Writing track level quantities as IP track state");
  edm4hep::TrackState ipState =
      writeTrackState(gctx, edm4hep::TrackState::AtIP, trackParams,
                      magneticField, fieldCache, logger);

  // Emit track states following the EDM4hep/LCIO convention used by
  // k4ActsTracking: the IP state first, then the measurement states ordered
  // inside-out (first hit ... last hit). Reverse the outside-in buffers so the
  // measurement states, and their tracker hits, come out in that order.
  std::reverse(outTrackStates.begin(), outTrackStates.end());
  std::reverse(outStates.begin(), outStates.end());

  to.addToTrackStates(ipState);
  for (std::size_t i = 0; i < outTrackStates.size(); ++i) {
    to.addToTrackStates(outTrackStates[i]);
    // Resolve the tracker hit here rather than in the loop above, so that the
    // callback is invoked in the same order the states are emitted.
    if (hitLookup) {
      if (std::optional<edm4hep::TrackerHit> hit = hitLookup(outStates[i]);
          hit.has_value()) {
        to.addToTrackerHits(hit.value());
      }
    }
  }
}

/// Write an Acts track to EDM4hep format, using a uniform magnetic field.
///
/// Convenience overload of @ref writeTrack for the common case of a constant
/// solenoidal field. Delegates to the @c Acts::MagneticFieldProvider form with
/// an @c Acts::ConstantBField.
///
/// @param gctx The geometry context
/// @param track The Acts track to convert
/// @param to The EDM4hep track to write to
/// @param Bz The (uniform) magnetic field z-component in Acts native units
/// @param hitLookup Optional callback mapping a measurement track state to its
///                  EDM4hep tracker hit (defaults to writing no tracker hits)
/// @param logger The logger instance
// NOLINTBEGIN(performance-unnecessary-value-param)
// `to` is a podio handle and is passed by value throughout this interface, to
// match the signature of the overload it forwards to. clang-tidy only sees it
// forwarded here, not mutated.
template <Acts::TrackProxyConcept track_proxy_t>
void writeTrack(const Acts::GeometryContext& gctx, track_proxy_t track,
                edm4hep::MutableTrack to, double Bz,
                const TrackerHitLookup& hitLookup = {},
                const Acts::Logger& logger = Acts::getDummyLogger()) {
  Acts::ConstantBField magneticField{Acts::Vector3{0, 0, Bz}};
  Acts::MagneticFieldContext mctx{};
  writeTrack(gctx, mctx, std::move(track), to, magneticField, hitLookup,
             logger);
}
// NOLINTEND(performance-unnecessary-value-param)

/// Read an EDM4hep track into Acts format, using a magnetic field provider.
///
/// This is the general form of @ref readTrack, and the inverse of the
/// corresponding @ref writeTrack overload. The perigee conversion of each track
/// state evaluates the local field via @p magneticField at that state's point
/// of closest approach, which supports spatially varying fields. That is the
/// same position @ref writeTrackState samples the field at, so the two
/// round-trip consistently. Note this is not in general the state's
/// @c referencePoint: for parameters expressed on a perigee the two differ by
/// d0/z0.
///
/// @param gctx The geometry context
/// @param mctx The magnetic field context
/// @param from The EDM4hep track to read
/// @param track The Acts track proxy to fill
/// @param magneticField The magnetic field provider, evaluated at each state
/// @param logger The logger instance
template <Acts::TrackProxyConcept track_proxy_t>
void readTrack(const Acts::GeometryContext& gctx,
               const Acts::MagneticFieldContext& mctx,
               const edm4hep::Track& from, track_proxy_t& track,
               const Acts::MagneticFieldProvider& magneticField,
               const Acts::Logger& logger = Acts::getDummyLogger()) {
  using namespace Acts;
  ACTS_VERBOSE("Reading track from EDM4hep");
  TrackStatePropMask mask = TrackStatePropMask::Smoothed;

  auto fieldCache = magneticField.makeCache(mctx);
  auto bzAtPosition = [&](const Acts::Vector3& position) {
    return detail::localBz(magneticField, fieldCache, position, logger);
  };

  std::optional<edm4hep::TrackState> ipState;

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

    auto params = detail::unpackTrackState(trackState);

    auto ts = track.appendTrackState(mask);
    ts.typeFlags().setIsMeasurement();

    // Evaluate the local field where the particle actually is, matching what
    // writeTrackState does
    auto converted = detail::convertTrackParametersFromEdm4hep(
        bzAtPosition(detail::pointOfClosestApproach(gctx, params)), params);

    ts.smoothed() = converted.parameters();
    ts.smoothedCovariance() =
        converted.covariance().value_or(BoundMatrix::Zero());
    ts.setReferenceSurface(params.surface);
  }

  if (!ipState.has_value()) {
    ACTS_ERROR("Did not find IP state in edm4hep input");
    throw std::runtime_error{"Did not find IP state in edm4hep input"};
  }

  detail::Parameters params = detail::unpackTrackState(ipState.value());

  // Evaluate the local field where the particle actually is, matching what
  // writeTrackState does
  auto converted = detail::convertTrackParametersFromEdm4hep(
      bzAtPosition(detail::pointOfClosestApproach(gctx, params)), params);

  ACTS_VERBOSE("IP state parameters: " << converted.parameters().transpose());
  ACTS_VERBOSE("-> covariance:\n"
               << converted.covariance().value_or(BoundMatrix::Zero()));

  track.parameters() = converted.parameters();
  track.covariance() = converted.covariance().value_or(BoundMatrix::Zero());
  track.setReferenceSurface(params.surface);

  track.chi2() = from.getChi2();
  track.nDoF() = from.getNdf();
  track.nHoles() = static_cast<unsigned int>(from.getNholes());
  track.nMeasurements() = track.nTrackStates();
}

/// Read an EDM4hep track into Acts format, using a uniform magnetic field.
///
/// Convenience overload of @ref readTrack for the common case of a constant
/// solenoidal field. Delegates to the @c Acts::MagneticFieldProvider form with
/// an @c Acts::ConstantBField.
///
/// @param gctx The geometry context
/// @param from The EDM4hep track to read
/// @param track The Acts track proxy to fill
/// @param Bz The (uniform) magnetic field z-component in Acts native units
/// @param logger The logger instance
template <Acts::TrackProxyConcept track_proxy_t>
void readTrack(const Acts::GeometryContext& gctx, const edm4hep::Track& from,
               track_proxy_t& track, double Bz,
               const Acts::Logger& logger = Acts::getDummyLogger()) {
  Acts::ConstantBField magneticField{Acts::Vector3{0, 0, Bz}};
  Acts::MagneticFieldContext mctx{};
  readTrack(gctx, mctx, from, track, magneticField, logger);
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
