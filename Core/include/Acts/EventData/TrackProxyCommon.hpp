// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Common.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/TrackParameterHelpers.hpp"
#include "Acts/EventData/Types.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/HashedString.hpp"
#include "Acts/Utilities/UnitVectors.hpp"

namespace Acts {

namespace detail_tp {
inline constexpr HashedString kTipIndexKey = hashString("tipIndex");
inline constexpr HashedString kStemIndexKey = hashString("stemIndex");
inline constexpr HashedString kMeasurementsKey = hashString("nMeasurements");
inline constexpr HashedString kHolesKey = hashString("nHoles");
inline constexpr HashedString kOutliersKey = hashString("nOutliers");
inline constexpr HashedString kSharedHitsKey = hashString("nSharedHits");
inline constexpr HashedString kSplitHitsKey = hashString("nSplitHits");
inline constexpr HashedString kChi2Key = hashString("chi2");
inline constexpr HashedString kNdfKey = hashString("ndf");
inline constexpr HashedString kNextKey = hashString("next");
}  // namespace detail_tp

/// Common CRTP implementation shared by the various track proxy front-ends.
/// The derived proxy only needs to expose the `component` helpers, while this
/// base class wires up the commonly used convenience methods.
///
/// @tparam Derived The proxy implementation inheriting from this base
/// @tparam index_t The index type used by the proxy
/// @tparam read_only Whether the proxy provides mutable access
template <typename Derived, typename index_t, bool read_only>
class TrackProxyCommon {
 protected:
  /// Cast to derived type
  /// @return Reference to the derived proxy
  constexpr Derived& derived() { return static_cast<Derived&>(*this); }
  /// Cast to derived type
  /// @return Const reference to the derived proxy
  constexpr const Derived& derived() const {
    return static_cast<const Derived&>(*this);
  }

 public:
  /// Index type used for referencing track states.
  using IndexType = index_t;

  /// Get the tip index, i.e. the entry point into the track state container.
  /// @return The tip index
  IndexType tipIndex() const {
    return derived().template component<IndexType, detail_tp::kTipIndexKey>();
  }

  /// Get a mutable reference to the tip index.
  /// @return Mutable reference to the tip index
  IndexType& tipIndex()
    requires(!read_only)
  {
    return derived().template component<IndexType, detail_tp::kTipIndexKey>();
  }

  /// Get the stem index, i.e. the innermost track state in the track.
  /// @return The stem index
  IndexType stemIndex() const {
    return derived().template component<IndexType, detail_tp::kStemIndexKey>();
  }

  /// Get a mutable reference to the stem index.
  /// @return Mutable reference to the stem index
  IndexType& stemIndex()
    requires(!read_only)
  {
    return derived().template component<IndexType, detail_tp::kStemIndexKey>();
  }

  /// Return whether the track is forward-linked, i.e. whether each track state
  /// has a valid next-state index. Forward-linking is required to use the
  /// inside-out @c trackStates() range.
  /// @return True if the track is forward-linked
  bool isForwardLinked() const { return stemIndex() != kTrackIndexInvalid; }

  /// Return the number of measurements assigned to this track.
  /// @return The number of measurements
  unsigned int nMeasurements() const {
    return derived()
        .template component<unsigned int, detail_tp::kMeasurementsKey>();
  }

  /// Return a mutable reference to the number of measurements.
  /// @return Mutable reference to the number of measurements
  unsigned int& nMeasurements()
    requires(!read_only)
  {
    return derived()
        .template component<unsigned int, detail_tp::kMeasurementsKey>();
  }

  /// Return the number of holes on this track.
  /// @return The number of holes
  unsigned int nHoles() const {
    return derived().template component<unsigned int, detail_tp::kHolesKey>();
  }

  /// Return a mutable reference to the number of holes.
  /// @return Mutable reference to the number of holes
  unsigned int& nHoles()
    requires(!read_only)
  {
    return derived().template component<unsigned int, detail_tp::kHolesKey>();
  }

  /// Return the number of outliers for this track.
  /// @return The number of outliers
  unsigned int nOutliers() const {
    return derived()
        .template component<unsigned int, detail_tp::kOutliersKey>();
  }

  /// Return a mutable reference to the number of outliers.
  /// @return Mutable reference to the number of outliers
  unsigned int& nOutliers()
    requires(!read_only)
  {
    return derived()
        .template component<unsigned int, detail_tp::kOutliersKey>();
  }

  /// Return the number of shared hits for this track.
  /// @return The number of shared hits
  unsigned int nSharedHits() const {
    return derived()
        .template component<unsigned int, detail_tp::kSharedHitsKey>();
  }

  /// Return a mutable reference to the number of shared hits.
  /// @return Mutable reference to the number of shared hits
  unsigned int& nSharedHits()
    requires(!read_only)
  {
    return derived()
        .template component<unsigned int, detail_tp::kSharedHitsKey>();
  }

  /// Return the number of split hits for this track.
  /// @return The number of split hits
  unsigned int nSplitHits() const {
    return derived()
        .template component<unsigned int, detail_tp::kSplitHitsKey>();
  }

  /// Return a mutable reference to the number of split hits.
  /// @return Mutable reference to the number of split hits
  unsigned int& nSplitHits()
    requires(!read_only)
  {
    return derived()
        .template component<unsigned int, detail_tp::kSplitHitsKey>();
  }

  /// Return the local chi-squared contribution.
  /// @return The chi-squared value
  float chi2() const {
    return derived().template component<float, detail_tp::kChi2Key>();
  }

  /// Return a mutable reference to the local chi-squared contribution.
  /// @return Mutable reference to the chi-squared value
  float& chi2()
    requires(!read_only)
  {
    return derived().template component<float, detail_tp::kChi2Key>();
  }

  /// Return the number of degrees of freedom.
  /// @return The number of degrees of freedom
  unsigned int nDoF() const {
    return derived().template component<unsigned int, detail_tp::kNdfKey>();
  }

  /// Return a mutable reference to the number of degrees of freedom.
  /// @return Mutable reference to the number of degrees of freedom
  unsigned int& nDoF()
    requires(!read_only)
  {
    return derived().template component<unsigned int, detail_tp::kNdfKey>();
  }

  /// Access the theta parameter of the track at the reference surface
  /// @return The theta parameter
  [[deprecated("Use theta(gctx) instead")]]
  double theta() const {
    return derived().parameters()[eBoundDir1];
  }

  /// Access the phi parameter of the track at the reference surface
  /// @return The phi parameter
  [[deprecated("Use phi(gctx) instead")]]
  double phi() const {
    return derived().parameters()[eBoundDir0];
  }

  /// Access the loc0 parameter of the track at the reference surface
  /// @return The loc0 parameter
  double loc0() const { return derived().parameters()[eBoundLoc0]; }

  /// Access the loc1 parameter of the track at the reference surface
  /// @return The loc1 parameter
  double loc1() const { return derived().parameters()[eBoundLoc1]; }

  /// Access the time parameter of the track at the reference surface
  /// @return The time parameter
  double time() const { return derived().parameters()[eBoundTime]; }

  /// Access the q/p (curvature) parameter of the track at the reference surface
  /// @return The q/p parameter
  double qOverP() const { return derived().parameters()[eBoundQOverP]; }

  /// Get the charge
  /// @return The charge value
  double charge() const {
    return derived().particleHypothesis().extractCharge(qOverP());
  }

  /// Get the absolute momentum
  /// @return The absolute momentum value
  double absoluteMomentum() const {
    return derived().particleHypothesis().extractMomentum(qOverP());
  }

  /// Get the transverse momentum
  /// @return The transverse momentum value
  [[deprecated("Use transverseMomentum(gctx) instead")]]
  double transverseMomentum() const {
    return std::sin(derived().theta()) * derived().absoluteMomentum();
  }

  /// Get a unit vector along the track direction at the reference surface
  /// @return The direction unit vector
  [[deprecated("Use direction(gctx) instead")]]
  Vector3 direction() const {
    return makeDirectionFromPhiTheta(derived().phi(), derived().theta());
  }

  /// Get the global momentum vector
  /// @return the global momentum vector
  [[deprecated("Use momentum(gctx) instead")]]
  Vector3 momentum() const {
    return derived().absoluteMomentum() * derived().direction();
  }

  /// Global azimuthal direction angle in radians.
  /// @param geoCtx Geometry context (unused by the canonical parametrization)
  /// @return Global phi
  double phi(const GeometryContext& /*geoCtx*/) const {
    return derived().parameters()[eBoundDir0];
  }
  /// Global polar direction angle in radians.
  /// @param geoCtx Geometry context (unused by the canonical parametrization)
  /// @return Global theta
  double theta(const GeometryContext& /*geoCtx*/) const {
    return derived().parameters()[eBoundDir1];
  }
  /// Global unit direction vector.
  /// @param geoCtx Geometry context (unused by the canonical parametrization)
  /// @return Normalized global direction
  Vector3 direction(const GeometryContext& /*geoCtx*/) const {
    return BoundParametrization::globalSpherical().direction(
        derived().parameters());
  }
  /// Global momentum vector.
  /// @param geoCtx Geometry context for direction interpretation
  /// @return Global momentum
  Vector3 momentum(const GeometryContext& geoCtx) const {
    return derived().absoluteMomentum() * direction(geoCtx);
  }
  /// Transverse momentum relative to the global z axis.
  /// @param geoCtx Geometry context for direction interpretation
  /// @return Transverse momentum magnitude
  double transverseMomentum(const GeometryContext& geoCtx) const {
    return std::sin(theta(geoCtx)) * derived().absoluteMomentum();
  }
  /// Global four-momentum vector: (px, py, pz, e).
  /// @param geoCtx Geometry context for direction interpretation
  /// @return Global four-momentum
  Vector4 fourMomentum(const GeometryContext& geoCtx) const {
    Vector4 p4;
    const Vector3 p3 = momentum(geoCtx);
    p4.template head<3>() = p3;
    const float m = derived().particleHypothesis().mass();
    p4[eEnergy] = std::sqrt(m * m + p3.squaredNorm());
    return p4;
  }

  /// Get the four-momentum vector: (px, py, pz, e)
  /// @return the four-momentum vector
  [[deprecated("Use fourMomentum(gctx) instead")]]
  Vector4 fourMomentum() const {
    Vector4 p4 = Vector4::Zero();

    Vector3 p3 = derived().momentum();
    p4[eMom0] = p3[eMom0];
    p4[eMom1] = p3[eMom1];
    p4[eMom2] = p3[eMom2];

    float m = derived().particleHypothesis().mass();
    p4[eEnergy] = std::sqrt(m * m + p3.squaredNorm());

    return p4;
  }
};

}  // namespace Acts
