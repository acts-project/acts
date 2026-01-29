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
#include "Acts/Utilities/HashedString.hpp"

namespace Acts {

namespace detail_tp {
inline constexpr HashedString kTipIndexKey = hashString("tipIndex");
inline constexpr HashedString kStemIndexKey = hashString("stemIndex");
inline constexpr HashedString kMeasurementsKey = hashString("nMeasurements");
inline constexpr HashedString kHolesKey = hashString("nHoles");
inline constexpr HashedString kOutliersKey = hashString("nOutliers");
inline constexpr HashedString kSharedHitsKey = hashString("nSharedHits");
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
  constexpr Derived& derived() { return static_cast<Derived&>(*this); }
  constexpr const Derived& derived() const {
    return static_cast<const Derived&>(*this);
  }

 public:
  /// Index type used for referencing track states.
  using IndexType = index_t;

  /// Get the tip index, i.e. the entry point into the track state container.
  IndexType tipIndex() const {
    return derived().template component<IndexType, detail_tp::kTipIndexKey>();
  }

  /// Get a mutable reference to the tip index.
  IndexType& tipIndex()
    requires(!read_only)
  {
    return derived().template component<IndexType, detail_tp::kTipIndexKey>();
  }

  /// Get the stem index, i.e. the innermost track state in the track.
  IndexType stemIndex() const {
    return derived().template component<IndexType, detail_tp::kStemIndexKey>();
  }

  /// Get a mutable reference to the stem index.
  IndexType& stemIndex()
    requires(!read_only)
  {
    return derived().template component<IndexType, detail_tp::kStemIndexKey>();
  }

  /// Return the number of measurements assigned to this track.
  unsigned int nMeasurements() const {
    return derived()
        .template component<unsigned int, detail_tp::kMeasurementsKey>();
  }

  /// Return a mutable reference to the number of measurements.
  unsigned int& nMeasurements()
    requires(!read_only)
  {
    return derived()
        .template component<unsigned int, detail_tp::kMeasurementsKey>();
  }

  /// Return the number of holes on this track.
  unsigned int nHoles() const {
    return derived().template component<unsigned int, detail_tp::kHolesKey>();
  }

  /// Return a mutable reference to the number of holes.
  unsigned int& nHoles()
    requires(!read_only)
  {
    return derived().template component<unsigned int, detail_tp::kHolesKey>();
  }

  /// Return the number of outliers for this track.
  unsigned int nOutliers() const {
    return derived()
        .template component<unsigned int, detail_tp::kOutliersKey>();
  }

  /// Return a mutable reference to the number of outliers.
  unsigned int& nOutliers()
    requires(!read_only)
  {
    return derived()
        .template component<unsigned int, detail_tp::kOutliersKey>();
  }

  /// Return the number of shared hits for this track.
  unsigned int nSharedHits() const {
    return derived()
        .template component<unsigned int, detail_tp::kSharedHitsKey>();
  }

  /// Return a mutable reference to the number of shared hits.
  unsigned int& nSharedHits()
    requires(!read_only)
  {
    return derived()
        .template component<unsigned int, detail_tp::kSharedHitsKey>();
  }

  /// Return the local chi-squared contribution.
  float chi2() const {
    return derived().template component<float, detail_tp::kChi2Key>();
  }

  /// Return a mutable reference to the local chi-squared contribution.
  float& chi2()
    requires(!read_only)
  {
    return derived().template component<float, detail_tp::kChi2Key>();
  }

  /// Return the number of degrees of freedom.
  unsigned int nDoF() const {
    return derived().template component<unsigned int, detail_tp::kNdfKey>();
  }

  /// Return a mutable reference to the number of degrees of freedom.
  unsigned int& nDoF()
    requires(!read_only)
  {
    return derived().template component<unsigned int, detail_tp::kNdfKey>();
  }

  /// Access the theta parameter of the track at the reference surface
  /// @return The theta parameter
  double theta() const { return derived().parameters()[eBoundTheta]; }

  /// Access the phi parameter of the track at the reference surface
  /// @return The phi parameter
  double phi() const { return derived().parameters()[eBoundPhi]; }

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
  double transverseMomentum() const {
    return std::sin(derived().theta()) * derived().absoluteMomentum();
  }

  /// Get a unit vector along the track direction at the reference surface
  /// @return The direction unit vector
  Vector3 direction() const {
    return makeDirectionFromPhiTheta(derived().phi(), derived().theta());
  }

  /// Get the global momentum vector
  /// @return the global momentum vector
  Vector3 momentum() const {
    return derived().absoluteMomentum() * derived().direction();
  }

  /// Get the four-momentum vector: (px, py, pz, e)
  /// @return the four-momentum vector
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
