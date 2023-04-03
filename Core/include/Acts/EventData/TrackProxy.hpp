// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/Charge.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/TrackStatePropMask.hpp"
#include "Acts/Utilities/UnitVectors.hpp"

#include <iterator>
#include <type_traits>

namespace Acts {

template <typename track_container_t, typename traj_t,
          template <typename> class holder_t>
class TrackContainer;

namespace detail_tc {
template <typename T, bool select>
using ConstIf = std::conditional_t<select, const T, T>;

/// Helper iterator to allow iteration over tracks via track proxies.
template <typename container_t, typename proxy_t, bool ReadOnly>
class TrackProxyIterator {
  using ProxyType = proxy_t;
  using IndexType = typename ProxyType::IndexType;
  using ContainerType = container_t;

 public:
  using iterator_category = std::random_access_iterator_tag;
  using value_type = ProxyType;
  using difference_type = std::ptrdiff_t;
  using pointer = void;
  using reference = void;

  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  TrackProxyIterator(container_t& container, IndexType itrack)
      : m_container(&container), m_itrack(itrack) {}

  template <bool RO = ReadOnly, typename = std::enable_if_t<RO>>
  TrackProxyIterator(const container_t& container, IndexType itrack)
      : m_container(&container), m_itrack(itrack) {}

  TrackProxyIterator& operator++() {
    m_itrack++;
    return *this;
  }
  TrackProxyIterator& operator--() {
    m_itrack--;
    return *this;
  }

  bool operator==(const TrackProxyIterator& other) const {
    return m_container == other.m_container && m_itrack == other.m_itrack;
  }

  bool operator!=(const TrackProxyIterator& other) const {
    return !(*this == other);
  }

  bool operator<(const TrackProxyIterator& other) const {
    return m_itrack < other.m_itrack;
  }

  bool operator>(const TrackProxyIterator& other) const {
    return m_itrack > other.m_itrack;
  }

  bool operator<=(const TrackProxyIterator& other) const {
    return m_itrack <= other.m_itrack;
  }

  bool operator>=(const TrackProxyIterator& other) const {
    return m_itrack >= other.m_itrack;
  }

  ProxyType operator*() const { return m_container->getTrack(m_itrack); }

  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  ProxyType operator*() {
    return m_container->getTrack(m_itrack);
  }

  TrackProxyIterator operator[](difference_type n) const {
    TrackProxyIterator copy = *this;
    copy += n;
    return copy;
  };

  TrackProxyIterator& operator+=(difference_type n) {
    m_itrack += n;
    return *this;
  }

  TrackProxyIterator operator-=(difference_type n) {
    m_itrack -= n;
    return *this;
  }

  friend difference_type operator-(const TrackProxyIterator& lhs,
                                   const TrackProxyIterator& rhs) {
    return lhs.m_itrack - rhs.m_itrack;
  }

  friend TrackProxyIterator operator+(const TrackProxyIterator& lhs,
                                      difference_type rhs) {
    TrackProxyIterator copy = lhs;
    copy += rhs;
    return copy;
  }

  friend TrackProxyIterator operator+(difference_type lhs,
                                      const TrackProxyIterator& rhs) {
    return rhs + lhs;
  }

  friend TrackProxyIterator operator-(const TrackProxyIterator& lhs,
                                      difference_type rhs) {
    return lhs + (-rhs);
  }

  friend TrackProxyIterator operator-(difference_type lhs,
                                      const TrackProxyIterator& rhs) {
    return rhs + (-lhs);
  }

 private:
  detail_lt::TransitiveConstPointer<ConstIf<ContainerType, ReadOnly>>
      m_container;
  IndexType m_itrack;
};

}  // namespace detail_tc

/// Proxy class representing a single track
/// @tparam track_container_t the container backend
/// @tparam trajectory_t the track state container backend
/// @tparam holder_t ownership management class for the backend
/// @tparam read_only true if this track container is not mutable
template <typename track_container_t, typename trajectory_t,
          template <typename> class holder_t, bool read_only = true>
class TrackProxy {
 public:
  static constexpr bool ReadOnly = read_only;
  using Container = track_container_t;
  using Trajectory = trajectory_t;

  using MutableTrackProxy =
      TrackProxy<track_container_t, trajectory_t, holder_t, false>;
  using ConstTrackProxy =
      TrackProxy<track_container_t, trajectory_t, holder_t, true>;

  using TrackStateProxy = typename Trajectory::TrackStateProxy;
  using ConstTrackStateProxy = typename Trajectory::ConstTrackStateProxy;

  using Parameters =
      typename detail_lt::Types<eBoundSize, false>::CoefficientsMap;
  using ConstParameters =
      typename detail_lt::Types<eBoundSize, true>::CoefficientsMap;

  using Covariance =
      typename detail_lt::Types<eBoundSize, false>::CovarianceMap;
  using ConstCovariance =
      typename detail_lt::Types<eBoundSize, true>::CovarianceMap;

  using IndexType = typename Container::IndexType;
  static constexpr IndexType kInvalid = Container::kInvalid;

#ifndef DOXYGEN
  friend TrackContainer<Container, Trajectory, holder_t>;
  friend MutableTrackProxy;
  friend ConstTrackProxy;
#endif

  /// Copy constructor from a mutable track proxy. This is always valid, either
  /// mutable to mutable or mutable to const
  /// @param other the other track state proxy
  TrackProxy(const MutableTrackProxy& other)
      : m_container{other.m_container}, m_index{other.m_index} {}

  /// Copy assignment operator from mutable track proxy. This is always valid,
  /// either mutable to mutable or mutable to const
  /// @param other the other track state proxy
  TrackProxy& operator=(const MutableTrackProxy& other) {
    m_container = other.m_container;
    m_index = other.m_index;
    return *this;
  }

  /// Retrieve a mutable reference to a component
  /// @tparam T The type of the component to access
  /// @tparam key String key for the component to access
  /// @return Mutable reference to the component given by @p key
  template <typename T, HashedString key, bool RO = ReadOnly,
            typename = std::enable_if_t<!RO>>
  constexpr T& component() {
    return m_container->template component<T, key>(m_index);
  }

  /// Retrieve a mutable reference to a component
  /// @tparam T The type of the component to access
  /// @param key String key for the component to access
  /// @return Mutable reference to the component given by @p key
  template <typename T, bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  constexpr T& component(HashedString key) {
    return m_container->template component<T>(key, m_index);
  }

  /// Retrieve a mutable reference to a component
  /// @tparam T The type of the component to access
  /// @param key String key for the component to access
  /// @note This might hash the @p key at runtime instead of compile-time
  /// @return Mutable reference to the component given by @p key
  template <typename T, bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  constexpr T& component(std::string_view key) {
    return m_container->template component<T>(hashString(key), m_index);
  }

  /// Retrieve a const reference to a component
  /// @tparam T The type of the component to access
  /// @tparam key String key for the component to access
  /// @return Const reference to the component given by @p key
  template <typename T, HashedString key>
  constexpr const T& component() const {
    return m_container->template component<T, key>(m_index);
  }

  /// Retrieve a const reference to a component
  /// @tparam T The type of the component to access
  /// @param key String key for the component to access
  /// @return Const reference to the component given by @p key
  template <typename T>
  constexpr const T& component(HashedString key) const {
    return m_container->template component<T>(key, m_index);
  }

  /// Retrieve a const reference to a component
  /// @tparam T The type of the component to access
  /// @param key String key for the component to access
  /// @note This might hash the @p key at runtime instead of compile-time
  /// @return Const reference to the component given by @p key
  template <typename T>
  constexpr const T& component(std::string_view key) const {
    return m_container->template component<T>(hashString(key), m_index);
  }

  /// Get the tip index, i.e. the entry point into the track state container
  /// @return the tip index by value
  IndexType tipIndex() const {
    return component<IndexType>(hashString("tipIndex"));
  }

  /// Get a mutable reference to the tip index, i.e. the entry point into the
  /// track container
  /// @return mutable reference to the tip index
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  IndexType& tipIndex() {
    return component<IndexType>(hashString("tipIndex"));
  }

  /// Get the reference surface of the track (e.g. the perigee)
  /// @return the reference surface
  const Surface& referenceSurface() const {
    return *component<std::shared_ptr<const Surface>,
                      hashString("referenceSurface")>();
  }

  // NOLINTBEGIN(performance-unnecessary-value-param)
  // looks like a false-positive. clang-tidy believes `srf` is not movable.
  /// Set a new reference surface for this track
  /// @param srf The surface to set
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  void setReferenceSurface(std::shared_ptr<const Surface> srf) {
    component<std::shared_ptr<const Surface>,
              hashString("referenceSurface")>() = std::move(srf);
  }
  // NOLINTEND(performance-unnecessary-value-param)

  /// Return whether a reference surface is associated to this track
  /// @return whether a surface exists or not
  bool hasReferenceSurface() const {
    return !!component<std::shared_ptr<const Surface>,
                       hashString("referenceSurface")>();
  }

  /// Get the parameters of the track at the reference surface (e.g. perigee).
  /// Const version
  /// @return Proxy vector for the parameters
  ConstParameters parameters() const {
    return m_container->parameters(m_index);
  }

  /// Get the covariance of the track at the reference surface (e.g. perigee).
  /// Const version
  /// @return Proxy matrix for the covariance
  ConstCovariance covariance() const {
    return m_container->covariance(m_index);
  }

  /// Get the parameters of the track at the reference surface (e.g. perigee).
  /// Mutable version
  /// @return Proxy vector for the parameters
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  Parameters parameters() {
    return m_container->parameters(m_index);
  }

  /// Get the covariance of the track at the reference surface (e.g. perigee).
  /// Mutable version
  /// @return Proxy matrix for the covariance
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  Covariance covariance() {
    return m_container->covariance(m_index);
  }

  ActsScalar charge() const {
    // Currently, neutral tracks are not supported here
    // @TODO: Evaluate if/how neutral 'tracks' should be accounted for
    return SinglyCharged{}.extractCharge(parameters()[eBoundQOverP]);
  }

  /// Access the theta parameter of the track at the reference surface
  /// @return The theta parameter
  ActsScalar theta() const { return parameters()[eBoundTheta]; }

  /// Access the phi parameter of the track at the reference surface
  /// @return The phi parameter
  ActsScalar phi() const { return parameters()[eBoundPhi]; }

  /// Access the loc0 parameter of the track at the reference surface
  /// @return The loc0 parameter
  ActsScalar loc0() const { return parameters()[eBoundLoc0]; }

  /// Access the loc1 parameter of the track at the reference surface
  /// @return The loc1 parameter
  ActsScalar loc1() const { return parameters()[eBoundLoc1]; }

  /// Access the time parameter of the track at the reference surface
  /// @return The time parameter
  ActsScalar time() const { return parameters()[eBoundTime]; }

  /// Access the q/p (curvature) parameter of the track at the reference surface
  /// @return The q/p parameter
  ActsScalar qOverP() const { return parameters()[eBoundQOverP]; }

  /// Get the absolute momentum of the tack
  /// @return The absolute track momentum
  ActsScalar absoluteMomentum() const {
    return SinglyCharged{}.extractMomentum(qOverP());
  }

  /// Get the transverse momentum of the track
  /// @return The track transverse momentum value
  ActsScalar transverseMomentum() const {
    return std::sin(theta()) * absoluteMomentum();
  }

  /// Get a unit vector along the track direction at the reference surface
  /// @return The direction unit vector
  Vector3 unitDirection() const {
    return makeDirectionUnitFromPhiTheta(phi(), theta());
  }

  /// Get the global momentum vector
  /// @return the global momentum vector
  Vector3 momentum() const { return absoluteMomentum() * unitDirection(); }

  /// Get a range over the track states of this track. Return value is
  /// compatible with range based for loop. Const version
  /// @return Track state range to iterate over
  auto trackStates() const { return m_container->trackStateRange(m_index); }

  /// Get a range over the track states of this track. Return value is
  /// compatible with range based for loop. Mutable version
  /// @return Track state range to iterate over
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  auto trackStates() {
    return m_container->trackStateRange(m_index);
  }

  /// Append a track state to this track. This will modify the tip index to
  /// point at the newly created track state, which will be directly after the
  /// previous track state at tip index.
  /// @param mask The allocation prop mask for the new track state
  /// @return The newly added track state
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  auto appendTrackState(TrackStatePropMask mask = TrackStatePropMask::All) {
    auto& tsc = m_container->trackStateContainer();
    auto ts = tsc.getTrackState(tsc.addTrackState(mask, tipIndex()));
    tipIndex() = ts.index();
    return ts;
  }

  /// Return the number of track states associated to this track
  /// @note This is calculated by iterating over the track states which is
  ///       somewhat expensive. Consider caching this value if you need It
  ///       more than once.
  /// @return The number of track states
  unsigned int nTrackStates() const {
    // @TODO: This should probably be cached, distance is expensive
    //        without random access
    if (tipIndex() == kInvalid) {
      // no tip index -> no track states
      return 0;
    }
    auto tsRange = trackStates();
    return std::distance(tsRange.begin(), tsRange.end());
  }

  /// Return the number of measurements for the track. Const version
  /// @return The number of measurements
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  unsigned int& nMeasurements() {
    return component<unsigned int>(hashString("nMeasurements"));
  }

  /// Return a mutable reference to the number of measurements for the track.
  /// Mutable version
  /// @return The number of measurements
  unsigned int nMeasurements() const {
    return component<unsigned int>(hashString("nMeasurements"));
  }

  /// Return a mutable reference to the number of holes for the track.
  /// Mutable version
  /// @return The number of holes
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  unsigned int& nHoles() {
    return component<unsigned int>(hashString("nHoles"));
  }

  /// Return the number of measurements for the track. Const version
  /// @return The number of measurements
  unsigned int nHoles() const {
    return component<unsigned int>(hashString("nHoles"));
  }

  /// Return a mutable reference to the number of outliers for the track.
  /// Mutable version
  /// @return The number of outliers
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  unsigned int& nOutliers() {
    return component<unsigned int>(hashString("nOutliers"));
  }

  /// Return the number of outliers for the track. Const version
  /// @return The number of outliers
  unsigned int nOutliers() const {
    return component<unsigned int>(hashString("nOutliers"));
  }

  /// Return a mutable reference to the number of shared hits for the track.
  /// Mutable version
  /// @return The number of shared hits
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  unsigned int& nSharedHits() {
    return component<unsigned int>(hashString("nSharedHits"));
  }

  /// Return the number of shared hits for the track. Const version
  /// @return The number of shared hits
  unsigned int nSharedHits() const {
    return component<unsigned int>(hashString("nSharedHits"));
  }

  /// Return a mutable reference to the chi squared
  /// Mutable version
  /// @return The chi squared
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  float& chi2() {
    return component<float>(hashString("chi2"));
  }

  /// Return the chi squared for the track. Const version
  /// @return The chi squared
  float chi2() const { return component<float>(hashString("chi2")); }

  /// Return a mutable reference to the number of degrees of freedom for the
  /// track. Mutable version
  /// @return The the number of degrees of freedom
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  unsigned int& nDoF() {
    return component<unsigned int>(hashString("ndf"));
  }

  /// Return the number of degrees of freedom for the track. Const version
  /// @return The number of degrees of freedom
  unsigned int nDoF() const {
    return component<unsigned int>(hashString("ndf"));
  }

  /// Return the index of this track in the track container
  /// @note This is separate from the tip index
  /// @return the track index
  IndexType index() const { return m_index; }

  /// Return a reference to the track container backend, mutable version.
  /// @return reference to the track container backend
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  auto& container() {
    return *m_container;
  }

  /// Copy the content of another track proxy into this one
  /// @tparam track_proxy_t the other track proxy's type
  /// @param other The the track proxy
  template <typename track_proxy_t, bool RO = ReadOnly,
            typename = std::enable_if_t<!RO>>
  void copyFrom(const track_proxy_t& other) {
    // @TODO: Add constraint on which track proxies are allowed,
    // this is only implicit right now

    tipIndex() = other.tipIndex();
    parameters() = other.parameters();
    covariance() = other.covariance();
    if (other.hasReferenceSurface()) {
      setReferenceSurface(other.referenceSurface().getSharedPtr());
    }
    nMeasurements() = other.nMeasurements();
    nHoles() = other.nHoles();

    // This will only be valid if the backends match and support this operation
    m_container->copyDynamicFrom(m_index, other.m_container->container(),
                                 other.m_index);
  }

  /// Return a reference to the track container backend, const version.
  /// @return reference to the track container backend
  const auto& container() const { return *m_container; }

  /// Equality operator with another track proxy
  /// Checks the container identity and the track index
  /// @return True if the track proxies refer to the same track
  bool operator==(const TrackProxy& other) const {
    return &(*m_container) == &(*other.m_container) && m_index == other.m_index;
  }

  // Track proxies are friends, not food!
  template <typename A, typename B, template <typename> class H, bool R>
  friend class TrackProxy;

 private:
  TrackProxy(detail_tc::ConstIf<TrackContainer<Container, Trajectory, holder_t>,
                                ReadOnly>& container,
             IndexType itrack)
      : m_container{&container}, m_index{itrack} {}

  detail_lt::TransitiveConstPointer<detail_tc::ConstIf<
      TrackContainer<Container, Trajectory, holder_t>, ReadOnly>>
      m_container;
  IndexType m_index;
};
}  // namespace Acts
