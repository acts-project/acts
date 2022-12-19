// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/Utilities/HashedString.hpp"

#include <any>
#include <cstddef>
#include <iterator>

namespace Acts {

template <typename track_container_t, typename traj_t,
          template <typename> class holder_t>
class TrackContainer;

namespace detail_tc {
template <typename T, bool select>
using ConstIf = std::conditional_t<select, const T, T>;

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

  /// Return the number of track states associated to this track
  /// @note This is calculated by iterating over the track states which is
  ///       somewhat expensive. Consider caching this value if you need It
  ///       more than once.
  /// @return The number of track states
  unsigned int nTrackStates() const {
    // @TODO: This should probably be cached, distance is expensive
    //        without random access
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

  /// Return a reference to the track container backend, const version.
  /// @return reference to the track container backend
  const auto& container() const { return *m_container; }

  /// Equality operator with another track proxy
  /// Checks the container identity and the track index
  /// @return True if the track proxies refer to the same track
  bool operator==(const TrackProxy& other) const {
    return &(*m_container) == &(*other.m_container) && m_index == other.m_index;
  }

 private:
  TrackProxy(ConstIf<TrackContainer<Container, Trajectory, holder_t>, ReadOnly>&
                 container,
             IndexType itrack)
      : m_container{&container}, m_index{itrack} {}

  detail_lt::TransitiveConstPointer<
      ConstIf<TrackContainer<Container, Trajectory, holder_t>, ReadOnly>>
      m_container;
  IndexType m_index;
};

/// Internal holder type for referencing a backend without ownership
template <typename T>
struct RefHolder {
  T* ptr;

  RefHolder(T* _ptr) : ptr{_ptr} {}
  RefHolder(T& ref) : ptr{&ref} {}

  const T& operator*() const { return *ptr; }
  T& operator*() { return *ptr; }

  const T* operator->() const { return ptr; }
  T* operator->() { return ptr; }
};

/// Internal holder type holding a backend container by value
template <typename T>
struct ValueHolder {
  T val;

  ValueHolder(T& _val) : val{_val} {}
  ValueHolder(T&& _val) : val{std::move(_val)} {}

  const T& operator*() const { return val; }
  T& operator*() { return val; }

  const T* operator->() const { return &val; }
  T* operator->() { return &val; }
};

template <template <typename...> class, template <typename...> class>
struct is_same_template : std::false_type {};

template <template <typename...> class T>
struct is_same_template<T, T> : std::true_type {};

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

template <typename T>
struct IsReadOnlyTrackContainer;

/// Track container interface class. This type represents a collections of
/// tracks. It uses a backend to store bothe the actual tracks and the
/// associated track states.
/// @tparam track_container_t the track container backend
/// @tparam traj_t the track state container backend
/// @tparam holder_t ownership management class for the backend
template <typename track_container_t, typename traj_t,
          template <typename> class holder_t = detail_tc::RefHolder>
class TrackContainer {
 public:
  static constexpr bool ReadOnly =
      IsReadOnlyTrackContainer<track_container_t>::value;
  static constexpr bool TrackStateReadOnly =
      IsReadOnlyMultiTrajectory<traj_t>::value;

  static_assert(ReadOnly == TrackStateReadOnly,
                "Either both track container and track state container need to "
                "be readonly or both have to be readwrite");

  using IndexType = MultiTrajectoryTraits::IndexType;
  static constexpr IndexType kInvalid = MultiTrajectoryTraits::kInvalid;

  using TrackProxy =
      detail_tc::TrackProxy<track_container_t, traj_t, holder_t, false>;
  using ConstTrackProxy =
      detail_tc::TrackProxy<track_container_t, traj_t, holder_t, true>;

#ifndef DOXYGEN
  friend TrackProxy;
  friend ConstTrackProxy;
#endif

  /// Constructor from a track container backend and a track state container
  /// backend
  /// @param container the track container backend
  /// @param traj the track state container backend
  TrackContainer(holder_t<track_container_t> container, holder_t<traj_t> traj)
      : m_container{std::move(container)}, m_traj{std::move(traj)} {}

  /// Constructor from references to a track container backend and to a track
  /// state container backend
  /// @note The track container will not assume ownership over the backends in this case.
  ///       You need to ensure suitable lifetime
  /// @param container the track container backend
  /// @param traj the track state container backend
  template <template <typename> class H = holder_t,
            typename = std::enable_if_t<
                detail_tc::is_same_template<H, detail_tc::RefHolder>::value>>
  TrackContainer(track_container_t& container, traj_t& traj)
      : m_container{&container}, m_traj{&traj} {}

  /// Get a const track proxy for a track index
  /// @param itrack the track index in the container
  /// @return A const track proxy for the index
  ConstTrackProxy getTrack(IndexType itrack) const { return {*this, itrack}; }

  /// Get a mutable track proxy for a track index
  /// @param itrack the track index in the container
  /// @return A mutable track proxy for the index
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  TrackProxy getTrack(IndexType itrack) {
    return {*this, itrack};
  }

  /// Get the size of the track container
  /// @return the sixe
  constexpr IndexType size() const { return m_container->size_impl(); }

  /// Add a track to the container. Note this only creates the logical track and
  /// allocates memory. You can combine this with @c getTrack to obtain a track proxy
  /// @return the index to the newly added track
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  IndexType addTrack() {
    return m_container->addTrack_impl();
  }

  /// Add a dymanic column to the track container
  /// @param key the name of the column to be added
  template <typename T, bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  constexpr void addColumn(const std::string& key) {
    m_container->template addColumn_impl<T>(key);
  }

  /// Check if this track container has a specific dynamic column
  /// @param key the key to check for
  constexpr bool hasColumn(const std::string& key) const {
    return m_container->hasColumn_impl(hashString(key));
  }

  /// Check if a this track container has a specific dynamic column
  /// @param key the key to check for
  constexpr bool hasColumn(HashedString key) const {
    return m_container->hasColumn_impl(key);
  }

  /// Get a mutable reference to the track container backend
  /// @return a mutable reference to the backend
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  auto& container() {
    return *m_container;
  }

  /// Get a const reference to the track container backend
  /// @return a const reference to the backend
  const auto& container() const { return *m_container; }

  /// Get a mutable reference to the track state container backend
  /// @return a mutable reference to the backend
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  auto& trackStateContainer() {
    return *m_traj;
  }

  /// Get a const reference to the track state container backend
  /// @return a const reference to the backend
  const auto& trackStateContainer() const { return *m_traj; }

  /// Get a mutable iterator to the first track in the container
  /// @return a mutable iterator to the first track
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  auto begin() {
    return detail_tc::TrackProxyIterator<std::decay_t<decltype(*this)>,
                                         TrackProxy, false>{*this, 0};
  }

  /// Get a past-the-end iterator for this container
  /// @return a past-the-end iterator
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  auto end() {
    return detail_tc::TrackProxyIterator<std::decay_t<decltype(*this)>,
                                         TrackProxy, false>{*this, size()};
  }

  /// Get an const iterator to the first track in the container
  /// @return a const iterator to the first track
  auto begin() const {
    return detail_tc::TrackProxyIterator<std::decay_t<decltype(*this)>,
                                         ConstTrackProxy, true>{*this, 0};
  }

  /// Get a past-the-end iterator for this container
  /// @return a past-the-end iterator
  auto end() const {
    return detail_tc::TrackProxyIterator<std::decay_t<decltype(*this)>,
                                         ConstTrackProxy, true>{*this, size()};
  }

 protected:
  template <typename T, HashedString key, bool RO = ReadOnly,
            typename = std::enable_if_t<!RO>>
  constexpr T& component(IndexType itrack) {
    return *std::any_cast<T*>(container().component_impl(key, itrack));
  }

  template <typename T, bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  constexpr T& component(HashedString key, IndexType itrack) {
    return *std::any_cast<T*>(container().component_impl(key, itrack));
  }

  template <typename T, HashedString key>
  constexpr const T& component(IndexType itrack) const {
    return *std::any_cast<const T*>(container().component_impl(key, itrack));
  }

  template <typename T>
  constexpr const T& component(HashedString key, IndexType itrack) const {
    return *std::any_cast<const T*>(container().component_impl(key, itrack));
  }

  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  constexpr typename TrackProxy::Parameters parameters(IndexType itrack) {
    return container().parameters(itrack);
  }

  constexpr typename ConstTrackProxy::ConstParameters parameters(
      IndexType itrack) const {
    return container().parameters(itrack);
  }

  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  constexpr typename TrackProxy::Covariance covariance(IndexType itrack) {
    return container().covariance(itrack);
  }

  constexpr typename ConstTrackProxy::ConstCovariance covariance(
      IndexType itrack) const {
    return container().covariance(itrack);
  }

  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  auto trackStateRange(IndexType itrack) {
    auto tip = component<IndexType, hashString("tipIndex")>(itrack);
    return m_traj->trackStateRange(tip);
  }

  auto trackStateRange(IndexType itrack) const {
    auto tip = component<IndexType, hashString("tipIndex")>(itrack);
    return m_traj->trackStateRange(tip);
  }

 private:
  detail_tc::ConstIf<holder_t<track_container_t>, ReadOnly> m_container;
  detail_tc::ConstIf<holder_t<traj_t>, ReadOnly> m_traj;
};

template <typename track_container_t, typename traj_t>
TrackContainer(track_container_t& container, traj_t& traj)
    -> TrackContainer<track_container_t, traj_t, detail_tc::RefHolder>;

template <typename track_container_t, typename traj_t>
TrackContainer(track_container_t&& container, traj_t&& traj)
    -> TrackContainer<track_container_t, traj_t, detail_tc::ValueHolder>;

/// Utility class that eases accessing dynamic columns in track containers
/// @tparam T the type of the value to access
/// @tparam ReadOnly true if this is a const accessor
template <typename T, bool ReadOnly>
struct TrackAccessorBase {
  HashedString key;

  /// Create the accessor from an already-hashed string key
  /// @param _key the key
  TrackAccessorBase(HashedString _key) : key{_key} {}
  /// Create the accessor from a string key
  /// @param _key the key
  TrackAccessorBase(const std::string& _key) : key{hashString(_key)} {}

  /// Access the stored key on the track given as an argument. Mutable version
  /// @tparam track_proxy_t the type of the track proxy
  /// @param track the track to access
  /// @return mutable reference to the column behind the key
  template <typename track_proxy_t, bool RO = ReadOnly,
            typename = std::enable_if_t<!RO>>
  T& operator()(track_proxy_t track) {
    static_assert(!track_proxy_t::ReadOnly,
                  "Cannot get mutable ref for const track proxy");
    return track.template component<T>(key);
  }

  /// Access the stored key on the track given as an argument. COnst version
  /// @tparam track_proxy_t the type of the track proxy
  /// @param track the track to access
  /// @return const reference to the column behind the key
  template <typename track_proxy_t, bool RO = ReadOnly,
            typename = std::enable_if_t<RO>>
  const T& operator()(track_proxy_t track) {
    if constexpr (track_proxy_t::ReadOnly) {
      return track.template component<T>(key);
    } else {
      typename track_proxy_t::ConstTrackProxy ctrack{track};
      return ctrack.template component<T>(key);
    }
  }
};

template <typename T>
using TrackAccessor = TrackAccessorBase<T, false>;
template <typename T>
using ConstTrackAccessor = TrackAccessorBase<T, true>;

}  // namespace Acts
