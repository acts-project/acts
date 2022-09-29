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

namespace Acts {

template <typename track_container_t, typename traj_t>
class TrackContainer;

namespace detail_tc {
template <typename T, bool select>
using ConstIf = std::conditional_t<select, const T, T>;

template <typename track_container_t, typename trajectory_t,
          bool ReadOnly = true>
class TrackProxy {
 public:
  using Container = track_container_t;
  using Trajectory = trajectory_t;

  using TrackStateProxy = typename Trajectory::TrackStateProxy;
  using ConstTrackStateProxy = typename Trajectory::ConstTrackStateProxy;

  using Parameters =
      typename detail_lt::Types<eBoundSize, ReadOnly>::CoefficientsMap;
  using Covariance =
      typename detail_lt::Types<eBoundSize, ReadOnly>::CovarianceMap;

  using IndexType = typename Container::IndexType;
  static constexpr IndexType kInvalid = Container::kInvalid;

  friend TrackContainer<Container, Trajectory>;

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

  constexpr Parameters parameters() const {
    return m_container->parameters(m_index);
  }

  constexpr Covariance covariance() const {
    return m_container->covariance(m_index);
  }

 private:
  TrackProxy(
      ConstIf<TrackContainer<Container, Trajectory>, ReadOnly>& container,
      IndexType itrack)
      : m_container{&container}, m_index{itrack} {}

  ConstIf<TrackContainer<Container, Trajectory>, ReadOnly>* m_container;
  IndexType m_index;
};
}  // namespace detail_tc

template <typename T>
struct isReadOnlyTrackContainer;

// template <typename derived_t>
// class TrackContainerBackend {
// public:
// using Derived = derived_t;

// static constexpr bool ReadOnly = isReadOnlyTrackContainer<Derived>::value;

// using IndexType = MultiTrajectoryTraits::IndexType;
// static constexpr IndexType kInvalid = MultiTrajectoryTraits::kInvalid;

// protected:
// TrackContainerBackend() = default;  // pseudo abstract base class

// private:
// /// Helper to static cast this to the Derived class for CRTP
// constexpr Derived& self() { return static_cast<Derived&>(*this); }
// /// Helper to static cast this to the Derived class for CRTP. Const version.
// constexpr const Derived& self() const {
// return static_cast<const Derived&>(*this);
// }

// public:
// };

template <typename track_container_t, typename traj_t>
class TrackContainer {
 public:
  static constexpr bool ReadOnly =
      isReadOnlyTrackContainer<track_container_t>::value;

  using IndexType = MultiTrajectoryTraits::IndexType;
  static constexpr IndexType kInvalid = MultiTrajectoryTraits::kInvalid;

  using TrackProxy = detail_tc::TrackProxy<track_container_t, traj_t, false>;
  using ConstTrackProxy =
      detail_tc::TrackProxy<track_container_t, traj_t, true>;

  friend TrackProxy;
  friend ConstTrackProxy;

  TrackContainer(track_container_t container, traj_t traj)
      : m_container{std::move(container)}, m_traj{std::move(traj)} {}

  ConstTrackProxy getTrack(IndexType itrack) const { return {*this, itrack}; }

  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  TrackProxy getTrack(IndexType itrack) {
    return {*this, itrack};
  }

  constexpr IndexType size() const { return m_container.size_impl(); }

  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  IndexType addTrack() {
    return m_container.addTrack_impl();
  }

  template <typename T, bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  constexpr void addColumn(const std::string& key) {
    m_container.template addColumn_impl<T>(key);
  }

  constexpr bool hasColumn(HashedString key) const {
    return m_container.hasColumn_impl(key);
  }

 protected:
  template <typename T, HashedString key, bool RO = ReadOnly,
            typename = std::enable_if_t<!RO>>
  constexpr T& component(IndexType itrack) {
    return *std::any_cast<T*>(m_container.component_impl(key, itrack));
  }

  template <typename T, bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  constexpr T& component(HashedString key, IndexType istate) {
    return *std::any_cast<T*>(m_container.component_impl(key, istate));
  }

  template <typename T, HashedString key>
  constexpr const T& component(IndexType itrack) const {
    return *std::any_cast<const T*>(m_container.component_impl(key, itrack));
  }

  template <typename T>
  constexpr const T& component(HashedString key, IndexType itrack) const {
    return *std::any_cast<const T*>(m_container.component_impl(key, itrack));
  }

  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  constexpr typename TrackProxy::Parameters parameters(IndexType itrack) {
    return m_container.parameters(itrack);
  }

  constexpr typename TrackProxy::Parameters parameters(IndexType itrack) {
    return m_container.parameters(itrack);
  }

  constexpr typename ConstTrackProxy::Parameters parameters(
      IndexType itrack) const {
    return m_container.parameters(itrack);
  }

  constexpr typename TrackProxy::Covariance covariance(IndexType itrack) {
    return m_container.covariance(itrack);
  }

  constexpr typename ConstTrackProxy::Covariance covariance(
      IndexType itrack) const {
    return m_container.covariance(itrack);
  }

 private:
  track_container_t m_container;
  traj_t m_traj;
};

}  // namespace Acts
