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

template <typename derived_t, typename traj_t>
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

  using IndexType = typename Container::IndexType;
  static constexpr IndexType kInvalid = Container::kInvalid;

  friend TrackContainer<Container, trajectory_t>;

 private:
  TrackProxy(
      ConstIf<TrackContainer<Container, trajectory_t>, ReadOnly>& container,
      IndexType itrack)
      : m_container{&container}, m_index{itrack} {}

  ConstIf<TrackContainer<Container, trajectory_t>, ReadOnly>* m_container;
  IndexType m_index;
};
}  // namespace detail_tc

template <typename T>
struct isReadOnlyTrackContainer;

template <typename derived_t>
class TrackContainerBackend {
 public:
  using Derived = derived_t;

  static constexpr bool ReadOnly = isReadOnlyTrackContainer<Derived>::value;

  using IndexType = MultiTrajectoryTraits::IndexType;
  static constexpr IndexType kInvalid = MultiTrajectoryTraits::kInvalid;

  // using TrackProxy = detail_tc::TrackProxy<Derived, traj_t, false>;
  // using ConstTrackProxy = detail_tc::TrackProxy<Derived, traj_t, true>;

 protected:
  TrackContainerBackend() = default;  // pseudo abstract base class

 private:
  /// Helper to static cast this to the Derived class for CRTP
  constexpr Derived& self() { return static_cast<Derived&>(*this); }
  /// Helper to static cast this to the Derived class for CRTP. Const version.
  constexpr const Derived& self() const {
    return static_cast<const Derived&>(*this);
  }

 public:
  // ConstTrackProxy getTrack(IndexType itrack) const { return {*this, itrack};
  // }

  // template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  // TrackProxy getTrack(IndexType itrack) {
  // return {*this, itrack};
  // }

  constexpr IndexType size() const { return self().size_impl(); }

  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  IndexType addTrack() {
    return self().addTrack_impl();
  }

  template <typename T, bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  constexpr void addColumn(const std::string& key) {
    self().template addColumn_impl<T>(key);
  }

  constexpr bool hasColumn(HashedString key) const {
    return self().hasColumn_impl(key);
  }

 protected:
  template <typename T, HashedString key, bool RO = ReadOnly,
            typename = std::enable_if_t<!RO>>
  constexpr T& component(IndexType itrack) {
    return *std::any_cast<T*>(self().component_impl(key, itrack));
  }

  template <typename T, bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  constexpr T& component(HashedString key, IndexType istate) {
    return *std::any_cast<T*>(self().component_impl(key, istate));
  }

  template <typename T, HashedString key>
  constexpr const T& component(IndexType itrack) const {
    return *std::any_cast<const T*>(self().component_impl(key, itrack));
  }

  template <typename T>
  constexpr const T& component(HashedString key, IndexType itrack) const {
    return *std::any_cast<const T*>(self().component_impl(key, itrack));
  }

  // private:
  // traj_t m_traj;
};

}  // namespace Acts
