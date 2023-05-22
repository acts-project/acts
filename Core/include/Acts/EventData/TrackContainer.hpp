// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/Charge.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/TrackProxy.hpp"
#include "Acts/EventData/Utils.hpp"
#include "Acts/Utilities/HashedString.hpp"
#include "Acts/Utilities/Holders.hpp"
#include "Acts/Utilities/UnitVectors.hpp"

#include <any>
#include <cstddef>
#include <iterator>

namespace Acts {

template <typename T>
struct IsReadOnlyTrackContainer;

/// Track container interface class. This type represents a collections of
/// tracks. It uses a backend to store bothe the actual tracks and the
/// associated track states.
/// @tparam track_container_t the track container backend
/// @tparam traj_t the track state container backend
/// @tparam holder_t ownership management class for the backend
template <typename track_container_t, typename traj_t,
          template <typename> class holder_t = detail::RefHolder>
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
      Acts::TrackProxy<track_container_t, traj_t, holder_t, false>;
  using ConstTrackProxy =
      Acts::TrackProxy<track_container_t, traj_t, holder_t, true>;

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
                detail::is_same_template<H, detail::RefHolder>::value>>
  TrackContainer(track_container_t& container, traj_t& traj)
      : m_container{&container}, m_traj{&traj} {}

  /// Constructor from const references to a track container backend and to a
  /// track state container backend
  /// @note The track container will not assume ownership over the backends in this case.
  ///       You need to ensure suitable lifetime
  /// @param container the track container backend
  /// @param traj the track state container backend
  template <
      template <typename> class H = holder_t,
      bool RO = (IsReadOnlyTrackContainer<track_container_t>::value &&
                 IsReadOnlyMultiTrajectory<traj_t>::value),
      typename = std::enable_if_t<
          detail::is_same_template<H, detail::ConstRefHolder>::value && RO>>
  TrackContainer(const track_container_t& container, const traj_t& traj)
      : m_container{&container}, m_traj{&traj} {}

  /// Get a const track proxy for a track index
  /// @param itrack the track index in the container
  /// @return A const track proxy for the index
  ConstTrackProxy getTrack(IndexType itrack) const {
    return {*this, itrack};
  }

  /// Get a mutable track proxy for a track index
  /// @param itrack the track index in the container
  /// @return A mutable track proxy for the index
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  TrackProxy getTrack(IndexType itrack) {
    return {*this, itrack};
  }

  /// Get the size of the track container
  /// @return the sixe
  constexpr IndexType size() const {
    return m_container->size_impl();
  }

  /// Add a track to the container. Note this only creates the logical track and
  /// allocates memory. You can combine this with @c getTrack to obtain a track proxy
  /// @return the index to the newly added track
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  IndexType addTrack() {
    auto track = getTrack(m_container->addTrack_impl());
    track.tipIndex() = kInvalid;
    return track.index();
  }

  /// Remove a track at index @p itrack from the container
  /// @note This invalidates all track proxies!
  /// @param itrack The index of the track to remmove
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  void removeTrack(IndexType itrack) {
    m_container->removeTrack_impl(itrack);
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
  const auto& container() const {
    return *m_container;
  }

  /// Get a mutable reference to the track state container backend
  /// @return a mutable reference to the backend
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  auto& trackStateContainer() {
    return *m_traj;
  }

  /// Retrieve the holder of the track state container
  /// @return The track state container including it's holder
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  auto& trackStateContainerHolder() {
    return m_traj;
  }

  /// Get a const reference to the track state container backend
  /// @return a const reference to the backend
  const auto& trackStateContainer() const {
    return *m_traj;
  }

  /// Retrieve the holder of the track state container
  /// @return The track state container including it's holder
  const auto& trackStateContainerHolder() const {
    return m_traj;
  }

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

  /// Helper function to make this track container match the dynamic columns of
  /// another one. This will only work if the track container supports this
  /// source, and depends on the implementation details of the dynamic columns
  /// of the container
  /// @tparam other_track_container_t Type of the other track container
  /// @param other The other track container
  template <typename other_track_container_t, bool RO = ReadOnly,
            typename = std::enable_if_t<!RO>>
  void ensureDynamicColumns(const other_track_container_t& other) {
    container().ensureDynamicColumns_impl(other.container());
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
  template <typename T, bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  void copyDynamicFrom(IndexType dstIdx, const T& src, IndexType srcIdx) {
    container().copyDynamicFrom_impl(dstIdx, src, srcIdx);
  }

  detail_tc::ConstIf<holder_t<track_container_t>, ReadOnly> m_container;
  detail_tc::ConstIf<holder_t<traj_t>, ReadOnly> m_traj;
};

template <typename track_container_t, typename traj_t>
TrackContainer(track_container_t& container, traj_t& traj)
    -> TrackContainer<track_container_t, traj_t, detail::RefHolder>;

template <typename track_container_t, typename traj_t>
TrackContainer(const track_container_t& container, const traj_t& traj)
    -> TrackContainer<track_container_t, traj_t, detail::ConstRefHolder>;

template <typename track_container_t, typename traj_t>
TrackContainer(track_container_t&& container, traj_t&& traj)
    -> TrackContainer<track_container_t, traj_t, detail::ValueHolder>;

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
  T& operator()(track_proxy_t track) const {
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
  const T& operator()(track_proxy_t track) const {
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
