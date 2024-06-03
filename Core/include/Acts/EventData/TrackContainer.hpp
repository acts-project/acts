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
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/TrackContainerBackendConcept.hpp"
#include "Acts/EventData/TrackProxy.hpp"
#include "Acts/EventData/Types.hpp"
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
/// tracks. It uses a backend to store both the actual tracks and the
/// associated track states.
/// @tparam track_container_t the track container backend
/// @tparam traj_t the track state container backend
/// @tparam holder_t ownership management class for the backend
template <ACTS_CONCEPT(TrackContainerBackend) track_container_t,
          typename traj_t,
          template <typename> class holder_t = detail::RefHolder>
class TrackContainer {
 public:
  /// Indicates if this track container is read-only, or if it can be modified
  static constexpr bool ReadOnly =
      IsReadOnlyTrackContainer<track_container_t>::value;

  /// Indicates if the track state container is read-only, or if it can be
  /// modified
  static constexpr bool TrackStateReadOnly =
      IsReadOnlyMultiTrajectory<traj_t>::value;

  static_assert(ReadOnly == TrackStateReadOnly,
                "Either both track container and track state container need to "
                "be readonly or both have to be readwrite");

  /// The index type of the track container, taken from the container backend
  using IndexType = TrackIndexType;

  /// Sentinel value that indicates an invalid index
  static constexpr IndexType kInvalid = kTrackIndexInvalid;

  /// Alias for the mutable version of a track proxy, with the same backends as
  /// this container
  using TrackProxy =
      Acts::TrackProxy<track_container_t, traj_t, holder_t, false>;

  /// Alias for the const version of a track proxy, with the same backends as
  /// this container
  using ConstTrackProxy =
      Acts::TrackProxy<track_container_t, traj_t, holder_t, true>;

#ifndef DOXYGEN
  friend TrackProxy;
  friend ConstTrackProxy;
#endif

  /// @anchor track_container_construction
  /// @name TrackContainer construction
  ///
  /// Constructors for the track container by using a set of backends
  /// (track + track state). The container can either take ownership of the
  /// backends or just hold references to them. This is driven by the @c
  /// holder_t template parameter, where you can also supply a custom holder
  /// type. Template deduction is used to try to guess the correct holder type.
  ///
  /// @{

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

  /// @}

  /// @anchor track_container_track_access
  /// @name TrackContainer track (proxy) access and manipulation
  ///
  /// These methods allow accessing tracks, i.e. adding or retrieving a track
  /// proxy that points at a specific track in the container.
  ///
  /// @{

  /// Get a const track proxy for a track index
  /// @param itrack the track index in the container
  /// @return A const track proxy for the index
  ConstTrackProxy getTrack(IndexType itrack) const {
    return {*this, itrack};
  }

  /// Get a mutable track proxy for a track index
  /// @note Only available if the track container is not read-only
  /// @param itrack the track index in the container
  /// @return A mutable track proxy for the index
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  TrackProxy getTrack(IndexType itrack) {
    return {*this, itrack};
  }

  /// Add a track to the container. Note this only creates the logical track and
  /// allocates memory. You can combine this with @c getTrack to obtain a track proxy
  /// @note Only available if the track container is not read-only
  /// @return the index to the newly added track
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  IndexType addTrack() {
    auto track = getTrack(m_container->addTrack_impl());
    track.tipIndex() = kInvalid;
    return track.index();
  }

  /// Add a track to the container and return a track proxy to it
  /// This effectively calls @c addTrack and @c getTrack
  /// @note Only available if the track container is not read-only
  /// @return a track proxy to the newly added track
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  TrackProxy makeTrack() {
    return getTrack(addTrack());
  }

  /// Remove a track at index @p itrack from the container
  /// @note Only available if the track container is not read-only
  /// @note This invalidates track proxies that point to tracks with larger
  ///       indices than @p itrack!
  /// @param itrack The index of the track to remove
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  void removeTrack(IndexType itrack) {
    m_container->removeTrack_impl(itrack);
  }

  /// Get a mutable iterator to the first track in the container
  /// @note Only available if the track container is not read-only
  /// @return a mutable iterator to the first track
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  auto begin() {
    return detail_tc::TrackProxyIterator<std::decay_t<decltype(*this)>,
                                         TrackProxy, false>{*this, 0};
  }

  /// Get a past-the-end iterator for this container
  /// @note Only available if the track container is not read-only
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

  /// @}

  /// @anchor track_container_columns
  /// @name TrackContainer column management
  /// TrackContainer can manage a set of common static columns, and dynamic
  /// columns that can be added at runtime. This set of methods allows you to
  /// manage the dynamic columns.
  /// @{

  /// Add a dymanic column to the track container
  /// @note Only available if the track container is not read-only
  /// @param key the name of the column to be added
  template <typename T, bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  constexpr void addColumn(std::string_view key) {
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

  /// Helper function to make this track container match the dynamic columns of
  /// another one. This will only work if the track container supports this
  /// source, and depends on the implementation details of the dynamic columns
  /// of the container
  /// @note Only available if the track container is not read-only
  /// @tparam other_track_container_t Type of the other track container
  /// @param other The other track container
  template <typename other_track_container_t, bool RO = ReadOnly,
            typename = std::enable_if_t<!RO>>
  void ensureDynamicColumns(const other_track_container_t& other) {
    container().ensureDynamicColumns_impl(other.container());
  }

  /// @}

  /// @anchor track_congtainer_backend_access
  /// @name TrackContainer backend access
  /// These methods allow accessing the backend of the track container. In most
  /// cases, this is not necessary for interacting with the container.
  /// @{

  /// Get a mutable reference to the track container backend
  /// @note Only available if the track container is not read-only
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
  /// @note Only available if the track container is not read-only
  /// @return a mutable reference to the backend
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  auto& trackStateContainer() {
    return *m_traj;
  }

  /// Retrieve the holder of the track state container
  /// @return The track state container including it's holder
  /// @note Only available if the track container is not read-only
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

  /// @}

  /// Get the size (number of tracks) of the track container
  /// @return the sixe
  constexpr IndexType size() const {
    return m_container->size_impl();
  }

  /// Clear the content of the track container
  /// @note Only available if the track container is not read-only
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  void clear() {
    m_container->clear();
    m_traj->clear();
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
  auto reverseTrackStateRange(IndexType itrack) {
    auto tip = component<IndexType, hashString("tipIndex")>(itrack);
    return m_traj->reverseTrackStateRange(tip);
  }

  auto reverseTrackStateRange(IndexType itrack) const {
    auto tip = component<IndexType, hashString("tipIndex")>(itrack);
    return m_traj->reverseTrackStateRange(tip);
  }

  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  auto forwardTrackStateRange(IndexType itrack) {
    auto stem = component<IndexType, hashString("stemIndex")>(itrack);
    if (stem == kInvalid) {
      throw std::invalid_argument{"Track has no stem index"};
    }
    return m_traj->forwardTrackStateRange(stem);
  }

  auto forwardTrackStateRange(IndexType itrack) const {
    auto stem = component<IndexType, hashString("stemIndex")>(itrack);
    if (stem == kInvalid) {
      throw std::invalid_argument{"Track has no stem index"};
    }
    return m_traj->forwardTrackStateRange(stem);
  }

 private:
  template <typename T, bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  void copyDynamicFrom(IndexType dstIdx, const T& src, IndexType srcIdx) {
    const auto& dynamicKeys = src.dynamicKeys_impl();
    for (const auto key : dynamicKeys) {
      std::any srcPtr = src.component_impl(key, srcIdx);
      container().copyDynamicFrom_impl(dstIdx, key, srcPtr);
    }
  }

  detail_tc::ConstIf<holder_t<track_container_t>, ReadOnly> m_container;
  detail_tc::ConstIf<holder_t<traj_t>, ReadOnly> m_traj;
};

template <ACTS_CONCEPT(TrackContainerBackend) track_container_t,
          typename traj_t>
TrackContainer(track_container_t& container, traj_t& traj)
    -> TrackContainer<track_container_t, traj_t, detail::RefHolder>;

template <ACTS_CONCEPT(TrackContainerBackend) track_container_t,
          typename traj_t>
TrackContainer(const track_container_t& container, const traj_t& traj)
    -> TrackContainer<track_container_t, traj_t, detail::ConstRefHolder>;

template <ACTS_CONCEPT(TrackContainerBackend) track_container_t,
          typename traj_t>
TrackContainer(track_container_t&& container, traj_t&& traj)
    -> TrackContainer<track_container_t, traj_t, detail::ValueHolder>;

}  // namespace Acts
