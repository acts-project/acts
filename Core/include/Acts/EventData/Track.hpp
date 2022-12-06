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

template <typename track_container_t, typename traj_t,
          template <typename> class holder_t>
class TrackContainer;

namespace detail_tc {
template <typename T, bool select>
using ConstIf = std::conditional_t<select, const T, T>;

template <typename track_container_t, typename trajectory_t,
          template <typename> class holder_t, bool ReadOnly = true>
class TrackProxy {
 public:
  using Container = track_container_t;
  using Trajectory = trajectory_t;

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

  friend TrackContainer<Container, Trajectory, holder_t>;

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

  IndexType tipIndex() const {
    return component<IndexType>(hashString("tipIndex"));
  }

  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  IndexType& tipIndex() {
    return component<IndexType>(hashString("tipIndex"));
  }

  const Surface& referenceSurface() const {
    return *component<std::shared_ptr<const Surface>,
                      hashString("referenceSurface")>();
  }

  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  void setReferenceSurface(std::shared_ptr<const Surface> srf) {
    component<std::shared_ptr<const Surface>,
              hashString("referenceSurface")>() = std::move(srf);
  }

  ConstParameters parameters() const {
    return m_container->parameters(m_index);
  }

  ConstCovariance covariance() const {
    return m_container->covariance(m_index);
  }

  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  Parameters parameters() {
    return m_container->parameters(m_index);
  }

  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  Covariance covariance() {
    return m_container->covariance(m_index);
  }

  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  auto trackStates() {
    return m_container->trackStateRange(m_index);
  }

  auto trackStates() const { return m_container->trackStateRange(m_index); }

  IndexType index() const { return m_index; }

  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  auto& container() {
    return *m_container;
  }

  const auto& container() const { return *m_container; }

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

template <typename container_t, typename proxy_t, bool ReadOnly>
class TrackProxyIterator {
  using ProxyType = proxy_t;
  using IndexType = typename ProxyType::IndexType;
  using ContainerType = container_t;

 public:
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

  bool operator==(const TrackProxyIterator& other) const {
    return m_container == other.m_container && m_itrack == other.m_itrack;
  }

  bool operator!=(const TrackProxyIterator& other) const {
    return !(*this == other);
  }

  ProxyType operator*() const { return m_container->getTrack(m_itrack); }

  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  ProxyType operator*() {
    return m_container->getTrack(m_itrack);
  }

 private:
  detail_lt::TransitiveConstPointer<ConstIf<ContainerType, ReadOnly>>
      m_container;
  IndexType m_itrack;
};

}  // namespace detail_tc

template <typename T>
struct isReadOnlyTrackContainer;

template <typename track_container_t, typename traj_t,
          template <typename> class holder_t = detail_tc::RefHolder>
class TrackContainer {
 public:
  static constexpr bool ReadOnly =
      isReadOnlyTrackContainer<track_container_t>::value;
  static constexpr bool TrackStateReadOnly =
      isReadOnlyMultiTrajectory<traj_t>::value;

  static_assert(ReadOnly == TrackStateReadOnly,
                "Either both track container and track state container need to "
                "be readonly or both have to be readwrite");

  using IndexType = MultiTrajectoryTraits::IndexType;
  static constexpr IndexType kInvalid = MultiTrajectoryTraits::kInvalid;

  using TrackProxy =
      detail_tc::TrackProxy<track_container_t, traj_t, holder_t, false>;
  using ConstTrackProxy =
      detail_tc::TrackProxy<track_container_t, traj_t, holder_t, true>;

  friend TrackProxy;
  friend ConstTrackProxy;

  TrackContainer(holder_t<track_container_t> container, holder_t<traj_t> traj)
      : m_container{std::move(container)}, m_traj{std::move(traj)} {}

  template <template <typename> class H = holder_t,
            typename = std::enable_if_t<
                detail_tc::is_same_template<H, detail_tc::RefHolder>::value>>
  TrackContainer(track_container_t& container, traj_t& traj)
      : m_container{&container}, m_traj{&traj} {}

  ConstTrackProxy getTrack(IndexType itrack) const { return {*this, itrack}; }

  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  TrackProxy getTrack(IndexType itrack) {
    return {*this, itrack};
  }

  constexpr IndexType size() const { return m_container->size_impl(); }

  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  IndexType addTrack() {
    return m_container->addTrack_impl();
  }

  template <typename T, bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  constexpr void addColumn(const std::string& key) {
    m_container->template addColumn_impl<T>(key);
  }

  constexpr bool hasColumn(HashedString key) const {
    return m_container->hasColumn_impl(key);
  }

  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  auto& container() {
    return *m_container;
  }

  const auto& container() const { return *m_container; }

  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  auto& trackStateContainer() {
    return *m_traj;
  }

  const auto& trackStateContainer() const { return *m_traj; }

  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  auto begin() {
    return detail_tc::TrackProxyIterator<std::decay_t<decltype(*this)>,
                                         TrackProxy, false>{*this, 0};
  }
  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  auto end() {
    return detail_tc::TrackProxyIterator<std::decay_t<decltype(*this)>,
                                         TrackProxy, false>{*this, size()};
  }

  auto begin() const {
    return detail_tc::TrackProxyIterator<std::decay_t<decltype(*this)>,
                                         ConstTrackProxy, true>{*this, 0};
  }

  auto end() const {
    return detail_tc::TrackProxyIterator<std::decay_t<decltype(*this)>,
                                         ConstTrackProxy, true>{*this, size()};
  }

 protected:
  template <typename T, HashedString key, bool RO = ReadOnly,
            typename = std::enable_if_t<!RO>>
  constexpr T& component(IndexType itrack) {
    return *std::any_cast<T*>(m_container->component_impl(key, itrack));
  }

  template <typename T, bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  constexpr T& component(HashedString key, IndexType istate) {
    return *std::any_cast<T*>(m_container->component_impl(key, istate));
  }

  template <typename T, HashedString key>
  constexpr const T& component(IndexType itrack) const {
    return *std::any_cast<const T*>(m_container->component_impl(key, itrack));
  }

  template <typename T>
  constexpr const T& component(HashedString key, IndexType itrack) const {
    return *std::any_cast<const T*>(m_container->component_impl(key, itrack));
  }

  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  constexpr typename TrackProxy::Parameters parameters(IndexType itrack) {
    return m_container->parameters(itrack);
  }

  constexpr typename TrackProxy::Parameters parameters(IndexType itrack) {
    return m_container->parameters(itrack);
  }

  constexpr typename ConstTrackProxy::Parameters parameters(
      IndexType itrack) const {
    return m_container->parameters(itrack);
  }

  constexpr typename TrackProxy::Covariance covariance(IndexType itrack) {
    return m_container->covariance(itrack);
  }

  constexpr typename ConstTrackProxy::Covariance covariance(
      IndexType itrack) const {
    return m_container->covariance(itrack);
  }

  template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
  auto trackStateRange(IndexType itrack) {
    auto tip = component<IndexType>(hashString("tipIndex"), itrack);
    return m_traj->trackStateRange(tip);
  }

  auto trackStateRange(IndexType itrack) const {
    auto tip = component<IndexType>(hashString("tipIndex"), itrack);
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

}  // namespace Acts
