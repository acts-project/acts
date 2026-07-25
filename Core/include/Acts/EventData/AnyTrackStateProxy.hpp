// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/TrackStatePropMask.hpp"
#include "Acts/EventData/TrackStateProxy.hpp"
#include "Acts/EventData/TrackStateProxyCommon.hpp"
#include "Acts/EventData/TrackStateProxyConcept.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/HashedString.hpp"

#include <any>
#include <cassert>
#include <memory>
#include <string_view>
#include <type_traits>
#include <utility>

namespace Acts {

/// Type-erased track state proxy for any trajectory backend.
/// @tparam read_only True for const access, false for mutable access.
template <bool read_only>
class AnyTrackStateProxy;

namespace detail_anytstate {

/// Type-erased range over the track states of a track (defined below, after
/// AnyTrackStateProxy).
template <bool reverse, bool read_only>
class AnyTrackStateRange;

class TrackStateHandlerConstBase {
 public:
  using ParametersMap =
      typename detail_tsp::FixedSizeTypes<eBoundSize, false>::CoefficientsMap;
  using ConstParametersMap =
      typename detail_tsp::FixedSizeTypes<eBoundSize, true>::CoefficientsMap;
  using CovarianceMap =
      typename detail_tsp::FixedSizeTypes<eBoundSize, false>::CovarianceMap;
  using ConstCovarianceMap =
      typename detail_tsp::FixedSizeTypes<eBoundSize, true>::CovarianceMap;
  using EffectiveCalibratedMap =
      typename detail_tsp::DynamicSizeTypes<false>::CoefficientsMap;
  using ConstEffectiveCalibratedMap =
      typename detail_tsp::DynamicSizeTypes<true>::CoefficientsMap;
  using EffectiveCalibratedCovarianceMap =
      typename detail_tsp::DynamicSizeTypes<false>::CovarianceMap;
  using ConstEffectiveCalibratedCovarianceMap =
      typename detail_tsp::DynamicSizeTypes<true>::CovarianceMap;

  virtual ~TrackStateHandlerConstBase() = default;

  virtual TrackIndexType index(TrackIndexType state) const { return state; }

  virtual TrackIndexType calibratedSize(const void* container,
                                        TrackIndexType index) const = 0;

  virtual ConstParametersMap parameters(
      const void* container, TrackIndexType parametersIndex) const = 0;

  virtual ConstCovarianceMap covariance(
      const void* container, TrackIndexType covarianceIndex) const = 0;

  virtual const double* calibratedData(const void* container,
                                       TrackIndexType index) const = 0;

  virtual const double* calibratedCovarianceData(
      const void* container, TrackIndexType index) const = 0;

  virtual bool has(const void* container, TrackIndexType index,
                   HashedString key) const = 0;

  virtual std::any component(const void* container, TrackIndexType index,
                             HashedString key) const = 0;

  virtual bool hasColumn(const void* container, HashedString key) const = 0;

  virtual const Surface* referenceSurface(const void* container,
                                          TrackIndexType index) const = 0;

  virtual bool hasReferenceSurface(const void* container,
                                   TrackIndexType index) const = 0;

  virtual bool hasUncalibratedSourceLink(const void* container,
                                         TrackIndexType index) const = 0;

  virtual SourceLink getUncalibratedSourceLink(const void* container,
                                               TrackIndexType index) const = 0;

  virtual ConstCovarianceMap jacobian(const void* container,
                                      TrackIndexType index) const = 0;
};

class TrackStateHandlerMutableBase : public TrackStateHandlerConstBase {
 public:
  using TrackStateHandlerConstBase::component;
  using TrackStateHandlerConstBase::covariance;
  using TrackStateHandlerConstBase::jacobian;
  using TrackStateHandlerConstBase::parameters;

  virtual ParametersMap parameters(void* container,
                                   TrackIndexType parametersIndex) const = 0;

  virtual CovarianceMap covariance(void* container,
                                   TrackIndexType covarianceIndex) const = 0;

  virtual double* calibratedDataMutable(void* container,
                                        TrackIndexType index) const = 0;

  virtual double* calibratedCovarianceDataMutable(
      void* container, TrackIndexType index) const = 0;

  virtual std::any component(void* container, TrackIndexType index,
                             HashedString key) const = 0;

  virtual CovarianceMap jacobian(void* container,
                                 TrackIndexType index) const = 0;

  virtual void unset(void* container, TrackIndexType index,
                     TrackStatePropMask target) const = 0;

  virtual void allocateCalibrated(void* container, TrackIndexType index,
                                  std::size_t measdim) const = 0;

  virtual void setUncalibratedSourceLink(void* container, TrackIndexType index,
                                         SourceLink&& sourceLink) const = 0;

  virtual void setReferenceSurface(
      void* container, TrackIndexType index,
      std::shared_ptr<const Surface> surface) const = 0;

  virtual void addTrackStateComponents(void* container, TrackIndexType index,
                                       TrackStatePropMask mask) const = 0;
};

template <typename trajectory_t>
struct TrackStateHandlerTraits {
  using Trajectory = std::remove_const_t<trajectory_t>;
  using MultiTrajectoryType = MultiTrajectory<Trajectory>;
  static constexpr bool ReadOnly =
      std::is_const_v<trajectory_t> || MultiTrajectoryType::ReadOnly;
};

template <typename trajectory_t,
          bool read_only = TrackStateHandlerTraits<trajectory_t>::ReadOnly>
class TrackStateHandler;

template <typename trajectory_t>
class TrackStateHandler<trajectory_t, true> /*final*/
    : public TrackStateHandlerConstBase {
  using MultiTrajectoryType =
      typename TrackStateHandlerTraits<trajectory_t>::MultiTrajectoryType;

 public:
  static const TrackStateHandler& instance() {
    static const TrackStateHandler s_instance;
    return s_instance;
  }

  TrackIndexType calibratedSize(const void* container,
                                TrackIndexType index) const final {
    assert(container != nullptr);
    const auto* traj = static_cast<const MultiTrajectoryType*>(container);
    return traj->calibratedSize(index);
  }

  ConstParametersMap parameters(const void* container,
                                TrackIndexType index) const final {
    assert(container != nullptr);
    const auto* traj = static_cast<const MultiTrajectoryType*>(container);
    return traj->parameters(index);
  }

  ConstCovarianceMap covariance(const void* container,
                                TrackIndexType index) const final {
    assert(container != nullptr);
    const auto* traj = static_cast<const MultiTrajectoryType*>(container);
    return traj->covariance(index);
  }

  const double* calibratedData(const void* container,
                               TrackIndexType index) const final {
    assert(container != nullptr);
    const auto* traj = static_cast<const MultiTrajectoryType*>(container);
    return traj->template calibrated<eBoundSize>(index).data();
  }

  const double* calibratedCovarianceData(const void* container,
                                         TrackIndexType index) const final {
    assert(container != nullptr);
    const auto* traj = static_cast<const MultiTrajectoryType*>(container);
    return traj->template calibratedCovariance<eBoundSize>(index).data();
  }

  const Surface* referenceSurface(const void* container,
                                  TrackIndexType index) const final {
    assert(container != nullptr);
    const auto* traj = static_cast<const MultiTrajectoryType*>(container);
    return traj->referenceSurface(index);
  }

  bool hasReferenceSurface(const void* container,
                           TrackIndexType index) const final {
    return referenceSurface(container, index) != nullptr;
  }

  bool hasUncalibratedSourceLink(const void* container,
                                 TrackIndexType index) const final {
    assert(container != nullptr);
    const auto* traj = static_cast<const MultiTrajectoryType*>(container);
    return traj->has(hashString("uncalibratedSourceLink"), index);
  }

  SourceLink getUncalibratedSourceLink(const void* container,
                                       TrackIndexType index) const final {
    assert(container != nullptr);
    const auto* traj = static_cast<const MultiTrajectoryType*>(container);
    return traj->getUncalibratedSourceLink(index);
  }

  ConstCovarianceMap jacobian(const void* container,
                              TrackIndexType index) const final {
    assert(container != nullptr);
    const auto* traj = static_cast<const MultiTrajectoryType*>(container);
    return traj->jacobian(index);
  }

  bool has(const void* container, TrackIndexType index,
           HashedString key) const final {
    assert(container != nullptr);
    const auto* traj = static_cast<const MultiTrajectoryType*>(container);
    return traj->has(key, index);
  }

  std::any component(const void* container, TrackIndexType index,
                     HashedString key) const final {
    assert(container != nullptr);
    const auto* traj = static_cast<const MultiTrajectoryType*>(container);
    return traj->self().component_impl(key, index);
  }

  bool hasColumn(const void* container, HashedString key) const final {
    assert(container != nullptr);
    const auto* traj = static_cast<const MultiTrajectoryType*>(container);
    return traj->hasColumn(key);
  }

 private:
  TrackStateHandler() = default;
};

template <typename trajectory_t>
class TrackStateHandler<trajectory_t, false> /*final*/
    : public TrackStateHandlerMutableBase {
  using MultiTrajectoryType =
      typename TrackStateHandlerTraits<trajectory_t>::MultiTrajectoryType;

 public:
  static const TrackStateHandler& instance() {
    static const TrackStateHandler s_instance;
    return s_instance;
  }

  TrackIndexType calibratedSize(const void* container,
                                TrackIndexType index) const final {
    assert(container != nullptr);
    const auto* traj = static_cast<const MultiTrajectoryType*>(container);
    return traj->calibratedSize(index);
  }

  ConstParametersMap parameters(const void* container,
                                TrackIndexType index) const final {
    assert(container != nullptr);
    const auto* traj = static_cast<const MultiTrajectoryType*>(container);
    return traj->parameters(index);
  }

  ParametersMap parameters(void* container, TrackIndexType index) const final {
    assert(container != nullptr);
    auto* traj = static_cast<MultiTrajectoryType*>(container);
    return traj->parameters(index);
  }

  ConstCovarianceMap covariance(const void* container,
                                TrackIndexType index) const final {
    assert(container != nullptr);
    const auto* traj = static_cast<const MultiTrajectoryType*>(container);
    return traj->covariance(index);
  }

  CovarianceMap covariance(void* container, TrackIndexType index) const final {
    assert(container != nullptr);
    auto* traj = static_cast<MultiTrajectoryType*>(container);
    return traj->covariance(index);
  }

  const double* calibratedData(const void* container,
                               TrackIndexType index) const final {
    assert(container != nullptr);
    const auto* traj = static_cast<const MultiTrajectoryType*>(container);
    return traj->template calibrated<eBoundSize>(index).data();
  }

  const double* calibratedCovarianceData(const void* container,
                                         TrackIndexType index) const final {
    assert(container != nullptr);
    const auto* traj = static_cast<const MultiTrajectoryType*>(container);
    return traj->template calibratedCovariance<eBoundSize>(index).data();
  }

  double* calibratedDataMutable(void* container,
                                TrackIndexType index) const final {
    assert(container != nullptr);
    auto* traj = static_cast<MultiTrajectoryType*>(container);
    return traj->template calibrated<eBoundSize>(index).data();
  }

  double* calibratedCovarianceDataMutable(void* container,
                                          TrackIndexType index) const final {
    assert(container != nullptr);
    auto* traj = static_cast<MultiTrajectoryType*>(container);
    return traj->template calibratedCovariance<eBoundSize>(index).data();
  }

  const Surface* referenceSurface(const void* container,
                                  TrackIndexType index) const final {
    assert(container != nullptr);
    const auto* traj = static_cast<const MultiTrajectoryType*>(container);
    return traj->referenceSurface(index);
  }

  bool hasReferenceSurface(const void* container,
                           TrackIndexType index) const final {
    return referenceSurface(container, index) != nullptr;
  }

  bool hasUncalibratedSourceLink(const void* container,
                                 TrackIndexType index) const final {
    assert(container != nullptr);
    const auto* traj = static_cast<const MultiTrajectoryType*>(container);
    return traj->has(hashString("uncalibratedSourceLink"), index);
  }

  SourceLink getUncalibratedSourceLink(const void* container,
                                       TrackIndexType index) const final {
    assert(container != nullptr);
    const auto* traj = static_cast<const MultiTrajectoryType*>(container);
    return traj->getUncalibratedSourceLink(index);
  }

  ConstCovarianceMap jacobian(const void* container,
                              TrackIndexType index) const final {
    assert(container != nullptr);
    const auto* traj = static_cast<const MultiTrajectoryType*>(container);
    return traj->jacobian(index);
  }

  CovarianceMap jacobian(void* container, TrackIndexType index) const final {
    assert(container != nullptr);
    auto* traj = static_cast<MultiTrajectoryType*>(container);
    return traj->jacobian(index);
  }

  bool has(const void* container, TrackIndexType index,
           HashedString key) const final {
    assert(container != nullptr);
    const auto* traj = static_cast<const MultiTrajectoryType*>(container);
    return traj->has(key, index);
  }

  std::any component(const void* container, TrackIndexType index,
                     HashedString key) const final {
    assert(container != nullptr);
    const auto* traj = static_cast<const MultiTrajectoryType*>(container);
    return traj->self().component_impl(key, index);
  }

  bool hasColumn(const void* container, HashedString key) const final {
    assert(container != nullptr);
    const auto* traj = static_cast<const MultiTrajectoryType*>(container);
    return traj->hasColumn(key);
  }

  std::any component(void* container, TrackIndexType index,
                     HashedString key) const final {
    assert(container != nullptr);
    auto* traj = static_cast<MultiTrajectoryType*>(container);
    return traj->self().component_impl(key, index);
  }

  void unset(void* container, TrackIndexType index,
             TrackStatePropMask target) const final {
    assert(container != nullptr);
    auto* traj = static_cast<MultiTrajectoryType*>(container);
    traj->unset(target, index);
  }

  void allocateCalibrated(void* container, TrackIndexType index,
                          std::size_t measdim) const final {
    assert(container != nullptr);
    auto* traj = static_cast<MultiTrajectoryType*>(container);
    traj->allocateCalibrated(index, measdim);
  }

  void setUncalibratedSourceLink(void* container, TrackIndexType index,
                                 SourceLink&& sourceLink) const final {
    assert(container != nullptr);
    auto* traj = static_cast<MultiTrajectoryType*>(container);
    traj->setUncalibratedSourceLink(index, std::move(sourceLink));
  }

  void setReferenceSurface(void* container, TrackIndexType index,
                           std::shared_ptr<const Surface> surface) const final {
    assert(container != nullptr);
    auto* traj = static_cast<MultiTrajectoryType*>(container);
    traj->setReferenceSurface(index, std::move(surface));
  }

  void addTrackStateComponents(void* container, TrackIndexType index,
                               TrackStatePropMask mask) const final {
    assert(container != nullptr);
    auto* traj = static_cast<MultiTrajectoryType*>(container);
    traj->addTrackStateComponents(index, mask);
  }

 private:
  TrackStateHandler() = default;
};

}  // namespace detail_anytstate

/// Type-erased track state proxy for any trajectory backend.
///
/// Provides the TrackStateProxy interface by delegating to a runtime handler.
/// @tparam read_only True for const access, false for mutable access.
template <bool read_only>
class AnyTrackStateProxy
    : public TrackStateProxyCommon<AnyTrackStateProxy<read_only>, read_only> {
  using Base = TrackStateProxyCommon<AnyTrackStateProxy<read_only>, read_only>;

  friend class TrackStateProxyCommon<AnyTrackStateProxy<read_only>, read_only>;

  using IndexType = Acts::TrackIndexType;

 public:
  /// Whether this proxy provides read-only or mutable access
  static constexpr bool ReadOnly = read_only;

  /// Mutable track state proxy type
  using MutableTrackState = AnyTrackStateProxy<false>;
  /// Const track state proxy type
  using ConstTrackState = AnyTrackStateProxy<true>;
  /// Const proxy type alias
  using ConstProxyType = AnyTrackStateProxy<true>;

  /// Mutable parameters map type
  using ParametersMap =
      detail_anytstate::TrackStateHandlerConstBase::ParametersMap;
  /// Const parameters map type
  using ConstParametersMap =
      detail_anytstate::TrackStateHandlerConstBase::ConstParametersMap;
  /// Mutable covariance map type
  using CovarianceMap =
      detail_anytstate::TrackStateHandlerConstBase::CovarianceMap;
  /// Const covariance map type
  using ConstCovarianceMap =
      detail_anytstate::TrackStateHandlerConstBase::ConstCovarianceMap;
  /// Mutable effective calibrated map type
  using MutableEffectiveCalibratedMap =
      detail_anytstate::TrackStateHandlerConstBase::EffectiveCalibratedMap;
  /// Const effective calibrated map type
  using ConstEffectiveCalibratedMap =
      detail_anytstate::TrackStateHandlerConstBase::ConstEffectiveCalibratedMap;
  /// Mutable effective calibrated covariance map type
  using MutableEffectiveCalibratedCovarianceMap = detail_anytstate::
      TrackStateHandlerConstBase::EffectiveCalibratedCovarianceMap;
  /// Const effective calibrated covariance map type
  using ConstEffectiveCalibratedCovarianceMap = detail_anytstate::
      TrackStateHandlerConstBase::ConstEffectiveCalibratedCovarianceMap;

  /// Container pointer type
  using ContainerPointer = std::conditional_t<ReadOnly, const void*, void*>;

  using Base::allocateCalibrated;
  using Base::calibrated;
  using Base::calibratedCovariance;
  using Base::chi2;
  using Base::covariance;
  using Base::effectiveCalibrated;
  using Base::effectiveCalibratedCovariance;
  using Base::filtered;
  using Base::filteredCovariance;
  using Base::getMask;
  using Base::hasCalibrated;
  using Base::hasFiltered;
  using Base::hasJacobian;
  using Base::hasPredicted;
  using Base::hasPrevious;
  using Base::hasProjector;
  using Base::hasSmoothed;
  using Base::hasUncalibratedSourceLink;
  using Base::parameters;
  using Base::pathLength;
  using Base::predicted;
  using Base::predictedCovariance;
  using Base::previous;
  using Base::projectorSubspaceHelper;
  using Base::projectorSubspaceIndices;
  using Base::setProjectorSubspaceIndices;
  using Base::smoothed;
  using Base::smoothedCovariance;
  using Base::typeFlags;

  /// Construct an `AnyTrackStateProxy` from a concrete track-state proxy.
  /// @tparam track_state_proxy_t Proxy type satisfying the concept.
  /// @param ts Proxy that supplies the trajectory backend and index.
  template <TrackStateProxyConcept track_state_proxy_t>
    requires((ReadOnly || !track_state_proxy_t::ReadOnly) &&
             !std::is_same_v<std::remove_cv_t<track_state_proxy_t>,
                             AnyTrackStateProxy<true>> &&
             !std::is_same_v<std::remove_cv_t<track_state_proxy_t>,
                             AnyTrackStateProxy<false>>)
  explicit AnyTrackStateProxy(track_state_proxy_t& ts)
      : m_container(nullptr), m_index(ts.m_istate) {
    using trajectory_t = typename track_state_proxy_t::Trajectory;
    auto* containerPtr = ts.rawTrajectoryPtr();
    if constexpr (ReadOnly) {
      m_container = static_cast<const void*>(containerPtr);
    } else {
      m_container = static_cast<void*>(containerPtr);
    }
    m_handler = &detail_anytstate::TrackStateHandler<trajectory_t>::instance();
  }

  /// Get the index of the underlying track state.
  /// @return Track state index within its container.
  TrackIndexType index() const { return m_index; }

  /// Check if a compile-time keyed component exists on this track state.
  /// @tparam key Component key encoded as a hashed string literal.
  /// @return True if the component is available.
  template <HashedString key>
  bool has() const {
    return constHandler()->has(containerPtr(), m_index, key);
  }

  /// Check if a hashed component exists on this track state.
  /// @param key Component identifier.
  /// @return True if the component is available.
  bool has(HashedString key) const {
    return constHandler()->has(containerPtr(), m_index, key);
  }

  /// Check if a string-named component exists on this track state.
  /// @param key Component identifier as a string.
  /// @return True if the component is available.
  bool has(std::string_view key) const { return has(hashStringDynamic(key)); }

  /// Check if the trajectory container exposes a column.
  /// @param key Column identifier.
  /// @return True if the column exists.
  bool hasColumn(HashedString key) const {
    return constHandler()->hasColumn(containerPtr(), key);
  }

  /// Check if the trajectory container exposes a column.
  /// @param key Column identifier as a string.
  /// @return True if the column exists.
  bool hasColumn(std::string_view key) const {
    return hasColumn(hashStringDynamic(key));
  }

  /// Access a const component through a compile-time key.
  /// @tparam T Component type.
  /// @tparam key Component key encoded as a hashed string literal.
  /// @return Const reference to the requested component.
  template <typename T, HashedString key>
  const T& component() const {
    std::any result = constHandler()->component(containerPtr(), m_index, key);
    return *std::any_cast<const T*>(result);
  }

  /// Access a const component by hashed key.
  /// @tparam T Component type.
  /// @param key Component identifier.
  /// @return Const reference to the requested component.
  template <typename T>
  const T& component(HashedString key) const {
    std::any result = constHandler()->component(containerPtr(), m_index, key);
    return *std::any_cast<const T*>(result);
  }

  /// Access a const component by string key.
  /// @tparam T Component type.
  /// @param key Component identifier as a string.
  /// @return Const reference to the requested component.
  template <typename T>
  const T& component(std::string_view key) const {
    return component<T>(hashStringDynamic(key));
  }

  /// Access a mutable component through a compile-time key.
  /// @tparam T Component type.
  /// @tparam key Component key encoded as a hashed string literal.
  /// @return Mutable reference to the requested component.
  template <typename T, HashedString key>
  T& component()
    requires(!ReadOnly)
  {
    std::any result =
        mutableHandler()->component(mutableContainerPtr(), m_index, key);
    return *std::any_cast<T*>(result);
  }

  /// Access a mutable component by hashed key.
  /// @tparam T Component type.
  /// @param key Component identifier.
  /// @return Mutable reference to the requested component.
  template <typename T>
  T& component(HashedString key)
    requires(!ReadOnly)
  {
    std::any result =
        mutableHandler()->component(mutableContainerPtr(), m_index, key);
    return *std::any_cast<T*>(result);
  }

  /// Access a mutable component by string key.
  /// @tparam T Component type.
  /// @param key Component identifier as a string.
  /// @return Mutable reference to the requested component.
  template <typename T>
  T& component(std::string_view key)
    requires(!ReadOnly)
  {
    return component<T>(hashStringDynamic(key));
  }

  /// Access the surface the state is referenced to.
  /// @return Reference surface of the track state.
  const Surface& referenceSurface() const {
    assert(hasReferenceSurface());
    const Surface* surface =
        constHandler()->referenceSurface(containerPtr(), m_index);
    assert(surface != nullptr);
    return *surface;
  }

  /// Check whether a reference surface is attached.
  /// @return True if a valid reference surface exists.
  bool hasReferenceSurface() const {
    return constHandler()->hasReferenceSurface(containerPtr(), m_index);
  }

  /// Assign a new reference surface to the track state.
  /// @param surface Surface that should be referenced.
  void setReferenceSurface(std::shared_ptr<const Surface> surface)
    requires(!ReadOnly)
  {
    mutableHandler()->setReferenceSurface(mutableContainerPtr(), m_index,
                                          std::move(surface));
  }

  /// Retrieve the original, uncalibrated source link.
  /// @return Copy of the stored source link.
  SourceLink getUncalibratedSourceLink() const {
    assert(hasUncalibratedSourceLink());
    return constHandler()->getUncalibratedSourceLink(containerPtr(), m_index);
  }

  /// Store an uncalibrated source link on this state.
  /// @param sourceLink Source link to copy into the track state.
  void setUncalibratedSourceLink(SourceLink sourceLink)
    requires(!ReadOnly)
  {
    mutableHandler()->setUncalibratedSourceLink(mutableContainerPtr(), m_index,
                                                std::move(sourceLink));
  }

  /// Retrieve the measurement dimension of the calibrated data.
  /// @return Number of calibrated measurement entries.
  TrackIndexType calibratedSize() const {
    return constHandler()->calibratedSize(containerPtr(), m_index);
  }

  /// Allocate memory for runtime-dimension calibrated data.
  /// @param measdim Number of measurement rows to reserve.
  void allocateCalibrated(std::size_t measdim)
    requires(!ReadOnly)
  {
    mutableHandler()->allocateCalibrated(mutableContainerPtr(), m_index,
                                         measdim);
    component<TrackIndexType, detail_tsp::kMeasDimKey>() =
        static_cast<TrackIndexType>(measdim);
  }

  /// Access the transport Jacobian.
  /// @return Const map referencing the Jacobian matrix.
  ConstCovarianceMap jacobian() const {
    assert(hasJacobian());
    return constHandler()->jacobian(containerPtr(), m_index);
  }

  /// Access the transport Jacobian.
  /// @return Mutable map referencing the Jacobian matrix.
  CovarianceMap jacobian()
    requires(!ReadOnly)
  {
    assert(hasJacobian());
    return mutableHandler()->jacobian(mutableContainerPtr(), m_index);
  }

  /// Remove dynamic components according to a mask.
  /// @param target Property mask describing which components to drop.
  void unset(TrackStatePropMask target)
    requires(!ReadOnly)
  {
    mutableHandler()->unset(mutableContainerPtr(), m_index, target);
  }

 protected:
  /// Access const parameters at specific index
  /// @param parIndex Parameter index
  /// @return Const parameters map
  ConstParametersMap parametersAtIndex(IndexType parIndex) const {
    return constHandler()->parameters(containerPtr(), parIndex);
  }

  /// Access mutable parameters at specific index
  /// @param parIndex Parameter index
  /// @return Mutable parameters map
  ParametersMap parametersAtIndexMutable(IndexType parIndex) const
    requires(!ReadOnly)
  {
    return mutableHandler()->parameters(mutableContainerPtr(), parIndex);
  }

  /// Access const covariance at specific index
  /// @param covIndex Covariance index
  /// @return Const covariance map
  ConstCovarianceMap covarianceAtIndex(IndexType covIndex) const {
    return constHandler()->covariance(containerPtr(), covIndex);
  }

  /// Access mutable covariance at specific index
  /// @param covIndex Covariance index
  /// @return Mutable covariance map
  CovarianceMap covarianceAtIndexMutable(IndexType covIndex) const
    requires(!ReadOnly)
  {
    return mutableHandler()->covariance(mutableContainerPtr(), covIndex);
  }

  /// Access mutable calibrated measurement data
  /// @return Pointer to mutable calibrated data
  double* calibratedDataMutable()
    requires(!ReadOnly)
  {
    return mutableHandler()->calibratedDataMutable(mutableContainerPtr(),
                                                   m_index);
  }

  /// Access const calibrated measurement data
  /// @return Pointer to const calibrated data
  const double* calibratedData() const {
    return constHandler()->calibratedData(containerPtr(), m_index);
  }

  /// Access mutable calibrated measurement covariance data
  /// @return Pointer to mutable calibrated covariance data
  double* calibratedCovarianceDataMutable()
    requires(!ReadOnly)
  {
    return mutableHandler()->calibratedCovarianceDataMutable(
        mutableContainerPtr(), m_index);
  }

  /// Access const calibrated measurement covariance data
  /// @return Pointer to const calibrated covariance data
  const double* calibratedCovarianceData() const {
    return constHandler()->calibratedCovarianceData(containerPtr(), m_index);
  }

 private:
  template <bool>
  friend class AnyTrackStateProxy;
  template <bool, bool>
  friend class detail_anytstate::AnyTrackStateRange;

  /// Create a sibling proxy pointing at a different state index within the
  /// same trajectory container. Used by the track-state ranges to walk the
  /// `previous()` / `next()` links.
  /// @param idx The state index to point at.
  /// @return A copy of this proxy rebound to @p idx.
  AnyTrackStateProxy withIndex(TrackIndexType idx) const {
    AnyTrackStateProxy copy(*this);
    copy.m_index = idx;
    return copy;
  }

  const detail_anytstate::TrackStateHandlerConstBase* constHandler() const {
    return m_handler;
  }

  const detail_anytstate::TrackStateHandlerMutableBase* mutableHandler() const
    requires(!ReadOnly)
  {
    return static_cast<const detail_anytstate::TrackStateHandlerMutableBase*>(
        m_handler);
  }

  const void* containerPtr() const { return m_container; }

  void* mutableContainerPtr() const
    requires(!ReadOnly)
  {
    return m_container;
  }

  ContainerPointer m_container{};
  TrackIndexType m_index{};
  const detail_anytstate::TrackStateHandlerConstBase* m_handler{};
};

namespace detail_anytstate {

/// Type-erased range over the track states of a track.
///
/// Mirrors @c detail_lt::TrackStateRange but yields type-erased
/// @c AnyTrackStateProxy objects. The range walks the `previous()` links when
/// @p reverse is true (outside-in) and the `next()` links otherwise
/// (inside-out). A default-constructed range is empty.
///
/// @tparam reverse True to iterate from the tip inwards, false to iterate from
///                 the stem outwards.
/// @tparam read_only True for const access, false for mutable access.
template <bool reverse, bool read_only>
class AnyTrackStateRange {
  using ProxyType = AnyTrackStateProxy<read_only>;

 public:
  /// Forward iterator over the track states. A `nullopt` proxy marks the
  /// past-the-end iterator.
  struct Iterator {
    /// The track state the iterator currently points at, or `nullopt` for end.
    std::optional<ProxyType> proxy;

    using iterator_category = std::forward_iterator_tag;
    using value_type = ProxyType;
    using difference_type = std::ptrdiff_t;
    using pointer = void;
    using reference = void;

    /// Advance to the next track state along the link direction.
    /// @return Reference to this iterator.
    Iterator& operator++() {
      if (!proxy) {
        return *this;
      }
      if constexpr (reverse) {
        if (proxy->hasPrevious()) {
          proxy = proxy->withIndex(proxy->previous());
        } else {
          proxy = std::nullopt;
        }
      } else {
        TrackIndexType next =
            proxy->template component<TrackIndexType, hashString("next")>();
        if (next != kTrackIndexInvalid) {
          proxy = proxy->withIndex(next);
        } else {
          proxy = std::nullopt;
        }
      }
      return *this;
    }

    /// Post-increment.
    /// @return Copy of the iterator before advancing.
    Iterator operator++(int) {
      Iterator tmp(*this);
      operator++();
      return tmp;
    }

    /// Compare two iterators. Two iterators are equal if both are past-the-end
    /// or both point at the same track state index.
    /// @param other The iterator to compare against.
    /// @return True if the iterators are equal.
    bool operator==(const Iterator& other) const {
      if (!proxy && !other.proxy) {
        return true;
      }
      if (proxy && other.proxy) {
        return proxy->index() == other.proxy->index();
      }
      return false;
    }

    /// Dereference to the current track state proxy.
    /// @return The current track state proxy.
    ProxyType operator*() const { return *proxy; }
  };

  /// Construct a range starting at @p begin.
  /// @param begin The first track state of the range.
  explicit AnyTrackStateRange(ProxyType begin) : m_begin{begin} {}
  /// Construct an empty range.
  AnyTrackStateRange() : m_begin{std::nullopt} {}

  /// @return Iterator to the first track state.
  Iterator begin() { return Iterator{m_begin}; }
  /// @return Past-the-end iterator.
  Iterator end() { return Iterator{std::nullopt}; }

 private:
  Iterator m_begin;
};

/// Const range iterating from the tip inwards (reverse).
using AnyConstReverseTrackStateRange = AnyTrackStateRange<true, true>;
/// Mutable range iterating from the tip inwards (reverse).
using AnyMutableReverseTrackStateRange = AnyTrackStateRange<true, false>;
/// Const range iterating from the stem outwards (forward).
using AnyConstTrackStateRange = AnyTrackStateRange<false, true>;
/// Mutable range iterating from the stem outwards (forward).
using AnyMutableTrackStateRange = AnyTrackStateRange<false, false>;

}  // namespace detail_anytstate

/// Type alias for a const track state proxy
using AnyConstTrackStateProxy = AnyTrackStateProxy<true>;
/// Type alias for a mutable track state proxy
using AnyMutableTrackStateProxy = AnyTrackStateProxy<false>;

static_assert(ConstTrackStateProxyConcept<AnyConstTrackStateProxy>);
static_assert(MutableTrackStateProxyConcept<AnyMutableTrackStateProxy>);

}  // namespace Acts
