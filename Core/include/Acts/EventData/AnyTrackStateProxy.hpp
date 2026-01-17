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
#include "Acts/EventData/SubspaceHelpers.hpp"
#include "Acts/EventData/TrackStatePropMask.hpp"
#include "Acts/EventData/TrackStateProxy.hpp"
#include "Acts/EventData/TrackStateProxyCommon.hpp"
#include "Acts/EventData/TrackStateProxyConcept.hpp"
#include "Acts/EventData/TrackStateType.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/HashedString.hpp"

#include <algorithm>
#include <any>
#include <cassert>
#include <memory>
#include <ranges>
#include <string_view>
#include <type_traits>
#include <utility>

namespace Acts {

template <bool read_only>
class AnyTrackStateProxy;

namespace detail_anytstate {

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
class TrackStateHandler<trajectory_t, true> final
    : public TrackStateHandlerConstBase {
  using MultiTrajectoryType =
      typename TrackStateHandlerTraits<trajectory_t>::MultiTrajectoryType;

 public:
  static const TrackStateHandler& instance() {
    static const TrackStateHandler s_instance;
    return s_instance;
  }

  TrackIndexType calibratedSize(const void* container,
                                TrackIndexType index) const override {
    assert(container != nullptr);
    const auto* traj = static_cast<const MultiTrajectoryType*>(container);
    return traj->calibratedSize(index);
  }

  ConstParametersMap parameters(const void* container,
                                TrackIndexType index) const override {
    assert(container != nullptr);
    const auto* traj = static_cast<const MultiTrajectoryType*>(container);
    return traj->parameters(index);
  }

  ConstCovarianceMap covariance(const void* container,
                                TrackIndexType index) const override {
    assert(container != nullptr);
    const auto* traj = static_cast<const MultiTrajectoryType*>(container);
    return traj->covariance(index);
  }

  const double* calibratedData(const void* container,
                               TrackIndexType index) const override {
    assert(container != nullptr);
    const auto* traj = static_cast<const MultiTrajectoryType*>(container);
    return traj->template calibrated<eBoundSize>(index).data();
  }

  const double* calibratedCovarianceData(const void* container,
                                         TrackIndexType index) const override {
    assert(container != nullptr);
    const auto* traj = static_cast<const MultiTrajectoryType*>(container);
    return traj->template calibratedCovariance<eBoundSize>(index).data();
  }

  const Surface* referenceSurface(const void* container,
                                  TrackIndexType index) const override {
    assert(container != nullptr);
    const auto* traj = static_cast<const MultiTrajectoryType*>(container);
    return traj->referenceSurface(index);
  }

  bool hasReferenceSurface(const void* container,
                           TrackIndexType index) const override {
    return referenceSurface(container, index) != nullptr;
  }

  SourceLink getUncalibratedSourceLink(const void* container,
                                       TrackIndexType index) const override {
    assert(container != nullptr);
    const auto* traj = static_cast<const MultiTrajectoryType*>(container);
    return traj->getUncalibratedSourceLink(index);
  }

  ConstCovarianceMap jacobian(const void* container,
                              TrackIndexType index) const override {
    assert(container != nullptr);
    const auto* traj = static_cast<const MultiTrajectoryType*>(container);
    return traj->jacobian(index);
  }

  bool has(const void* container, TrackIndexType index,
           HashedString key) const override {
    assert(container != nullptr);
    const auto* traj = static_cast<const MultiTrajectoryType*>(container);
    return traj->has(key, index);
  }

  std::any component(const void* container, TrackIndexType index,
                     HashedString key) const override {
    assert(container != nullptr);
    const auto* traj = static_cast<const MultiTrajectoryType*>(container);
    return traj->self().component_impl(key, index);
  }

  bool hasColumn(const void* container, HashedString key) const override {
    assert(container != nullptr);
    const auto* traj = static_cast<const MultiTrajectoryType*>(container);
    return traj->hasColumn(key);
  }

 private:
  TrackStateHandler() = default;
};

template <typename trajectory_t>
class TrackStateHandler<trajectory_t, false> final
    : public TrackStateHandlerMutableBase {
  using MultiTrajectoryType =
      typename TrackStateHandlerTraits<trajectory_t>::MultiTrajectoryType;

 public:
  static const TrackStateHandler& instance() {
    static const TrackStateHandler s_instance;
    return s_instance;
  }

  TrackIndexType calibratedSize(const void* container,
                                TrackIndexType index) const override {
    assert(container != nullptr);
    const auto* traj = static_cast<const MultiTrajectoryType*>(container);
    return traj->calibratedSize(index);
  }

  ConstParametersMap parameters(const void* container,
                                TrackIndexType index) const override {
    assert(container != nullptr);
    const auto* traj = static_cast<const MultiTrajectoryType*>(container);
    return traj->parameters(index);
  }

  ParametersMap parameters(void* container,
                           TrackIndexType index) const override {
    assert(container != nullptr);
    auto* traj = static_cast<MultiTrajectoryType*>(container);
    return traj->parameters(index);
  }

  ConstCovarianceMap covariance(const void* container,
                                TrackIndexType index) const override {
    assert(container != nullptr);
    const auto* traj = static_cast<const MultiTrajectoryType*>(container);
    return traj->covariance(index);
  }

  CovarianceMap covariance(void* container,
                           TrackIndexType index) const override {
    assert(container != nullptr);
    auto* traj = static_cast<MultiTrajectoryType*>(container);
    return traj->covariance(index);
  }

  const double* calibratedData(const void* container,
                               TrackIndexType index) const override {
    assert(container != nullptr);
    const auto* traj = static_cast<const MultiTrajectoryType*>(container);
    return traj->template calibrated<eBoundSize>(index).data();
  }

  const double* calibratedCovarianceData(const void* container,
                                         TrackIndexType index) const override {
    assert(container != nullptr);
    const auto* traj = static_cast<const MultiTrajectoryType*>(container);
    return traj->template calibratedCovariance<eBoundSize>(index).data();
  }

  double* calibratedDataMutable(void* container,
                                TrackIndexType index) const override {
    assert(container != nullptr);
    auto* traj = static_cast<MultiTrajectoryType*>(container);
    return traj->template calibrated<eBoundSize>(index).data();
  }

  double* calibratedCovarianceDataMutable(void* container,
                                          TrackIndexType index) const override {
    assert(container != nullptr);
    auto* traj = static_cast<MultiTrajectoryType*>(container);
    return traj->template calibratedCovariance<eBoundSize>(index).data();
  }

  const Surface* referenceSurface(const void* container,
                                  TrackIndexType index) const override {
    assert(container != nullptr);
    const auto* traj = static_cast<const MultiTrajectoryType*>(container);
    return traj->referenceSurface(index);
  }

  bool hasReferenceSurface(const void* container,
                           TrackIndexType index) const override {
    return referenceSurface(container, index) != nullptr;
  }

  SourceLink getUncalibratedSourceLink(const void* container,
                                       TrackIndexType index) const override {
    assert(container != nullptr);
    const auto* traj = static_cast<const MultiTrajectoryType*>(container);
    return traj->getUncalibratedSourceLink(index);
  }

  ConstCovarianceMap jacobian(const void* container,
                              TrackIndexType index) const override {
    assert(container != nullptr);
    const auto* traj = static_cast<const MultiTrajectoryType*>(container);
    return traj->jacobian(index);
  }

  CovarianceMap jacobian(void* container, TrackIndexType index) const override {
    assert(container != nullptr);
    auto* traj = static_cast<MultiTrajectoryType*>(container);
    return traj->jacobian(index);
  }

  bool has(const void* container, TrackIndexType index,
           HashedString key) const override {
    assert(container != nullptr);
    const auto* traj = static_cast<const MultiTrajectoryType*>(container);
    return traj->has(key, index);
  }

  std::any component(const void* container, TrackIndexType index,
                     HashedString key) const override {
    assert(container != nullptr);
    const auto* traj = static_cast<const MultiTrajectoryType*>(container);
    return traj->self().component_impl(key, index);
  }

  bool hasColumn(const void* container, HashedString key) const override {
    assert(container != nullptr);
    const auto* traj = static_cast<const MultiTrajectoryType*>(container);
    return traj->hasColumn(key);
  }

  std::any component(void* container, TrackIndexType index,
                     HashedString key) const override {
    assert(container != nullptr);
    auto* traj = static_cast<MultiTrajectoryType*>(container);
    return traj->self().component_impl(key, index);
  }

  void unset(void* container, TrackIndexType index,
             TrackStatePropMask target) const override {
    assert(container != nullptr);
    auto* traj = static_cast<MultiTrajectoryType*>(container);
    traj->unset(target, index);
  }

  void allocateCalibrated(void* container, TrackIndexType index,
                          std::size_t measdim) const override {
    assert(container != nullptr);
    auto* traj = static_cast<MultiTrajectoryType*>(container);
    traj->allocateCalibrated(index, measdim);
  }

  void setUncalibratedSourceLink(void* container, TrackIndexType index,
                                 SourceLink&& sourceLink) const override {
    assert(container != nullptr);
    auto* traj = static_cast<MultiTrajectoryType*>(container);
    traj->setUncalibratedSourceLink(index, std::move(sourceLink));
  }

  void setReferenceSurface(
      void* container, TrackIndexType index,
      std::shared_ptr<const Surface> surface) const override {
    assert(container != nullptr);
    auto* traj = static_cast<MultiTrajectoryType*>(container);
    traj->setReferenceSurface(index, std::move(surface));
  }

  void addTrackStateComponents(void* container, TrackIndexType index,
                               TrackStatePropMask mask) const override {
    assert(container != nullptr);
    auto* traj = static_cast<MultiTrajectoryType*>(container);
    traj->addTrackStateComponents(index, mask);
  }

 private:
  TrackStateHandler() = default;
};

}  // namespace detail_anytstate

template <bool read_only>
class AnyTrackStateProxy
    : public TrackStateProxyCommon<AnyTrackStateProxy<read_only>, read_only> {
  using Base = TrackStateProxyCommon<AnyTrackStateProxy<read_only>, read_only>;

  friend class TrackStateProxyCommon<AnyTrackStateProxy<read_only>, read_only>;

  using IndexType = Acts::TrackIndexType;

 public:
  static constexpr bool ReadOnly = read_only;

  using MutableTrackState = AnyTrackStateProxy<false>;
  using ConstTrackState = AnyTrackStateProxy<true>;
  using ConstProxyType = AnyTrackStateProxy<true>;

  using ParametersMap =
      detail_anytstate::TrackStateHandlerConstBase::ParametersMap;
  using ConstParametersMap =
      detail_anytstate::TrackStateHandlerConstBase::ConstParametersMap;
  using CovarianceMap =
      detail_anytstate::TrackStateHandlerConstBase::CovarianceMap;
  using ConstCovarianceMap =
      detail_anytstate::TrackStateHandlerConstBase::ConstCovarianceMap;
  using MutableEffectiveCalibratedMap =
      detail_anytstate::TrackStateHandlerConstBase::EffectiveCalibratedMap;
  using ConstEffectiveCalibratedMap =
      detail_anytstate::TrackStateHandlerConstBase::ConstEffectiveCalibratedMap;
  using MutableEffectiveCalibratedCovarianceMap = detail_anytstate::
      TrackStateHandlerConstBase::EffectiveCalibratedCovarianceMap;
  using ConstEffectiveCalibratedCovarianceMap = detail_anytstate::
      TrackStateHandlerConstBase::ConstEffectiveCalibratedCovarianceMap;

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
    requires(ReadOnly || !track_state_proxy_t::ReadOnly)
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
    assert(has(detail_tsp::kUncalibratedKey));
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
  ConstParametersMap parametersAtIndex(IndexType parIndex) const {
    return constHandler()->parameters(containerPtr(), parIndex);
  }

  ParametersMap parametersAtIndexMutable(IndexType parIndex) const
    requires(!ReadOnly)
  {
    return mutableHandler()->parameters(mutableContainerPtr(), parIndex);
  }

  ConstCovarianceMap covarianceAtIndex(IndexType covIndex) const {
    return constHandler()->covariance(containerPtr(), covIndex);
  }

  CovarianceMap covarianceAtIndexMutable(IndexType covIndex) const
    requires(!ReadOnly)
  {
    return mutableHandler()->covariance(mutableContainerPtr(), covIndex);
  }

  double* calibratedDataMutable()
    requires(!ReadOnly)
  {
    return mutableHandler()->calibratedDataMutable(mutableContainerPtr(),
                                                   m_index);
  }

  const double* calibratedData() const {
    return constHandler()->calibratedData(containerPtr(), m_index);
  }

  double* calibratedCovarianceDataMutable()
    requires(!ReadOnly)
  {
    return mutableHandler()->calibratedCovarianceDataMutable(
        mutableContainerPtr(), m_index);
  }

  const double* calibratedCovarianceData() const {
    return constHandler()->calibratedCovarianceData(containerPtr(), m_index);
  }

 private:
  template <bool>
  friend class AnyTrackStateProxy;

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

using AnyConstTrackStateProxy = AnyTrackStateProxy<true>;
using AnyMutableTrackStateProxy = AnyTrackStateProxy<false>;

static_assert(ConstTrackStateProxyConcept<AnyConstTrackStateProxy>);
static_assert(MutableTrackStateProxyConcept<AnyMutableTrackStateProxy>);

}  // namespace Acts
