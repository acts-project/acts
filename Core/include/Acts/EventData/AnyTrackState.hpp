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
#include "Acts/EventData/TrackStateProxyConcept.hpp"
#include "Acts/EventData/TrackStateType.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/HashedString.hpp"
#include "Acts/Utilities/Helpers.hpp"

#include <algorithm>
#include <any>
#include <bitset>
#include <cassert>
#include <memory>
#include <optional>
#include <ranges>
#include <stdexcept>
#include <string_view>
#include <type_traits>
#include <utility>

namespace Acts {

template <bool read_only>
class AnyTrackState;

namespace detail_anytstate {

using TrackStateParametersMap =
    typename detail_lt::FixedSizeTypes<eBoundSize, false>::CoefficientsMap;
using TrackStateConstParametersMap =
    typename detail_lt::FixedSizeTypes<eBoundSize, true>::CoefficientsMap;
using TrackStateCovarianceMap =
    typename detail_lt::FixedSizeTypes<eBoundSize, false>::CovarianceMap;
using TrackStateConstCovarianceMap =
    typename detail_lt::FixedSizeTypes<eBoundSize, true>::CovarianceMap;
using TrackStateEffectiveCalibratedMap =
    typename detail_lt::DynamicSizeTypes<false>::CoefficientsMap;
using TrackStateConstEffectiveCalibratedMap =
    typename detail_lt::DynamicSizeTypes<true>::CoefficientsMap;
using TrackStateEffectiveCalibratedCovarianceMap =
    typename detail_lt::DynamicSizeTypes<false>::CovarianceMap;
using TrackStateConstEffectiveCalibratedCovarianceMap =
    typename detail_lt::DynamicSizeTypes<true>::CovarianceMap;

class TrackStateHandlerConstBase {
 public:
  virtual ~TrackStateHandlerConstBase() = default;

  virtual TrackIndexType index(TrackIndexType state) const { return state; }

  virtual TrackIndexType calibratedSize(const void* container,
                                        TrackIndexType index) const = 0;

  virtual TrackStateConstParametersMap parameters(
      const void* container, TrackIndexType parametersIndex) const = 0;

  virtual TrackStateConstCovarianceMap covariance(
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

  virtual TrackStateConstCovarianceMap jacobian(const void* container,
                                                TrackIndexType index) const = 0;
};

class TrackStateHandlerMutableBase : public TrackStateHandlerConstBase {
 public:
  using TrackStateHandlerConstBase::component;
  using TrackStateHandlerConstBase::covariance;
  using TrackStateHandlerConstBase::jacobian;
  using TrackStateHandlerConstBase::parameters;

  virtual TrackStateParametersMap parameters(
      void* container, TrackIndexType parametersIndex) const = 0;

  virtual TrackStateCovarianceMap covariance(
      void* container, TrackIndexType covarianceIndex) const = 0;

  virtual double* calibratedDataMutable(void* container,
                                        TrackIndexType index) const = 0;

  virtual double* calibratedCovarianceDataMutable(
      void* container, TrackIndexType index) const = 0;

  virtual std::any component(void* container, TrackIndexType index,
                             HashedString key) const = 0;

  virtual TrackStateCovarianceMap jacobian(void* container,
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

  virtual void copyDynamicFrom(void* container, TrackIndexType self,
                               const void* otherContainer,
                               TrackIndexType other) const = 0;
};

template <typename trajectory_t>
class TrackStateHandler final : public TrackStateHandlerMutableBase {
  using MultiTrajectoryType = MultiTrajectory<trajectory_t>;

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

  TrackStateConstParametersMap parameters(const void* container,
                                          TrackIndexType index) const override {
    assert(container != nullptr);
    const auto* traj = static_cast<const MultiTrajectoryType*>(container);
    return traj->parameters(index);
  }

  TrackStateParametersMap parameters(void* container,
                                     TrackIndexType index) const override {
    assert(container != nullptr);
    auto* traj = static_cast<MultiTrajectoryType*>(container);
    return traj->parameters(index);
  }

  TrackStateConstCovarianceMap covariance(const void* container,
                                          TrackIndexType index) const override {
    assert(container != nullptr);
    const auto* traj = static_cast<const MultiTrajectoryType*>(container);
    return traj->covariance(index);
  }

  TrackStateCovarianceMap covariance(void* container,
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

  TrackStateConstCovarianceMap jacobian(const void* container,
                                        TrackIndexType index) const override {
    assert(container != nullptr);
    const auto* traj = static_cast<const MultiTrajectoryType*>(container);
    return traj->jacobian(index);
  }

  TrackStateCovarianceMap jacobian(void* container,
                                   TrackIndexType index) const override {
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

  void copyDynamicFrom(void* container, TrackIndexType self,
                       const void* otherContainer,
                       TrackIndexType other) const override {
    assert(container != nullptr);
    assert(otherContainer != nullptr);
    auto* traj = static_cast<MultiTrajectoryType*>(container);
    const auto* otherTraj =
        static_cast<const MultiTrajectoryType*>(otherContainer);
    traj->copyDynamicFrom(self, *otherTraj, other);
  }

 private:
  TrackStateHandler() = default;
};

}  // namespace detail_anytstate

template <bool read_only>
class AnyTrackState {
 public:
  static constexpr bool ReadOnly = read_only;

  using MutableTrackState = AnyTrackState<false>;
  using ConstTrackState = AnyTrackState<true>;
  using ConstProxyType = AnyTrackState<true>;

  using ParametersMap = detail_anytstate::TrackStateParametersMap;
  using ConstParametersMap = detail_anytstate::TrackStateConstParametersMap;
  using CovarianceMap = detail_anytstate::TrackStateCovarianceMap;
  using ConstCovarianceMap = detail_anytstate::TrackStateConstCovarianceMap;
  using MutableEffectiveCalibratedMap =
      detail_anytstate::TrackStateEffectiveCalibratedMap;
  using ConstEffectiveCalibratedMap =
      detail_anytstate::TrackStateConstEffectiveCalibratedMap;
  using MutableEffectiveCalibratedCovarianceMap =
      detail_anytstate::TrackStateEffectiveCalibratedCovarianceMap;
  using ConstEffectiveCalibratedCovarianceMap =
      detail_anytstate::TrackStateConstEffectiveCalibratedCovarianceMap;

  using ContainerPointer = std::conditional_t<ReadOnly, const void*, void*>;

  template <TrackStateProxyConcept track_state_proxy_t>
    requires(ReadOnly || !track_state_proxy_t::ReadOnly)
  explicit AnyTrackState(const track_state_proxy_t& ts)
      : m_container(nullptr), m_index(ts.m_istate), m_handler(nullptr) {
    using trajectory_t = typename track_state_proxy_t::Trajectory;
    auto* containerPtr = ts.rawTrajectoryPtr();
    if constexpr (ReadOnly) {
      m_container = static_cast<const void*>(containerPtr);
    } else {
      m_container = static_cast<void*>(containerPtr);
    }
    m_handler = &detail_anytstate::TrackStateHandler<trajectory_t>::instance();
  }

  TrackIndexType index() const { return m_index; }

  template <HashedString key>
  bool has() const {
    return constHandler()->has(containerPtr(), m_index, key);
  }

  bool has(HashedString key) const {
    return constHandler()->has(containerPtr(), m_index, key);
  }

  bool has(std::string_view key) const { return has(hashStringDynamic(key)); }

  bool hasColumn(HashedString key) const {
    return constHandler()->hasColumn(containerPtr(), key);
  }

  bool hasColumn(std::string_view key) const {
    return hasColumn(hashStringDynamic(key));
  }

  template <typename T, HashedString key>
  const T& component() const {
    std::any result = constHandler()->component(containerPtr(), m_index, key);
    return *std::any_cast<const T*>(result);
  }

  template <typename T>
  const T& component(HashedString key) const {
    std::any result = constHandler()->component(containerPtr(), m_index, key);
    return *std::any_cast<const T*>(result);
  }

  template <typename T>
  const T& component(std::string_view key) const {
    return component<T>(hashStringDynamic(key));
  }

  template <typename T, HashedString key>
  T& component()
    requires(!ReadOnly)
  {
    std::any result =
        mutableHandler()->component(mutableContainerPtr(), m_index, key);
    return *std::any_cast<T*>(result);
  }

  template <typename T>
  T& component(HashedString key)
    requires(!ReadOnly)
  {
    std::any result =
        mutableHandler()->component(mutableContainerPtr(), m_index, key);
    return *std::any_cast<T*>(result);
  }

  template <typename T>
  T& component(std::string_view key)
    requires(!ReadOnly)
  {
    return component<T>(hashStringDynamic(key));
  }

  TrackIndexType previous() const {
    return component<TrackIndexType, detail_tsp::kPreviousKey>();
  }

  TrackIndexType& previous()
    requires(!ReadOnly)
  {
    return component<TrackIndexType, detail_tsp::kPreviousKey>();
  }

  bool hasPrevious() const { return previous() != kTrackIndexInvalid; }

  TrackStatePropMask getMask() const {
    using PM = TrackStatePropMask;
    PM mask = PM::None;
    if (hasPredicted()) {
      mask |= PM::Predicted;
    }
    if (hasFiltered()) {
      mask |= PM::Filtered;
    }
    if (hasSmoothed()) {
      mask |= PM::Smoothed;
    }
    if (hasJacobian()) {
      mask |= PM::Jacobian;
    }
    if (hasCalibrated()) {
      mask |= PM::Calibrated;
    }
    return mask;
  }

  const Surface& referenceSurface() const {
    assert(hasReferenceSurface());
    const Surface* surface =
        constHandler()->referenceSurface(containerPtr(), m_index);
    assert(surface != nullptr);
    return *surface;
  }

  bool hasReferenceSurface() const {
    return constHandler()->hasReferenceSurface(containerPtr(), m_index);
  }

  void setReferenceSurface(std::shared_ptr<const Surface> surface)
    requires(!ReadOnly)
  {
    mutableHandler()->setReferenceSurface(mutableContainerPtr(), m_index,
                                          std::move(surface));
  }

  SourceLink getUncalibratedSourceLink() const {
    assert(has(detail_tsp::kUncalibratedKey));
    return constHandler()->getUncalibratedSourceLink(containerPtr(), m_index);
  }

  void setUncalibratedSourceLink(SourceLink sourceLink)
    requires(!ReadOnly)
  {
    mutableHandler()->setUncalibratedSourceLink(mutableContainerPtr(), m_index,
                                                std::move(sourceLink));
  }

  bool hasPredicted() const { return has(detail_tsp::kPredictedKey); }

  bool hasFiltered() const { return has(detail_tsp::kFilteredKey); }

  bool hasSmoothed() const { return has(detail_tsp::kSmoothedKey); }

  bool hasJacobian() const { return has(detail_tsp::kJacobianKey); }

  bool hasProjector() const { return has(detail_tsp::kProjectorKey); }

  bool hasCalibrated() const { return has(detail_tsp::kCalibratedKey); }

  TrackIndexType calibratedSize() const {
    return constHandler()->calibratedSize(containerPtr(), m_index);
  }

  void allocateCalibrated(std::size_t measdim)
    requires(!ReadOnly)
  {
    mutableHandler()->allocateCalibrated(mutableContainerPtr(), m_index,
                                         measdim);
    component<TrackIndexType, detail_tsp::kMeasDimKey>() =
        static_cast<TrackIndexType>(measdim);
  }

  template <typename val_t, typename cov_t>
  void allocateCalibrated(const Eigen::DenseBase<val_t>& val,
                          const Eigen::DenseBase<cov_t>& cov)
    requires(!ReadOnly && Concepts::eigen_base_is_fixed_size<val_t> &&
             Concepts::eigen_bases_have_same_num_rows<val_t, cov_t> &&
             Concepts::eigen_base_is_square<cov_t> &&
             Eigen::PlainObjectBase<val_t>::RowsAtCompileTime <=
                 static_cast<std::underlying_type_t<BoundIndices>>(eBoundSize))
  {
    constexpr std::size_t measdim =
        static_cast<std::size_t>(val_t::RowsAtCompileTime);
    allocateCalibrated(measdim);
    calibrated<measdim>() = val;
    calibratedCovariance<measdim>() = cov;
  }

  ConstEffectiveCalibratedMap effectiveCalibrated() const {
    const double* data =
        constHandler()->calibratedData(containerPtr(), m_index);
    const auto size = calibratedSize();
    return ConstEffectiveCalibratedMap(data, size);
  }

  MutableEffectiveCalibratedMap effectiveCalibrated()
    requires(!ReadOnly)
  {
    double* data =
        mutableHandler()->calibratedDataMutable(mutableContainerPtr(), m_index);
    const auto size = calibratedSize();
    return MutableEffectiveCalibratedMap(data, size);
  }

  ConstEffectiveCalibratedCovarianceMap effectiveCalibratedCovariance() const {
    const double* data =
        constHandler()->calibratedCovarianceData(containerPtr(), m_index);
    const auto size = calibratedSize();
    return ConstEffectiveCalibratedCovarianceMap(data, size, size);
  }

  MutableEffectiveCalibratedCovarianceMap effectiveCalibratedCovariance()
    requires(!ReadOnly)
  {
    double* data = mutableHandler()->calibratedCovarianceDataMutable(
        mutableContainerPtr(), m_index);
    const auto size = calibratedSize();
    return MutableEffectiveCalibratedCovarianceMap(data, size, size);
  }

  BoundSubspaceIndices projectorSubspaceIndices() const {
    assert(hasProjector());
    const auto& serialized =
        component<SerializedSubspaceIndices, detail_tsp::kProjectorKey>();
    return deserializeSubspaceIndices<eBoundSize>(serialized);
  }

  template <std::ranges::sized_range index_range_t>
    requires(!ReadOnly)
  void setProjectorSubspaceIndices(const index_range_t& indices) {
    component<SerializedSubspaceIndices, detail_tsp::kProjectorKey>() =
        serializeSubspaceIndices<eBoundSize>(indices);
    component<TrackIndexType, detail_tsp::kMeasDimKey>() =
        static_cast<TrackIndexType>(indices.size());
  }

  ConstParametersMap parameters() const {
    if (hasSmoothed()) {
      return smoothed();
    } else if (hasFiltered()) {
      return filtered();
    }
    return predicted();
  }

  ConstCovarianceMap covariance() const {
    if (hasSmoothed()) {
      return smoothedCovariance();
    } else if (hasFiltered()) {
      return filteredCovariance();
    }
    return predictedCovariance();
  }

  ConstCovarianceMap jacobian() const {
    assert(hasJacobian());
    return constHandler()->jacobian(containerPtr(), m_index);
  }

  CovarianceMap jacobian()
    requires(!ReadOnly)
  {
    assert(hasJacobian());
    return mutableHandler()->jacobian(mutableContainerPtr(), m_index);
  }

  double pathLength() const {
    return component<double, detail_tsp::kPathLengthKey>();
  }

  double& pathLength()
    requires(!ReadOnly)
  {
    return component<double, detail_tsp::kPathLengthKey>();
  }

  ConstParametersMap predicted() const {
    assert(hasPredicted());
    return parametersFromComponent(detail_tsp::kPredictedKey);
  }

  ParametersMap predicted()
    requires(!ReadOnly)
  {
    return mutableParametersFromComponent(detail_tsp::kPredictedKey);
  }

  ConstCovarianceMap predictedCovariance() const {
    assert(hasPredicted());
    return covarianceFromComponent(detail_tsp::kPredictedKey);
  }

  CovarianceMap predictedCovariance()
    requires(!ReadOnly)
  {
    return mutableCovarianceFromComponent(detail_tsp::kPredictedKey);
  }

  ConstParametersMap filtered() const {
    assert(hasFiltered());
    return parametersFromComponent(detail_tsp::kFilteredKey);
  }

  ParametersMap filtered()
    requires(!ReadOnly)
  {
    return mutableParametersFromComponent(detail_tsp::kFilteredKey);
  }

  ConstCovarianceMap filteredCovariance() const {
    assert(hasFiltered());
    return covarianceFromComponent(detail_tsp::kFilteredKey);
  }

  CovarianceMap filteredCovariance()
    requires(!ReadOnly)
  {
    return mutableCovarianceFromComponent(detail_tsp::kFilteredKey);
  }

  ConstParametersMap smoothed() const {
    assert(hasSmoothed());
    return parametersFromComponent(detail_tsp::kSmoothedKey);
  }

  ParametersMap smoothed()
    requires(!ReadOnly)
  {
    return mutableParametersFromComponent(detail_tsp::kSmoothedKey);
  }

  ConstCovarianceMap smoothedCovariance() const {
    assert(hasSmoothed());
    return covarianceFromComponent(detail_tsp::kSmoothedKey);
  }

  CovarianceMap smoothedCovariance()
    requires(!ReadOnly)
  {
    return mutableCovarianceFromComponent(detail_tsp::kSmoothedKey);
  }

  ConstTrackStateType typeFlags() const {
    const auto raw =
        component<TrackStateType::raw_type, detail_tsp::kTypeFlagsKey>();
    return ConstTrackStateType{raw};
  }

  TrackStateType typeFlags()
    requires(!ReadOnly)
  {
    auto& raw =
        component<TrackStateType::raw_type, detail_tsp::kTypeFlagsKey>();
    return TrackStateType{raw};
  }

  float chi2() const { return component<float, detail_tsp::kChi2Key>(); }

  float& chi2()
    requires(!ReadOnly)
  {
    return component<float, detail_tsp::kChi2Key>();
  }

  template <std::size_t measdim>
  typename TrackStateTraits<measdim, true>::Calibrated calibrated() const {
    assert(calibratedSize() == static_cast<TrackIndexType>(measdim));
    const double* data =
        constHandler()->calibratedData(containerPtr(), m_index);
    return typename TrackStateTraits<measdim, true>::Calibrated(data);
  }

  template <std::size_t measdim>
  typename TrackStateTraits<measdim, false>::Calibrated calibrated()
    requires(!ReadOnly)
  {
    assert(calibratedSize() == static_cast<TrackIndexType>(measdim));
    double* data =
        mutableHandler()->calibratedDataMutable(mutableContainerPtr(), m_index);
    return typename TrackStateTraits<measdim, false>::Calibrated(data);
  }

  template <std::size_t measdim>
  typename TrackStateTraits<measdim, true>::CalibratedCovariance
  calibratedCovariance() const {
    assert(calibratedSize() == static_cast<TrackIndexType>(measdim));
    const double* data =
        constHandler()->calibratedCovarianceData(containerPtr(), m_index);
    return typename TrackStateTraits<measdim, true>::CalibratedCovariance(data);
  }

  template <std::size_t measdim>
  typename TrackStateTraits<measdim, false>::CalibratedCovariance
  calibratedCovariance()
    requires(!ReadOnly)
  {
    assert(calibratedSize() == static_cast<TrackIndexType>(measdim));
    double* data = mutableHandler()->calibratedCovarianceDataMutable(
        mutableContainerPtr(), m_index);
    return
        typename TrackStateTraits<measdim, false>::CalibratedCovariance(data);
  }

  void unset(TrackStatePropMask target)
    requires(!ReadOnly)
  {
    mutableHandler()->unset(mutableContainerPtr(), m_index, target);
  }

 private:
  template <bool>
  friend class AnyTrackState;

  TrackIndexType componentIndexValue(HashedString key) const {
    assert(has(key));
    const auto& ref = component<TrackIndexType>(key);
    return ref;
  }

  ConstParametersMap parametersFromComponent(HashedString key) const {
    const auto parIndex = componentIndexValue(key);
    return constHandler()->parameters(containerPtr(), parIndex);
  }

  ParametersMap mutableParametersFromComponent(HashedString key) const
    requires(!ReadOnly)
  {
    const auto parIndex = componentIndexValue(key);
    return mutableHandler()->parameters(mutableContainerPtr(), parIndex);
  }

  ConstCovarianceMap covarianceFromComponent(HashedString key) const {
    const auto covIndex = componentIndexValue(key);
    return constHandler()->covariance(containerPtr(), covIndex);
  }

  CovarianceMap mutableCovarianceFromComponent(HashedString key) const
    requires(!ReadOnly)
  {
    const auto covIndex = componentIndexValue(key);
    return mutableHandler()->covariance(mutableContainerPtr(), covIndex);
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

  const void* containerPtr() const {
    if constexpr (ReadOnly) {
      return m_container;
    } else {
      return static_cast<const void*>(m_container);
    }
  }

  void* mutableContainerPtr() const
    requires(!ReadOnly)
  {
    return m_container;
  }

  ContainerPointer m_container;
  TrackIndexType m_index;
  const detail_anytstate::TrackStateHandlerConstBase* m_handler;
};

using AnyConstTrackState = AnyTrackState<true>;
using AnyMutableTrackState = AnyTrackState<false>;

static_assert(ConstTrackStateProxyConcept<AnyConstTrackState>);
static_assert(MutableTrackStateProxyConcept<AnyMutableTrackState>);

}  // namespace Acts
