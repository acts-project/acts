// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/SubspaceHelpers.hpp"
#include "Acts/EventData/SourceLink.hpp"
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

struct AnyTrackStateTrajectoryTag {
  using ConstTrackStateProxy = AnyTrackState<true>;
  using TrackStateProxy = AnyTrackState<false>;
};

namespace detail {

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

  virtual TrackStateConstEffectiveCalibratedMap effectiveCalibrated(
      const void* container, TrackIndexType index) const = 0;

  virtual TrackStateConstEffectiveCalibratedCovarianceMap
  effectiveCalibratedCovariance(const void* container,
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
  using TrackStateHandlerConstBase::jacobian;
  using TrackStateHandlerConstBase::covariance;
  using TrackStateHandlerConstBase::parameters;

  virtual TrackStateParametersMap parameters(void* container,
                                             TrackIndexType parametersIndex) const = 0;

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

  virtual void shareFrom(void* container, TrackIndexType self,
                         TrackIndexType other, TrackStatePropMask shareSource,
                         TrackStatePropMask shareTarget) const = 0;

  virtual void unset(void* container, TrackIndexType index,
                     TrackStatePropMask target) const = 0;

  virtual void allocateCalibrated(void* container, TrackIndexType index,
                                  std::size_t measdim) const = 0;

  virtual void setUncalibratedSourceLink(void* container, TrackIndexType index,
                                         SourceLink&& sourceLink) const = 0;

  virtual void setReferenceSurface(void* container, TrackIndexType index,
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

  TrackStateConstEffectiveCalibratedMap effectiveCalibrated(
      const void* container, TrackIndexType index) const override {
    assert(container != nullptr);
    const auto* traj = static_cast<const MultiTrajectoryType*>(container);
    return traj->effectiveCalibrated(index);
  }

  TrackStateConstEffectiveCalibratedCovarianceMap effectiveCalibratedCovariance(
      const void* container, TrackIndexType index) const override {
    assert(container != nullptr);
    const auto* traj = static_cast<const MultiTrajectoryType*>(container);
    return traj->effectiveCalibratedCovariance(index);
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

  void shareFrom(void* container, TrackIndexType self, TrackIndexType other,
                 TrackStatePropMask shareSource,
                 TrackStatePropMask shareTarget) const override {
    assert(container != nullptr);
    auto* traj = static_cast<MultiTrajectoryType*>(container);
    traj->shareFrom(self, other, shareSource, shareTarget);
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

  void setReferenceSurface(void* container, TrackIndexType index,
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

}  // namespace detail

template <bool read_only>
class AnyTrackState {
 public:
  static constexpr bool ReadOnly = read_only;

  using MutableTrackState = AnyTrackState<false>;
  using ConstTrackState = AnyTrackState<true>;
  using ConstProxyType = AnyTrackState<true>;
  using Trajectory = AnyTrackStateTrajectoryTag;

  using ParametersMap = detail::TrackStateParametersMap;
  using ConstParametersMap = detail::TrackStateConstParametersMap;
  using CovarianceMap = detail::TrackStateCovarianceMap;
  using ConstCovarianceMap = detail::TrackStateConstCovarianceMap;
  using MutableEffectiveCalibratedMap =
      detail::TrackStateEffectiveCalibratedMap;
  using ConstEffectiveCalibratedMap =
      detail::TrackStateConstEffectiveCalibratedMap;
  using MutableEffectiveCalibratedCovarianceMap =
      detail::TrackStateEffectiveCalibratedCovarianceMap;
  using ConstEffectiveCalibratedCovarianceMap =
      detail::TrackStateConstEffectiveCalibratedCovarianceMap;

  using ContainerPointer = std::conditional_t<ReadOnly, const void*, void*>;

  class ContainerView {
   public:
    ContainerView() = default;
    ContainerView(const detail::TrackStateHandlerConstBase* handler,
                  const void* container)
        : m_handler(handler), m_container(container) {}

    bool hasColumn(HashedString key) const {
      assert(m_handler != nullptr && m_container != nullptr);
      return m_handler->hasColumn(m_container, key);
    }

   private:
    friend class AnyTrackState;
    const detail::TrackStateHandlerConstBase* m_handler = nullptr;
    const void* m_container = nullptr;
  };

  AnyTrackState()
      : m_container(nullptr), m_index(kTrackIndexInvalid), m_handler(nullptr) {}

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
    m_handler = &detail::TrackStateHandler<trajectory_t>::instance();
  }

  bool isValid() const {
    return m_container != nullptr && m_handler != nullptr &&
           m_index != kTrackIndexInvalid;
  }

  explicit operator bool() const { return isValid(); }

  TrackIndexType index() const {
    assert(isValid());
    return m_index;
  }

  template <HashedString key>
  bool has() const {
    assert(isValid());
    return constHandler()->has(containerPtr(), m_index, key);
  }

  bool has(HashedString key) const {
    assert(isValid());
    return constHandler()->has(containerPtr(), m_index, key);
  }

  bool has(std::string_view key) const {
    return has(hashStringDynamic(key));
  }

  bool hasColumn(HashedString key) const {
    assert(isValid());
    return constHandler()->hasColumn(containerPtr(), key);
  }

  bool hasColumn(std::string_view key) const {
    return hasColumn(hashStringDynamic(key));
  }

  template <typename T, HashedString key>
  const T& component() const {
    assert(isValid());
    std::any result = constHandler()->component(containerPtr(), m_index, key);
    return *std::any_cast<const T*>(result);
  }

  template <typename T>
  const T& component(HashedString key) const {
    assert(isValid());
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
    assert(isValid());
    std::any result =
        mutableHandler()->component(mutableContainerPtr(), m_index, key);
    return *std::any_cast<T*>(result);
  }

  template <typename T>
  T& component(HashedString key)
    requires(!ReadOnly)
  {
    assert(isValid());
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

  ContainerView container() const {
    assert(isValid());
    return ContainerView{constHandler(), containerPtr()};
  }

  TrackIndexType previous() const {
    assert(isValid());
    return component<TrackIndexType, detail_tsp::kPreviousKey>();
  }

  TrackIndexType& previous()
    requires(!ReadOnly)
  {
    assert(isValid());
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
    assert(isValid());
    assert(hasReferenceSurface());
    const Surface* surface =
        constHandler()->referenceSurface(containerPtr(), m_index);
    assert(surface != nullptr);
    return *surface;
  }

  bool hasReferenceSurface() const {
    assert(isValid());
    return constHandler()->hasReferenceSurface(containerPtr(), m_index);
  }

  void setReferenceSurface(std::shared_ptr<const Surface> surface)
    requires(!ReadOnly)
  {
    assert(isValid());
    mutableHandler()->setReferenceSurface(mutableContainerPtr(), m_index,
                                          std::move(surface));
  }

  SourceLink getUncalibratedSourceLink() const {
    assert(isValid());
    assert(has(detail_tsp::kUncalibratedKey));
    return constHandler()->getUncalibratedSourceLink(containerPtr(), m_index);
  }

  void setUncalibratedSourceLink(SourceLink sourceLink)
    requires(!ReadOnly)
  {
    assert(isValid());
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
    assert(isValid());
    return constHandler()->calibratedSize(containerPtr(), m_index);
  }

  void allocateCalibrated(std::size_t measdim)
    requires(!ReadOnly)
  {
    assert(isValid());
    mutableHandler()->allocateCalibrated(mutableContainerPtr(), m_index,
                                         measdim);
    component<TrackIndexType, detail_tsp::kMeasDimKey>() =
        static_cast<TrackIndexType>(measdim);
  }

  template <typename val_t, typename cov_t>
  void allocateCalibrated(const Eigen::DenseBase<val_t>& val,
                          const Eigen::DenseBase<cov_t>& cov)
    requires(!ReadOnly)
  {
    constexpr std::size_t measdim = static_cast<std::size_t>(val_t::RowsAtCompileTime);
    static_assert(measdim > 0, "Measurement dimension must be static");
    static_assert(measdim <= eBoundSize,
                  "Measurement dimension exceeds bound parameters");
    static_assert(cov_t::RowsAtCompileTime == val_t::RowsAtCompileTime);
    static_assert(cov_t::ColsAtCompileTime == val_t::RowsAtCompileTime);
    allocateCalibrated(measdim);
    calibrated<measdim>() = val;
    calibratedCovariance<measdim>() = cov;
  }

  ConstEffectiveCalibratedMap effectiveCalibrated() const {
    assert(isValid());
    const double* data =
        constHandler()->calibratedData(containerPtr(), m_index);
    const auto size = calibratedSize();
    return ConstEffectiveCalibratedMap(data, size);
  }

  MutableEffectiveCalibratedMap effectiveCalibrated()
    requires(!ReadOnly)
  {
    assert(isValid());
    double* data =
        mutableHandler()->calibratedDataMutable(mutableContainerPtr(), m_index);
    const auto size = calibratedSize();
    return MutableEffectiveCalibratedMap(data, size);
  }

  ConstEffectiveCalibratedCovarianceMap effectiveCalibratedCovariance() const {
    assert(isValid());
    const double* data =
        constHandler()->calibratedCovarianceData(containerPtr(), m_index);
    const auto size = calibratedSize();
    return ConstEffectiveCalibratedCovarianceMap(data, size, size);
  }

  MutableEffectiveCalibratedCovarianceMap effectiveCalibratedCovariance()
    requires(!ReadOnly)
  {
    assert(isValid());
    double* data = mutableHandler()->calibratedCovarianceDataMutable(
        mutableContainerPtr(), m_index);
    const auto size = calibratedSize();
    return MutableEffectiveCalibratedCovarianceMap(data, size, size);
  }

  BoundSubspaceIndices projectorSubspaceIndices() const {
    assert(isValid());
    assert(hasProjector());
    const auto& serialized =
        component<SerializedSubspaceIndices, detail_tsp::kProjectorKey>();
    return deserializeSubspaceIndices<eBoundSize>(serialized);
  }

  template <std::ranges::sized_range index_range_t>
    requires(!ReadOnly)
  void setProjectorSubspaceIndices(const index_range_t& indices) {
    assert(isValid());
    storeProjectorIndices(encodeSubspaceIndices(indices),
                          static_cast<std::size_t>(indices.size()));
  }

  template <typename Derived>
  void setProjector(const Eigen::MatrixBase<Derived>& projector)
    requires(!ReadOnly)
  {
    static_assert(Derived::ColsAtCompileTime == eBoundSize,
                  "Projector must span the bound parameter dimension");
    std::size_t measdim = static_cast<std::size_t>(projector.rows());
    auto indices = projectorIndicesFromMatrix(projector.derived(), measdim);
    storeProjectorIndices(indices, measdim);
  }

  void setProjectorBitset(ProjectorBitset bitset)
    requires(!ReadOnly)
  {
    assert(isValid());
    std::bitset<eBoundSize * eBoundSize> bs(bitset);
    const auto matrix =
        bitsetToMatrix<ActsMatrix<eBoundSize, eBoundSize>>(bs);
    std::size_t measdim = eBoundSize;
    auto indices = projectorIndicesFromMatrix(matrix, measdim);
    storeProjectorIndices(indices, measdim);
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
    assert(isValid());
    assert(hasJacobian());
    return constHandler()->jacobian(containerPtr(), m_index);
  }

  CovarianceMap jacobian()
    requires(!ReadOnly)
  {
    assert(isValid());
    assert(hasJacobian());
    return mutableHandler()->jacobian(mutableContainerPtr(), m_index);
  }

  double pathLength() const {
    assert(isValid());
    return component<double, detail_tsp::kPathLengthKey>();
  }

  double& pathLength()
    requires(!ReadOnly)
  {
    assert(isValid());
    return component<double, detail_tsp::kPathLengthKey>();
  }

  ConstParametersMap predicted() const {
    assert(isValid());
    assert(hasPredicted());
    return parametersFromComponent(detail_tsp::kPredictedKey);
  }

  ParametersMap predicted()
    requires(!ReadOnly)
  {
    assert(isValid());
    return mutableParametersFromComponent(detail_tsp::kPredictedKey);
  }

  ConstCovarianceMap predictedCovariance() const {
    assert(isValid());
    assert(hasPredicted());
    return covarianceFromComponent(detail_tsp::kPredictedKey);
  }

  CovarianceMap predictedCovariance()
    requires(!ReadOnly)
  {
    assert(isValid());
    return mutableCovarianceFromComponent(detail_tsp::kPredictedKey);
  }

  ConstParametersMap filtered() const {
    assert(isValid());
    assert(hasFiltered());
    return parametersFromComponent(detail_tsp::kFilteredKey);
  }

  ParametersMap filtered()
    requires(!ReadOnly)
  {
    assert(isValid());
    return mutableParametersFromComponent(detail_tsp::kFilteredKey);
  }

  ConstCovarianceMap filteredCovariance() const {
    assert(isValid());
    assert(hasFiltered());
    return covarianceFromComponent(detail_tsp::kFilteredKey);
  }

  CovarianceMap filteredCovariance()
    requires(!ReadOnly)
  {
    assert(isValid());
    return mutableCovarianceFromComponent(detail_tsp::kFilteredKey);
  }

  ConstParametersMap smoothed() const {
    assert(isValid());
    assert(hasSmoothed());
    return parametersFromComponent(detail_tsp::kSmoothedKey);
  }

  ParametersMap smoothed()
    requires(!ReadOnly)
  {
    assert(isValid());
    return mutableParametersFromComponent(detail_tsp::kSmoothedKey);
  }

  ConstCovarianceMap smoothedCovariance() const {
    assert(isValid());
    assert(hasSmoothed());
    return covarianceFromComponent(detail_tsp::kSmoothedKey);
  }

  CovarianceMap smoothedCovariance()
    requires(!ReadOnly)
  {
    assert(isValid());
    return mutableCovarianceFromComponent(detail_tsp::kSmoothedKey);
  }

  ConstTrackStateType typeFlags() const {
    assert(isValid());
    const auto raw =
        component<TrackStateType::raw_type, detail_tsp::kTypeFlagsKey>();
    return ConstTrackStateType{raw};
  }

  TrackStateType typeFlags()
    requires(!ReadOnly)
  {
    assert(isValid());
    auto& raw =
        component<TrackStateType::raw_type, detail_tsp::kTypeFlagsKey>();
    return TrackStateType{raw};
  }

  float chi2() const {
    assert(isValid());
    return component<float, detail_tsp::kChi2Key>();
  }

  float& chi2()
    requires(!ReadOnly)
  {
    assert(isValid());
    return component<float, detail_tsp::kChi2Key>();
  }

  template <std::size_t measdim>
  typename TrackStateTraits<measdim, true>::Calibrated calibrated() const {
    assert(isValid());
    [[maybe_unused]] const auto size = calibratedSize();
    assert(size == static_cast<TrackIndexType>(measdim));
    const double* data =
        constHandler()->calibratedData(containerPtr(), m_index);
    return typename TrackStateTraits<measdim, true>::Calibrated(data);
  }

  template <std::size_t measdim>
  typename TrackStateTraits<measdim, false>::Calibrated calibrated()
    requires(!ReadOnly)
  {
    assert(isValid());
    [[maybe_unused]] const auto size = calibratedSize();
    assert(size == static_cast<TrackIndexType>(measdim));
    double* data =
        mutableHandler()->calibratedDataMutable(mutableContainerPtr(), m_index);
    return typename TrackStateTraits<measdim, false>::Calibrated(data);
  }

  template <std::size_t measdim>
  typename TrackStateTraits<measdim, true>::CalibratedCovariance
  calibratedCovariance() const {
    assert(isValid());
    [[maybe_unused]] const auto size = calibratedSize();
    assert(size == static_cast<TrackIndexType>(measdim));
    const double* data =
        constHandler()->calibratedCovarianceData(containerPtr(), m_index);
    return typename TrackStateTraits<measdim, true>::CalibratedCovariance(data);
  }

  template <std::size_t measdim>
  typename TrackStateTraits<measdim, false>::CalibratedCovariance
  calibratedCovariance()
    requires(!ReadOnly)
  {
    assert(isValid());
    [[maybe_unused]] const auto size = calibratedSize();
    assert(size == static_cast<TrackIndexType>(measdim));
    double* data = mutableHandler()->calibratedCovarianceDataMutable(
        mutableContainerPtr(), m_index);
    return
        typename TrackStateTraits<measdim, false>::CalibratedCovariance(data);
  }

  void unset(TrackStatePropMask target)
    requires(!ReadOnly)
  {
    assert(isValid());
    mutableHandler()->unset(mutableContainerPtr(), m_index, target);
  }

  void shareFrom(TrackStatePropMask shareSource, TrackStatePropMask shareTarget)
    requires(!ReadOnly)
  {
    shareFrom(*this, shareSource, shareTarget);
  }

  void shareFrom(TrackStatePropMask component)
    requires(!ReadOnly)
  {
    shareFrom(component, component);
  }

  template <bool otherReadOnly>
  void shareFrom(const AnyTrackState<otherReadOnly>& other,
                 TrackStatePropMask shareSource, TrackStatePropMask shareTarget)
    requires(!ReadOnly)
  {
    assert(isValid());
    assert(other.isValid());
    if (m_handler != other.m_handler) {
      throw std::invalid_argument(
          "Cannot share components between incompatible handlers");
    }
    assert(containerPtr() == other.containerPtr() &&
           "Cannot share components across MultiTrajectories");
    mutableHandler()->shareFrom(mutableContainerPtr(), m_index, other.m_index,
                                shareSource, shareTarget);
  }

  template <bool otherReadOnly>
  void shareFrom(const AnyTrackState<otherReadOnly>& other,
                 TrackStatePropMask component)
    requires(!ReadOnly)
  {
    shareFrom(other, component, component);
  }

  template <TrackStateProxyConcept proxy_t>
  void shareFrom(const proxy_t& other, TrackStatePropMask shareSource,
                 TrackStatePropMask shareTarget)
    requires(!ReadOnly)
  {
    AnyTrackState<true> otherAny(other);
    shareFrom(otherAny, shareSource, shareTarget);
  }

  template <TrackStateProxyConcept proxy_t>
  void shareFrom(const proxy_t& other, TrackStatePropMask component)
    requires(!ReadOnly)
  {
    shareFrom(other, component, component);
  }

  template <bool otherReadOnly>
  void copyFrom(const AnyTrackState<otherReadOnly>& other,
                TrackStatePropMask mask = TrackStatePropMask::All,
                bool onlyAllocated = true)
    requires(!ReadOnly)
  {
    using PM = TrackStatePropMask;
    assert(isValid());
    assert(other.isValid());
    if (m_handler != other.m_handler) {
      throw std::invalid_argument(
          "Cannot copy components between incompatible handlers");
    }

    auto src = other.getMask() & mask;
    if (src == TrackStatePropMask::None || mask == TrackStatePropMask::None) {
      // still propagate scalar quantities below
    }

    auto dest = getMask();
    if (onlyAllocated) {
      if (ACTS_CHECK_BIT(src, PM::Calibrated) && !hasCalibrated() &&
          other.hasCalibrated()) {
        allocateCalibrated(other.calibratedSize());
        dest |= PM::Calibrated;
      }

      auto missing =
          static_cast<std::underlying_type_t<PM>>((src ^ dest) & src);
      if ((missing != 0 || dest == TrackStatePropMask::None ||
           src == TrackStatePropMask::None) &&
          mask != TrackStatePropMask::None) {
        throw std::runtime_error(
            "Attempt track state copy with incompatible allocations");
      }
    } else {
      mutableHandler()->addTrackStateComponents(mutableContainerPtr(), m_index,
                                                mask);
    }

    if (ACTS_CHECK_BIT(src, PM::Predicted)) {
      predicted() = other.predicted();
      predictedCovariance() = other.predictedCovariance();
    }

    if (ACTS_CHECK_BIT(src, PM::Filtered)) {
      filtered() = other.filtered();
      filteredCovariance() = other.filteredCovariance();
    }

    if (ACTS_CHECK_BIT(src, PM::Smoothed)) {
      smoothed() = other.smoothed();
      smoothedCovariance() = other.smoothedCovariance();
    }

    if (other.hasUncalibratedSourceLink()) {
      setUncalibratedSourceLink(other.getUncalibratedSourceLink());
    }

    if (ACTS_CHECK_BIT(src, PM::Jacobian) && other.hasJacobian()) {
      jacobian() = other.jacobian();
    }

    if (ACTS_CHECK_BIT(src, PM::Calibrated) && other.hasCalibrated()) {
      auto* self = this;
      visit_measurement(other.calibratedSize(), [&](auto N) {
        constexpr std::size_t measdim = decltype(N)::value;
        self->allocateCalibrated(
            other.template calibrated<measdim>().eval(),
            other.template calibratedCovariance<measdim>().eval());
      });
      setProjectorSubspaceIndices(other.projectorSubspaceIndices());
    }

    chi2() = other.chi2();
    pathLength() = other.pathLength();
    typeFlags() = other.typeFlags();

    if (other.hasReferenceSurface()) {
      setReferenceSurface(other.referenceSurface().getSharedPtr());
    }

    mutableHandler()->copyDynamicFrom(mutableContainerPtr(), m_index,
                                      other.containerPtr(), other.m_index);
  }

  template <TrackStateProxyConcept proxy_t>
  void copyFrom(const proxy_t& other, TrackStatePropMask mask = TrackStatePropMask::All,
                bool onlyAllocated = true)
    requires(!ReadOnly)
  {
    AnyTrackState<true> otherAny(other);
    copyFrom(otherAny, mask, onlyAllocated);
  }

 private:
  template <bool>
  friend class AnyTrackState;

  template <std::ranges::sized_range index_range_t>
  BoundSubspaceIndices encodeSubspaceIndices(
      const index_range_t& indices) const {
    if (indices.size() > eBoundSize) {
      throw std::out_of_range("Projector dimension exceeds bound parameters");
    }
    BoundSubspaceIndices bound = kBoundSubspaceIndicesInvalid;
    std::size_t pos = 0;
    for (auto index : indices) {
      auto converted = static_cast<std::uint8_t>(index);
      if (converted >= eBoundSize) {
        throw std::out_of_range("Projector index out of range");
      }
      bound[pos++] = converted;
    }
    for (; pos < eBoundSize; ++pos) {
      bound[pos] = static_cast<std::uint8_t>(eBoundSize);
    }
    return bound;
  }

  template <typename Derived>
  BoundSubspaceIndices projectorIndicesFromMatrix(
      const Eigen::MatrixBase<Derived>& projector, std::size_t& measdim) const {
    BoundSubspaceIndices indices = kBoundSubspaceIndicesInvalid;
    measdim = 0;
    const auto rows = static_cast<std::size_t>(projector.rows());
    for (std::size_t row = 0; row < rows; ++row) {
      std::optional<std::uint8_t> chosen;
      for (std::size_t col = 0; col < static_cast<std::size_t>(projector.cols());
           ++col) {
        if (projector(row, col) != 0.) {
          chosen = static_cast<std::uint8_t>(col);
          break;
        }
      }
      if (!chosen) {
        throw std::invalid_argument(
            "Projector matrix row does not define a subspace index");
      }
      indices[row] = *chosen;
      ++measdim;
    }
    for (std::size_t row = rows; row < eBoundSize; ++row) {
      indices[row] = static_cast<std::uint8_t>(eBoundSize);
    }
    return indices;
  }

  void storeProjectorIndices(const BoundSubspaceIndices& indices,
                             std::size_t measdim)
    requires(!ReadOnly)
  {
    component<SerializedSubspaceIndices, detail_tsp::kProjectorKey>() =
        serializeSubspaceIndices<eBoundSize>(indices);
    component<TrackIndexType, detail_tsp::kMeasDimKey>() =
        static_cast<TrackIndexType>(measdim);
  }

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

  const detail::TrackStateHandlerConstBase* constHandler() const {
    return m_handler;
  }

  const detail::TrackStateHandlerMutableBase* mutableHandler() const
    requires(!ReadOnly)
  {
    return static_cast<const detail::TrackStateHandlerMutableBase*>(m_handler);
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
  const detail::TrackStateHandlerConstBase* m_handler;
};

using AnyConstTrackState = AnyTrackState<true>;
using AnyMutableTrackState = AnyTrackState<false>;

static_assert(ConstTrackStateProxyConcept<AnyConstTrackState>);
static_assert(MutableTrackStateProxyConcept<AnyMutableTrackState>);

}  // namespace Acts
