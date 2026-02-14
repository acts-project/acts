// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/MultiTrajectoryBackendConcept.hpp"
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/TrackStatePropMask.hpp"
#include "Acts/EventData/Types.hpp"
#include "Acts/EventData/detail/DynamicKeyIterator.hpp"
#include "Acts/Utilities/HashedString.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "ActsPlugins/EDM4hep/PodioDynamicColumns.hpp"
#include "ActsPlugins/EDM4hep/PodioTrackContainer.hpp"
#include "ActsPlugins/EDM4hep/PodioUtil.hpp"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#include "ActsPodioEdm/BoundParametersCollection.h"
#include "ActsPodioEdm/JacobianCollection.h"
#include "ActsPodioEdm/TrackStateCollection.h"
#include "ActsPodioEdm/TrackStateInfo.h"
#pragma GCC diagnostic pop

#include <any>
#include <memory>
#include <stdexcept>
#include <tuple>
#include <type_traits>

#include <podio/CollectionBase.h>
#include <podio/Frame.h>

namespace ActsPlugins {
/// @addtogroup edm4hep_plugin
/// @{

class MutablePodioTrackStateContainer;
class ConstPodioTrackStateContainer;

}  // namespace ActsPlugins

namespace Acts {
template <>
struct IsReadOnlyMultiTrajectory<ActsPlugins::ConstPodioTrackStateContainer>
    : std::true_type {};

template <>
struct IsReadOnlyMultiTrajectory<ActsPlugins::MutablePodioTrackStateContainer>
    : std::false_type {};
}  // namespace Acts

namespace ActsPlugins {

/// Base class for PODIO track state containers
class PodioTrackStateContainerBase {
 public:
  /// Mutable parameters map type
  using Parameters =
      typename Acts::detail_tsp::FixedSizeTypes<Acts::eBoundSize,
                                                false>::CoefficientsMap;
  /// Mutable covariance map type
  using Covariance =
      typename Acts::detail_tsp::FixedSizeTypes<Acts::eBoundSize,
                                                false>::CovarianceMap;

  /// Const parameters map type
  using ConstParameters =
      typename Acts::detail_tsp::FixedSizeTypes<Acts::eBoundSize,
                                                true>::CoefficientsMap;
  /// Const covariance map type
  using ConstCovariance =
      typename Acts::detail_tsp::FixedSizeTypes<Acts::eBoundSize,
                                                true>::CovarianceMap;

 protected:
  /// Check if a component exists for a track state
  /// @param instance Container instance
  /// @param key Component key
  /// @param istate Track state index
  /// @return True if component exists
  template <typename T>
  static constexpr bool has_impl(T& instance, Acts::HashedString key,
                                 Acts::TrackIndexType istate) {
    constexpr auto kInvalid = Acts::kTrackIndexInvalid;
    using namespace Acts::HashedStringLiteral;
    auto trackState = instance.m_collection->at(istate);
    const auto& data = trackState.getData();
    switch (key) {
      case "predicted"_hash:
        return data.ipredicted != kInvalid;
      case "filtered"_hash:
        return data.ifiltered != kInvalid;
      case "smoothed"_hash:
        return data.ismoothed != kInvalid;
      case "calibrated"_hash:
        return data.measdim != kInvalid;
      case "calibratedCov"_hash:
        return data.measdim != kInvalid;
      case "jacobian"_hash:
        return data.ijacobian != kInvalid;
      case "projector"_hash:
        return data.hasProjector;
      case "uncalibratedSourceLink"_hash:
        return data.uncalibratedIdentifier != PodioUtil::kNoIdentifier;
      case "previous"_hash:
      case "next"_hash:
      case "measdim"_hash:
      case "referenceSurface"_hash:
      case "chi2"_hash:
      case "pathLength"_hash:
      case "typeFlags"_hash:
        return true;
      default:
        return instance.m_dynamic.contains(key);
    }

    return false;
  }

  /// Get a component from a track state
  /// @param instance Container instance
  /// @param key Component key
  /// @param istate Track state index
  /// @return Component value as std::any
  template <bool EnsureConst, typename T>
  static std::any component_impl(T& instance, Acts::HashedString key,
                                 Acts::TrackIndexType istate) {
    if constexpr (EnsureConst) {
      static_assert(std::is_const_v<std::remove_reference_t<T>>,
                    "Is not const");
    }
    using namespace Acts::HashedStringLiteral;
    auto trackState = instance.m_collection->at(istate);
    std::conditional_t<EnsureConst, const ActsPodioEdm::TrackStateInfo*,
                       ActsPodioEdm::TrackStateInfo*>
        dataPtr;
    if constexpr (EnsureConst) {
      dataPtr = &trackState.getData();
    } else {
      dataPtr = &PodioUtil::getDataMutable(trackState);
    }
    auto& data = *dataPtr;
    switch (key) {
      case "previous"_hash:
        return &data.previous;
      case "next"_hash:
        return &data.next;
      case "predicted"_hash:
        return &data.ipredicted;
      case "filtered"_hash:
        return &data.ifiltered;
      case "smoothed"_hash:
        return &data.ismoothed;
      case "projector"_hash:
        return &data.projector;
      case "measdim"_hash:
        return &data.measdim;
      case "chi2"_hash:
        return &data.chi2;
      case "pathLength"_hash:
        return &data.pathLength;
      case "typeFlags"_hash:
        return &data.typeFlags;
      default:
        auto it = instance.m_dynamic.find(key);
        if (it == instance.m_dynamic.end()) {
          throw std::runtime_error("Unable to handle this component");
        }
        std::conditional_t<EnsureConst,
                           const podio_detail::ConstDynamicColumnBase*,
                           podio_detail::DynamicColumnBase*>
            col = it->second.get();
        assert(col && "Dynamic column is null");
        return col->get(istate);
    }
  }

  /// Check if a column exists
  /// @param instance Container instance
  /// @param key Column key
  /// @return True if column exists
  template <typename T>
  static constexpr bool hasColumn_impl(T& instance, Acts::HashedString key) {
    using namespace Acts::HashedStringLiteral;
    switch (key) {
      case "predicted"_hash:
      case "filtered"_hash:
      case "smoothed"_hash:
      case "jacobian"_hash:
      case "projector"_hash:
      case "previous"_hash:
      case "next"_hash:
      case "uncalibratedSourceLink"_hash:
      case "referenceSurface"_hash:
      case "measdim"_hash:
      case "chi2"_hash:
      case "pathLength"_hash:
      case "typeFlags"_hash:
        return true;
      default:
        return instance.m_dynamic.contains(key);
    }
  }

  /// Populate surface buffer from track state collection
  /// @param helper Conversion helper
  /// @param collection Track state collection
  /// @param surfaces Output surface buffer
  static void populateSurfaceBuffer(
      const PodioUtil::ConversionHelper& helper,
      const ActsPodioEdm::TrackStateCollection& collection,
      std::vector<std::shared_ptr<const Acts::Surface>>& surfaces) noexcept {
    surfaces.reserve(collection.size());
    for (ActsPodioEdm::TrackState trackState : collection) {
      surfaces.push_back(PodioUtil::convertSurfaceFromPodio(
          helper, trackState.getReferenceSurface()));
    }
  }
};

/// Read-only track state container backend using podio for storage
class ConstPodioTrackStateContainer final
    : public PodioTrackStateContainerBase,
      public Acts::MultiTrajectory<ConstPodioTrackStateContainer> {
 public:
  /// Constructor from collections
  /// @param helper Conversion helper
  /// @param trackStates Track state collection
  /// @param params Parameters collection
  /// @param jacs Jacobian collection
  ConstPodioTrackStateContainer(
      const PodioUtil::ConversionHelper& helper,
      const ActsPodioEdm::TrackStateCollection& trackStates,
      const ActsPodioEdm::BoundParametersCollection& params,
      const ActsPodioEdm::JacobianCollection& jacs)
      : m_helper{helper},
        m_collection{&trackStates},
        m_params{&params},
        m_jacs{&jacs} {
    // Not much we can do to recover dynamic columns here
    populateSurfaceBuffer(m_helper, *m_collection, m_surfaces);
  }

  /// Construct a const track state container from a mutable
  /// @warning If the source mutable container is modified, this container
  ///          will be corrupted, as surface buffer and dynamic column state can
  ///          not be synchronized!
  /// @param other The mutable container to construct from
  explicit ConstPodioTrackStateContainer(
      const MutablePodioTrackStateContainer& other);

  /// Constructor from frame
  /// @param helper Conversion helper
  /// @param frame Frame containing track state data
  /// @param suffix Optional collection name suffix
  ConstPodioTrackStateContainer(const PodioUtil::ConversionHelper& helper,
                                const podio::Frame& frame,
                                const std::string& suffix = "")
      : m_helper{helper},
        m_collection{nullptr},
        m_params{nullptr},
        m_jacs{nullptr} {
    std::string s = suffix.empty() ? suffix : "_" + suffix;

    std::vector<std::string> available = frame.getAvailableCollections();

    std::string trackStatesKey = "trackStates" + s;
    std::string paramsKey = "trackStateParameters" + s;
    std::string jacsKey = "trackStateJacobians" + s;

    if (!Acts::rangeContainsValue(available, trackStatesKey)) {
      throw std::runtime_error{"Track state collection '" + trackStatesKey +
                               "' not found in frame"};
    }

    if (!Acts::rangeContainsValue(available, paramsKey)) {
      throw std::runtime_error{"Track state parameters collection '" +
                               paramsKey + "' not found in frame"};
    }

    if (!Acts::rangeContainsValue(available, jacsKey)) {
      throw std::runtime_error{"Track state jacobian collection '" + jacsKey +
                               "' not found in frame"};
    }

    loadCollection<ActsPodioEdm::TrackStateCollection>(m_collection, frame,
                                                       trackStatesKey);
    loadCollection<ActsPodioEdm::BoundParametersCollection>(m_params, frame,
                                                            paramsKey);
    loadCollection<ActsPodioEdm::JacobianCollection>(m_jacs, frame, jacsKey);

    populateSurfaceBuffer(m_helper, *m_collection, m_surfaces);

    podio_detail::recoverDynamicColumns(frame, trackStatesKey, m_dynamic);
  }

  /// @copydoc Acts::MultiTrajectory::parameters(IndexType) const
  using Acts::MultiTrajectory<ConstPodioTrackStateContainer>::parameters;
  /// @copydoc Acts::MultiTrajectory::covariance(IndexType) const
  using Acts::MultiTrajectory<ConstPodioTrackStateContainer>::covariance;

  /// Get dynamic column keys
  /// @return Range of dynamic column keys
  Acts::detail::DynamicKeyRange<podio_detail::ConstDynamicColumnBase>
  dynamicKeys_impl() const {
    return {m_dynamic.begin(), m_dynamic.end()};
  }

 private:
  template <typename collection_t>
  static void loadCollection(collection_t const*& dest,
                             const podio::Frame& frame,
                             const std::string& key) {
    const auto* collection = frame.get(key);

    if (const auto* d = dynamic_cast<const collection_t*>(collection);
        d != nullptr) {
      dest = d;
    } else {
      throw std::runtime_error{"Unable to get collection " + key};
    }
  }

 public:
  /// Get parameters for a track state
  /// @param istate Track state index
  /// @return Track parameters
  ConstParameters parameters_impl(IndexType istate) const {
    return ConstParameters{m_params->at(istate).getData().values.data()};
  }

  /// Get covariance for a track state
  /// @param istate Track state index
  /// @return Track covariance matrix
  ConstCovariance covariance_impl(IndexType istate) const {
    return ConstCovariance{m_params->at(istate).getData().covariance.data()};
  }

  /// Get jacobian for a track state
  /// @param istate Track state index
  /// @return Jacobian matrix
  ConstCovariance jacobian_impl(IndexType istate) const {
    IndexType ijacobian = m_collection->at(istate).getData().ijacobian;
    return ConstCovariance{m_jacs->at(ijacobian).getData().values.data()};
  }

  /// Get calibrated measurement
  /// @param index Track state index
  /// @return Calibrated measurement
  template <std::size_t measdim>
  ConstTrackStateProxy::ConstCalibrated<measdim> calibrated_impl(
      IndexType index) const {
    return ConstTrackStateProxy::ConstCalibrated<measdim>{
        m_collection->at(index).getData().measurement.data()};
  }

  /// Get calibrated measurement covariance
  /// @param index Track state index
  /// @return Calibrated measurement covariance
  template <std::size_t measdim>
  ConstTrackStateProxy::ConstCalibratedCovariance<measdim>
  calibratedCovariance_impl(IndexType index) const {
    return ConstTrackStateProxy::ConstCalibratedCovariance<measdim>{
        m_collection->at(index).getData().measurementCovariance.data()};
  }

  /// Get number of track states
  /// @return Number of track states
  IndexType size_impl() const { return m_collection->size(); }

  /// Get a component from a track state
  /// @param key Component key
  /// @param istate Track state index
  /// @return Component value
  std::any component_impl(Acts::HashedString key, IndexType istate) const {
    return PodioTrackStateContainerBase::component_impl<true>(*this, key,
                                                              istate);
  }

  /// Check if a column exists
  /// @param key Column key
  /// @return True if column exists
  constexpr bool hasColumn_impl(Acts::HashedString key) const {
    return PodioTrackStateContainerBase::hasColumn_impl(*this, key);
  }

  /// Check if a component exists for a track state
  /// @param key Component key
  /// @param istate Track state index
  /// @return True if component exists
  constexpr bool has_impl(Acts::HashedString key, IndexType istate) const {
    return PodioTrackStateContainerBase::has_impl(*this, key, istate);
  }

  /// Get calibrated measurement size
  /// @param istate Track state index
  /// @return Measurement size
  Acts::TrackIndexType calibratedSize_impl(IndexType istate) const {
    return m_collection->at(istate).getData().measdim;
  }

  /// Get uncalibrated source link for a track state
  /// @param istate Track state index
  /// @return Uncalibrated source link
  Acts::SourceLink getUncalibratedSourceLink_impl(IndexType istate) const {
    return m_helper.get().identifierToSourceLink(
        m_collection->at(istate).getData().uncalibratedIdentifier);
  }

  /// Get reference surface for a track state
  /// @param istate Track state index
  /// @return Reference surface pointer
  const Acts::Surface* referenceSurface_impl(IndexType istate) const {
    return m_surfaces.at(istate).get();
  }

 private:
  friend class PodioTrackStateContainerBase;

  std::reference_wrapper<const PodioUtil::ConversionHelper> m_helper;
  const ActsPodioEdm::TrackStateCollection* m_collection;
  const ActsPodioEdm::BoundParametersCollection* m_params;
  const ActsPodioEdm::JacobianCollection* m_jacs;
  std::vector<std::shared_ptr<const Acts::Surface>> m_surfaces;

  std::unordered_map<Acts::HashedString,
                     std::unique_ptr<podio_detail::ConstDynamicColumnBase>>
      m_dynamic;
  std::vector<Acts::HashedString> m_dynamicKeys;
};

static_assert(
    Acts::IsReadOnlyMultiTrajectory<ConstPodioTrackStateContainer>::value,
    "MutablePodioTrackStateContainer should not be read-only");

static_assert(
    Acts::ConstMultiTrajectoryBackend<ConstPodioTrackStateContainer>,
    "ConstPodioTrackStateContainer does not fulfill TrackContainerBackend");

/// Mutable Podio-based track state container implementation
class MutablePodioTrackStateContainer final
    : public PodioTrackStateContainerBase,
      public Acts::MultiTrajectory<MutablePodioTrackStateContainer> {
 public:
  /// Constructor
  /// @param helper Conversion helper
  explicit MutablePodioTrackStateContainer(PodioUtil::ConversionHelper& helper)
      : m_helper{helper} {
    m_collection = std::make_unique<ActsPodioEdm::TrackStateCollection>();
    m_jacs = std::make_unique<ActsPodioEdm::JacobianCollection>();
    m_params = std::make_unique<ActsPodioEdm::BoundParametersCollection>();

    populateSurfaceBuffer(m_helper, *m_collection, m_surfaces);
  }

  /// Get const track parameters for a track state
  /// @param istate Track state index
  /// @return Const parameters
  ConstParameters parameters_impl(IndexType istate) const {
    return ConstParameters{m_params->at(istate).getData().values.data()};
  }

  /// Get mutable track parameters for a track state
  /// @param istate Track state index
  /// @return Mutable parameters
  Parameters parameters_impl(IndexType istate) {
    return Parameters{
        PodioUtil::getDataMutable(m_params->at(istate)).values.data()};
  }

  /// Get const covariance matrix for a track state
  /// @param istate Track state index
  /// @return Const covariance matrix
  ConstCovariance covariance_impl(IndexType istate) const {
    return ConstCovariance{m_params->at(istate).getData().covariance.data()};
  }

  /// Get mutable covariance matrix for a track state
  /// @param istate Track state index
  /// @return Mutable covariance matrix
  Covariance covariance_impl(IndexType istate) {
    return Covariance{
        PodioUtil::getDataMutable(m_params->at(istate)).covariance.data()};
  }

  /// Get the jacobian matrix for a track state (const version)
  /// @param istate Track state index
  /// @return Const jacobian matrix
  ConstCovariance jacobian_impl(IndexType istate) const {
    IndexType ijacobian = m_collection->at(istate).getData().ijacobian;
    return ConstCovariance{m_jacs->at(ijacobian).getData().values.data()};
  }

  /// Get the jacobian matrix for a track state (mutable version)
  /// @param istate Track state index
  /// @return Mutable jacobian matrix
  Covariance jacobian_impl(IndexType istate) {
    IndexType ijacobian = m_collection->at(istate).getData().ijacobian;
    return Covariance{
        PodioUtil::getDataMutable(m_jacs->at(ijacobian)).values.data()};
  }

  /// Get calibrated measurement vector (const version)
  /// @tparam measdim Dimension of the measurement
  /// @param index Track state index
  /// @return Const calibrated measurement vector
  template <std::size_t measdim>
  ConstTrackStateProxy::ConstCalibrated<measdim> calibrated_impl(
      IndexType index) const {
    return ConstTrackStateProxy::ConstCalibrated<measdim>{
        m_collection->at(index).getData().measurement.data()};
  }

  /// Get calibrated measurement vector (mutable version)
  /// @tparam measdim Dimension of the measurement
  /// @param index Track state index
  /// @return Mutable calibrated measurement vector
  template <std::size_t measdim>
  TrackStateProxy::Calibrated<measdim> calibrated_impl(IndexType index) {
    return TrackStateProxy::Calibrated<measdim>{
        PodioUtil::getDataMutable(m_collection->at(index)).measurement.data()};
  }

  /// Get calibrated measurement covariance (const version)
  /// @tparam measdim Dimension of the measurement
  /// @param index Track state index
  /// @return Const calibrated covariance matrix
  template <std::size_t measdim>
  ConstTrackStateProxy::ConstCalibratedCovariance<measdim>
  calibratedCovariance_impl(IndexType index) const {
    return ConstTrackStateProxy::ConstCalibratedCovariance<measdim>{
        m_collection->at(index).getData().measurementCovariance.data()};
  }

  /// Get calibrated measurement covariance (mutable version)
  /// @tparam measdim Dimension of the measurement
  /// @param index Track state index
  /// @return Mutable calibrated covariance matrix
  template <std::size_t measdim>
  TrackStateProxy::CalibratedCovariance<measdim> calibratedCovariance_impl(
      IndexType index) {
    return TrackStateProxy::CalibratedCovariance<measdim>{
        PodioUtil::getDataMutable(m_collection->at(index))
            .measurementCovariance.data()};
  }

  /// Get the number of track states
  /// @return Number of track states in the container
  IndexType size_impl() const { return m_collection->size(); }

  /// Get a component by key from the track state (const version)
  /// @param key Component key
  /// @param istate Track state index
  /// @return Component value
  std::any component_impl(Acts::HashedString key, IndexType istate) const {
    return PodioTrackStateContainerBase::component_impl<true>(*this, key,
                                                              istate);
  }

  /// Get a component by key from the track state (mutable version)
  /// @param key Component key
  /// @param istate Track state index
  /// @return Component value
  std::any component_impl(Acts::HashedString key, IndexType istate) {
    return PodioTrackStateContainerBase::component_impl<false>(*this, key,
                                                               istate);
  }

  /// Check if a dynamic column exists
  /// @param key Column key
  /// @return True if column exists
  constexpr bool hasColumn_impl(Acts::HashedString key) const {
    return PodioTrackStateContainerBase::hasColumn_impl(*this, key);
  }

  /// Check if a track state component is present
  /// @param key Component key
  /// @param istate Track state index
  /// @return True if component is present
  constexpr bool has_impl(Acts::HashedString key, IndexType istate) const {
    return PodioTrackStateContainerBase::has_impl(*this, key, istate);
  }

  /// Add a new track state and return its index
  /// @param mask Track state component mask
  /// @param iprevious Index of previous track state
  /// @return Index of the new track state
  IndexType addTrackState_impl(
      Acts::TrackStatePropMask mask = Acts::TrackStatePropMask::All,
      Acts::TrackIndexType iprevious = Acts::kTrackIndexInvalid) {
    auto trackState = m_collection->create();
    auto& data = PodioUtil::getDataMutable(trackState);
    data.previous = iprevious;
    data.ipredicted = kInvalid;
    data.ifiltered = kInvalid;
    data.ismoothed = kInvalid;
    data.ijacobian = kInvalid;

    PodioUtil::getReferenceSurfaceMutable(trackState).surfaceType =
        PodioUtil::kNoSurface;

    using enum Acts::TrackStatePropMask;

    if (ACTS_CHECK_BIT(mask, Predicted)) {
      m_params->create();
      data.ipredicted = m_params->size() - 1;
    }
    if (ACTS_CHECK_BIT(mask, Filtered)) {
      m_params->create();
      data.ifiltered = m_params->size() - 1;
    }
    if (ACTS_CHECK_BIT(mask, Smoothed)) {
      m_params->create();
      data.ismoothed = m_params->size() - 1;
    }
    if (ACTS_CHECK_BIT(mask, Jacobian)) {
      m_jacs->create();
      data.ijacobian = m_jacs->size() - 1;
    }
    data.measdim = kInvalid;
    data.hasProjector = false;
    if (ACTS_CHECK_BIT(mask, Calibrated)) {
      data.hasProjector = true;
    }
    m_surfaces.emplace_back();

    data.uncalibratedIdentifier = PodioUtil::kNoIdentifier;
    assert(m_collection->size() == m_surfaces.size() &&
           "Inconsistent surface buffer");

    for (const auto& [key, vec] : m_dynamic) {
      vec->add();
    }

    return m_collection->size() - 1;
  }

  /// Add track state components based on mask
  /// @param istate Track state index
  /// @param mask Track state component mask
  void addTrackStateComponents_impl(IndexType istate,
                                    Acts::TrackStatePropMask mask) {
    auto& data = PodioUtil::getDataMutable(m_collection->at(istate));

    using enum Acts::TrackStatePropMask;

    if (ACTS_CHECK_BIT(mask, Predicted) && data.ipredicted == kInvalid) {
      m_params->create();
      data.ipredicted = m_params->size() - 1;
    }

    if (ACTS_CHECK_BIT(mask, Filtered) && data.ifiltered == kInvalid) {
      m_params->create();
      data.ifiltered = m_params->size() - 1;
    }

    if (ACTS_CHECK_BIT(mask, Smoothed) && data.ismoothed == kInvalid) {
      m_params->create();
      data.ismoothed = m_params->size() - 1;
    }

    if (ACTS_CHECK_BIT(mask, Jacobian) && data.ijacobian == kInvalid) {
      m_jacs->create();
      data.ijacobian = m_jacs->size() - 1;
    }

    if (ACTS_CHECK_BIT(mask, Calibrated) && !data.hasProjector) {
      data.hasProjector = true;
    }
  }

  /// Share data from another track state
  /// @param iself Index of the destination track state
  /// @param iother Index of the source track state
  /// @param shareSource Component to share from source
  /// @param shareTarget Target component to share to
  void shareFrom_impl(Acts::TrackIndexType iself, Acts::TrackIndexType iother,
                      Acts::TrackStatePropMask shareSource,
                      Acts::TrackStatePropMask shareTarget) {
    auto& self = PodioUtil::getDataMutable(m_collection->at(iself));
    auto& other = PodioUtil::getDataMutable(m_collection->at(iother));

    assert(ACTS_CHECK_BIT(getTrackState(iother).getMask(), shareSource) &&
           "Source has incompatible allocation");

    using enum Acts::TrackStatePropMask;

    IndexType sourceIndex{kInvalid};
    switch (shareSource) {
      case Predicted:
        sourceIndex = other.ipredicted;
        break;
      case Filtered:
        sourceIndex = other.ifiltered;
        break;
      case Smoothed:
        sourceIndex = other.ismoothed;
        break;
      case Jacobian:
        sourceIndex = other.ijacobian;
        break;
      default:
        throw std::domain_error{"Unable to share this component"};
    }

    assert(sourceIndex != kInvalid);

    switch (shareTarget) {
      case Predicted:
        assert(shareSource != Jacobian);
        self.ipredicted = sourceIndex;
        break;
      case Filtered:
        assert(shareSource != Jacobian);
        self.ifiltered = sourceIndex;
        break;
      case Smoothed:
        assert(shareSource != Jacobian);
        self.ismoothed = sourceIndex;
        break;
      case Jacobian:
        assert(shareSource == Jacobian);
        self.ijacobian = sourceIndex;
        break;
      default:
        throw std::domain_error{"Unable to share this component"};
    }
  }

  /// Unset a track state component
  /// @param target Component to unset
  /// @param istate Track state index
  void unset_impl(Acts::TrackStatePropMask target,
                  Acts::TrackIndexType istate) {
    auto& data = PodioUtil::getDataMutable(m_collection->at(istate));

    switch (target) {
      using enum Acts::TrackStatePropMask;

      case Predicted:
        data.ipredicted = kInvalid;
        break;
      case Filtered:
        data.ifiltered = kInvalid;
        break;
      case Smoothed:
        data.ismoothed = kInvalid;
        break;
      case Jacobian:
        data.ijacobian = kInvalid;
        break;
      case Calibrated:
        data.measdim = kInvalid;
        break;
      default:
        throw std::domain_error{"Unable to unset this component"};
    }
  }

  /// Clear all track states
  void clear_impl() {
    m_collection->clear();
    m_params->clear();
    m_surfaces.clear();
    for (const auto& [key, vec] : m_dynamic) {
      vec->clear();
    }
  }

  /// Add a dynamic column
  /// @tparam T Column value type
  /// @param key Column key
  template <typename T>
  constexpr void addColumn_impl(std::string_view key) {
    Acts::HashedString hashedKey = Acts::hashStringDynamic(key);
    m_dynamic.insert(
        {hashedKey, std::make_unique<podio_detail::DynamicColumn<T>>(key)});
  }

  /// Allocate calibrated measurement and covariance
  /// @tparam val_t Type of the measurement vector
  /// @tparam cov_t Type of the covariance matrix
  /// @param istate Track state index
  /// @param val Measurement vector to store
  /// @param cov Covariance matrix to store
  template <typename val_t, typename cov_t>
  void allocateCalibrated_impl(IndexType istate,
                               const Eigen::DenseBase<val_t>& val,
                               const Eigen::DenseBase<cov_t>& cov)
    requires(Acts::Concepts::eigen_base_is_fixed_size<val_t> &&
             Eigen::PlainObjectBase<val_t>::RowsAtCompileTime <=
                 toUnderlying(Acts::eBoundSize) &&
             Acts::Concepts::eigen_bases_have_same_num_rows<val_t, cov_t> &&
             Acts::Concepts::eigen_base_is_square<cov_t>)
  {
    constexpr std::size_t measdim = val_t::RowsAtCompileTime;

    auto& data = PodioUtil::getDataMutable(m_collection->at(istate));

    if (data.measdim != kInvalid && data.measdim != measdim) {
      throw std::invalid_argument{
          "Measurement dimension does not match the allocated dimension"};
    }

    data.measdim = measdim;

    Eigen::Map<Acts::Vector<measdim>> valMap(data.measurement.data());
    valMap = val;

    Eigen::Map<Acts::SquareMatrix<measdim>> covMap(
        data.measurementCovariance.data());
    covMap = cov;
  }

  /// Set the uncalibrated source link for a track state
  /// @param istate Track state index
  /// @param sourceLink Source link to set
  void setUncalibratedSourceLink_impl(IndexType istate,
                                      const Acts::SourceLink& sourceLink) {
    PodioUtil::Identifier id =
        m_helper.get().sourceLinkToIdentifier(sourceLink);
    auto& data = PodioUtil::getDataMutable(m_collection->at(istate));
    data.uncalibratedIdentifier = id;
  }

  /// Set the reference surface for a track state
  /// @param istate Track state index
  /// @param surface Reference surface to set
  void setReferenceSurface_impl(IndexType istate,
                                std::shared_ptr<const Acts::Surface> surface) {
    auto trackState = m_collection->at(istate);
    trackState.setReferenceSurface(
        PodioUtil::convertSurfaceToPodio(m_helper, *surface));
    m_surfaces.at(istate) = std::move(surface);
  }

  /// Get the size of the calibrated measurement
  /// @param istate Track state index
  /// @return Size of the calibrated measurement
  Acts::TrackIndexType calibratedSize_impl(IndexType istate) const {
    return m_collection->at(istate).getData().measdim;
  }

  /// Get the uncalibrated source link for a track state
  /// @param istate Track state index
  /// @return Uncalibrated source link
  Acts::SourceLink getUncalibratedSourceLink_impl(IndexType istate) const {
    return m_helper.get().identifierToSourceLink(
        m_collection->at(istate).getData().uncalibratedIdentifier);
  }

  /// Get the reference surface for a track state
  /// @param istate Track state index
  /// @return Pointer to the reference surface
  const Acts::Surface* referenceSurface_impl(IndexType istate) const {
    return m_surfaces.at(istate).get();
  }

  /// Release collections into a podio frame
  /// @param frame The podio frame to release into
  /// @param suffix Optional suffix for collection names
  void releaseInto(podio::Frame& frame, const std::string& suffix = "") {
    std::string s = suffix;
    if (!s.empty()) {
      s = "_" + s;
    }
    frame.put(std::move(m_collection), "trackStates" + s);
    frame.put(std::move(m_params), "trackStateParameters" + s);
    frame.put(std::move(m_jacs), "trackStateJacobians" + s);
    m_surfaces.clear();

    for (const auto& [key, col] : m_dynamic) {
      col->releaseInto(frame, "trackStates" + s + "_extra__");
    }

    m_dynamic.clear();
  }

  /// Get dynamic column keys
  /// @return Range of dynamic column keys
  Acts::detail::DynamicKeyRange<podio_detail::DynamicColumnBase>
  dynamicKeys_impl() const {
    return {m_dynamic.begin(), m_dynamic.end()};
  }

  /// Copy dynamic column data from another track state
  /// @param dstIdx Destination track state index
  /// @param key Column key
  /// @param srcPtr Source pointer
  void copyDynamicFrom_impl(IndexType dstIdx, Acts::HashedString key,
                            const std::any& srcPtr) {
    auto it = m_dynamic.find(key);
    if (it == m_dynamic.end()) {
      throw std::invalid_argument{
          "Destination container does not have matching dynamic column"};
    }

    it->second->copyFrom(dstIdx, srcPtr);
  }

 private:
  friend class PodioTrackStateContainerBase;
  friend class ConstPodioTrackStateContainer;

  std::reference_wrapper<PodioUtil::ConversionHelper> m_helper;
  std::unique_ptr<ActsPodioEdm::TrackStateCollection> m_collection;
  std::unique_ptr<ActsPodioEdm::BoundParametersCollection> m_params;
  std::unique_ptr<ActsPodioEdm::JacobianCollection> m_jacs;
  std::vector<std::shared_ptr<const Acts::Surface>> m_surfaces;

  std::unordered_map<Acts::HashedString,
                     std::unique_ptr<podio_detail::DynamicColumnBase>>
      m_dynamic;
  std::vector<Acts::HashedString> m_dynamicKeys;
};

static_assert(
    !Acts::IsReadOnlyMultiTrajectory<MutablePodioTrackStateContainer>::value,
    "MutablePodioTrackStateContainer should not be read-only");

static_assert(
    Acts::MutableMultiTrajectoryBackend<MutablePodioTrackStateContainer>,
    "MutablePodioTrackStateContainer does not fulfill "
    "TrackStateContainerBackend");

inline ConstPodioTrackStateContainer::ConstPodioTrackStateContainer(
    const MutablePodioTrackStateContainer& other)
    : m_helper{other.m_helper},
      m_collection{other.m_collection.get()},
      m_params{other.m_params.get()},
      m_jacs{other.m_jacs.get()},
      m_surfaces{other.m_surfaces} {
  for (const auto& [key, col] : other.m_dynamic) {
    m_dynamic.insert({key, col->asConst()});
  }
}

/// @}
}  // namespace ActsPlugins
