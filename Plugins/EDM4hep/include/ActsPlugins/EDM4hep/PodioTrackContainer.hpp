// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/detail/DynamicKeyIterator.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "ActsPlugins/EDM4hep/PodioDynamicColumns.hpp"
#include "ActsPlugins/EDM4hep/PodioUtil.hpp"
#include "ActsPodioEdm/Surface.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#include "ActsPodioEdm/ParticleHypothesis.h"
#include "ActsPodioEdm/Track.h"
#include "ActsPodioEdm/TrackCollection.h"
#include "ActsPodioEdm/TrackInfo.h"
#pragma GCC diagnostic pop

#include <stdexcept>
#include <type_traits>

#include <podio/Frame.h>

namespace ActsPlugins {
/// @addtogroup edm4hep_plugin
/// @{

class MutablePodioTrackContainer;
class ConstPodioTrackContainer;

}  // namespace ActsPlugins

namespace Acts {

template <>
struct IsReadOnlyTrackContainer<ActsPlugins::MutablePodioTrackContainer>
    : std::false_type {};

template <>
struct IsReadOnlyTrackContainer<ActsPlugins::ConstPodioTrackContainer>
    : std::true_type {};
}  // namespace Acts

namespace ActsPlugins {

/// Base class for PODIO track containers
class PodioTrackContainerBase {
 public:
  /// Track index type
  using IndexType = Acts::TrackIndexType;
  /// Invalid track index
  static constexpr auto kInvalid = Acts::kTrackIndexInvalid;
  /// Maximum measurement size
  static constexpr auto MeasurementSizeMax = Acts::kMeasurementSizeMax;

  /// Mutable track parameters type
  using Parameters =
      typename Acts::detail_tsp::FixedSizeTypes<Acts::eBoundSize,
                                                false>::CoefficientsMap;
  /// Mutable track covariance type
  using Covariance =
      typename Acts::detail_tsp::FixedSizeTypes<Acts::eBoundSize,
                                                false>::CovarianceMap;

  /// Const track parameters type
  using ConstParameters =
      typename Acts::detail_tsp::FixedSizeTypes<Acts::eBoundSize,
                                                true>::CoefficientsMap;
  /// Const track covariance type
  using ConstCovariance =
      typename Acts::detail_tsp::FixedSizeTypes<Acts::eBoundSize,
                                                true>::CovarianceMap;

 protected:
  /// Constructor
  /// @param helper Conversion helper
  explicit PodioTrackContainerBase(const PodioUtil::ConversionHelper& helper)
      : m_helper{helper} {}

  /// Get component implementation
  /// @param instance Container instance
  /// @param key Component key
  /// @param itrack Track index
  /// @return Component value
  template <bool EnsureConst, typename T>
  static std::any component_impl(T& instance, Acts::HashedString key,
                                 IndexType itrack) {
    using namespace Acts::HashedStringLiteral;
    if constexpr (EnsureConst) {
      static_assert(std::is_const_v<std::remove_reference_t<T>>,
                    "Is not const");
    }

    using namespace Acts::HashedStringLiteral;
    auto track = instance.m_collection->at(itrack);
    std::conditional_t<EnsureConst, const ActsPodioEdm::TrackInfo*,
                       ActsPodioEdm::TrackInfo*>
        dataPtr;
    if constexpr (EnsureConst) {
      dataPtr = &track.getData();
    } else {
      dataPtr = &PodioUtil::getDataMutable(track);
    }
    auto& data = *dataPtr;
    switch (key) {
      case "tipIndex"_hash:
        return &data.tipIndex;
      case "stemIndex"_hash:
        return &data.stemIndex;
      case "params"_hash:
        return data.parameters.data();
      case "cov"_hash:
        return data.covariance.data();
      case "nMeasurements"_hash:
        return &data.nMeasurements;
      case "nHoles"_hash:
        return &data.nHoles;
      case "chi2"_hash:
        return &data.chi2;
      case "ndf"_hash:
        return &data.ndf;
      case "nOutliers"_hash:
        return &data.nOutliers;
      case "nSharedHits"_hash:
        return &data.nSharedHits;
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
        return col->get(itrack);
    }
  }

  /// Get dynamic keys implementation
  /// @param instance Container instance
  /// @return Dynamic key range
  template <typename T>
  static auto dynamicKeys_impl(T& instance) {
    using column_type =
        typename decltype(instance.m_dynamic)::mapped_type::element_type;
    return Acts::detail::DynamicKeyRange<column_type>{
        instance.m_dynamic.begin(), instance.m_dynamic.end()};
  }

  /// Get particle hypothesis implementation
  /// @param instance Container instance
  /// @param itrack Track index
  /// @return Particle hypothesis
  template <typename T>
  static Acts::ParticleHypothesis particleHypothesis_impl(T& instance,
                                                          IndexType itrack) {
    auto track = instance.m_collection->at(itrack);
    const auto& src = track.getParticleHypothesis();
    return Acts::ParticleHypothesis{static_cast<Acts::PdgParticle>(src.absPdg),
                                    src.mass, Acts::ChargeHypothesis{src.absQ}};
  }

  /// Populate surface buffer from track collection
  /// @param helper Conversion helper
  /// @param collection Track collection
  /// @param surfaces Surface buffer to populate
  static void populateSurfaceBuffer(
      const PodioUtil::ConversionHelper& helper,
      const ActsPodioEdm::TrackCollection& collection,
      std::vector<std::shared_ptr<const Acts::Surface>>& surfaces) noexcept {
    surfaces.reserve(collection.size());
    for (ActsPodioEdm::Track track : collection) {
      surfaces.push_back(PodioUtil::convertSurfaceFromPodio(
          helper, track.getReferenceSurface()));
    }
  }

  /// Conversion helper
  std::reference_wrapper<const PodioUtil::ConversionHelper> m_helper;
  /// Surface buffer
  std::vector<std::shared_ptr<const Acts::Surface>> m_surfaces;
};

/// Mutable Podio-based track container implementation
class MutablePodioTrackContainer : public PodioTrackContainerBase {
 public:
  /// Constructor
  /// @param helper Conversion helper
  explicit MutablePodioTrackContainer(const PodioUtil::ConversionHelper& helper)
      : PodioTrackContainerBase{helper},
        m_collection{std::make_unique<ActsPodioEdm::TrackCollection>()} {
    populateSurfaceBuffer(m_helper, *m_collection, m_surfaces);
  }

  /// Copy constructor
  /// @param other Source container
  MutablePodioTrackContainer(const MutablePodioTrackContainer& other);
  /// Move constructor
  /// @param other Source container
  MutablePodioTrackContainer(MutablePodioTrackContainer&& other) = default;

  /// Constructor from const container
  /// @param other Source container
  explicit MutablePodioTrackContainer(const ConstPodioTrackContainer& other);

  // BEGIN INTERFACE HELPER

 private:
  std::shared_ptr<const Acts::Surface> getOrCreateSurface(IndexType itrack) {
    std::shared_ptr<const Acts::Surface>& ptr = m_surfaces.at(itrack);
    if (!ptr) {
      ActsPodioEdm::Track track = m_collection->at(itrack);
      ptr = PodioUtil::convertSurfaceFromPodio(m_helper,
                                               track.getReferenceSurface());
    }
    return ptr;
  }

 public:
  /// Access component
  /// @param key Column key
  /// @param itrack Track index
  /// @return Component value
  std::any component_impl(Acts::HashedString key, IndexType itrack) {
    return PodioTrackContainerBase::component_impl<false>(*this, key, itrack);
  }

  /// Access component (const)
  /// @param key Column key
  /// @param itrack Track index
  /// @return Component value
  std::any component_impl(Acts::HashedString key, IndexType itrack) const {
    return PodioTrackContainerBase::component_impl<true>(*this, key, itrack);
  }

  /// Check if column exists
  /// @param key Column key
  /// @return True if column exists
  bool hasColumn_impl(Acts::HashedString key) const {
    return m_dynamic.contains(key);
  }

  /// Get container size
  /// @return Number of tracks
  std::size_t size_impl() const { return m_collection->size(); }

  /// Clear all tracks
  void clear() { m_collection->clear(); }

  // END INTERFACE HELPER

  /// Get reference surface
  /// @param itrack Track index
  /// @return Reference surface pointer
  const Acts::Surface* referenceSurface_impl(IndexType itrack) const {
    return m_surfaces.at(itrack).get();
  }

  /// Get particle hypothesis
  /// @param itrack Track index
  /// @return Particle hypothesis
  Acts::ParticleHypothesis particleHypothesis_impl(IndexType itrack) const {
    return PodioTrackContainerBase::particleHypothesis_impl(*this, itrack);
  }

  /// Set reference surface
  /// @param itrack Track index
  /// @param surface Reference surface
  void setReferenceSurface_impl(IndexType itrack,
                                std::shared_ptr<const Acts::Surface> surface) {
    auto track = m_collection->at(itrack);
    if (surface == nullptr) {
      track.setReferenceSurface({.surfaceType = PodioUtil::kNoSurface,
                                 .identifier = PodioUtil::kNoIdentifier});
      m_surfaces.at(itrack) = nullptr;
    } else {
      track.setReferenceSurface(
          PodioUtil::convertSurfaceToPodio(m_helper, *surface));
      m_surfaces.at(itrack) = std::move(surface);
    }
  }

 public:
  // BEGIN INTERFACE

  /// Add a new track
  /// @return Index of the added track
  IndexType addTrack_impl() {
    auto track = m_collection->create();
    PodioUtil::getReferenceSurfaceMutable(track).surfaceType =
        PodioUtil::kNoSurface;
    m_surfaces.emplace_back();
    for (const auto& [key, vec] : m_dynamic) {
      vec->add();
    }
    return m_collection->size() - 1;
  };

  /// Remove a track
  /// @param itrack Track index
  void removeTrack_impl(IndexType itrack);

  /// Add a column
  /// @tparam T Column value type
  /// @param key Column key
  template <typename T>
  constexpr void addColumn_impl(std::string_view key) {
    Acts::HashedString hashedKey = Acts::hashStringDynamic(key);
    m_dynamic.insert(
        {hashedKey, std::make_unique<podio_detail::DynamicColumn<T>>(key)});
  }

  /// Get track parameters
  /// @param itrack Track index
  /// @return Track parameters
  Parameters parameters(IndexType itrack) {
    return Parameters{
        PodioUtil::getDataMutable(m_collection->at(itrack)).parameters.data()};
  }

  /// Get track parameters (const)
  /// @param itrack Track index
  /// @return Track parameters
  ConstParameters parameters(IndexType itrack) const {
    return ConstParameters{
        m_collection->at(itrack).getData().parameters.data()};
  }

  /// Get track covariance
  /// @param itrack Track index
  /// @return Track covariance
  Covariance covariance(IndexType itrack) {
    return Covariance{
        PodioUtil::getDataMutable(m_collection->at(itrack)).covariance.data()};
  }

  /// Get track covariance (const)
  /// @param itrack Track index
  /// @return Track covariance
  ConstCovariance covariance(IndexType itrack) const {
    return ConstCovariance{
        m_collection->at(itrack).getData().covariance.data()};
  }

  /// Copy dynamic column from source
  /// @param dstIdx Destination track index
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

  /// Ensure dynamic columns match other container
  /// @param other Other container
  void ensureDynamicColumns_impl(const MutablePodioTrackContainer& other);

  /// Reserve storage
  void reserve(IndexType /*size*/) {}

  /// Get track collection
  /// @return Track collection
  ActsPodioEdm::TrackCollection& trackCollection() { return *m_collection; }

  /// Release into frame
  /// @param frame Destination frame
  /// @param suffix Collection name suffix
  void releaseInto(podio::Frame& frame, const std::string& suffix = "") {
    std::string s = suffix;
    if (!s.empty()) {
      s = "_" + s;
    }
    frame.put(std::move(m_collection), "tracks" + s);
    m_surfaces.clear();

    for (const auto& [key, col] : m_dynamic) {
      col->releaseInto(frame, "tracks" + s + "_extra__");
    }
  }

  /// Get dynamic keys
  /// @return Dynamic key range
  Acts::detail::DynamicKeyRange<podio_detail::DynamicColumnBase>
  dynamicKeys_impl() const {
    return PodioTrackContainerBase::dynamicKeys_impl(*this);
  }

  /// Set particle hypothesis
  /// @param itrack Track index
  /// @param particleHypothesis Particle hypothesis
  void setParticleHypothesis_impl(
      IndexType itrack, const Acts::ParticleHypothesis& particleHypothesis) {
    ActsPodioEdm::ParticleHypothesis pHypo;
    pHypo.absPdg = particleHypothesis.absolutePdg();
    pHypo.mass = particleHypothesis.mass();
    pHypo.absQ = particleHypothesis.absoluteCharge();
    m_collection->at(itrack).setParticleHypothesis(pHypo);
  }

  // END INTERFACE

 private:
  friend PodioTrackContainerBase;

  std::unique_ptr<ActsPodioEdm::TrackCollection> m_collection;
  std::vector<Acts::HashedString> m_dynamicKeys;
  std::unordered_map<Acts::HashedString,
                     std::unique_ptr<podio_detail::DynamicColumnBase>>
      m_dynamic;
};

static_assert(
    Acts::TrackContainerBackend<MutablePodioTrackContainer>,
    "MutablePodioTrackContainer does not fulfill TrackContainerBackend");

/// Read-only track container backend using podio for storage
class ConstPodioTrackContainer : public PodioTrackContainerBase {
 public:
  /// Constructor from collection
  /// @param helper Conversion helper
  /// @param collection Track collection
  ConstPodioTrackContainer(const PodioUtil::ConversionHelper& helper,
                           const ActsPodioEdm::TrackCollection& collection)
      : PodioTrackContainerBase{helper}, m_collection{&collection} {
    // Not much we can do to recover dynamic columns here
    populateSurfaceBuffer(m_helper, *m_collection, m_surfaces);
  }

  /// Constructor from frame
  /// @param helper Conversion helper
  /// @param frame Podio frame
  /// @param suffix Collection name suffix
  ConstPodioTrackContainer(const PodioUtil::ConversionHelper& helper,
                           const podio::Frame& frame,
                           const std::string& suffix = "")
      : PodioTrackContainerBase{helper} {
    std::string s = suffix.empty() ? suffix : "_" + suffix;
    std::string tracksKey = "tracks" + s;

    std::vector<std::string> available = frame.getAvailableCollections();
    if (!Acts::rangeContainsValue(available, tracksKey)) {
      throw std::runtime_error{"Track collection '" + tracksKey +
                               "' not found in frame"};
    }

    const auto* collection = frame.get(tracksKey);

    if (const auto* d =
            dynamic_cast<const ActsPodioEdm::TrackCollection*>(collection);
        d != nullptr) {
      m_collection = d;
    } else {
      throw std::runtime_error{"Unable to get collection " + tracksKey};
    }

    populateSurfaceBuffer(m_helper, *m_collection, m_surfaces);

    podio_detail::recoverDynamicColumns(frame, tracksKey, m_dynamic);
  }

  /// Get component implementation
  /// @param key Component key
  /// @param itrack Track index
  /// @return Component value
  std::any component_impl(Acts::HashedString key, IndexType itrack) const {
    return PodioTrackContainerBase::component_impl<true>(*this, key, itrack);
  }

  /// Check if column exists implementation
  /// @param key Column key
  /// @return True if column exists
  bool hasColumn_impl(Acts::HashedString key) const {
    return m_dynamic.contains(key);
  }

  /// Get size implementation
  /// @return Number of tracks
  std::size_t size_impl() const { return m_collection->size(); }

  /// Get reference surface implementation
  /// @param itrack Track index
  /// @return Reference surface pointer
  const Acts::Surface* referenceSurface_impl(IndexType itrack) const {
    return m_surfaces.at(itrack).get();
  }

  /// Get particle hypothesis implementation
  /// @param itrack Track index
  /// @return Particle hypothesis
  Acts::ParticleHypothesis particleHypothesis_impl(IndexType itrack) const {
    return PodioTrackContainerBase::particleHypothesis_impl(*this, itrack);
  }

  /// Get track parameters
  /// @param itrack Track index
  /// @return Track parameters
  ConstParameters parameters(IndexType itrack) const {
    return ConstParameters{
        m_collection->at(itrack).getData().parameters.data()};
  }

  /// Get track covariance
  /// @param itrack Track index
  /// @return Track covariance matrix
  ConstCovariance covariance(IndexType itrack) const {
    return ConstCovariance{
        m_collection->at(itrack).getData().covariance.data()};
  }

  /// Get the underlying track collection
  /// @return Track collection
  const ActsPodioEdm::TrackCollection& trackCollection() {
    return *m_collection;
  }

  /// Get dynamic keys implementation
  /// @return Range of dynamic column keys
  Acts::detail::DynamicKeyRange<podio_detail::ConstDynamicColumnBase>
  dynamicKeys_impl() const {
    return PodioTrackContainerBase::dynamicKeys_impl(*this);
  }

 private:
  friend PodioTrackContainerBase;

  const ActsPodioEdm::TrackCollection* m_collection;
  std::unordered_map<Acts::HashedString,
                     std::unique_ptr<podio_detail::ConstDynamicColumnBase>>
      m_dynamic;
  std::vector<Acts::HashedString> m_dynamicKeys;
};

static_assert(
    Acts::ConstTrackContainerBackend<ConstPodioTrackContainer>,
    "ConstPodioTrackContainer does not fulfill ConstTrackContainerBackend");

/// @}
}  // namespace ActsPlugins
