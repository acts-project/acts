// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/detail/DynamicColumn.hpp"
#include "Acts/Plugins/Podio/PodioDynamicColumns.hpp"
#include "Acts/Plugins/Podio/PodioUtil.hpp"
#include "ActsPodioEdm/Track.h"
#include "ActsPodioEdm/TrackCollection.h"
#include "ActsPodioEdm/TrackInfo.h"

#include <mutex>
#include <stdexcept>
#include <type_traits>

#include <podio/Frame.h>

namespace Acts {

class MutablePodioTrackContainer;
class ConstPodioTrackContainer;

template <>
struct IsReadOnlyTrackContainer<MutablePodioTrackContainer> : std::false_type {
};

template <>
struct IsReadOnlyTrackContainer<ConstPodioTrackContainer> : std::true_type {};

class PodioTrackContainerBase {
 public:
  using IndexType = MultiTrajectoryTraits::IndexType;
  static constexpr auto kInvalid = MultiTrajectoryTraits::kInvalid;
  static constexpr auto MeasurementSizeMax =
      MultiTrajectoryTraits::MeasurementSizeMax;

  using Parameters =
      typename detail_lt::Types<eBoundSize, false>::CoefficientsMap;
  using Covariance =
      typename detail_lt::Types<eBoundSize, false>::CovarianceMap;

  using ConstParameters =
      typename detail_lt::Types<eBoundSize, true>::CoefficientsMap;
  using ConstCovariance =
      typename detail_lt::Types<eBoundSize, true>::CovarianceMap;

 protected:
  PodioTrackContainerBase(const PodioUtil::ConversionHelper& helper)
      : m_helper{helper} {}

  template <bool EnsureConst, typename T>
  static std::any component_impl(T& instance, HashedString key,
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
      dataPtr = &track.data();
    }
    auto& data = *dataPtr;
    switch (key) {
      case "tipIndex"_hash:
        return &data.tipIndex;
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

  static void populateSurfaceBuffer(
      const PodioUtil::ConversionHelper& helper,
      const ActsPodioEdm::TrackCollection& collection,
      std::vector<std::shared_ptr<const Surface>>& surfaces) noexcept {
    surfaces.reserve(collection.size());
    for (ActsPodioEdm::Track track : collection) {
      surfaces.push_back(PodioUtil::convertSurfaceFromPodio(
          helper, track.getReferenceSurface()));
    }
  }

  std::reference_wrapper<const PodioUtil::ConversionHelper> m_helper;
  std::vector<std::shared_ptr<const Surface>> m_surfaces;
};

class MutablePodioTrackContainer : public PodioTrackContainerBase {
 public:
  MutablePodioTrackContainer(const PodioUtil::ConversionHelper& helper)
      : PodioTrackContainerBase{helper},
        m_collection{std::make_unique<ActsPodioEdm::TrackCollection>()} {
    populateSurfaceBuffer(m_helper, *m_collection, m_surfaces);
  }

  MutablePodioTrackContainer(const MutablePodioTrackContainer& other);
  MutablePodioTrackContainer(MutablePodioTrackContainer&& other) = default;

  MutablePodioTrackContainer(const ConstPodioTrackContainer& other);

  // BEGIN INTERFACE HELPER

 private:
  std::shared_ptr<const Surface> getOrCreateSurface(IndexType itrack) {
    std::shared_ptr<const Surface>& ptr = m_surfaces.at(itrack);
    if (!ptr) {
      ActsPodioEdm::Track track = m_collection->at(itrack);
      ptr = PodioUtil::convertSurfaceFromPodio(m_helper,
                                               track.getReferenceSurface());
    }
    return ptr;
  }

 public:
  std::any component_impl(HashedString key, IndexType itrack) {
    return PodioTrackContainerBase::component_impl<false>(*this, key, itrack);
  }

  std::any component_impl(HashedString key, IndexType itrack) const {
    return PodioTrackContainerBase::component_impl<true>(*this, key, itrack);
  }

  bool hasColumn_impl(HashedString key) const {
    return m_dynamic.find(key) != m_dynamic.end();
  }

  size_t size_impl() const { return m_collection->size(); }
  // END INTERFACE HELPER

  const Surface* referenceSurface_impl(IndexType itrack) const {
    return m_surfaces.at(itrack).get();
  }

  void setReferenceSurface_impl(IndexType itrack,
                                std::shared_ptr<const Surface> surface) {
    auto track = m_collection->at(itrack);
    track.setReferenceSurface(
        PodioUtil::convertSurfaceToPodio(m_helper, *surface));
    m_surfaces.at(itrack) = std::move(surface);
  }

 public:
  // BEGIN INTERFACE

  IndexType addTrack_impl() {
    auto track = m_collection->create();
    track.referenceSurface().surfaceType = PodioUtil::kNoSurface;
    m_surfaces.emplace_back();
    for (const auto& [key, vec] : m_dynamic) {
      vec->add();
    }
    return m_collection->size() - 1;
  };

  void removeTrack_impl(IndexType itrack);

  template <typename T>
  constexpr void addColumn_impl(const std::string& key) {
    m_dynamic.insert({hashString(key),
                      std::make_unique<podio_detail::DynamicColumn<T>>(key)});
  }

  Parameters parameters(IndexType itrack) {
    return Parameters{m_collection->at(itrack).data().parameters.data()};
  }

  ConstParameters parameters(IndexType itrack) const {
    return ConstParameters{
        m_collection->at(itrack).getData().parameters.data()};
  }

  Covariance covariance(IndexType itrack) {
    return Covariance{m_collection->at(itrack).data().covariance.data()};
  }

  ConstCovariance covariance(IndexType itrack) const {
    return ConstCovariance{
        m_collection->at(itrack).getData().covariance.data()};
  }

  // @TODO What's the equivalent of this?
  void copyDynamicFrom_impl(IndexType dstIdx,
                            const MutablePodioTrackContainer& src,
                            IndexType srcIdx);

  void ensureDynamicColumns_impl(const MutablePodioTrackContainer& other);

  void reserve(IndexType /*size*/) {}

  ActsPodioEdm::TrackCollection& trackCollection() { return *m_collection; }

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

  // END INTERFACE

 private:
  friend PodioTrackContainerBase;

  std::unique_ptr<ActsPodioEdm::TrackCollection> m_collection;
  std::unordered_map<HashedString,
                     std::unique_ptr<podio_detail::DynamicColumnBase>>
      m_dynamic;
};

ACTS_STATIC_CHECK_CONCEPT(TrackContainerBackend, MutablePodioTrackContainer);

class ConstPodioTrackContainer : public PodioTrackContainerBase {
 public:
  ConstPodioTrackContainer(const PodioUtil::ConversionHelper& helper,
                           const ActsPodioEdm::TrackCollection& collection)
      : PodioTrackContainerBase{helper}, m_collection{&collection} {
    populateSurfaceBuffer(m_helper, *m_collection, m_surfaces);
  }

  ConstPodioTrackContainer(const PodioUtil::ConversionHelper& helper,
                           const podio::Frame& frame,
                           const std::string& suffix = "")
      : PodioTrackContainerBase{helper} {
    std::string s = suffix.empty() ? suffix : "_" + suffix;
    std::string tracksKey = "tracks" + s;

    std::vector<std::string> available = frame.getAvailableCollections();
    if (std::find(available.begin(), available.end(), tracksKey) ==
        available.end()) {
      throw std::runtime_error{"Track collection '" + tracksKey +
                               "'not found in frame"};
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

    // let's find dynamic columns
    using types =
        std::tuple<int32_t, int64_t, uint32_t, uint64_t, float, double>;

    for (const auto& col : available) {
      std::string prefix = tracksKey + "_extra__";
      size_t p = col.find(prefix);
      if (p == std::string::npos) {
        continue;
      }
      std::string dynName = col.substr(prefix.size());
      const podio::CollectionBase* coll = frame.get(col);

      std::unique_ptr<podio_detail::ConstDynamicColumnBase> up;

      std::apply(
          [&](auto... args) {
            auto inner = [&](auto arg) {
              if (up) {
                return;
              }
              using T = decltype(arg);
              const auto* dyn =
                  dynamic_cast<const podio::UserDataCollection<T>*>(coll);
              if (dyn == nullptr) {
                return;
              }
              up = std::make_unique<podio_detail::ConstDynamicColumn<T>>(
                  dynName, *dyn);
            };

            ((inner(args)), ...);
          },
          types{});

      if (!up) {
        throw std::runtime_error{"Dynamic column '" + dynName +
                                 "' is not of allowed type"};
      }

      m_dynamic.insert({hashString(dynName), std::move(up)});
    }
  }

  std::any component_impl(HashedString key, IndexType itrack) const {
    return PodioTrackContainerBase::component_impl<true>(*this, key, itrack);
  }

  bool hasColumn_impl(HashedString key) const {
    return m_dynamic.find(key) != m_dynamic.end();
  }

  size_t size_impl() const { return m_collection->size(); }

  const Surface* referenceSurface_impl(IndexType itrack) const {
    return m_surfaces.at(itrack).get();
  }

  ConstParameters parameters(IndexType itrack) const {
    return ConstParameters{
        m_collection->at(itrack).getData().parameters.data()};
  }

  ConstCovariance covariance(IndexType itrack) const {
    return ConstCovariance{
        m_collection->at(itrack).getData().covariance.data()};
  }

  const ActsPodioEdm::TrackCollection& trackCollection() {
    return *m_collection;
  }

 private:
  friend PodioTrackContainerBase;

  const ActsPodioEdm::TrackCollection* m_collection;
  std::unordered_map<HashedString,
                     std::unique_ptr<podio_detail::ConstDynamicColumnBase>>
      m_dynamic;
};

ACTS_STATIC_CHECK_CONCEPT(ConstTrackContainerBackend, ConstPodioTrackContainer);

}  //  namespace Acts
