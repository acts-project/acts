// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/Track.hpp"
#include "Acts/EventData/detail/DynamicColumn.hpp"

#include <unordered_map>

namespace Acts {

namespace detail_vtc {

class VectorTrackContainerBase {
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
  VectorTrackContainerBase() = default;

  VectorTrackContainerBase(const VectorTrackContainerBase& other)
      : m_tipIndex{other.m_tipIndex},
        m_params{other.m_params},
        m_cov{other.m_cov},
        m_referenceSurfaces{other.m_referenceSurfaces} {
    for (const auto& [key, value] : other.m_dynamic) {
      m_dynamic.insert({key, value->clone()});
    }
  };

  VectorTrackContainerBase(VectorTrackContainerBase&& other) = default;

  // BEGIN INTERFACE HELPER

  template <bool EnsureConst, typename T>
  static std::any component_impl(T& instance, HashedString key,
                                 IndexType itrack) {
    using namespace Acts::HashedStringLiteral;
    if constexpr (EnsureConst) {
      static_assert(std::is_const_v<std::remove_reference_t<T>>,
                    "Is not const");
    }

    using namespace Acts::HashedStringLiteral;
    switch (key) {
      case "tipIndex"_hash:
        return &instance.m_tipIndex[itrack];
      case "params"_hash:
        return &instance.m_params[itrack];
      case "cov"_hash:
        return &instance.m_cov[itrack];
      case "referenceSurface"_hash:
        return &instance.m_referenceSurfaces[itrack];
      case "nMeasurements"_hash:
        return &instance.m_nMeasurements[itrack];
      case "nHoles"_hash:
        return &instance.m_nHoles[itrack];
      default:
        auto it = instance.m_dynamic.find(key);
        if (it == instance.m_dynamic.end()) {
          throw std::runtime_error("Unable to handle this component");
        }

        std::conditional_t<EnsureConst, const detail::DynamicColumnBase*,
                           detail::DynamicColumnBase*>
            col = it->second.get();
        assert(col && "Dynamic column is null");
        return col->get(itrack);
    }
  }

 public:
  constexpr bool hasColumn_impl(HashedString key) const {
    using namespace Acts::HashedStringLiteral;
    switch (key) {
      default:
        return m_dynamic.find(key) != m_dynamic.end();
    }
  }

  std::size_t size_impl() const { return m_tipIndex.size(); }
  // END INTERFACE HELPER

  std::vector<IndexType> m_tipIndex;
  std::vector<typename detail_lt::Types<eBoundSize>::Coefficients> m_params;
  std::vector<typename detail_lt::Types<eBoundSize>::Covariance> m_cov;
  std::vector<std::shared_ptr<const Surface>> m_referenceSurfaces;

  std::vector<unsigned int> m_nMeasurements;
  std::vector<unsigned int> m_nHoles;

  std::unordered_map<HashedString, std::unique_ptr<detail::DynamicColumnBase>>
      m_dynamic;
};

}  // namespace detail_vtc

class VectorTrackContainer;
template <>
struct IsReadOnlyTrackContainer<VectorTrackContainer> : std::false_type {};

class VectorTrackContainer final : public detail_vtc::VectorTrackContainerBase {
 public:
  VectorTrackContainer() : VectorTrackContainerBase{} {}
  VectorTrackContainer(const VectorTrackContainer& other) = default;
  VectorTrackContainer(VectorTrackContainer&&) = default;

 public:
  // BEGIN INTERFACE

  std::any component_impl(HashedString key, IndexType itrack) {
    return detail_vtc::VectorTrackContainerBase::component_impl<false>(
        *this, key, itrack);
  }

  std::any component_impl(HashedString key, IndexType itrack) const {
    return detail_vtc::VectorTrackContainerBase::component_impl<true>(
        *this, key, itrack);
  }

  IndexType addTrack_impl() {
    m_tipIndex.emplace_back();

    m_params.emplace_back();
    m_cov.emplace_back();
    m_referenceSurfaces.emplace_back();

    m_nMeasurements.emplace_back();
    m_nHoles.emplace_back();

    // dynamic columns
    for (auto& [key, vec] : m_dynamic) {
      vec->add();
    }

    return m_tipIndex.size() - 1;
  }

  template <typename T>
  constexpr void addColumn_impl(const std::string& key) {
    m_dynamic.insert(
        {hashString(key), std::make_unique<detail::DynamicColumn<T>>()});
  }

  Parameters parameters(IndexType itrack) {
    return Parameters{m_params[itrack].data()};
  }

  ConstParameters parameters(IndexType itrack) const {
    return ConstParameters{m_params[itrack].data()};
  }

  Covariance covariance(IndexType itrack) {
    return Covariance{m_cov[itrack].data()};
  }

  ConstCovariance covariance(IndexType itrack) const {
    return ConstCovariance{m_cov[itrack].data()};
  }

  // END INTERFACE
};

class ConstVectorTrackContainer;
template <>
struct IsReadOnlyTrackContainer<ConstVectorTrackContainer> : std::true_type {};

class ConstVectorTrackContainer final
    : public detail_vtc::VectorTrackContainerBase {
 public:
  ConstVectorTrackContainer() : VectorTrackContainerBase{} {}

  ConstVectorTrackContainer(const ConstVectorTrackContainer& other) = default;
  ConstVectorTrackContainer(const VectorTrackContainer& other)
      : VectorTrackContainerBase{other} {}

  ConstVectorTrackContainer(ConstVectorTrackContainer&&) = default;
  ConstVectorTrackContainer(VectorTrackContainer&& other)
      : VectorTrackContainerBase{std::move(other)} {}

 public:
  // BEGIN INTERFACE

  std::any component_impl(HashedString key, IndexType itrack) const {
    return detail_vtc::VectorTrackContainerBase::component_impl<true>(
        *this, key, itrack);
  }

  ConstParameters parameters(IndexType itrack) const {
    return ConstParameters{m_params[itrack].data()};
  }

  // END INTERFACE
};

}  // namespace Acts
