// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/Track.hpp"

#include <unordered_map>

namespace Acts {

namespace detail_vtc {

class VectorTrackContainerBase {
 public:
  using IndexType = MultiTrajectoryTraits::IndexType;
  static constexpr auto kInvalid = MultiTrajectoryTraits::kInvalid;
  static constexpr auto MeasurementSizeMax =
      MultiTrajectoryTraits::MeasurementSizeMax;

 protected:
  VectorTrackContainerBase() = default;

  VectorTrackContainerBase(const VectorTrackContainerBase& other) {
    for (const auto& [key, value] : other.m_dynamic) {
      m_dynamic.insert({key, value->clone()});
    }
  };

  struct DynamicColumnBase {
    virtual ~DynamicColumnBase() = 0;

    virtual std::any get(size_t i) = 0;
    virtual std::any get(size_t i) const = 0;

    virtual void add() = 0;
    virtual void clear() = 0;

    virtual std::unique_ptr<DynamicColumnBase> clone() const = 0;
  };

  template <typename T>
  struct DynamicColumn : public DynamicColumnBase {
    ~DynamicColumn() override = default;

    std::any get(size_t i) override {
      assert(i < m_vector.size() && "DynamicColumn out of bounds");
      return &m_vector[i];
    }

    std::any get(size_t i) const override {
      assert(i < m_vector.size() && "DynamicColumn out of bounds");
      return &m_vector[i];
    }

    void add() override { m_vector.emplace_back(); }
    void clear() override { m_vector.clear(); }

    std::unique_ptr<DynamicColumnBase> clone() const override {
      return std::make_unique<DynamicColumn<T>>(*this);
    }

    std::vector<T> m_vector;
  };

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
        return instance.m_tipIndex[itrack];
      default:
        auto it = instance.m_dynamic.find(key);
        if (it == instance.m_dynamic.end()) {
          throw std::runtime_error("Unable to handle this component");
        }
        auto& col = it->second;
        assert(col && "Dynamic column is null");
        return col->get(itrack);
    }
  }

  template <typename T>
  static constexpr bool hasColumn_impl(T& instance, HashedString key) {
    using namespace Acts::HashedStringLiteral;
    switch (key) {
      default:
        return instance.m_dynamic.find(key) != instance.m_dynamic.end();
    }
  }
  // END INTERFACE HELPER

  std::vector<IndexType> m_tipIndex;

  std::unordered_map<HashedString, std::unique_ptr<DynamicColumnBase>>
      m_dynamic;
};

}  // namespace detail_vtc

class VectorTrackContainer;
template <>
struct isReadOnlyTrackContainer<VectorTrackContainer> : std::false_type {};

class VectorTrackContainer final : public detail_vtc::VectorTrackContainerBase {
 public:
  VectorTrackContainer() = default;
  VectorTrackContainer(const VectorTrackContainer& other)
      : VectorTrackContainerBase{other} {}

  VectorTrackContainer(VectorTrackContainer&&) = default;
  VectorTrackContainer& operator=(const VectorTrackContainer&) = default;
  VectorTrackContainer& operator=(VectorTrackContainer&&) = default;

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
    return m_tipIndex.size() - 1;
  }

  // END INTERFACE
};

class ConstVectorTrackContainer;
template <>
struct isReadOnlyTrackContainer<ConstVectorTrackContainer> : std::true_type {};

class ConstVectorTrackContainer final
    : public detail_vtc::VectorTrackContainerBase {
 public:
  ConstVectorTrackContainer() = default;

  ConstVectorTrackContainer(const ConstVectorTrackContainer& other)
      : VectorTrackContainerBase{other} {}

  ConstVectorTrackContainer(const VectorTrackContainer& other)
      : VectorTrackContainerBase{other} {}

  ConstVectorTrackContainer(ConstVectorTrackContainer&&) = default;
  ConstVectorTrackContainer& operator=(const ConstVectorTrackContainer&) =
      default;
  ConstVectorTrackContainer& operator=(ConstVectorTrackContainer&&) = default;

 public:
  // BEGIN INTERFACE

  std::any component_impl(HashedString key, IndexType itrack) const {
    return detail_vtc::VectorTrackContainerBase::component_impl<true>(
        *this, key, itrack);
  }

  // END INTERFACE
};

}  // namespace Acts
