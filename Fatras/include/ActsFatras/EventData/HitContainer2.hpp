// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "ActsFatras/EventData/Hit.hpp"
#include "ActsFatras/EventData/ParticleContainer2.hpp"

#include <unordered_map>

namespace ActsFatras {

using HitIndex = std::uint32_t;
using HitIndexSubset = std::span<const HitIndex>;

class HitContainer2 {
 public:
  /// Type alias for particle index in container
  using Index = HitIndex;
  /// Type alias for subset of particle indices
  using IndexSubset = HitIndexSubset;

  /// Returns the number of hits in the container.
  /// @return The number of hits in the container.
  [[nodiscard]] std::uint32_t size() const noexcept { return m_hits.size(); }
  /// Checks if the container is empty.
  /// @return True if the container is empty, false otherwise.
  [[nodiscard]] bool empty() const noexcept { return m_hits.empty(); }

  /// Reserves space for the given number of hits.
  /// @param size The number of hits to reserve space for.
  void reserve(std::uint32_t size) noexcept;

  /// Clears the container, removing all hits and columns.
  void clear() noexcept;

  const std::vector<Hit>& hits() const noexcept { return m_hits; }

  const std::unordered_map<ParticleIndex, std::vector<HitIndex>>&
  hitsByParticles() const noexcept {
    return m_hitsByParticles;
  }

  const std::unordered_map<Acts::GeometryIdentifier, std::vector<HitIndex>>&
  hitsBySurfaces() const noexcept {
    return m_hitsBySurfaces;
  }

  Hit& push_back(const Hit& hit);

  Hit& emplace_back(Acts::GeometryIdentifier geometryId,
                    ParticleIndex particleId, const Acts::Vector4& pos4,
                    const Acts::Vector4& before4, const Acts::Vector4& after4,
                    std::int32_t index = -1);

  Hit& operator[](std::size_t index) { return m_hits[index]; }
  const Hit& operator[](std::size_t index) const { return m_hits[index]; }

  Hit& at(std::size_t index) { return m_hits.at(index); }
  const Hit& at(std::size_t index) const { return m_hits.at(index); }

  auto begin() noexcept { return m_hits.begin(); }
  auto end() noexcept { return m_hits.end(); }
  auto begin() const noexcept { return m_hits.begin(); }
  auto end() const noexcept { return m_hits.end(); }

  std::span<const HitIndex> hitIndicesByParticle(
      ParticleIndex particleId) const;
  std::span<const HitIndex> hitIndicesBySurface(
      Acts::GeometryIdentifier geometryId) const;

  /// Subset facade over arbitrary index sets.
  template <bool read_only>
  class Subset : public Acts::detail::ContainerSubset<
                     Subset<read_only>, Subset<true>, HitContainer2,
                     std::conditional_t<read_only, const Hit&, Hit&>,
                     std::span<const Index>, read_only> {
   public:
    /// Base class type
    using Base = Acts::detail::ContainerSubset<
        Subset<read_only>, Subset<true>, HitContainer2,
        std::conditional_t<read_only, const Hit&, Hit&>, std::span<const Index>,
        read_only>;

    using Base::Base;
  };

  /// Type alias for mutable subset of hits
  using MutableSubset = Subset<false>;
  /// Type alias for const subset of hits
  using ConstSubset = Subset<true>;

  /// Creates a mutable subset of hits from the given index subset.
  /// @param subset The index subset to create the subset from.
  /// @return A mutable subset of hits.
  MutableSubset subset(const IndexSubset& subset) noexcept {
    return MutableSubset(*this, subset);
  }
  /// Creates a const subset of hits from the given index subset.
  /// @param subset The index subset to create the subset from.
  /// @return A const subset of hits.
  ConstSubset subset(const IndexSubset& subset) const noexcept {
    return ConstSubset(*this, subset);
  }

  ConstSubset hitsByParticle(ParticleIndex particleId) const {
    return subset(hitIndicesByParticle(particleId));
  }
  ConstSubset hitsBySurface(Acts::GeometryIdentifier geometryId) const {
    return subset(hitIndicesBySurface(geometryId));
  }

 private:
  std::vector<Hit> m_hits;
  std::unordered_map<ParticleIndex, std::vector<HitIndex>> m_hitsByParticles;
  std::unordered_map<Acts::GeometryIdentifier, std::vector<HitIndex>>
      m_hitsBySurfaces;
};

}  // namespace ActsFatras
