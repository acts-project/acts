// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsFatras/EventData/ProcessType.hpp"

#include <boost/container/flat_set.hpp>

namespace ActsExamples {

class SimVertexBarcode {
 public:
  using PrimaryVertexId = SimBarcode::PrimaryVertexId;
  using SecondaryVertexId = SimBarcode::SecondaryVertexId;
  using ParticleId = SimBarcode::ParticleId;
  using GenerationId = SimBarcode::GenerationId;
  using SubParticleId = SimBarcode::SubParticleId;

  explicit constexpr SimVertexBarcode(SimBarcode barcode)
      : m_id(barcode.vertexId()) {}

  constexpr SimVertexBarcode() = default;

  /// Return the barcode.
  constexpr SimBarcode barcode() const { return m_id; }

  /// Return the primary vertex identifier.
  constexpr PrimaryVertexId vertexPrimary() const {
    return m_id.vertexPrimary();
  }
  /// Return the secondary vertex identifier.
  constexpr SecondaryVertexId vertexSecondary() const {
    return m_id.vertexSecondary();
  }
  /// Return the generation identifier.
  constexpr GenerationId generation() const { return m_id.generation(); }

  /// Create a new barcode with a different primary vertex identifier.
  [[nodiscard]]
  constexpr SimVertexBarcode withVertexPrimary(PrimaryVertexId id) const {
    return SimVertexBarcode(m_id.withVertexPrimary(id));
  }
  /// Create a new barcode with a different secondary vertex identifier.
  [[nodiscard]]
  constexpr SimVertexBarcode withVertexSecondary(SecondaryVertexId id) const {
    return SimVertexBarcode(m_id.withVertexSecondary(id));
  }
  /// Create a new barcode with a different generation identifier.
  [[nodiscard]]
  constexpr SimVertexBarcode withGeneration(GenerationId id) const {
    return SimVertexBarcode(m_id.withGeneration(id));
  }

  std::size_t hash() const { return m_id.hash(); }

 private:
  /// The vertex ID
  /// Note that only primary, secondary and generation should be set
  SimBarcode m_id;

  friend constexpr bool operator<(SimVertexBarcode lhs, SimVertexBarcode rhs) {
    return lhs.m_id < rhs.m_id;
  }

  friend constexpr bool operator==(SimVertexBarcode lhs, SimVertexBarcode rhs) {
    return lhs.m_id == rhs.m_id;
  }

  friend inline std::ostream& operator<<(std::ostream& os,
                                         SimVertexBarcode idx) {
    return os << idx.m_id;
  }
};

/// A simulated vertex e.g. from a physics process.
struct SimVertex {
  /// The vertex ID
  SimVertexBarcode id = SimVertexBarcode(SimBarcode::Invalid());
  /// The vertex four-position
  Acts::Vector4 position4 = Acts::Vector4::Zero();
  /// The vertex process type
  ActsFatras::ProcessType process = ActsFatras::ProcessType::eUndefined;
  /// The incoming particles into the vertex
  SimBarcodeContainer incoming;
  /// The outgoing particles from the vertex
  SimBarcodeContainer outgoing;

  /// Construct the vertex from a position and optional process type.
  ///
  /// @param position4_ the vertex four-position
  /// @param process_ the process type that generated this vertex
  ///
  /// Associated particles are left empty by default and must be filled by the
  /// user after construction.
  SimVertex(
      SimVertexBarcode id_, const Acts::Vector4& position4_,
      ActsFatras::ProcessType process_ = ActsFatras::ProcessType::eUndefined)
      : id(id_), position4(position4_), process(process_) {}
  // explicitly default rule-of-five.
  SimVertex() = default;
  SimVertex(const SimVertex&) = default;
  SimVertex(SimVertex&&) = default;
  SimVertex& operator=(const SimVertex&) = default;
  SimVertex& operator=(SimVertex&&) = default;

  constexpr SimVertexBarcode vertexId() const { return id; }
  /// The vertex three-position.
  auto position() const { return position4.head<3>(); }
  /// The vertex time.
  double time() const { return position4[3]; }
};

namespace detail {
struct CompareVertexId {
  using is_transparent = void;
  constexpr bool operator()(const SimVertex& lhs, const SimVertex& rhs) const {
    return lhs.vertexId() < rhs.vertexId();
  }
  constexpr bool operator()(SimVertexBarcode lhs, const SimVertex& rhs) const {
    return lhs < rhs.vertexId();
  }
  constexpr bool operator()(const SimVertex& lhs, SimVertexBarcode rhs) const {
    return lhs.vertexId() < rhs;
  }
};
}  // namespace detail

/// Store vertices ordered by vertex identifier.
using SimVertexContainer =
    ::boost::container::flat_set<SimVertex, detail::CompareVertexId>;

}  // namespace ActsExamples

// specialize std::hash so Barcode can be used e.g. in an unordered_map
namespace std {
template <>
struct hash<ActsExamples::SimVertexBarcode> {
  auto operator()(ActsExamples::SimVertexBarcode barcode) const noexcept {
    return barcode.hash();
  }
};
}  // namespace std
