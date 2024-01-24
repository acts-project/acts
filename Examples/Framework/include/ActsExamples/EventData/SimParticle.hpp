// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Utilities/GroupBy.hpp"
#include "ActsFatras/EventData/Particle.hpp"

#include <boost/container/flat_set.hpp>

namespace ActsExamples {

class SimParticle : public ::ActsFatras::Particle {
 public:
  using ::ActsFatras::Particle::Particle;

  SimParticle(Particle&& particle) : Particle(std::move(particle)) {}
  SimParticle(const Particle& particle) : Particle(particle) {}

  SimParticle* parent() { return m_parent; }
  const SimParticle* parent() const { return m_parent; }
  void setParent(SimParticle* parent) { m_parent = parent; }

  auto childrenBegin() { return m_children.begin(); }
  auto childrenBegin() const { return m_children.begin(); }
  auto childrenEnd() { return m_children.end(); }
  auto childrenEnd() const { return m_children.end(); }

  auto children() const { return makeRange(childrenBegin(), childrenEnd()); }

  void addChild(SimParticle& child) { m_children.push_back(&child); }

 private:
  SimParticle* m_parent = nullptr;
  std::vector<SimParticle*> m_children = {};
};

namespace detail {
struct CompareParticleId {
  using is_transparent = void;
  constexpr bool operator()(const SimParticle& lhs,
                            const SimParticle& rhs) const {
    return lhs.particleId() < rhs.particleId();
  }
  constexpr bool operator()(ActsFatras::Barcode lhs,
                            const SimParticle& rhs) const {
    return lhs < rhs.particleId();
  }
  constexpr bool operator()(const SimParticle& lhs,
                            ActsFatras::Barcode rhs) const {
    return lhs.particleId() < rhs;
  }
};
struct PrimaryVertexIdGetter {
  constexpr ActsFatras::Barcode operator()(const SimParticle& particle) const {
    return ActsFatras::Barcode().setVertexPrimary(
        particle.particleId().vertexPrimary());
  }
};
struct SecondaryVertexIdGetter {
  constexpr ActsFatras::Barcode operator()(const SimParticle& particle) const {
    return ActsFatras::Barcode()
        .setVertexPrimary(particle.particleId().vertexPrimary())
        .setVertexSecondary(particle.particleId().vertexSecondary());
  }
};
}  // namespace detail

// using SimParticle = ::ActsFatras::Particle;
/// Store particles ordered by particle identifier.
using SimParticleContainer =
    ::boost::container::flat_set<SimParticle, detail::CompareParticleId>;

/// Iterate over groups of particles belonging to the same primary vertex.
inline GroupBy<SimParticleContainer::const_iterator,
               detail::PrimaryVertexIdGetter>
groupByPrimaryVertex(const SimParticleContainer& container) {
  return makeGroupBy(container, detail::PrimaryVertexIdGetter());
}

/// Iterate over groups of particles belonging to the same secondary vertex.
///
/// For each primary vertex, this yields one group of particles belonging
/// directly to the primary vertex (secondary vertex id 0) and a group for
/// each secondary vertex.
inline GroupBy<SimParticleContainer::const_iterator,
               detail::SecondaryVertexIdGetter>
groupBySecondaryVertex(const SimParticleContainer& container) {
  return makeGroupBy(container, detail::SecondaryVertexIdGetter());
}

void graphvizSimParticleContainer(std::ostream& os,
                                  const SimParticleContainer& container);

}  // namespace ActsExamples
