// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/EventData/SimParticle.hpp"

#include <set>

namespace ActsExamples {
namespace {
std::string vid(unsigned int vtx) {
  return "V" + std::to_string(vtx);
}
std::string pid(const SimParticle& particle) {
  return "P" + std::to_string(particle.particleId().value());
};

std::string plabel(const SimParticle& particle) {
  std::stringstream ss;
  ss << particle.pdg() << "\\n(" << particle.particleId() << ")";
  return ss.str();
};

void printParticle(std::ostream& os, const SimParticle& parent) {
  for (const auto* child : parent.children()) {
    os << pid(parent) << " -> " << pid(*child);

    if (parent.particleId().vertexSecondary() ==
        child->particleId().vertexSecondary()) {
      os << " [style=dashed]";
    }
    os << ";\n";

    os << pid(*child) << " [label=\"" << plabel(*child) << "\"];\n";

    std::cout << parent.particleId() << " -> " << child->particleId()
              << std::endl;

    printParticle(os, *child);
  }
}

}  // namespace

void graphvizSimParticleContainer(std::ostream& os,
                                  const SimParticleContainer& container) {
  os << "digraph Event {\n";

  std::set<unsigned int> primaryVertices;

  // find root particles and group by PV
  for (const auto& particle : container) {
    if (particle.particleId().generation() > 0) {
      continue;
    }
    primaryVertices.insert(particle.particleId().vertexPrimary());
    os << vid(particle.particleId().vertexPrimary()) << " -> " << pid(particle)
       << ";\n";
    os << pid(particle) << " [label=\"" << plabel(particle) << "\"];\n";

    printParticle(os, particle);
  }

  for (unsigned int vtx : primaryVertices) {
    os << vid(vtx) << " [shape=box,label=\"PV" << vtx << "\"];\n";
  }

  os << "}";
}
}  // namespace ActsExamples
