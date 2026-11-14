// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsFatras/Synthetic/EventCsvWriter.hpp"

#include <cstdint>
#include <fstream>

namespace ActsFatras {

void Synthetic::writeEventCsv(const Event& event, const DetectorLayout& layout,
                              const std::string& prefix) {
  std::ofstream spFile(prefix + "_spacepoints.csv");
  spFile << "x,y,z,r,phi,layer,cylinder,particle,primary\n";
  const auto particleColumn =
      event.spacePoints.column<std::uint32_t>("particleId");
  const auto layerColumn = event.spacePoints.column<std::uint32_t>("layerId");
  for (const auto sp : event.spacePoints) {
    const std::uint32_t layer = sp.extra(layerColumn);
    const std::uint32_t particle = sp.extra(particleColumn);
    spFile << sp.x() << ',' << sp.y() << ',' << sp.z() << ',' << sp.r() << ','
           << sp.phi() << ',' << layer << ','
           << (layout.layers[layer].shape == SurfaceShape::Cylinder ? 1 : 0)
           << ',' << particle << ','
           << (event.particles[particle].primary() ? 1 : 0) << '\n';
  }

  std::ofstream particleFile(prefix + "_particles.csv");
  particleFile << "pt,eta,phi,d0,z0,charge,productionRadius,"
                  "productionZ,primary,numHits,generation,productionSurface\n";
  for (const GeneratedParticle& particle : event.particles) {
    particleFile << particle.pt << ',' << particle.eta << ',' << particle.phi
                 << ',' << particle.d0 << ',' << particle.z0 << ','
                 << particle.charge << ',' << particle.productionRadius << ','
                 << particle.productionZ << ',' << (particle.primary() ? 1 : 0)
                 << ',' << particle.numHits << ','
                 << static_cast<std::uint32_t>(particle.generation)
                 << ','
                 // -1 rather than the sentinel, which a reader would plot
                 << (particle.productionSurface == kNoSurface
                         ? -1
                         : static_cast<std::int64_t>(
                               particle.productionSurface))
                 << '\n';
  }
}

}  // namespace ActsFatras
