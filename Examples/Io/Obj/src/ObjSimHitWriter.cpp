// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Obj/ObjSimHitWriter.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Utilities/Paths.hpp"
#include "ActsFatras/EventData/Barcode.hpp"
#include "ActsFatras/EventData/Hit.hpp"

#include <fstream>
#include <stdexcept>
#include <unordered_map>
#include <vector>

ActsExamples::ObjSimHitWriter::ObjSimHitWriter(
    const ActsExamples::ObjSimHitWriter::Config& config,
    Acts::Logging::Level level)
    : WriterT(config.inputSimHits, "ObjSimHitWriter", level), m_cfg(config) {
  // inputSimHits is already checked by base constructor
  if (m_cfg.outputStem.empty()) {
    throw std::invalid_argument("Missing output filename stem");
  }
}

ActsExamples::ProcessCode ActsExamples::ObjSimHitWriter::writeT(
    const AlgorithmContext& ctx, const ActsExamples::SimHitContainer& simHits) {
  // ensure exclusive access to tree/file while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);

  // open per-event file for all simhit components
  std::string pathSimHit = perEventFilepath(
      m_cfg.outputDir, m_cfg.outputStem + ".obj", ctx.eventNumber);

  std::ofstream os(pathSimHit, std::ofstream::out | std::ofstream::trunc);
  if (!os) {
    throw std::ios_base::failure("Could not open '" + pathSimHit +
                                 "' to write");
  }
  // Initialize the vertex counter
  unsigned int vCounter = 0;

  if (!m_cfg.drawConnections) {
    // Write data from internal immplementation
    for (const auto& simHit : simHits) {
      // local simhit information in global coord.
      const Acts::Vector4& globalPos4 = simHit.fourPosition();

      os << "v " << globalPos4[Acts::ePos0] / Acts::UnitConstants::mm << " "
         << globalPos4[Acts::ePos1] / Acts::UnitConstants::mm << " "
         << globalPos4[Acts::ePos2] / Acts::UnitConstants::mm << std::endl;

    }  // end simHit loop
  } else {
    // We need to associate first
    std::unordered_map<uint64_t, std::vector<Acts::Vector4>> particleHits;
    // pre-loop over hits
    for (const auto& simHit : simHits) {
      if (particleHits.find(simHit.particleId().value()) ==
          particleHits.end()) {
        particleHits[simHit.particleId().value()] = {};
      }
      particleHits[simHit.particleId().value()].push_back(
          simHit.fourPosition());
    }
    // Draw loop
    unsigned int lOffset = 1;
    for (auto [pId, pHits] : particleHits) {
      // Draw the particle hits
      std::sort(pHits.begin(), pHits.end(),
                [](const Acts::Vector4& a, const Acts::Vector4& b) {
                  return a[Acts::eTime] < b[Acts::eTime];
                });
      os << "o particle_" << pId << std::endl;
      for (const auto& hit : pHits) {
        os << "v " << hit[Acts::ePos0] / Acts::UnitConstants::mm << " "
           << hit[Acts::ePos1] / Acts::UnitConstants::mm << " "
           << hit[Acts::ePos2] / Acts::UnitConstants::mm << std::endl;
      }
      // Draw the line
      for (unsigned int iv = lOffset + 1; iv < lOffset + pHits.size(); ++iv) {
        os << "l " << iv - 1 << " " << iv << '\n';
      }
      os << '\n';
      // Increase the offset count
      lOffset += pHits.size();
    }
  }
  os.close();

  return ActsExamples::ProcessCode::SUCCESS;
}

ActsExamples::ProcessCode ActsExamples::ObjSimHitWriter::finalize() {
  return ActsExamples::ProcessCode::SUCCESS;
}
