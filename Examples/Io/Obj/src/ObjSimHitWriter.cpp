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

#include <unsupported/Eigen/Splines>

namespace {

template <typename input_vector_type>
std::vector<Acts::Vector3> interpolatedPoints(
    const std::vector<input_vector_type>& inputs, unsigned int nPoints) {
  Eigen::MatrixXd points(3, inputs.size());
  for (unsigned int i = 0; i < inputs.size(); ++i) {
    points.col(i) = inputs[i].template head<3>().transpose();
  }
  Eigen::Spline<Acts::ActsScalar, 3> spline3D =
      Eigen::SplineFitting<Eigen::Spline<Acts::ActsScalar, 3>>::Interpolate(
          points, 2);

  std::vector<Acts::Vector3> output;
  Acts::ActsScalar step = 1. / (nPoints - 1);
  for (unsigned int i = 0; i < nPoints; ++i) {
    Acts::ActsScalar t = i * step;
    Eigen::Vector3d point = spline3D(t);
    output.push_back(Acts::Vector3(point[0], point[1], point[2]));
  }
  return output;
}

}  // namespace

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
  std::string pathSimTrajectory = perEventFilepath(
      m_cfg.outputDir, m_cfg.outputStem + "_trajectory.obj", ctx.eventNumber);

  std::ofstream osHits(pathSimHit, std::ofstream::out | std::ofstream::trunc);
  std::ofstream osTrajectory(pathSimTrajectory,
                             std::ofstream::out | std::ofstream::trunc);

  if (!osHits || !osTrajectory) {
    throw std::ios_base::failure("Could not open '" + pathSimHit + "' or '" +
                                 pathSimTrajectory + "' to write");
  }

  // Only hit plotting
  if (!m_cfg.drawConnections) {
    // Write data from internal immplementation
    for (const auto& simHit : simHits) {
      // local simhit information in global coord.
      const Acts::Vector4& globalPos4 = simHit.fourPosition();

      osHits << "v " << globalPos4[Acts::ePos0] / Acts::UnitConstants::mm << " "
             << globalPos4[Acts::ePos1] / Acts::UnitConstants::mm << " "
             << globalPos4[Acts::ePos2] / Acts::UnitConstants::mm << std::endl;

    }  // end simHit loop
  } else {
    // We need to associate first
    std::unordered_map<uint64_t, std::vector<Acts::Vector4>> particleHits;
    // Pre-loop over hits ... write those below threshold
    for (const auto& simHit : simHits) {
      Acts::ActsScalar momentum = simHit.momentum4Before().head<3>().norm();
      if (momentum < m_cfg.momentumThreshold) {
        ACTS_VERBOSE("Skipping : Hit below threshold: " << momentum);
        continue;
      } else if (momentum < m_cfg.momentumThresholdTraj) {
        ACTS_VERBOSE(
            "Skipping (trajectory): Hit below threshold: " << momentum);
        osHits << "v "
               << simHit.fourPosition()[Acts::ePos0] / Acts::UnitConstants::mm
               << " "
               << simHit.fourPosition()[Acts::ePos1] / Acts::UnitConstants::mm
               << " "
               << simHit.fourPosition()[Acts::ePos2] / Acts::UnitConstants::mm
               << std::endl;
        continue;
      }
      ACTS_VERBOSE("Accepting : Hit above threshold: " << momentum);

      if (particleHits.find(simHit.particleId().value()) ==
          particleHits.end()) {
        particleHits[simHit.particleId().value()] = {};
      }
      particleHits[simHit.particleId().value()].push_back(
          simHit.fourPosition());
    }
    // Draw loop
    unsigned int lOffset = 1;
    for (auto& [pId, pHits] : particleHits) {
      // Draw the particle hits
      std::sort(pHits.begin(), pHits.end(),
                [](const Acts::Vector4& a, const Acts::Vector4& b) {
                  return a[Acts::eTime] < b[Acts::eTime];
                });

      osHits << "o particle_" << pId << std::endl;
      for (const auto& hit : pHits) {
        osHits << "v " << hit[Acts::ePos0] / Acts::UnitConstants::mm << " "
               << hit[Acts::ePos1] / Acts::UnitConstants::mm << " "
               << hit[Acts::ePos2] / Acts::UnitConstants::mm << std::endl;
      }
      osHits << '\n';

      // Interpolate the points
      std::vector<Acts::Vector3> trajectory;
      if (pHits.size() < 3) {
        for (const auto& hit : pHits) {
          trajectory.push_back(hit.template head<3>());
        }
      } else {
        trajectory =
            interpolatedPoints(pHits, pHits.size() * m_cfg.nInterpolatedPoints);
      }

      osTrajectory << "o particle_trajectory_" << pId << std::endl;
      for (const auto& hit : trajectory) {
        osTrajectory << "v " << hit[Acts::ePos0] / Acts::UnitConstants::mm
                     << " " << hit[Acts::ePos1] / Acts::UnitConstants::mm << " "
                     << hit[Acts::ePos2] / Acts::UnitConstants::mm << std::endl;
      }
      // Draw the line
      for (unsigned int iv = lOffset + 1; iv < lOffset + trajectory.size();
           ++iv) {
        osTrajectory << "l " << iv - 1 << " " << iv << '\n';
      }
      osTrajectory << '\n';
      // Increase the offset count
      lOffset += trajectory.size();
    }
  }
  osHits.close();
  osTrajectory.close();

  return ActsExamples::ProcessCode::SUCCESS;
}

ActsExamples::ProcessCode ActsExamples::ObjSimHitWriter::finalize() {
  return ActsExamples::ProcessCode::SUCCESS;
}
