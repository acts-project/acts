// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/AlignmentMillePede/MillePedeSolverAlgorithm.hpp"

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsPlugins/Mille/MillePedeResultReader.hpp"
#include "ActsPlugins/Mille/MillePedeSolver.hpp"
#include "ActsPlugins/Mille/MillePedeSteering.hpp"

namespace ActsExamples {

ProcessCode MillePedeSolverAlgorithm::finalize() {
  ACTS_INFO("=== Generating alignment steering file for Millepede-II ===");
  ActsPlugins::ActsToMille::MillePedeSteering steer(logger().level());
  auto confPath = steer.generateSteeringFile(m_cfg.solverConfig.steeringFile,
                                             m_cfg.steeringConfig);
  if (confPath.empty()) {
    ACTS_ERROR(
        "Failed to generate the steering file, aborting alignment attempt");
    return ProcessCode::ABORT;
  }
  ACTS_INFO("=== Proceeding to run Millepede-II alignment fit ===");

  ActsPlugins::ActsToMille::MillePedeSolver solver(logger().level());

  const auto& res = solver.solve(m_cfg.solverConfig);

  if (!res.ok()) {
    ACTS_ERROR("=== Alignment FAILED! ===");
    return ProcessCode::ABORT;
  }

  ActsPlugins::ActsToMille::MillePedeResultReader reader(logger().level());

  auto readRes = reader.readParameters(res->resultsFile);

  if (!readRes.ok()) {
    ACTS_ERROR("Failed to read the MP result file!");
    return ProcessCode::ABORT;
  }
  ACTS_INFO("=== Finished alignment with exit code "
            << res->exitCode << " (" << res->exitMessage << "), found "
            << readRes->size() << " alignment parameters ===");
  // for now, do not further process the alignment results, as we are in
  // finalize()

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
