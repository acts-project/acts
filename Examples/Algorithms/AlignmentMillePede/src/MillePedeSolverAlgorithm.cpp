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

namespace ActsExamples {

ProcessCode MillePedeSolverAlgorithm::finalize() {
  ACTS_INFO("=== Proceeding to run Millepede-II alignment fit ===");

  ActsPlugins::MillePedeSolver solver(logger().clone());

  const auto& solverResult = solver.solve(m_cfg.solverConfig);

  if (!solverResult.ok()) {
    ACTS_ERROR("=== Alignment FAILED! ===");
    return ProcessCode::ABORT;
  }

  const auto& readRes =
      ActsPlugins::readMillePedeResult(solverResult->resultsFile, logger());

  if (!readRes.ok()) {
    ACTS_ERROR("Failed to read the MP result file!");
    return ProcessCode::ABORT;
  }
  ACTS_INFO("=== Finished alignment and readback with exit code "
            << solverResult->exitCode << " (" << solverResult->exitMessage
            << "), found " << readRes->size() << " alignment parameters ===");

  // for now, do not further process the alignment results, as we are in
  // finalize().

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
