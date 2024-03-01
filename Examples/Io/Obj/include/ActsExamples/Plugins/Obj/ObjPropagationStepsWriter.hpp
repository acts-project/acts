// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Propagator/detail/SteppingLogger.hpp"
#include "ActsExamples/Framework/WriterT.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <fstream>

namespace ActsExamples {

/// @class ObjPropagationStepsWriter
///
/// Write out the steps of test propgations for stepping validation
/// Writes one file per event with form:
///
///     event000000001-propagation-steps.obj
///     event000000002-propagation-steps.obj
///
/// One Thread per write call and hence thread safe
template <typename step_t>
class ObjPropagationStepsWriter
    : public WriterT<std::vector<std::vector<step_t>>> {
 public:
  struct Config {
    std::string collection;           ///< which collection to write
    std::string outputDir;            ///< where to place output files
    double outputScalor = 1.0;        ///< scale output values
    std::size_t outputPrecision = 6;  ///< floating point precision
  };

  /// Constructor with arguments
  ///
  /// @param cfg configuration struct
  /// @param level Output logging level
  ObjPropagationStepsWriter(const Config& cfg,
                            Acts::Logging::Level level = Acts::Logging::INFO)
      : WriterT<std::vector<std::vector<step_t>>>(cfg.collection,
                                                  "ObjSpacePointWriter", level),
        m_cfg(cfg) {
    if (m_cfg.collection.empty()) {
      throw std::invalid_argument("Missing input collection");
    }
  }

  /// Virtual destructor
  ~ObjPropagationStepsWriter() override = default;

  /// End-of-run hook
  ProcessCode finalize() override { return ActsExamples::ProcessCode::SUCCESS; }

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;  ///!< Internal configuration representation

 protected:
  /// This implementation holds the actual writing method
  /// and is called by the WriterT<>::write interface
  ProcessCode writeT(
      const AlgorithmContext& context,
      const std::vector<std::vector<step_t>>& stepCollection) override {
    // open per-event file
    std::string path = ActsExamples::perEventFilepath(
        m_cfg.outputDir, "propagation-steps.obj", context.eventNumber);
    std::ofstream os(path, std::ofstream::out | std::ofstream::trunc);
    if (!os) {
      throw std::ios_base::failure("Could not open '" + path + "' to write");
    }

    // Initialize the vertex counter
    unsigned int vCounter = 0;

    for (auto& steps : stepCollection) {
      // At least three points to draw
      if (steps.size() > 2) {
        // We start from one
        ++vCounter;
        for (auto& step : steps) {
          // Write the space point
          os << "v " << m_cfg.outputScalor * step.position.x() << " "
             << m_cfg.outputScalor * step.position.y() << " "
             << m_cfg.outputScalor * step.position.z() << '\n';
        }
        // Write out the line - only if we have at least two points created
        std::size_t vBreak = vCounter + steps.size() - 1;
        for (; vCounter < vBreak; ++vCounter) {
          os << "l " << vCounter << " " << vCounter + 1 << '\n';
        }
      }
    }
    return ActsExamples::ProcessCode::SUCCESS;
  }
};

}  // namespace ActsExamples
