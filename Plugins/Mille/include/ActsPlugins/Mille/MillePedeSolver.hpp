// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "ActsPlugins/Mille/MillePedeError.hpp"

#include <filesystem>
#include <vector>

namespace ActsPlugins {

/// class wrapping an external call to the 'pede' solver
/// program of the Millepede-II alignment toolkit.
class MillePedeSolver {
 public:
  /// @brief abstract summary of the exit codes returned by pede.
  enum class ExitStatus {
    NotFinishedOrCrashed,  /// job is still running or crashed before exiting
    NominalExit,           /// nominal exit
    TolerableWarnings,     /// exit with tolerable warnings, considered ok
    SeriousWarnings,       /// exit with serious warnings, should investigate
    NoSolution,            /// exit without solution (usually: rank deficit)
    Aborted                /// Aborted due to errors
  };

  /// @brief configuration for running pede.
  /// Currently, most details are delegated to the steering file syntax.
  struct Config {
    /// steering file with run options
    std::string steeringFile = "pedeSteerMaster.txt";

    /// directory to run in. Default: current work dir
    std::optional<std::filesystem::path> workDir;

    /// extra CLI options
    std::vector<std::string> extraOpts = {};

    /// destination for result file - default: keep original
    std::optional<std::string> resFileName;
    /// destination for the cout/cerr printout
    /// from pede - default: print to terminal
    std::optional<std::string> redirectStdout;
    /// destination for the log file - default: keep original
    std::optional<std::string> logFileName;
    /// destination for the histogram file - default: keep original
    std::optional<std::string> histoFileName;
    /// destination for the eigenvector file - default: keep original
    std::optional<std::string> evFileName;
  };

  /// @brief package the result of the alignment fit
  struct Result {
    int exitCode = -1;  /// raw pede exit code
    ExitStatus exitStatus =
        ExitStatus::NotFinishedOrCrashed;  /// summary exit status
    std::string exitMessage = "";          /// detailed exit message
    std::filesystem::path resultsFile;     /// file containing parameter results
    std::filesystem::path logFile;         /// log file
    std::filesystem::path histoFile;       /// file with validation histograms
    std::filesystem::path evFile;          /// file with eigenvectors
  };

  /// @brief constructor - nothing to do as minimal internal state carried
  explicit MillePedeSolver(std::unique_ptr<const Acts::Logger> _logger =
                               Acts::getDefaultLogger("MillePedeSolver",
                                                      Acts::Logging::INFO))
      : m_logger(std::move(_logger)) {}

  /// @brief Runs the solving.
  /// Can take some time for large fits.
  /// Will invoke pede, await the exit, and parse
  /// the output.
  /// @param cfg: The configuration to use
  Acts::Result<Result> solve(const Config& cfg) const;

 private:
  std::unique_ptr<const Acts::Logger> m_logger;

  /// Private access to the logger
  const Acts::Logger& logger() const { return *m_logger; }
};
}  // namespace ActsPlugins
