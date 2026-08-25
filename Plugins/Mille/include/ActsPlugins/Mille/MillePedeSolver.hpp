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

namespace ActsPlugins::ActsToMille {

/// class wrapping an external call to the 'pede' solver
/// program of the Millepede-II alignment toolkit.
class MillePedeSolver {
 public:
  /// @brief abstract summary of the exit codes returned by pede.
  enum class mpExitStatus {
    notFinishedOrCrashed,  /// job is still running or crashed before exiting
    nominalExit,           /// nominal exit
    tolerableWarnings,     /// exit with tolerable warnings, considered ok
    seriousWarnings,       /// exit with serious warnings, should investigate
    noSolution,            /// exit without solution (usually: rank deficit)
    aborted                /// aborted due to errors
  };

  /// @brief configuration for running pede.
  /// Currently, most details are delegated to the steering file syntax.
  struct Config {
    std::string steeringFile =
        "pedeSteerMaster.txt";                /// steering file with run options
    std::filesystem::path workDir = "";       /// directory to run in
    std::vector<std::string> extraOpts = {};  /// extra CLI options
    std::string resFileName =
        "";  /// destination for result file - if empty, keep original
    std::string logFileName =
        "";  /// destination for the log file - if empty, keep original
    std::string histoFileName =
        "";  /// destination for the histogram file - if empty, keep original
    std::string evFileName =
        "";  /// destination for the eigenvector file - if empty, keep original
  };

  /// @brief package the result of the alignment fit
  struct mpResult {
    int exitCode = -1;  /// raw pede exit code
    mpExitStatus exitStatus =
        mpExitStatus::notFinishedOrCrashed;  /// summary exit status
    std::string exitMessage = "";            /// detailed exit message
    std::filesystem::path resultsFile;  /// file containing parameter results
    std::filesystem::path logFile;      /// log file
    std::filesystem::path histoFile;    /// file with validation histograms
    std::filesystem::path evFile;       /// file with eigenvectors
  };

  /// @brief constructor - nothing to do as minimal internal state carried
  explicit MillePedeSolver(Acts::Logging::Level level = Acts::Logging::INFO)
      : m_logger(Acts::getDefaultLogger("MillePedeSolver", level)) {}

  /// @brief Runs the solving.
  /// Can take some time for large fits.
  /// Will invoke pede, await the exit, and parse
  /// the output.
  /// @param cfg: The configuration to use
  Acts::Result<mpResult> solve(Config cfg) const;

 private:
  /// @brief translation of the detailed pede code to a summary status
  /// Will translate the range of ~30 possible MP exit codes to a simplified
  /// enum value that gives a high-level summary of the status.
  /// @param theExitCode: The integer exit code found in the millepede.end file
  static mpExitStatus interpretExit(int theExitCode);

  /// @brief Reads the `millepede.end` file and parses its content
  /// @param mpend: The file to read, assumed that the user has checked for existence before
  /// @return a tuple containing the original integer exit code, a simplified exit status enum,
  /// and the additional status message emitted by pede.
  std::tuple<int, mpExitStatus, std::string> readDetailedExit(
      const std::filesystem::path& mpend) const;

  /// @brief Utility method to check if an output exists, copy it if the user requested a relocation,
  /// and returns the final resolved output location, which will be:
  /// - if the output does not exist: empty path
  /// - if the output does exist and was not relocated: original location
  /// - if the output does exist, and was relocated: user-specified location
  /// - if the output does exist, was requested for relocation but the copy
  /// failed: original location
  /// @param originalLoc: Original expected location of the file
  /// @param userLoc: Desired final location of the file. If empty, no copy will be made.
  std::filesystem::path copyIfRequested(
      const std::filesystem::path& originalLoc,
      const std::filesystem::path& userLoc) const;

  std::unique_ptr<const Acts::Logger> m_logger;  /// logger

  /// Private access to the logger
  const Acts::Logger& logger() const { return *m_logger; }
};
}  // namespace ActsPlugins::ActsToMille
