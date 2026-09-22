// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Mille/MillePedeSolver.hpp"

#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "ActsPlugins/Mille/MillePedeError.hpp"
#include "ActsPlugins/Mille/detail/RunSolverProcess.hpp"

#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <tuple>

namespace ActsPlugins::ActsToMille {

Acts::Result<MillePedeSolver::MpResult> MillePedeSolver::solve(
    const Config& cfg) const {
  ACTS_INFO("=== Proceeding to run Millepede-II alignment fit ===");

  // determine where the user wishes to run
  std::filesystem::path workDir = std::filesystem::current_path();

  if (!cfg.workDir.empty()) {
    workDir = cfg.workDir;
  }

  // create the work dir if it does not exist already
  if (!std::filesystem::exists(workDir)) {
    std::filesystem::create_directories(workDir);
  }

  ACTS_DEBUG("Will run the alignment fit in folder '"
             << std::filesystem::absolute(workDir) << "'");

  // make sure the steering file exists (or is empty)
  std::vector<std::string> mpArgs = cfg.extraOpts;

  if (std::filesystem::path steerPath = cfg.steeringFile; !steerPath.empty()) {
    if (steerPath.is_relative()) {
      steerPath = std::filesystem::current_path() / steerPath;
    }
    if (!std::filesystem::exists(steerPath)) {
      ACTS_ERROR("Steering file " << steerPath
                                  << " does not seem to exist - aborting!");
      return Acts::Result<MpResult>::failure(MillePedeError::SteeringNotFound);
    }
    mpArgs.push_back(steerPath);
  }

  ACTS_INFO(" Calling pede, this may take a while depending on problem size");
  // now run the fit
  ChildProcessStatus pedeProcessStatus =
      runSolverProcess("pede", mpArgs, workDir, cfg.redirectStdout);

  if (pedeProcessStatus == ChildProcessStatus::ProgNotFound) {
    ACTS_ERROR("Pede executable could not be found, aborting!");
    return Acts::Result<MpResult>::failure(
        MillePedeError::InstallationNotFound);
  }
  if (pedeProcessStatus == ChildProcessStatus::FailedWorkDir) {
    ACTS_ERROR("Could not navigate to the run directory `" << workDir
                                                           << "`, aborting!");
    return Acts::Result<MpResult>::failure(MillePedeError::SolverCrash);
  } else if (pedeProcessStatus != ChildProcessStatus::OK) {
    ACTS_ERROR("Pede invocation failed. Did not run alignment.");
    return Acts::Result<MpResult>::failure(MillePedeError::SolverCrash);
  }

  /// now check the detailed exit code - pede speaks fortranese, so we need to
  /// read this from a file

  // step 1: Ensure the file is actually there
  std::filesystem::path mpExit =
      workDir / std::filesystem::path("millepede.end");
  if (!std::filesystem::exists(mpExit)) {
    ACTS_ERROR(
        "Failed to find the Pede exit code file. The alignment has likely "
        "failed.");
    return Acts::Result<MpResult>::failure(MillePedeError::SolverCrash);
  }

  // step 2: Read and interpret the exit info file
  auto [exitCode, exitStatus, exitMessage] = readDetailedExit(mpExit);

  // step 3: Tell the user what happened

  if (exitStatus == MpExitStatus::Aborted) {
    ACTS_ERROR("Pede Aborted due to errors:\n   "
               << exitMessage << "\n"
               << "Please check the millepede log files in " << workDir
               << " for more information.");
    return Acts::Result<MpResult>::failure(MillePedeError::InvalidSolution);
  }

  else if (exitStatus == MpExitStatus::NoSolution) {
    ACTS_WARNING("Pede did not find a solution:\n   "
                 << exitMessage << "\n"
                 << "Please check the millepede log files in " << workDir
                 << " for more information.");
  }

  else if (exitStatus == MpExitStatus::SeriousWarnings) {
    ACTS_WARNING("Pede exited with severe warnings:\n   "
                 << exitMessage << "\n"
                 << "You should check the millepede log files in " << workDir
                 << " for more information.");
  }

  else if (exitStatus == MpExitStatus::TolerableWarnings) {
    ACTS_INFO("Pede exited with tolerable warnings:\n   "
              << exitMessage << "\n"
              << "You can check the millepede log files in " << workDir
              << " for more information.");
  } else {
    ACTS_INFO("Pede exited nominally:\n   "
              << exitMessage << "\n"
              << "You can check the millepede log files in " << workDir
              << " for more information.");
  }

  // step 4: Relocate output if configured by the user.

  auto outputResultFile =
      copyIfRequested(workDir / "millepede.res", cfg.resFileName);
  auto outputLogFile =
      copyIfRequested(workDir / "millepede.log", cfg.logFileName);
  auto outputHisFile =
      copyIfRequested(workDir / "millepede.his", cfg.histoFileName);
  auto outputEvFile =
      copyIfRequested(workDir / "millepede.eve", cfg.evFileName);

  /// the one file we **always** expect to exist is the result file.
  /// The others may not exist depending on the user steering options.
  if (outputResultFile.empty()) {
    ACTS_ERROR("Failed to find Millepede result file expected at "
               << workDir / "millepede.res");
    return Acts::Result<MpResult>::failure(MillePedeError::SolutionNotReadable);
  }

  ACTS_INFO("=== Congratulations, the alignment finished! ===");

  return MpResult{
      exitCode,         exitStatus, exitMessage,
      outputResultFile,  /// file containing parameter results
      outputLogFile,     /// log file
      outputHisFile,     /// file with validation histograms
      outputEvFile,      /// file with eigenvectors
  };
}

std::tuple<int, MillePedeSolver::MpExitStatus, std::string>
MillePedeSolver::readDetailedExit(const std::filesystem::path& mpend) const {
  std::ifstream exitCodeFile(mpend);  // open the Pede exit code file
  std::string statusMessage = "";
  int retCode = 0;
  // read the pede exit code from the file.
  exitCodeFile >> retCode;
  // some acrobatics to get a nice formatting of the remaining text
  std::stringstream remainingCont;
  remainingCont << exitCodeFile.rdbuf();
  statusMessage = remainingCont.str();
  statusMessage = statusMessage.substr(statusMessage.find_first_not_of(' '));
  statusMessage =
      statusMessage.substr(0, statusMessage.find_last_not_of(" \n") + 1);
  exitCodeFile.close();

  MpExitStatus status = interpretExit(retCode);
  return std::make_tuple(retCode, status, statusMessage);
}

MillePedeSolver::MpExitStatus MillePedeSolver::interpretExit(int theExitCode) {
  // implementation follows
  // https://millepede.pages.desy.de/millepede-ii/exit_code_page.html
  using enum MillePedeSolver::MpExitStatus;
  if (theExitCode < 0) {
    return NotFinishedOrCrashed;
  } else if (theExitCode == 0) {
    return NominalExit;
  } else if (theExitCode == 1) {
    return TolerableWarnings;
  } else if (theExitCode <= 4) {
    return SeriousWarnings;
  } else if (theExitCode == 5) {
    return NoSolution;
  } else {
    return Aborted;
  }
}

std::filesystem::path MillePedeSolver::copyIfRequested(
    const std::filesystem::path& originalLoc,
    const std::filesystem::path& userLoc) const {
  if (!std::filesystem::exists(originalLoc)) {
    return "";
  } else if (userLoc.empty()) {
    return originalLoc;
  } else {
    if (!std::filesystem::copy_file(
            originalLoc, userLoc,
            std::filesystem::copy_options::overwrite_existing)) {
      ACTS_WARNING(" Failed to copy output file " << originalLoc << " to "
                                                  << userLoc);
      return originalLoc;
    }
    return userLoc;
  }
}
}  // namespace ActsPlugins::ActsToMille
