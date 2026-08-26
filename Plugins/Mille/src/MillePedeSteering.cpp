// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Mille/MillePedeSteering.hpp"

#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "ActsPlugins/Mille/MillePedeError.hpp"

#include <filesystem>
#include <format>
#include <fstream>
#include <map>

using namespace ActsPlugins;

namespace {
using Strategy = ActsPlugins::MillePedeSolutionStrategy;
using EqualityConstraint = ActsPlugins::MillePedeEqualityConstraint;
using Config = ActsPlugins::MillePedeSteeringConfig;

void addInputBlock(std::ofstream& out, const Config& conf,
                   const Acts::Logger& logger) {
  if (conf.inputFiles.empty()) {
    ACTS_WARNING(
        "You are about to generate a steering file without any input! Suit "
        "yourself, but this will probably not align anything...");
  }
  out << "*** Input files *** " << std::endl;
  out << "Cfiles " << std::endl;
  for (const std::string& f : conf.inputFiles) {
    std::filesystem::path p = f;
    /// always write absolute paths, to avoid issues when
    /// we run pede in a sub-directory
    if (p.is_relative()) {
      p = std::filesystem::current_path() / p;
    }
    out << p.string() << std::endl;
  }
  out << std::endl;
}

void addSolutionMethod(std::ofstream& out, const Config& conf,
                       const Acts::Logger& logger) {
  static const std::map<Strategy, std::string> stratNames{
      {Strategy::Inversion, "inversion"},
      {Strategy::Diagonalization, "diagonalization"},
      {Strategy::Decomposition, "decomposition"},
      {Strategy::FullMINRES, "fullMINRES"},
      {Strategy::SparseMINRES, "sparseMINRES"},
      {Strategy::FullMINRES_QLP, "fullMINRES-QLP"},
      {Strategy::SparseMINRES_QLP, "sparseMINRES-QLP"},
      {Strategy::FullLAPACK, "fullLAPACK"},
      {Strategy::UnpackedLAPACK, "unpackedLAPACK"},
      {Strategy::SparsePARDISO, "sparsePARDISO"},
  };
  auto found = stratNames.find(conf.strategy);
  if (found == stratNames.end()) {
    ACTS_ERROR(
        "Strategy specified in the steering not yet supported in "
        "addSolutionMethod - this is a bug, please report to the devs");
    return;
  }
  out << "*** Solution method *** " << std::endl;
  out << std::format("method {}  {}  {}", found->second, conf.minIterations,
                     conf.convergenceLimit)
      << std::endl
      << std::endl;
}

void addConstraint(std::ofstream& out, const EqualityConstraint& constraint,
                   const Acts::Logger& /*logger*/) {
  out << std::format("Constraint {}", constraint.constraint);
  for (const auto& [label, weight] : constraint.labelsAndWeights) {
    out << std::format("    {}     {}", label, weight) << std::endl;
  }
  out << std::endl;
}

void addConstraintsBlock(std::ofstream& out, const Config& conf,
                         const Acts::Logger& logger) {
  out << "* Constraints to respect in the fit " << std::endl;
  for (const EqualityConstraint& con : conf.constraints) {
    addConstraint(out, con, logger);
  }
  out << std::endl;
}

void addOthers(std::ofstream& out, const Config& conf,
               const Acts::Logger& /*logger*/) {
  out << "*** Other steering flags " << std::endl;
  out << std::format("entries {} {} {}", std::get<0>(conf.entriesCut),
                     std::get<1>(conf.entriesCut), std::get<2>(conf.entriesCut))
      << std::endl;
  out << std::format("outlierdownweighting {} ", conf.outlierDownweighting)
      << std::endl;
  out << std::format("dwfractioncut {} ", conf.downweightFractionCut)
      << std::endl;
  out << std::format("threads {} {} ", conf.nOMPthreads, conf.nIOthreads)
      << std::endl;
  out << std::format("matiter {} ", conf.matIter) << std::endl;
  out << std::format("printcounts {} ", conf.printCounts) << std::endl;
  out << std::format("chisqcut {} {} ", conf.chi2Cut.first, conf.chi2Cut.second)
      << std::endl;
  if (conf.monitorResiduals) {
    out << " monitorresiduals" << std::endl;
  }
  if (conf.monitorPulls) {
    out << " monitorpulls" << std::endl;
  }
  if (conf.skipEmptyCons) {
    out << " skipemptycons" << std::endl;
  }
  if (conf.countRecords) {
    out << " countrecords" << std::endl;
  }
  out << std::endl;
  out << "*** User-specified config flags " << std::endl;
  for (const std::string& extraLine : conf.extraLines) {
    out << extraLine << std::endl;
  }
  out << std::endl;
}
}  // namespace

Acts::Result<std::filesystem::path> ActsPlugins::generateMillePedeSteeringFile(
    const std::filesystem::path& destination, const Config& config,
    const Acts::Logger& logger) {
  std::ofstream theSteer(destination);
  if (!theSteer.is_open()) {
    ACTS_ERROR("Unable to create the MP steering file " << destination);
    return Acts::Result<std::filesystem::path>::failure(
        MillePedeError::UnableToWriteSteering);
  }

  theSteer << "*** Auto-generated Millepede steering file from ACTS *** "
           << std::endl;

  std::filesystem::path consFile = "";
  if (!config.constraints.empty()) {
    consFile = destination;
    consFile.replace_extension("cons");
    theSteer << "*** Constraints file " << std::endl;
    theSteer << consFile << std::endl;
    std::ofstream theCons(consFile);
    if (!theSteer.is_open()) {
      ACTS_ERROR("Unable to create the MP constraints file " << consFile);
      return Acts::Result<std::filesystem::path>::failure(
          MillePedeError::UnableToWriteSteering);
    }
    addConstraintsBlock(theCons, config, logger());
    theCons.close();
  }

  addInputBlock(theSteer, config, logger());
  addSolutionMethod(theSteer, config, logger());
  addOthers(theSteer, config, logger());

  theSteer.close();
  return destination;
}
