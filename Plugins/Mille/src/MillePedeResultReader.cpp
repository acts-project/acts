// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Mille/MillePedeResultReader.hpp"

#include "Acts/Utilities/Result.hpp"
#include "ActsPlugins/Mille/MillePedeError.hpp"

#include <fstream>
#include <optional>
#include <system_error>

using namespace ActsPlugins::ActsToMille;

Acts::Result<std::vector<mpParameterResult>>
MillePedeResultReader::readParameters(
    const std::filesystem::path& theResultFile) const {
  std::vector<mpParameterResult> res;
  std::ifstream resFile(theResultFile);
  if (!resFile.is_open()) {
    ACTS_ERROR(" Failed to read the MP results file at " << theResultFile);
    return Acts::Result<std::vector<mpParameterResult>>::failure(
        MillePedeError::SolutionNotReadable);
  }
  while (!resFile.eof()) {
    std::string resLine;
    std::getline(resFile, resLine);
    auto parseOut = parseMpLine(resLine);
    if (parseOut.has_value()) {
      res.push_back(*parseOut);
    }
  }
  resFile.close();
  return res;
}
std::optional<mpParameterResult> MillePedeResultReader::parseMpLine(
    const std::string& resLine) const {
  mpParameterResult par;
  par.sigma = -1;
  par.nRecords = -1;
  std::stringstream sstr(resLine);
  // skip comment lines
  std::string firstWord = "";
  sstr >> firstWord;
  if (!firstWord.empty() && firstWord.starts_with("#")) {
    ACTS_DEBUG(" Skipping a commented line in the MP results file ");
    return std::nullopt;
  }
  sstr.seekg(0);
  sstr >> par.label >> par.val >> par.start >> par.delta;
  if (sstr.fail()) {
    ACTS_WARNING(" Failed to read a line of the MP results file ");
    return std::nullopt;
  }
  if (!sstr.eof()) {
    sstr >> par.sigma;
  }
  if (!sstr.eof()) {
    sstr >> par.nRecords;
  }
  return par;
}
