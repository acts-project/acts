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

Acts::Result<std::vector<MpParameterResult>>
MillePedeResultReader::readParameters(
    const std::filesystem::path& mpFile) const {
  std::vector<MpParameterResult> res;
  std::ifstream resFile(mpFile);
  if (!resFile.is_open()) {
    ACTS_ERROR(" Failed to read the MP results file '" << mpFile << "'");
    return Acts::Result<std::vector<MpParameterResult>>::failure(
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
std::optional<MpParameterResult> MillePedeResultReader::parseMpLine(
    const std::string& resLine) const {
  // skip blank lines
  if (resLine.empty()) {
    return std::nullopt;
  }
  MpParameterResult par;
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
  if (!firstWord.empty() && firstWord == "Parameter") {
    ACTS_DEBUG(" Skipping the header line of the MP results file ");
    return std::nullopt;
  }
  sstr.seekg(0);
  sstr >> par.label >> par.val >> par.start;
  if (sstr.fail()) {
    ACTS_WARNING(" Failed to read a line of the MP results file ");
    return std::nullopt;
  }
  // the following three elements are not guaranteed
  // to exist, depending on the Pede configuration.
  if (!sstr.eof()) {
    sstr >> par.delta;
  }
  if (!sstr.eof()) {
    sstr >> par.sigma;
  }
  if (!sstr.eof()) {
    sstr >> par.nRecords;
  }
  return par;
}
