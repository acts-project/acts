// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Json/BetheHeitlerApproxJsonReader.hpp"

#include "ActsPlugins/Json/BetheHeitlerApproxJsonConverter.hpp"

#include <fstream>
#include <stdexcept>

#include <nlohmann/json.hpp>

namespace Acts {

namespace {

void readJsonFile(const std::string& filepath, nlohmann::json& j) {
  std::ifstream in(filepath);
  if (!in) {
    throw std::invalid_argument("Could not open JSON file '" + filepath + "'");
  }
  in >> j;
  in.close();
}

}  // namespace

std::shared_ptr<const PolynomialBetheHeitlerApprox>
loadBetheHeitlerApproxFromJson(const std::string& filepath, bool clampToRange,
                               double noChangeLimit,
                               double singleGaussianLimit) {
  nlohmann::json j;
  readJsonFile(filepath, j);

  // Parse multiple x0 ranges
  std::vector<BetheHeitlerApproxJsonConverter::RangeData> ranges;

  if (j.contains("ranges")) {
    for (const auto& jrange : j["ranges"]) {
      BetheHeitlerApproxJsonConverter::RangeData range;
      from_json(jrange, range);
      ranges.push_back(range);
    }
  }

  if (ranges.empty()) {
    throw std::invalid_argument(
        "JSON file must contain 'ranges' array with at least one range");
  }

  // Sort ranges by lowX0 to ensure they are properly ordered
  std::sort(ranges.begin(), ranges.end(),
            [](const auto& a, const auto& b) { return a.lowX0 < b.lowX0; });

  // Build the PolynomialBetheHeitlerApprox with the ranges
  return std::make_shared<const PolynomialBetheHeitlerApprox>(
      ranges, clampToRange, noChangeLimit, singleGaussianLimit);
}

}  // namespace Acts
