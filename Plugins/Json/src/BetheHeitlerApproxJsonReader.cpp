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

std::shared_ptr<const AtlasBetheHeitlerApprox> loadBetheHeitlerApproxFromJson(
    const std::string& filepath, bool clampToRange, double noChangeLimit,
    double singleGaussianLimit) {
  nlohmann::json j;
  readJsonFile(filepath, j);

  // Parse limits from JSON
  double transform = true;

  if (j.contains("transform")) {
    transform = j["transform"].get<bool>();
  }

  // Parse multiple x0 ranges
  std::vector<BetheHeitlerApproxJsonConverter::RangeData> ranges;

  if (j.contains("ranges")) {
    for (const auto& jrange : j["ranges"]) {
      BetheHeitlerApproxJsonConverter::RangeData range;
      from_json(jrange, range);
      ranges.push_back(range);
    }
  }

  // Build low and high data from ranges
  // For now, just use the first range as low and the last as high
  // This will be expanded when we support multiple ranges
  AtlasBetheHeitlerApprox::Data lowData;
  AtlasBetheHeitlerApprox::Data highData;

  if (!ranges.empty()) {
    // Use the first range for low data
    for (const auto& r : ranges) {
      AtlasBetheHeitlerApprox::PolyData poly;
      poly.weightCoeffs = r.weightCoeffs;
      poly.meanCoeffs = r.meanCoeffs;
      poly.varCoeffs = r.varCoeffs;
      lowData.push_back(poly);
    }

    // For now, use the same as high data
    // In the future, we'll support truly distinct ranges
    highData = lowData;
  }

  // Use first range low_x0 as low limit, last range high_x0 as high limit
  double lowLimit = 0.10;   // Default
  double highLimit = 0.20;  // Default

  if (!ranges.empty()) {
    lowLimit = ranges.front().lowX0;
    highLimit = ranges.back().highX0;
  }

  // For now, use low data for both (same as current behavior)
  // When we expand to multiple ranges, this will need to be restructured
  return std::make_shared<const AtlasBetheHeitlerApprox>(
      lowData, lowData, transform, transform, lowLimit, highLimit, clampToRange,
      noChangeLimit, singleGaussianLimit);
}

}  // namespace Acts
