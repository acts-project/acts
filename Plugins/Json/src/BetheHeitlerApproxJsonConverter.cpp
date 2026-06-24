// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Json/BetheHeitlerApproxJsonConverter.hpp"

#include "Acts/Utilities/RangeXD.hpp"

#include <fstream>
#include <stdexcept>

void Acts::to_json(nlohmann::json& j,
                   const PolynomialBetheHeitlerApprox::PolyData& data) {
  j["weight_coeffs"] = data.weightCoeffs;
  j["mean_coeffs"] = data.meanCoeffs;
  j["var_coeffs"] = data.varCoeffs;
}

void from_json(const nlohmann::json& j,
               PolynomialBetheHeitlerApprox::PolyData& data) {
  data.weightCoeffs = j.at("weight_coeffs").get<std::vector<double>>();
  data.meanCoeffs = j.at("mean_coeffs").get<std::vector<double>>();
  data.varCoeffs = j.at("var_coeffs").get<std::vector<double>>();
}

void to_json(nlohmann::json& j,
             const PolynomialBetheHeitlerApprox::RangeData& data) {
  j["low_x0"] = data.range.min();
  j["high_x0"] = data.range.max();
  j["transform"] = data.transform;

  nlohmann::json components = nlohmann::json::array();
  for (const auto& cmp : data.data) {
    nlohmann::json jcmp;
    to_json(jcmp, cmp);
    components.push_back(jcmp);
  }
  j["components"] = components;
}

void from_json(const nlohmann::json& j,
               PolynomialBetheHeitlerApprox::RangeData& data) {
  data.range = Range1D<double>{j.at("low_x0").get<double>(),
                               j.at("high_x0").get<double>()};
  data.transform = j.value("transform", true);

  if (j.contains("components")) {
    data.data.clear();
    for (const auto& jcmp : j["components"]) {
      PolynomialBetheHeitlerApprox::PolyData component;
      from_json(jcmp, component);
      data.data.push_back(component);
    }
  } else if (j.contains("weight_coeffs")) {
    PolynomialBetheHeitlerApprox::PolyData component;
    from_json(j, component);
    data.data = {component};
  } else {
    throw std::runtime_error(
        "JSON range data must contain either 'components' array or "
        "'weight_coeffs' (flat format)");
  }
}

PolynomialBetheHeitlerApprox loadBetheHeitlerApproxFromJson(
    const std::string& filepath, bool clampToRange, double noChangeLimit,
    double singleGaussianLimit) {
  std::ifstream in(filepath);
  if (!in) {
    throw std::invalid_argument("Could not open JSON file '" + filepath + "'");
  }
  nlohmann::json j;
  in >> j;

  if (!j.contains("ranges")) {
    throw std::invalid_argument(
        "JSON file must contain 'ranges' array with at least one range");
  }

  std::vector<PolynomialBetheHeitlerApprox::RangeData> ranges;
  for (const auto& jrange : j["ranges"]) {
    PolynomialBetheHeitlerApprox::RangeData range;
    from_json(jrange, range);
    ranges.push_back(range);
  }

  if (ranges.empty()) {
    throw std::invalid_argument(
        "JSON file must contain 'ranges' array with at least one range");
  }

  return PolynomialBetheHeitlerApprox(std::move(ranges), clampToRange,
                                      noChangeLimit, singleGaussianLimit);
}

}  // namespace Acts
