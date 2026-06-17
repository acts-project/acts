// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Json/BetheHeitlerApproxJsonConverter.hpp"

#include <cstddef>

namespace Acts {

namespace BetheHeitlerApproxJsonConverter {

using Data = std::vector<PolynomialBetheHeitlerApprox::PolyData>;

void to_json(nlohmann::json& j,
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

void to_json(nlohmann::json& j, const Data& data) {
  j = nlohmann::json::array();
  for (const auto& component : data) {
    nlohmann::json jcmp;
    to_json(jcmp, component);
    j.push_back(jcmp);
  }
}

void from_json(const nlohmann::json& j, Data& data) {
  data.clear();
  for (const auto& jcmp : j) {
    PolynomialBetheHeitlerApprox::PolyData component;
    from_json(jcmp, component);
    data.push_back(component);
  }
}

void to_json(nlohmann::json& j, const RangeData& data) {
  j["low_x0"] = data.lowX0;
  j["high_x0"] = data.highX0;
  j["weight_coeffs"] = data.weightCoeffs;
  j["mean_coeffs"] = data.meanCoeffs;
  j["var_coeffs"] = data.varCoeffs;
}

void from_json(const nlohmann::json& j, RangeData& data) {
  data.lowX0 = j.at("low_x0").get<double>();
  data.highX0 = j.at("high_x0").get<double>();

  // Handle both flat coefficients and nested components format
  if (j.contains("components")) {
    // Nested components format - flatten all components into single coefficient
    // arrays For now, just use the first component
    if (!j["components"].empty()) {
      const auto& firstComp = j["components"][0];
      data.weightCoeffs = firstComp.at("weight").get<std::vector<double>>();
      data.meanCoeffs = firstComp.at("mean").get<std::vector<double>>();
      data.varCoeffs = firstComp.at("var").get<std::vector<double>>();
    }
  } else {
    // Flat coefficients format (original)
    data.weightCoeffs = j.at("weight_coeffs").get<std::vector<double>>();
    data.meanCoeffs = j.at("mean_coeffs").get<std::vector<double>>();
    data.varCoeffs = j.at("var_coeffs").get<std::vector<double>>();
  }
}

}  // namespace BetheHeitlerApproxJsonConverter

}  // namespace Acts
