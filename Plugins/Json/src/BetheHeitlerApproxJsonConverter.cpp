// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Json/BetheHeitlerApproxJsonConverter.hpp"

namespace Acts {

namespace BetheHeitlerApproxJsonConverter {

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

void to_json(nlohmann::json& j,
             const PolynomialBetheHeitlerApprox::RangeData& data) {
  j["low_x0"] = data.lowX0;
  j["high_x0"] = data.highX0;
  j["transform"] = data.transform;
  to_json(j, data.data);
}

void from_json(const nlohmann::json& j,
               PolynomialBetheHeitlerApprox::RangeData& data) {
  data.lowX0 = j.at("low_x0").get<double>();
  data.highX0 = j.at("high_x0").get<double>();
  data.transform = j.value("transform", true);

  if (j.contains("components")) {
    // Components array format
    data.data.clear();
    for (const auto& jcmp : j["components"]) {
      PolynomialBetheHeitlerApprox::PolyData component;
      from_json(jcmp, component);
      data.data.push_back(component);
    }
  } else if (j.contains("weight_coeffs")) {
    // Flat coefficients format (single component) - convert to Data
    PolynomialBetheHeitlerApprox::PolyData component;
    from_json(j, component);
    data.data = {component};
  } else {
    throw std::runtime_error(
        "JSON range data must contain either 'components' array or "
        "'weight_coeffs' (flat format)");
  }
}

}  // namespace BetheHeitlerApproxJsonConverter

}  // namespace Acts
