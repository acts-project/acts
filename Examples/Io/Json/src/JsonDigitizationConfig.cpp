// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Json/JsonDigitizationConfig.hpp"

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Plugins/Json/UtilitiesJsonConverter.hpp"
#include "ActsExamples/Digitization/Smearers.hpp"

#include <fstream>
#include <functional>

void ActsExamples::to_json(nlohmann::json& j,
                           const ActsExamples::ParameterSmearingConfig& psc) {
  j["index"] = psc.index;
  // Gauss:
  const Digitization::Gauss* gauss =
      psc.smearFunction.target<const Digitization::Gauss>();
  if (gauss != nullptr) {
    j["type"] = "Gauss";
    j["mean"] = 0;
    j["stddev"] = gauss->sigma;
  }
  // Truncated gauss:
  const Digitization::GaussTrunc* gaussT =
      psc.smearFunction.target<const Digitization::GaussTrunc>();
  if (gaussT != nullptr) {
    j["type"] = "GaussTrunc";
    j["mean"] = 0;
    j["stddev"] = gaussT->sigma;
    j["range"] = gaussT->range;
  }
  // Clipped gauss:
  const Digitization::GaussClipped* gaussC =
      psc.smearFunction.target<const Digitization::GaussClipped>();
  if (gaussC != nullptr) {
    j["type"] = "GaussClipped";
    j["mean"] = 0;
    j["stddev"] = gaussC->sigma;
    j["range"] = gaussC->range;
    j["max_attempts"] = gaussC->maxAttemps;
  }
  // Uniform
  const Digitization::Uniform* uniform =
      psc.smearFunction.target<const Digitization::Uniform>();
  if (uniform != nullptr) {
    j["type"] = "Uniform";
    j["bindata"] = nlohmann::json(uniform->binningData);
  }
  // Digital
  const Digitization::Digital* digital =
      psc.smearFunction.target<const Digitization::Digital>();
  if (uniform != nullptr) {
    j["type"] = "Digitial";
    j["bindata"] = nlohmann::json(digital->binningData);
  }
}

void ActsExamples::from_json(const nlohmann::json& j,
                             ActsExamples::ParameterSmearingConfig& psc) {
  std::string sType = j["type"];

  psc.index = static_cast<Acts::BoundIndices>(j["index"]);

  if (sType == "Gauss") {
    psc.smearFunction = Digitization::Gauss(j["stddev"]);
  } else if (sType == "GaussTrunc") {
    Acts::ActsScalar sigma = j["stddev"];
    std::pair<Acts::ActsScalar, Acts::ActsScalar> range = j["range"];
    psc.smearFunction = Digitization::GaussTrunc(sigma, range);
  } else if (sType == "GaussClipped") {
    Acts::ActsScalar sigma = j["stddev"];
    std::pair<Acts::ActsScalar, Acts::ActsScalar> range = j["range"];
    psc.smearFunction = Digitization::GaussClipped(sigma, range);
  } else if (sType == "Uniform") {
    Acts::BinningData bd;
    from_json(j["bindata"], bd);
    psc.smearFunction = Digitization::Uniform(std::move(bd));
  } else if (sType == "Digitial") {
    Acts::BinningData bd;
    from_json(j["bindata"], bd);
    psc.smearFunction = Digitization::Uniform(std::move(bd));
  }
}

void ActsExamples::to_json(nlohmann::json& j,
                           const ActsExamples::GeometricConfig& gdc) {
  std::vector<size_t> indices;
  for (const auto& idx : gdc.indices) {
    indices.push_back(static_cast<size_t>(idx));
  }
  j["indices"] = indices;
  j["segmentation"] = nlohmann::json(gdc.segmentation);
  j["thickness"] = gdc.thickness;
  j["threshold"] = gdc.threshold;
  j["digital"] = gdc.digital;
}

void ActsExamples::from_json(const nlohmann::json& j,
                             ActsExamples::GeometricConfig& gdc) {
  for (const auto& jidx : j["indices"]) {
    gdc.indices.push_back(static_cast<Acts::BoundIndices>(jidx));
  }
  from_json(j["segmentation"], gdc.segmentation);
  gdc.thickness = j["thickness"];
  gdc.threshold = j["threshold"];
  gdc.digital = j["digital"];
}

void ActsExamples::to_json(nlohmann::json& j,
                           const ActsExamples::SmearingConfig& sdc) {
  for (const auto& sc : sdc) {
    j.push_back(nlohmann::json(sc));
  }
}

void ActsExamples::from_json(const nlohmann::json& j,
                             ActsExamples::SmearingConfig& sdc) {
  for (const auto& jpsc : j) {
    ActsExamples::ParameterSmearingConfig psc;
    from_json(jpsc, psc);
    sdc.push_back(psc);
  }
}

void ActsExamples::to_json(nlohmann::json& j,
                           const ActsExamples::DigiComponentsConfig& dc) {
  if (not dc.geometricDigiConfig.indices.empty()) {
    j["geometric"] = nlohmann::json(dc.geometricDigiConfig);
  }
  if (not dc.smearingDigiConfig.empty()) {
    j["smearing"] = nlohmann::json(dc.smearingDigiConfig);
  }
}

void ActsExamples::from_json(const nlohmann::json& j,
                             ActsExamples::DigiComponentsConfig& dc) {
  if (j.find("geometric") != j.end()) {
    nlohmann::json jgdc = j["geometric"];
    from_json(jgdc, dc.geometricDigiConfig);
  }
  if (j.find("smearing") != j.end()) {
    nlohmann::json jsdc = j["smearing"];
    from_json(jsdc, dc.smearingDigiConfig);
  }
}

Acts::GeometryHierarchyMap<ActsExamples::DigiComponentsConfig>
ActsExamples::readDigiConfigFromJson(const std::string& path) {
  nlohmann::json djson;
  if (path.empty()) {
    return Acts::GeometryHierarchyMap<ActsExamples::DigiComponentsConfig>();
  }
  std::ifstream infile(path, std::ifstream::in | std::ifstream::binary);
  // rely on exception for error handling
  infile.exceptions(std::ofstream::failbit | std::ofstream::badbit);
  infile >> djson;
  return DigiConfigConverter("digitization-configuration").fromJson(djson);
}

void ActsExamples::writeDigiConfigToJson(
    const Acts::GeometryHierarchyMap<DigiComponentsConfig>& cfg,
    const std::string& path) {
  std::ofstream outfile(path, std::ofstream::out | std::ofstream::binary);
  // rely on exception for error handling
  outfile.exceptions(std::ofstream::failbit | std::ofstream::badbit);
  outfile << DigiConfigConverter("digitization-configuration")
                 .toJson(cfg, nullptr)
                 .dump(2);
}
