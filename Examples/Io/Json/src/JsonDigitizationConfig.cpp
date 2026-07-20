// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Json/JsonDigitizationConfig.hpp"

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/BinningData.hpp"
#include "Acts/Utilities/IAxis.hpp"
#include "Acts/Utilities/IMultiAxis.hpp"
#include "ActsExamples/Digitization/Smearers.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "ActsFatras/Digitization/UncorrelatedHitSmearer.hpp"
#include "ActsPlugins/Json/UtilitiesJsonConverter.hpp"

#include <cstddef>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

namespace ActsExamples {
namespace {
void to_json(nlohmann::json& j,
             const ActsFatras::SingleParameterSmearFunction<RandomEngine>& f) {
  // Gauss:
  auto gauss = f.target<const Digitization::Gauss>();
  if (gauss != nullptr) {
    j["type"] = "Gauss";
    j["mean"] = 0;
    j["stddev"] = gauss->sigma;
    return;
  }
  // Truncated gauss:
  auto gaussT = f.target<const Digitization::GaussTrunc>();
  if (gaussT != nullptr) {
    j["type"] = "GaussTrunc";
    j["mean"] = 0;
    j["stddev"] = gaussT->sigma;
    j["range"] = gaussT->range;
    return;
  }
  // Clipped gauss:
  auto gaussC = f.target<const Digitization::GaussClipped>();
  if (gaussC != nullptr) {
    j["type"] = "GaussClipped";
    j["mean"] = 0;
    j["stddev"] = gaussC->sigma;
    j["range"] = gaussC->range;
    j["max_attempts"] = gaussC->maxAttemps;
    return;
  }
  // Uniform
  auto uniform = f.target<const Digitization::Uniform>();
  if (uniform != nullptr) {
    j["type"] = "Uniform";
    j["bindata"] = nlohmann::json(uniform->binningData);
    return;
  }
  // Digital
  auto digital = f.target<const Digitization::Digital>();
  if (digital != nullptr) {
    j["type"] = "Digital";
    j["bindata"] = nlohmann::json(digital->binningData);
    return;
  }
  // Exact
  auto exact = f.target<const Digitization::Exact>();
  if (exact != nullptr) {
    j["type"] = "Exact";
    j["stddev"] = exact->sigma;
    return;
  }

  throw std::runtime_error("Unable to serialize smearer");
}

void from_json(const nlohmann::json& j,
               ActsFatras::SingleParameterSmearFunction<RandomEngine>& f) {
  std::string sType = j["type"];

  if (sType == "Gauss") {
    f = Digitization::Gauss(j["stddev"]);
  } else if (sType == "GaussTrunc") {
    double sigma = j["stddev"];
    std::pair<double, double> range = j["range"];
    f = Digitization::GaussTrunc(sigma, range);
  } else if (sType == "GaussClipped") {
    double sigma = j["stddev"];
    std::pair<double, double> range = j["range"];
    f = Digitization::GaussClipped(sigma, range);
  } else if (sType == "Uniform") {
    Acts::BinningData bd;
    from_json(j["bindata"], bd);
    f = Digitization::Uniform(bd);
  } else if (sType == "Digital") {
    Acts::BinningData bd;
    from_json(j["bindata"], bd);
    f = Digitization::Digital(bd);
  } else if (sType == "Exact") {
    f = Digitization::Exact(j["stddev"]);
  } else {
    throw std::invalid_argument("Unknown smearer type '" + sType + "'");
  }
}

}  // namespace
}  // namespace ActsExamples

void ActsExamples::to_json(nlohmann::json& j,
                           const ParameterSmearingConfig& psc) {
  j["index"] = psc.index;
  j["forcePositiveValues"] = psc.forcePositiveValues;
  to_json(j, psc.smearFunction);
}

void ActsExamples::from_json(const nlohmann::json& j,
                             ParameterSmearingConfig& psc) {
  psc.index = static_cast<Acts::BoundIndices>(j["index"]);
  if (j.find("forcePositiveValues") != j.end()) {
    psc.forcePositiveValues = j["forcePositiveValues"];
  }
  from_json(j, psc.smearFunction);
}

void ActsExamples::to_json(nlohmann::json& j, const GeometricConfig& gdc) {
  std::vector<std::size_t> indices;
  for (const auto& idx : gdc.indices) {
    indices.push_back(static_cast<std::size_t>(idx));
  }
  j["indices"] = indices;
  Acts::BinUtility segmentation;
  if (gdc.segmentation != nullptr) {
    for (const Acts::IAxis& axis : *gdc.segmentation) {
      const Acts::AxisDirection axisDir = axis.getDirection().value();
      const Acts::BinningOption bOption =
          axis.getBoundaryType() == Acts::AxisBoundaryType::Closed
              ? Acts::closed
              : Acts::open;
      if (axis.isEquidistant()) {
        segmentation += Acts::BinUtility(axis.getNBins(), axis.getMin(),
                                         axis.getMax(), bOption, axisDir);
      } else {
        const std::vector<double> edges = axis.getBinEdges();
        std::vector<float> floatEdges(edges.begin(), edges.end());
        segmentation += Acts::BinUtility(floatEdges, bOption, axisDir);
      }
    }
  }
  j["segmentation"] = segmentation;
  j["thickness"] = gdc.thickness;
  j["threshold"] = gdc.threshold;
  j["digital"] = gdc.digital;
  if (j.find("charge-smearing") != j.end()) {
    to_json(j["charge-smearing"], gdc.chargeSmearer);
  }
}

void ActsExamples::from_json(const nlohmann::json& j, GeometricConfig& gdc) {
  for (const auto& jidx : j["indices"]) {
    gdc.indices.push_back(static_cast<Acts::BoundIndices>(jidx));
  }
  Acts::BinUtility segmentation;
  from_json(j["segmentation"], segmentation);
  std::vector<std::unique_ptr<Acts::IAxis>> axes;
  for (const Acts::BinningData& binData : segmentation.binningData()) {
    const Acts::AxisDirection axisDir = binData.binvalue;
    const Acts::AxisType axisType = binData.type == Acts::equidistant
                                        ? Acts::AxisType::Equidistant
                                        : Acts::AxisType::Variable;
    const Acts::AxisBoundaryType abType = binData.option == Acts::open
                                              ? Acts::AxisBoundaryType::Open
                                              : Acts::AxisBoundaryType::Closed;
    if (axisType == Acts::AxisType::Equidistant) {
      axes.push_back(Acts::IAxis::createEquidistant(
          abType, binData.min, binData.max, binData.bins(), axisDir));
    } else {
      const std::vector<double> edges(binData.boundaries().begin(),
                                      binData.boundaries().end());
      axes.push_back(Acts::IAxis::createVariable(abType, edges, axisDir));
    }
  }
  if (axes.size() == 1) {
    gdc.segmentation = Acts::IMultiAxis::create(*axes[0]);
  } else if (axes.size() == 2) {
    gdc.segmentation = Acts::IMultiAxis::create(*axes[0], *axes[1]);
  } else if (!axes.empty()) {
    throw std::invalid_argument(
        "GeometricConfig: only 1D and 2D segmentations are supported");
  }
  gdc.thickness = j["thickness"];
  gdc.threshold = j["threshold"];
  gdc.digital = j["digital"];
  if (j.find("variances") != j.end()) {
    /// Read the variances from the json file
    auto jvariances = j["variances"];
    for (const auto& jvar : jvariances) {
      auto idx =
          static_cast<Acts::BoundIndices>(jvar["index"].get<std::size_t>());
      auto rms = jvar["rms"].get<std::vector<double>>();
      auto vars = rms;
      // Square the RMS values to get the variances
      std::transform(vars.begin(), vars.end(), vars.begin(),
                     [](double val) { return val * val; });
      gdc.varianceMap[idx] = vars;
    }
  }
  if (j.find("charge-smearing") != j.end()) {
    from_json(j["charge-smearing"], gdc.chargeSmearer);
  }
}

void ActsExamples::to_json(nlohmann::json& j, const SmearingConfig& sdc) {
  for (const auto& sc : sdc.params) {
    j.push_back(nlohmann::json(sc));
  }
}

void ActsExamples::from_json(const nlohmann::json& j, SmearingConfig& sdc) {
  for (const auto& jpsc : j) {
    ParameterSmearingConfig psc;
    from_json(jpsc, psc);
    sdc.params.push_back(psc);
  }
}

void ActsExamples::to_json(nlohmann::json& j, const DigiComponentsConfig& dc) {
  if (!dc.geometricDigiConfig.indices.empty()) {
    j["geometric"] = nlohmann::json(dc.geometricDigiConfig);
  }
  if (!dc.smearingDigiConfig.params.empty()) {
    j["smearing"] = nlohmann::json(dc.smearingDigiConfig);
  }
}

void ActsExamples::from_json(const nlohmann::json& j,
                             DigiComponentsConfig& dc) {
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
    return Acts::GeometryHierarchyMap<DigiComponentsConfig>();
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
