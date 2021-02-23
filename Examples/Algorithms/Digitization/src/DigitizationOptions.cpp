// This file is part of the Acts project.
//
// Copyright (C) 2020-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Digitization/DigitizationOptions.hpp"

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Digitization/JsonDigitizationConfig.hpp"
#include "ActsExamples/Digitization/Smearers.hpp"
#include "ActsExamples/Digitization/SmearingConfig.hpp"
#include "ActsExamples/Utilities/Options.hpp"

#include <fstream>
#include <numeric>
#include <string>

#include <boost/program_options.hpp>
#include <nlohmann/json.hpp>

void ActsExamples::Options::addDigitizationOptions(Description& desc) {
  using boost::program_options::bool_switch;
  using boost::program_options::value;

  // each volume configuration is one logical block
  //
  //   --digi-smear-volume-id=8
  //   --digi-smear-indices=0:1:5 # loc0, loc1, and time
  //   --digi-smear-types=0:0:3   # loc{0,1} uses gaussian, time uses uniform
  //   # parameter 0: loc0 gaussian width
  //   # parameter 1: loc1 gaussian width
  //   # parameter 2-4: time pitch,min,max
  //   --digi-smear-parameters=10:20:2.5:-25:25
  //
  // which can be repeated as often as needed
  //
  //   --digi-smear-volume-id=11
  //   --digi-smear-indices=1       # loc1
  //   --digi-smear-types=0         # loc1 uses gaussian
  //   --digi-smear-parameters=12.5 # loc1 gaussian width
  //
  auto opt = desc.add_options();
  opt("digi-config-file", value<std::string>()->default_value(""),
      "Configuration (.json) file for digitization description, overwrites "
      "smearing options input on command line.");
  opt("digi-geometric-3d", bool_switch(),
      "Geometric: Switching geometric digitisation in 3D on");
  opt("digi-smear", bool_switch(), "Smearing: Switching hit smearing on");
  opt("digi-smear-volume", value<std::vector<int>>(),
      "Smearing Input: sensitive volume identifiers.");
  opt("digi-smear-indices", value<std::vector<VariableIntegers>>(),
      "Smearing Input: smear parameter indices for this volume.");
  opt("digi-smear-types", value<std::vector<VariableIntegers>>(),
      "Smearing Input: smear function types as 0 (gauss), 1 (truncated gauss), "
      "2 (clipped gauss), 3 (uniform), 4 (digital).");
  opt("digi-smear-parameters", value<std::vector<VariableReals>>(),
      "Smearing Input: smear parameters depending on the smearing type, 1 "
      "parameter for simple gauss, 3 for all others (1 parameter, 2 range "
      "values.");
}

namespace {

enum SmearingTypes : int {
  eGauss = 0,
  eGaussTruncated = 1,
  eGaussClipped = 2,
  eUniform = 3,
  eDigital = 4,
};

constexpr size_t numConfigParametersForType(int type) {
  // Gaussian smearing requires only a width/standard deviation parameter
  // Everything else requires a width/pitch and a range
  return (static_cast<SmearingTypes>(type) == eGauss) ? 1u : 3u;
}

ActsFatras::SingleParameterSmearFunction<ActsExamples::RandomEngine>
makeSmearFunctionForType(SmearingTypes smearingType, const double* parameters) {
  using namespace ActsExamples::Digitization;

  switch (smearingType) {
    case eGauss:
      return Gauss(parameters[0]);
    case eGaussTruncated:
      return GaussTrunc(parameters[0], {parameters[1u], parameters[2u]});
    case eGaussClipped:
      return GaussClipped(parameters[0], {parameters[1u], parameters[2u]});
    case eUniform:
      return Uniform(parameters[0], {parameters[1u], parameters[2u]});
    case eDigital:
      return Digital(parameters[0], {parameters[1u], parameters[2u]});
  }
  return nullptr;
}

}  // namespace

ActsExamples::Digitization::AlgorithmConfig::AlgorithmConfig(const Options::Variables &vars, Acts::GeometryHierarchyMap<DigitizationConfig> &&digiCfgs)
  : isSimpleSmearer(vars["digi-smear"].as<bool>())
{
  digitizationConfigs = std::move(digiCfgs);
  if (isSimpleSmearer)
    smearingConfig(vars);
}

void ActsExamples::Digitization::AlgorithmConfig::smearingConfig(const Options::Variables &variables) {
  ACTS_LOCAL_LOGGER(
      Acts::getDefaultLogger("SmearingOptions", Acts::Logging::INFO));

  using namespace Acts::UnitLiterals;
  using namespace ActsExamples::Options;

  // Smear configuration with command line input
  // only limited smearing configuration possible:
  // complex smearing option configuration has to be done with a
  // digitization config file

  // in case of an error, we always return a configuration struct with empty
  // smearers configuration. this will be caught later on during the algorithm
  // construction.

  auto volumes = variables["digi-smear-volume"].as<std::vector<int>>();
  if (volumes.empty()) {
    // no configured volumes are not considered an error at this stage
    return;
  }

  if (not variables["digi-config-file"].as<std::string>().empty()) {
    ACTS_WARNING("Smearing configuration on command-line will override .json configuration!");
  }

  auto indices =
      variables["digi-smear-indices"].as<std::vector<VariableIntegers>>();
  auto types =
      variables["digi-smear-types"].as<std::vector<VariableIntegers>>();
  auto parameters =
      variables["digi-smear-parameters"].as<std::vector<VariableReals>>();
  if (indices.size() != volumes.size()) {
    ACTS_ERROR("Inconsistent digi-smear-indices options. Expected "
               << volumes.size() << ", but received " << indices.size());
    return;
  }
  if (types.size() != volumes.size()) {
    ACTS_ERROR("Inconsistent digi-smear-types options. Expected "
               << volumes.size() << ", but received " << types.size());
    return;
  }
  if (parameters.size() != volumes.size()) {
    ACTS_ERROR("Inconsistent digi-smear-parameters options. Expected "
               << volumes.size() << ", but received " << parameters.size());
    return;
  }

  // construct the input for the smearer configuation
  std::vector<std::pair<Acts::GeometryIdentifier, DigitizationConfig>>
      smearersInput;
  for (size_t ivol = 0; ivol < volumes.size(); ++ivol) {
    Acts::GeometryIdentifier geoId =
        Acts::GeometryIdentifier(0).setVolume(volumes[ivol]);
    const auto& volIndices = indices[ivol].values;
    const auto& volTypes = types[ivol].values;
    const auto& volParameters = parameters[ivol].values;

    if (volTypes.size() != volIndices.size()) {
      ACTS_ERROR("Inconsistent number of digi-smear-types values for volume "
                 << volumes[ivol] << ". Expected " << indices.size()
                 << ", but received " << types.size());
      return;
    }
    // count the expected number of smearing configuration parameters
    size_t expectedNumParameters = 0;
    for (auto smearingType : volTypes) {
      expectedNumParameters += numConfigParametersForType(smearingType);
    }
    if (volParameters.size() != expectedNumParameters) {
      ACTS_ERROR(
          "Inconsistent number of digi-smear-parameters values for volume "
          << volumes[ivol] << ". Expected " << expectedNumParameters
          << ", but received " << types.size());
      return;
    }

    // create the smearing configuration for this geometry identifier
    SmearingConfig geoCfg;
    geoCfg.reserve(volIndices.size());

    for (size_t iidx = 0, ipar = 0; iidx < volIndices.size(); ++iidx) {
      const auto paramIndex = static_cast<Acts::BoundIndices>(volIndices[iidx]);
      const auto smearingType = static_cast<SmearingTypes>(volTypes[iidx]);
      const double* smearingParameters = &volParameters[ipar];
      ipar += numConfigParametersForType(smearingType);

      ParameterSmearingConfig parCfg;
      parCfg.index = paramIndex;
      parCfg.smearFunction =
          makeSmearFunctionForType(smearingType, smearingParameters);
      if (not parCfg.smearFunction) {
        ACTS_ERROR("Invalid smearing type for entry "
                   << iidx << " for volume " << volumes[ivol] << ". Type "
                   << volTypes[iidx] << " is not valid");
        return;
      }
      geoCfg.emplace_back(std::move(parCfg));
    }
    DigitizationConfig dcfg;
    dcfg.smearingDigiConfig = geoCfg;
    smearersInput.emplace_back(geoId, std::move(dcfg));
  }
  // set the smearer configuration from the prepared input
  digitizationConfigs =
      Acts::GeometryHierarchyMap<DigitizationConfig>(std::move(smearersInput));
}

Acts::GeometryHierarchyMap<ActsExamples::DigitizationConfig>
ActsExamples::Options::readConfigFromJson(const std::string& path) {
  nlohmann::json djson;
  auto in = std::ifstream(path, std::ifstream::in | std::ifstream::binary);
  if (in.good()) {
    // Get the json file for the configuration
    in >> djson;
    in.close();
  }  // TODO handle errors
  return DigiConfigConverter("digitization-configuration").fromJson(djson);
}

std::shared_ptr<ActsExamples::IAlgorithm>
ActsExamples::Options::createDigitizationAlgorithm(
    ActsExamples::Digitization::AlgorithmConfig& cfg,
    Acts::Logging::Level lvl) {
  if (cfg.isSimpleSmearer)
    return std::make_shared<SmearingAlgorithm>(cfg, lvl);
  else
    return std::make_shared<DigitizationAlgorithm>(cfg, lvl);
}

std::vector<
    std::pair<Acts::GeometryIdentifier, std::vector<Acts::BoundIndices>>>
ActsExamples::Options::getBoundIndices(
    ActsExamples::Digitization::AlgorithmConfig& cfg) {
  std::vector<
      std::pair<Acts::GeometryIdentifier, std::vector<Acts::BoundIndices>>>
      bIndexInput;

  for (size_t ibi = 0; ibi < cfg.digitizationConfigs.size(); ++ibi) {
    Acts::GeometryIdentifier geoID = cfg.digitizationConfigs.idAt(ibi);
    const auto dCfg = cfg.digitizationConfigs.valueAt(ibi);
    std::vector<Acts::BoundIndices> boundIndices;
    if (cfg.isSimpleSmearer) {
      for (const auto& sConfig : dCfg.smearingDigiConfig) {
        boundIndices.push_back(sConfig.index);
      }
    } else {
      boundIndices.insert(boundIndices.end(),
                          dCfg.geometricDigiConfig.indices.begin(),
                          dCfg.geometricDigiConfig.indices.end());
    }
    bIndexInput.push_back({geoID, boundIndices});
  }
  return bIndexInput;
}
