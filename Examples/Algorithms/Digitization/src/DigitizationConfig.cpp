// This file is part of the Acts project.
//
// Copyright (C) 2020-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Digitization/DigitizationConfig.hpp"

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Digitization/DigitizationAlgorithm.hpp"
#include "ActsExamples/Digitization/Smearers.hpp"
#include "ActsExamples/Digitization/SmearingAlgorithm.hpp"
#include "ActsExamples/Digitization/SmearingConfig.hpp"
#include "ActsExamples/Utilities/Options.hpp"

#include <numeric>
#include <string>

#include <boost/program_options.hpp>

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

ActsExamples::DigitizationConfig::DigitizationConfig(
    const Options::Variables& vars,
    Acts::GeometryHierarchyMap<DigiComponentsConfig>&& digiCfgs)
    : isSimpleSmearer(vars["digi-smear"].as<bool>()),
      doMerge(vars["digi-merge"].as<bool>()),
      mergeNsigma(vars["digi-merge-nsigma"].as<double>()),
      mergeCommonCorner(vars["digi-merge-common-corner"].as<bool>()) {
  digitizationConfigs = std::move(digiCfgs);
  if (isSimpleSmearer)
    smearingConfig(vars);
}

ActsExamples::DigitizationConfig::DigitizationConfig(
    Acts::GeometryHierarchyMap<DigiComponentsConfig>&& digiCfgs)
    : isSimpleSmearer(false),
      doMerge(false),
      mergeNsigma(1.0),
      mergeCommonCorner(false) {
  digitizationConfigs = std::move(digiCfgs);
}

void ActsExamples::DigitizationConfig::smearingConfig(
    const Options::Variables& variables) {
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

  // no configured volumes are not considered an error at this stage
  if (not variables.count("digi-smear-volume"))
    return;
  auto volumes = variables["digi-smear-volume"].as<std::vector<int>>();
  if (volumes.empty())
    return;

  if (not variables["digi-config-file"].as<std::string>().empty()) {
    ACTS_WARNING(
        "Smearing configuration on command-line will override .json "
        "configuration!");
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
  std::vector<std::pair<Acts::GeometryIdentifier, DigiComponentsConfig>>
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
    DigiComponentsConfig dcfg;
    dcfg.smearingDigiConfig = geoCfg;
    smearersInput.emplace_back(geoId, std::move(dcfg));
  }
  // set the smearer configuration from the prepared input
  digitizationConfigs = Acts::GeometryHierarchyMap<DigiComponentsConfig>(
      std::move(smearersInput));
}

std::shared_ptr<ActsExamples::IAlgorithm>
ActsExamples::createDigitizationAlgorithm(ActsExamples::DigitizationConfig& cfg,
                                          Acts::Logging::Level lvl) {
  if (cfg.isSimpleSmearer)
    return std::make_shared<SmearingAlgorithm>(cfg, lvl);
  else
    return std::make_shared<DigitizationAlgorithm>(cfg, lvl);
}

std::vector<
    std::pair<Acts::GeometryIdentifier, std::vector<Acts::BoundIndices>>>
ActsExamples::DigitizationConfig::getBoundIndices() const {
  std::vector<
      std::pair<Acts::GeometryIdentifier, std::vector<Acts::BoundIndices>>>
      bIndexInput;

  for (size_t ibi = 0; ibi < digitizationConfigs.size(); ++ibi) {
    Acts::GeometryIdentifier geoID = digitizationConfigs.idAt(ibi);
    const auto dCfg = digitizationConfigs.valueAt(ibi);
    std::vector<Acts::BoundIndices> boundIndices;
    if (isSimpleSmearer) {
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
