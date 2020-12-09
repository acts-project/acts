// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Digitization/DigitizationOptions.hpp"

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Digitization/Smearers.hpp"
#include "ActsExamples/Utilities/Options.hpp"

#include <numeric>
#include <string>

#include <boost/program_options.hpp>

void ActsExamples::Options::addDigitizationOptions(
    ActsExamples::Options::Description& desc) {
  using boost::program_options::bool_switch;
  using boost::program_options::value;
  using namespace std;

  auto opt = desc.add_options();
  opt("digi-config-file", value<string>(),
      "Configuration (.json) file for digitization description, overwrites "
      "options input on command line.");
  opt("digi-geometric-3d", bool_switch(),
      "Geometric: Switching geometric digitisation in 3D on");
  opt("digi-smearing", bool_switch(), "Smearing: Switching hit smearing on");
  opt("digi-smear-volume-id",
      value<read_series>()->multitoken()->default_value({}),
      "Smearing Input: sensitive volume identifiers.");
  opt("digi-smear-dimensions",
      value<read_series>()->multitoken()->default_value({}),
      "Smearing Input: smear parameters for this volume.");
  opt("digi-smear-parIndices",
      value<read_series>()->multitoken()->default_value({}),
      "Smearing Input: smear parameter indices for this volume.");
  opt("digi-smear-types", value<read_series>()->multitoken()->default_value({}),
      "Smearing Input: smear function types as 0 (gauss), 1 (truncated gauss), "
      "2 (clipped gauss), 3 (uniform), 4 (digital).");
  opt("digi-smear-parameters",
      value<read_range>()->multitoken()->default_value({}),
      "Smearing Input: smear parameters depending on the smearing type, 1 "
      "parameter for simple gauss, 4 for all others (1 parameter, 2 absolute "
      "range "
      "values, 1 range indicator: both neg < 0, neg/pos = 0, both pos > 0).");
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
  // Everything else requires a width and a range
  return (static_cast<SmearingTypes>(type) == eGauss) ? 1u : 4u;
}

ActsFatras::SingleParameterSmearFunction<ActsExamples::RandomEngine>
makeSmearFunctionForType(SmearingTypes smearingType, const double* parameters) {
  using namespace ActsExamples::Digitization;

  // this is an artifact from the command line parsing. all numbers must be
  // positive, so we need to determine the signs of the range limits with an
  // additional parameter. see command line flag documentation above for
  // details.
  double signLow = 1;
  double signHigh = 1;
  if (smearingType != eGauss) {
    if (parameters[3u] < 0) {
      signLow = -1;
      signHigh = -1;
    } else if (0 < parameters[3u]) {
      signLow = 1;
      signHigh = 1;
    } else {
      signLow = -1;
      signHigh = 1;
    }
  }

  switch (smearingType) {
    case eGauss:
      return Gauss(parameters[0]);
    case eGaussTruncated:
      return GaussTrunc(parameters[0],
                        {signLow * parameters[1u], signHigh * parameters[2u]});
    case eGaussClipped:
      return GaussClipped(
          parameters[0], {signLow * parameters[1u], signHigh * parameters[2u]});
    case eUniform:
      return Uniform(parameters[0],
                     {signLow * parameters[1u], signHigh * parameters[2u]});
    case eDigital:
      return Digital(parameters[0],
                     {signLow * parameters[1u], signHigh * parameters[2u]});
  }
  return nullptr;
}

}  // namespace

ActsExamples::SmearingAlgorithm::Config
ActsExamples::Options::readSmearingConfig(
    const boost::program_options::variables_map& variables) {
  ACTS_LOCAL_LOGGER(
      Acts::getDefaultLogger("SmearingOptions", Acts::Logging::INFO));

  using namespace Acts::UnitLiterals;

  SmearingAlgorithm::Config smearCfg;

  // Smear configuration with command line input
  // only limited smearing configuration possible
  // TODO add configuration file for more complex smearing setups

  // in case of an error, we always return a configuration struct with empty
  // smearers configuration. this will be caught later on during the algorithm
  // construction.

  auto volumes = variables["digi-smear-volume-id"].as<read_series>();
  if (volumes.empty()) {
    // no configured volumes are not considered an error at this stage
    return smearCfg;
  }
  auto dimensions = variables["digi-smear-dimensions"].as<read_series>();
  if (dimensions.size() != volumes.size()) {
    ACTS_ERROR("Inconsistent digi-smear-dimensions number of entries. Expected "
               << volumes.size() << ", but received " << dimensions.size());
    return smearCfg;
  }
  // count the expected number of parameter indices
  auto expectedNumIndices =
      std::accumulate(dimensions.begin(), dimensions.end(), size_t(0u));
  auto indices = variables["digi-smear-parIndices"].as<read_series>();
  if (indices.size() != expectedNumIndices) {
    ACTS_ERROR("Inconsistent digi-smear-parIndices number of entries. Expected "
               << expectedNumIndices << ", but received " << indices.size());
    return smearCfg;
  }
  auto types = variables["digi-smear-types"].as<read_series>();
  if (types.size() != indices.size()) {
    ACTS_ERROR("Inconsistent digi-smear-types number of entries. Expected "
               << indices.size() << ", but received " << types.size());
    return smearCfg;
  }
  // count the expected number of smearing configuration parameters
  size_t expectedNumParameters = 0;
  for (auto smearingType : types) {
    expectedNumParameters += numConfigParametersForType(smearingType);
  }
  auto parameters = variables["digi-smear-parameters"].as<read_range>();
  if (parameters.size() != expectedNumParameters) {
    ACTS_ERROR("Inconsistent digi-smear-parameters number of entries. Expected "
               << expectedNumParameters << ", but received "
               << parameters.size());
    return smearCfg;
  }

  // construct the input for the smearer configuation
  std::vector<
      std::pair<Acts::GeometryIdentifier, SmearingAlgorithm::SmearerConfig>>
      smearersInput;
  // position of the current volume in the input vectors
  size_t ivol = 0;
  // position of the current dimension index/type in the input vectors
  size_t iidx = 0;
  // position of the current smear function parameter in the input vector
  size_t ipar = 0;
  for (; ivol < volumes.size(); ++ivol) {
    Acts::GeometryIdentifier geoId =
        Acts::GeometryIdentifier(0).setVolume(volumes[ivol]);
    // number of smeared track parameters
    size_t ndim = dimensions[ivol];

    // create the smearing configuration for this geometry identifier
    SmearingAlgorithm::SmearerConfig geoCfg;
    geoCfg.reserve(ndim);

    for (auto nidx = iidx + ndim; iidx < nidx; ++iidx) {
      const auto paramIndex = static_cast<Acts::BoundIndices>(indices[iidx]);
      const auto smearingType = static_cast<SmearingTypes>(types[iidx]);
      const double* smearingParameters = &parameters[ipar];
      ipar += numConfigParametersForType(smearingType);

      SmearingAlgorithm::ParameterSmearerConfig parCfg;
      parCfg.index = paramIndex;
      parCfg.smearFunction =
          makeSmearFunctionForType(smearingType, smearingParameters);
      if (not parCfg.smearFunction) {
        ACTS_ERROR("Invalid smearing type for entry "
                   << iidx << ". Type " << types[iidx] << " is not valid");
        return smearCfg;
      }
      geoCfg.emplace_back(std::move(parCfg));
    }
    smearersInput.emplace_back(geoId, std::move(geoCfg));
  }
  // set the smearer configuration from the prepared input
  smearCfg.smearers =
      Acts::GeometryHierarchyMap<SmearingAlgorithm::SmearerConfig>(
          std::move(smearersInput));

  return smearCfg;
}
