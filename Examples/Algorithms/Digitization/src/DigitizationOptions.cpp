// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Digitization/DigitizationOptions.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"
#include "ActsExamples/Digitization/HitSmearers.hpp"
#include "ActsExamples/Utilities/Options.hpp"
#include <boost/program_options.hpp>

#include <string>

void ActsExamples::Options::addDigitizationOptions(
    ActsExamples::Options::Description& desc) {
  using boost::program_options::value;
  using namespace std;

  auto opt = desc.add_options();
  opt("digi-input-hits", value<string>()->default_value("hits"),
      "Name of the input hit collection.");
  opt("digi-config-file", value<string>()->default_value(""),
      "Configuration (.json) file for digitization description, overwrites "
      "options input on command line.");
  opt("digi-geometric-3d", value<bool>()->default_value(true),
      "Geometric: Switching geometric digitisation in 3D on");
  opt("digi-smearing", value<bool>()->default_value(false),
      "Smearing: Switching geometric digitisation in 3D on");
  opt("digi-smear-output", value<string>()->default_value(""),
      "Smearing Output: Name of the output measurement collection.");
  opt("digi-smear-volume-id",
      value<read_series>()->multitoken()->default_value({}),
      "Smearing Input: sensitive volume identifiers.");
  opt("digi-smear-dimensions",
      value<read_series>()->multitoken()->default_value({}),
      "Smearing Input: smear parameters for this volume.");
  opt("digi-smear-binvalues",
      value<read_series>()->multitoken()->default_value({}),
      "Smearing Input: smear binning values for this volume.");
  opt("digi-smear-types", value<read_series>()->multitoken()->default_value({}),
      "Smearing Input: smear function types as 0 (gauss), 1 (truncated gauss), "
      "2 "
      "(clipped gauss), 3 (uniform).");
  opt("digi-smear-parameters",
      value<read_range>()->multitoken()->default_value({}),
      "Smearing Input: smear parameters depending on the smearing type.");
}

ActsExamples::SmearingAlgorithm::Config
ActsExamples::Options::readSmearingConfig(
    const boost::program_options::variables_map& variables) {
  ACTS_LOCAL_LOGGER(
      Acts::getDefaultLogger("SmearingOptions", Acts::Logging::INFO));

  using namespace Acts::UnitLiterals;

  SmearingAlgorithm::Config smearCfg;

  smearCfg.inputSimulatedHits = variables["digi-input-hits"].as<std::string>();
  smearCfg.outputMeasurements =
      variables["digi-smear-output"].as<std::string>();

  // Smear configuration with command line input
  // only limted smearing possible, use configuration for more
  // complex smearing setups
  auto volumes = variables["digi-smear-volume-id"].as<read_series>();
  if (not volumes.empty()) {
    HitSmearers::FunctionGenerator sFnc;

    auto vdims = variables["digi-smear-dimensions"].as<read_series>();
    auto bvalues = variables["digi-smear-binvalues"].as<read_series>();
    auto types = variables["digi-smear-types"].as<read_series>();
    auto parameters = variables["digi-smear-parameters"].as<read_range>();

    // Count the number of needed parameters
    // - type 0 needs 1 parameters
    // - type 1-3 needs 3 parameters
    size_t sumpars = 0;
    std::for_each(types.begin(), types.end(),
                  [&](int n) { sumpars += (n == 0) ? 1 : 3; });

    if (volumes.size() == vdims.size() and bvalues.size() == types.size() and
        parameters.size() == sumpars) {
      size_t vpos = 0;
      size_t ppos = 0;
      size_t padd = 0;

      ACTS_DEBUG("Volume parameters properly read from command line.")

      /// Extract the parameters from the global vector
      ///
      /// @param ftype the function type
      ///
      /// @return an extracted vector
      auto extract = [&](int ftype) -> std::vector<double> {
        padd = (ftype == 0) ? 1 : 3;
        std::vector<double> fpars = {parameters.begin() + ppos,
                                     parameters.begin() + ppos + padd};
        ppos += padd;
        return fpars;
      };

      for (unsigned int iv = 0; iv < volumes.size(); ++iv) {
        SmearingAlgorithm::SupportedSmearer smearer;

        int volID = volumes[iv];
        Acts::GeometryIdentifier volumeGeometryId =
            Acts::GeometryIdentifier(0).setVolume(volID);

        int volDim = vdims[iv];
        int ftype = 0;
        // 1 - dimensional, either loc 0 or loc 1
        if (volDim == 1) {
          using OneFnc = std::array<ActsFatras::SmearFunction<RandomEngine>, 1>;
          std::array<int, 1> ftypes;
          std::array<std::vector<double>, 1> fparameters;
          ftype = types[vpos++];
          ftypes[0] = ftype;
          fparameters[0] = extract(ftype);
          OneFnc smearFunctions = sFnc.generate<1>(ftypes, fparameters);
          // Create and fill
          if ((Acts::BoundIndices)bvalues[vpos - 1] == Acts::eBoundLoc0) {
            smearCfg.configured = true;
            smearer =
                std::pair<ActsFatras::BoundParametersSmearer<Acts::eBoundLoc0>,
                          OneFnc>(
                    ActsFatras::BoundParametersSmearer<Acts::eBoundLoc0>(),
                    std::move(smearFunctions));
          } else if ((Acts::BoundIndices)bvalues[vpos - 1] ==
                     Acts::eBoundLoc0) {
            smearCfg.configured = true;
            smearer =
                std::pair<ActsFatras::BoundParametersSmearer<Acts::eBoundLoc1>,
                          OneFnc>(
                    ActsFatras::BoundParametersSmearer<Acts::eBoundLoc1>(),
                    std::move(smearFunctions));
          }
        } else if (volDim == 2) {
          // 2 - dimensional,  loc0 x loc1, loc0 x time, loc1 x time
          using TwoFnc = std::array<ActsFatras::SmearFunction<RandomEngine>, 2>;
          std::array<int, 2> ftypes;
          std::array<std::vector<double>, 2> fparameters;
          for (size_t idim = 0; idim < (size_t)volDim; ++idim) {
            ftype = types[vpos++];
            ftypes[idim] = ftype;
            fparameters[idim] = extract(ftype);
          }
          TwoFnc smearFunctions = sFnc.generate<2>(ftypes, fparameters);
          if ((Acts::BoundIndices)bvalues[vpos - 2] == Acts::eBoundLoc0 and
              (Acts::BoundIndices) bvalues[vpos - 1] == Acts::eBoundTime) {
            smearCfg.configured = true;
            smearer =
                std::pair<ActsFatras::BoundParametersSmearer<Acts::eBoundLoc0,
                                                             Acts::eBoundTime>,
                          TwoFnc>(
                    ActsFatras::BoundParametersSmearer<Acts::eBoundLoc0,
                                                       Acts::eBoundTime>(),
                    std::move(smearFunctions));
          } else if ((Acts::BoundIndices)bvalues[vpos - 2] ==
                         Acts::eBoundLoc1 and
                     (Acts::BoundIndices)
                             bvalues[vpos - 1] == Acts::eBoundTime) {
            smearCfg.configured = true;
            smearer =
                std::pair<ActsFatras::BoundParametersSmearer<Acts::eBoundLoc1,
                                                             Acts::eBoundTime>,
                          TwoFnc>(
                    ActsFatras::BoundParametersSmearer<Acts::eBoundLoc1,
                                                       Acts::eBoundTime>(),
                    std::move(smearFunctions));
          } else if ((Acts::BoundIndices)bvalues[vpos - 2] ==
                         Acts::eBoundLoc0 and
                     (Acts::BoundIndices)
                             bvalues[vpos - 1] == Acts::eBoundLoc1) {
            smearCfg.configured = true;
            smearer =
                std::pair<ActsFatras::BoundParametersSmearer<Acts::eBoundLoc0,
                                                             Acts::eBoundLoc1>,
                          TwoFnc>(
                    ActsFatras::BoundParametersSmearer<Acts::eBoundLoc0,
                                                       Acts::eBoundLoc1>(),
                    std::move(smearFunctions));
          }

        } else if (volDim == 3) {
          // 3 - dimensional,  loc0 x loc1 x time
          using ThreeFnc =
              std::array<ActsFatras::SmearFunction<RandomEngine>, 3>;
          std::array<int, 3> ftypes;
          std::array<std::vector<double>, 3> fparameters;
          for (size_t idim = 0; idim < (size_t)volDim; ++idim) {
            ftype = types[vpos++];
            ftypes[idim] = ftype;
            fparameters[idim] = extract(ftype);
          }
          ThreeFnc smearFunctions = sFnc.generate<3>(ftypes, fparameters);
          if ((Acts::BoundIndices)bvalues[vpos - 3] == Acts::eBoundLoc0 and
              (Acts::BoundIndices) bvalues[vpos - 2] == Acts::eBoundLoc1 and
              (Acts::BoundIndices) bvalues[vpos - 1] == Acts::eBoundTime) {
            smearCfg.configured = true;
            smearer = std::pair<
                ActsFatras::BoundParametersSmearer<
                    Acts::eBoundLoc0, Acts::eBoundLoc1, Acts::eBoundTime>,
                ThreeFnc>(
                ActsFatras::BoundParametersSmearer<
                    Acts::eBoundLoc0, Acts::eBoundLoc1, Acts::eBoundTime>(),
                std::move(smearFunctions));
          }
        }
        // fill the smearer into the configuration map
        smearCfg.smearers.emplace_hint(smearCfg.smearers.end(),
                                       volumeGeometryId, std::move(smearer));
      }
    }
  }
  return smearCfg;
}
