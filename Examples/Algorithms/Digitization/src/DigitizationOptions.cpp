// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Digitization/DigitizationOptions.hpp"

#include "ActsExamples/Digitization/Smearers.hpp"
#include "ActsExamples/Utilities/Options.hpp"
#include <Acts/Utilities/Logger.hpp>
#include <Acts/Utilities/ParameterDefinitions.hpp>

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
  opt("digi-smearing", bool_switch(),
      "Smearing: Switching geometric digitisation in 3D on");
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
      "values, 1 range indicator: both neg < 0, pos/neg = 0, both pos > 0).");
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
    Digitization::SmearingFunctionGenerator sFnc;

    auto vdims = variables["digi-smear-dimensions"].as<read_series>();
    auto pindices = variables["digi-smear-parIndices"].as<read_series>();
    auto types = variables["digi-smear-types"].as<read_series>();
    auto parameters = variables["digi-smear-parameters"].as<read_range>();

    // Count the number of needed parameters
    // - type 0 needs 1 parameters
    // - type 1-4 needs 4 parameters:
    //      1 principle parameter, 2 range parameters, 1 range indicator
    size_t sumpars = 0;
    std::for_each(types.begin(), types.end(),
                  [&](int n) { sumpars += (n == 0) ? 1 : 4; });

    if (volumes.size() == vdims.size() and pindices.size() == types.size() and
        parameters.size() == sumpars) {
      size_t vpos = 0;
      size_t ppos = 0;

      ACTS_DEBUG("Volume parameters properly read from command line.")

      /// Extract the parameters from the global vector
      ///
      /// @param ftype the function type
      ///
      /// @return an extracted vector
      auto extract = [&](int ftype) -> std::vector<double> {
        size_t padd = (ftype == 0) ? 1 : 3;
        std::vector<double> fpars = {parameters.begin() + ppos,
                                     parameters.begin() + ppos + padd};
        // unfortunately boost::program_options can not deal with negative
        // values within a sequence, we thus have to flag the range at defined
        // with command line options:
        //    < 0 (both negative), 0 (symmetric), > 0 both positive
        if (padd == 3) {
          double rangeIndic = parameters[ppos + padd];
          fpars[1] *= (rangeIndic <= 0) ? -1. : 1.;
          fpars[2] *= (rangeIndic < 0) ? -1. : 1.;
          ++padd;
        }
        ppos += padd;
        return fpars;
      };

      std::vector<std::pair<Acts::GeometryIdentifier,
                            SmearingAlgorithm::SupportedSmearer> >
          smearers;
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
          if ((Acts::BoundIndices)pindices[vpos - 1] == Acts::eBoundLoc0) {
            smearCfg.configured = true;
            smearer =
                std::pair<ActsFatras::BoundParametersSmearer<Acts::eBoundLoc0>,
                          OneFnc>(
                    ActsFatras::BoundParametersSmearer<Acts::eBoundLoc0>(),
                    std::move(smearFunctions));
          } else if ((Acts::BoundIndices)pindices[vpos - 1] ==
                     Acts::eBoundLoc1) {
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
          if ((Acts::BoundIndices)pindices[vpos - 2] == Acts::eBoundLoc0 and
              (Acts::BoundIndices) pindices[vpos - 1] == Acts::eBoundTime) {
            smearCfg.configured = true;
            smearer =
                std::pair<ActsFatras::BoundParametersSmearer<Acts::eBoundLoc0,
                                                             Acts::eBoundTime>,
                          TwoFnc>(
                    ActsFatras::BoundParametersSmearer<Acts::eBoundLoc0,
                                                       Acts::eBoundTime>(),
                    std::move(smearFunctions));
          } else if ((Acts::BoundIndices)pindices[vpos - 2] ==
                         Acts::eBoundLoc1 and
                     (Acts::BoundIndices)
                             pindices[vpos - 1] == Acts::eBoundTime) {
            smearCfg.configured = true;
            smearer =
                std::pair<ActsFatras::BoundParametersSmearer<Acts::eBoundLoc1,
                                                             Acts::eBoundTime>,
                          TwoFnc>(
                    ActsFatras::BoundParametersSmearer<Acts::eBoundLoc1,
                                                       Acts::eBoundTime>(),
                    std::move(smearFunctions));
          } else if ((Acts::BoundIndices)pindices[vpos - 2] ==
                         Acts::eBoundLoc0 and
                     (Acts::BoundIndices)
                             pindices[vpos - 1] == Acts::eBoundLoc1) {
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
          if ((Acts::BoundIndices)pindices[vpos - 3] == Acts::eBoundLoc0 and
              (Acts::BoundIndices) pindices[vpos - 2] == Acts::eBoundLoc1 and
              (Acts::BoundIndices) pindices[vpos - 1] == Acts::eBoundTime) {
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
        smearers.push_back({volumeGeometryId, std::move(smearer)});
      }
      smearCfg.smearers =
          Acts::GeometryHierarchyMap<SmearingAlgorithm::SupportedSmearer>(
              std::move(smearers));

    } else if (parameters.size() != sumpars) {
      ACTS_ERROR("Expected " << sumpars << " parameters, but received "
                             << parameters.size());
    }
  }
  return smearCfg;
}
