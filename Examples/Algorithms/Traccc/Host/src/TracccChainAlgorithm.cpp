// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Acts Examples include(s)
#include "ActsExamples/Traccc/Host/TracccChainAlgorithm.hpp"
#include "ActsExamples/Traccc/Common/Conversion/SeedConversion.hpp"
#include "ActsExamples/Traccc/Common/Conversion/SpacePointConversion.hpp"

// System include(s).
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <type_traits>

ActsExamples::Traccc::Host::TracccChainAlgorithm::TracccChainAlgorithm(
    Config cfg, Acts::Logging::Level lvl)
    : ActsExamples::Traccc::Common::TracccChainAlgorithmBase(std::move(cfg),
                                                             std::move(lvl)),
      clusterizationAlgorithm(hostMemoryResource),
      spacepointFormationAlgorithm(hostMemoryResource),
      seedingAlgorithm(m_cfg.chainConfig->seedfinderConfig,
                       m_cfg.chainConfig->spacepointGridConfig,
                       m_cfg.chainConfig->seedfilterConfig, hostMemoryResource),
      trackParametersEstimationAlgorithm(hostMemoryResource),
      findingAlgorithm(m_cfg.chainConfig->findingConfig),
      fittingAlgorithm(m_cfg.chainConfig->fittingConfig),
      ambiguityResolutionAlgorithm(
          m_cfg.chainConfig->ambiguityResolutionConfig) {}

ActsExamples::ProcessCode
ActsExamples::Traccc::Host::TracccChainAlgorithm::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  vecmem::host_memory_resource mr;

  typename HostTypes::ClusterizationAlgorithmType::output_type measurements{
      &mr};
  typename HostTypes::SpacepointFormationAlgorithmType::output_type spacepoints{
      &mr};
  typename HostTypes::SeedingAlgorithmType::output_type seeds{&mr};
  typename HostTypes::TrackParametersEstimationAlgorithmType::output_type
      params{&mr};
  typename HostTypes::FindingAlgorithmType::output_type trackCandidates{&mr};
  typename HostTypes::FittingAlgorithmType::output_type trackStates{&mr};
  typename HostTypes::AmbiguityResolutionAlgorithmType::output_type
      resolvedTrackStates{&mr};

  const auto cellsMap = m_inputCells(ctx);

  auto [cells, modules] = converter.convertCells(cellsMap, &mr);

  measurements = clusterizationAlgorithm(vecmem::get_data(cells),
                                         vecmem::get_data(modules));

  ACTS_INFO("Ran the clusterization algorithm");

  ACTS_INFO("Checking if measurements match with Acts...");

  const auto actsMeasurements = m_inputMeasurements(ctx);
  auto measurementMap = converter.createMeasurementMap(measurements, actsMeasurements);

  ACTS_INFO("Measurements match (found a 1:1 mapping of indexes between traccc and Acts measurements)");

  spacepoints = spacepointFormationAlgorithm(vecmem::get_data(measurements),
                                             vecmem::get_data(modules));

  auto [convertedSpacePoints, spacePointMap] = ActsExamples::Traccc::Common::Conversion::convertSpacePoints(spacepoints, measurementMap);
  
  ACTS_INFO("Ran the spacepoint formation algorithm");

  seeds = seedingAlgorithm(spacepoints);

  auto [convertedSeeds, seedMap] = ActsExamples::Traccc::Common::Conversion::convertSeeds(seeds, vecmem::get_data(spacepoints), spacePointMap);

  ACTS_INFO("Ran the seeding algorithm");

  const typename FieldType::view_t fieldView(field);

  // Traccc expects a field vector of a constant field.
  params = trackParametersEstimationAlgorithm(spacepoints, seeds,
                                              fieldView.at(0.f, 0.f, 0.f));

  ACTS_INFO("Ran the parameters estimation algorithm");

  trackCandidates = findingAlgorithm(detector, field, measurements, params);

  ACTS_INFO("Ran the finding algorithm");

  trackStates = fittingAlgorithm(detector, field, trackCandidates);

  ACTS_INFO("Ran the fitting algorithm");

  if (m_cfg.enableAmbiguityResolution){
    trackStates = ambiguityResolutionAlgorithm(trackStates);

    ACTS_INFO("Ran the ambiguity resolution algorithm");
  }
  else{
    ACTS_INFO("Skipped the ambiguity resolution algorithm");
  }

  auto result = converter.convertTracks(trackStates, measurementMap);

  m_outputSeeds(ctx, std::move(convertedSeeds));
  m_outputSpacePoints(ctx, std::move(convertedSpacePoints));
  m_outputTracks(ctx, std::move(result));
  return ActsExamples::ProcessCode::SUCCESS;
}
