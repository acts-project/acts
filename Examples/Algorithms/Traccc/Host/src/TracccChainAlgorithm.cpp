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



std::tuple<vecmem::vector<traccc::measurement>, vecmem::vector<traccc::spacepoint>, vecmem::vector<traccc::seed>> ActsExamples::Traccc::Host::TracccChainAlgorithm::runDigitization(const vecmem::vector<traccc::cell>& cells, const vecmem::vector<traccc::cell_module>& modules, vecmem::host_memory_resource& mr) const {
  typename HostTypes::ClusterizationAlgorithmType::output_type measurements{
      &mr};
  typename HostTypes::SpacepointFormationAlgorithmType::output_type spacepoints{
      &mr};
  typename HostTypes::SeedingAlgorithmType::output_type seeds{&mr};

  measurements = clusterizationAlgorithm(vecmem::get_data(cells),
                                         vecmem::get_data(modules));
  
  ACTS_INFO("Ran the clusterization algorithm");

  spacepoints = spacepointFormationAlgorithm(vecmem::get_data(measurements),
                                             vecmem::get_data(modules));

  ACTS_INFO("Ran the spacepoint formation algorithm");

  seeds = seedingAlgorithm(spacepoints);

  ACTS_INFO("Ran the seeding algorithm");

  return std::make_tuple(std::move(measurements), std::move(spacepoints), std::move(seeds));
}

traccc::host_container<traccc::fitting_result<traccc::default_algebra>, traccc::track_state<traccc::default_algebra>> ActsExamples::Traccc::Host::TracccChainAlgorithm::runReconstruction(const vecmem::vector<traccc::measurement> measurements, const vecmem::vector<traccc::spacepoint> spacepoints,  const vecmem::vector<traccc::seed> seeds, vecmem::host_memory_resource& mr) const {
  typename HostTypes::TrackParametersEstimationAlgorithmType::output_type
    params{&mr};
  typename HostTypes::FindingAlgorithmType::output_type trackCandidates{&mr};
  typename HostTypes::FittingAlgorithmType::output_type tracks{&mr};
  
  const typename FieldType::view_t fieldView(field);

  // Traccc expects a field vector of a constant field.
  params = trackParametersEstimationAlgorithm(spacepoints, seeds,
                                              fieldView.at(0.f, 0.f, 0.f));

  ACTS_INFO("Ran the parameters estimation algorithm");

  trackCandidates = findingAlgorithm(detector, field, measurements, params);

  ACTS_INFO("Ran the finding algorithm");

  tracks = fittingAlgorithm(detector, field, trackCandidates);

  ACTS_INFO("Ran the fitting algorithm");

  if (m_cfg.enableAmbiguityResolution){
    tracks = ambiguityResolutionAlgorithm(tracks);

    ACTS_INFO("Ran the ambiguity resolution algorithm");
  }
  else{
    ACTS_INFO("Skipped the ambiguity resolution algorithm");
  }

  return tracks;
}

