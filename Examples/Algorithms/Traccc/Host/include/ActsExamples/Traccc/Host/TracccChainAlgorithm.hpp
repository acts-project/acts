// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// Traccc plugin include(s)
#include "ActsExamples/Traccc/Common/TracccChainAlgorithmBase.hpp"
#include "ActsExamples/Traccc/Host/Types.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

namespace ActsExamples::Traccc::Host {

class TracccChainAlgorithm final : public Common::TracccChainAlgorithmBase {
 public:
  /// Construct the traccc algorithm.
  ///
  /// @param cfg is the algorithm configuration
  /// @param lvl is the logging level
  TracccChainAlgorithm(Config cfg, Acts::Logging::Level lvl);

  std::tuple<vecmem::vector<traccc::measurement>, vecmem::vector<traccc::spacepoint>, vecmem::vector<traccc::seed>> runDigitization(const vecmem::vector<traccc::cell>& cells, const vecmem::vector<traccc::cell_module>& modules, vecmem::host_memory_resource& mr) const override;

  traccc::host_container<traccc::fitting_result<traccc::default_algebra>, traccc::track_state<traccc::default_algebra>> runReconstruction(const vecmem::vector<traccc::measurement> measurements, const vecmem::vector<traccc::spacepoint> spacepoints,  const vecmem::vector<traccc::seed> seeds, vecmem::host_memory_resource& mr) const override;

 private:
  using HostTypes =
      typename ActsExamples::Chain::Host::Types<typename FieldType::view_t>;

  typename HostTypes::ClusterizationAlgorithmType clusterizationAlgorithm;
  typename HostTypes::SpacepointFormationAlgorithmType
      spacepointFormationAlgorithm;
  typename HostTypes::SeedingAlgorithmType seedingAlgorithm;
  typename HostTypes::TrackParametersEstimationAlgorithmType
      trackParametersEstimationAlgorithm;
  typename HostTypes::FindingAlgorithmType findingAlgorithm;
  typename HostTypes::FittingAlgorithmType fittingAlgorithm;
  typename HostTypes::AmbiguityResolutionAlgorithmType
      ambiguityResolutionAlgorithm;
};

}  // namespace ActsExamples::Traccc::Host
