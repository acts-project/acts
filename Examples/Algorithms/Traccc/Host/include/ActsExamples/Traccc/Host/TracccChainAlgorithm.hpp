// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Plugins/Covfie/FieldConversion.hpp"
#include "Acts/Plugins/Traccc/DigitizationConfig.hpp"
#include "ActsExamples/EventData/Cluster.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/SimSeed.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Traccc/Host/Types.hpp"
#include "ActsExamples/Traccc/TracccChainConfig.hpp"

#include <covfie/core/backend/primitive/array.hpp>
#include <vecmem/memory/host_memory_resource.hpp>

#include "traccc/geometry/geometry.hpp"

namespace ActsExamples::Traccc::Host {

class TracccChainAlgorithm : public IAlgorithm {
 public:
  using Config = Common::TracccChainConfig;

  /// Construct the traccc algorithm.
  ///
  /// @param cfg is the algorithm configuration
  /// @param lvl is the logging level
  TracccChainAlgorithm(const Config& cfg, Acts::Logging::Level lvl);

 private:
  using field_t =
      covfie::field<covfie::backend::constant<covfie::vector::float3,
                                              covfie::vector::float3>>;
  using detector_t =
      detray::detector<detray::default_metadata, detray::host_container_types>;
  using cell_map_t =
      std::map<Acts::GeometryIdentifier, std::vector<Cluster::Cell>>;

  Config m_cfg;
  vecmem::host_memory_resource m_hostMemoryResource;
  const field_t m_field;
  const std::map<Acts::GeometryIdentifier, traccc::transform3>
      m_surfaceTransforms;
  const std::map<Acts::GeometryIdentifier, detray::geometry::barcode>
      m_barcodeMap;
  const Acts::TracccPlugin::DigitizationConfig m_digitizationConfig;

  using Types = ActsExamples::Traccc::Host::Types;

  typename Types::ClusterizationAlgorithmType clusterizationAlgorithm;
  typename Types::SpacepointFormationAlgorithmType spacepointFormationAlgorithm;
  typename Types::SeedingAlgorithmType seedingAlgorithm;
  typename Types::TrackParametersEstimationAlgorithmType
      trackParametersEstimationAlgorithm;
  typename Types::FindingAlgorithmType findingAlgorithm;
  typename Types::FittingAlgorithmType fittingAlgorithm;
  typename Types::AmbiguityResolutionAlgorithmType ambiguityResolutionAlgorithm;

  ReadDataHandle<cell_map_t> m_inputCells{this, "InputCells"};
  ReadDataHandle<MeasurementContainer> m_inputMeasurements{this,
                                                           "InputMeasurements"};
  ReadDataHandle<SimSpacePointContainer> m_inputSpacePoints{this,
                                                            "InputSpacePoints"};
  ReadDataHandle<SimSeedContainer> m_inputSeeds{this, "InputSeeds"};
  WriteDataHandle<std::vector<SimSpacePoint>> m_outputSpacePoints{
      this, "OutputSpacePoints"};
  WriteDataHandle<std::vector<SimSeed>> m_outputSeeds{this, "OutputSeeds"};
  WriteDataHandle<ConstTrackContainer> m_outputTracks{this, "OutputTracks"};

 public:
  ActsExamples::ProcessCode execute(
      const ActsExamples::AlgorithmContext& ctx) const override;

  const Config& config() const { return m_cfg; }
};

}  // namespace ActsExamples::Traccc::Host
