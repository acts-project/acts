// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/SpacePointFormation2/PixelSpacePointBuilder.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <string>
#include <unordered_map>

#include <traccc/edm/measurement_collection.hpp>
#include <traccc/edm/spacepoint_collection.hpp>
#include <vecmem/memory/host_memory_resource.hpp>

#include "traccc/geometry/detector.hpp"
#include "traccc/geometry/host_detector.hpp"
#include "traccc/io/read_detector.hpp"

namespace ActsExamples {

class ActsMeasToTracccAlg final : public IAlgorithm {
 public:
  struct Config {
    std::string detectorFile = "";
    std::string inputActsMeasurements = "measurements";
    std::string outputDetrayToActsMap = "detray-to-acts-map";
    std::string outputTracccMeasurements = "acts-to-traccc-measurements";
    std::vector<int> pixelVolumes;
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;
  };

  mutable vecmem::host_memory_resource m_mr;
  traccc::host_detector m_host_det;
  std::unordered_map<Acts::GeometryIdentifier, std::uint64_t> m_actsToDetrayMap;
  std::unordered_map<std::uint64_t, Acts::GeometryIdentifier> m_detrayToActsMap;

  explicit ActsMeasToTracccAlg(
      const Config& cfg, std::unique_ptr<const Acts::Logger> logger = nullptr);

  ProcessCode execute(const AlgorithmContext& ctx) const override;
  void buildSurfaceMap(
      const Acts::TrackingGeometry& trackingGeometry,
      const std::string& detectorFile,
      std::unordered_map<std::uint64_t, Acts::GeometryIdentifier>&
          m_detrayToActsMap,
      std::unordered_map<Acts::GeometryIdentifier, std::uint64_t>&
          m_actsToDetrayMap);
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;

  ReadDataHandle<MeasurementContainer> m_inputActsMeasurements{
      this, "inputActsMeasurements"};
  WriteDataHandle<std::unordered_map<std::uint64_t, Acts::GeometryIdentifier>>
      m_outputDetrayToActsMap{this, "outputDetrayToActsMap"};
  WriteDataHandle<traccc::edm::measurement_collection::host>
      m_outputTracccMeasurements{this, "outputTracccMeasurements"};
};

}  // namespace ActsExamples
