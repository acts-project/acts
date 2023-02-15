// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsFatras/EventData/Barcode.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "ActsExamples/Framework/BareAlgorithm.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"

#include <string>
#include <vector>

namespace ActsExamples {
  
class MeasurementMapSelectorAlgorithm final : public BareAlgorithm {
 public:
  struct Config {
    /// Input spacepoints collection.
    std::string inputMeasurementParticleMap;
    
    /// Input source links
    std::string inputSourceLinks;

    /// Output protoTracks collection.
    std::string outputMeasurementParticleMap;

    /// What spacepoints to keep
    std::vector<Acts::GeometryIdentifier> geometrySelection;
  };

  /// Constructor of the track finding algorithm
  ///
  /// @param cfg is the config struct to configure the algorithm
  /// @param level is the logging level
  MeasurementMapSelectorAlgorithm(Config cfg, Acts::Logging::Level lvl) :
    BareAlgorithm("SourceLinkSelection", lvl), m_cfg(cfg) {}

  virtual ~MeasurementMapSelectorAlgorithm() {}

  /// Filter the measurements
  ///
  /// @param ctx is the algorithm context that holds event-wise information
  /// @return a process code to steer the algorithm flow
  ActsExamples::ProcessCode execute(
      const ActsExamples::AlgorithmContext& ctx) const final {
    using HitParticlesMap = IndexMultimap<ActsFatras::Barcode>;
    const auto& inputSourceLinks = ctx.eventStore.get<IndexSourceLinkContainer>(m_cfg.inputSourceLinks);
    const auto& inputMap =
        ctx.eventStore.get<HitParticlesMap>(m_cfg.inputMeasurementParticleMap);
        
    HitParticlesMap outputMap;
    
    for (const auto geoId : m_cfg.geometrySelection) {
      auto range = selectLowestNonZeroGeometryObject(inputSourceLinks, geoId);
      auto groupedByModule = makeGroupBy(range, detail::GeometryIdGetter());

      for (const auto [moduleGeoId, moduleSourceLinks] : groupedByModule) {
        for (const auto& sourceLink : moduleSourceLinks) {
          const auto [begin, end] = inputMap.equal_range(sourceLink.index());
          
          outputMap.insert(begin, end);
        }
      }
    }
    
    ctx.eventStore.add(m_cfg.outputMeasurementParticleMap, std::move(outputMap));
    
    return ProcessCode::SUCCESS;
  }

  const Config& config() const { return m_cfg; }

 private:
  // configuration
  Config m_cfg;
};

}  // namespace ActsExamples
 

