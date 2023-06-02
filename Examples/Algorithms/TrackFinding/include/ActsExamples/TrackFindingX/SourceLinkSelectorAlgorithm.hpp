// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/Framework/BareAlgorithm.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <string>
#include <vector>

namespace ActsExamples {

class SourceLinkSelectorAlgorithm final : public BareAlgorithm {
 public:
  struct Config {
    /// Input spacepoints collection.
    std::string inputSourceLinks;

    /// Output protoTracks collection.
    std::string outputSourceLinks;

    /// What spacepoints to keep
    std::vector<Acts::GeometryIdentifier> geometrySelection;
  };

  /// Constructor of the track finding algorithm
  ///
  /// @param cfg is the config struct to configure the algorithm
  /// @param level is the logging level
  SourceLinkSelectorAlgorithm(Config cfg, Acts::Logging::Level lvl)
      : BareAlgorithm("SourceLinkSelection", lvl), m_cfg(cfg) {}

  virtual ~SourceLinkSelectorAlgorithm() {}

  /// Filter the measurements
  ///
  /// @param ctx is the algorithm context that holds event-wise information
  /// @return a process code to steer the algorithm flow
  ActsExamples::ProcessCode execute(
      const ActsExamples::AlgorithmContext& ctx) const final {
    const auto& inputSourceLinks =
        ctx.eventStore.get<IndexSourceLinkContainer>(m_cfg.inputSourceLinks);

    IndexSourceLinkContainer outputSourceLinks;

    for (const auto geoId : m_cfg.geometrySelection) {
      auto range = selectLowestNonZeroGeometryObject(inputSourceLinks, geoId);
      auto groupedByModule = makeGroupBy(range, detail::GeometryIdGetter());

      for (auto [moduleGeoId, moduleSourceLinks] : groupedByModule) {
        for (auto& sourceLink : moduleSourceLinks) {
          outputSourceLinks.insert(sourceLink);
        }
      }
    }

    ctx.eventStore.add<IndexSourceLinkContainer>(m_cfg.outputSourceLinks,
                                                 std::move(outputSourceLinks));

    return ProcessCode::SUCCESS;
  }

  const Config& config() const { return m_cfg; }

 private:
  // configuration
  Config m_cfg;
};

}  // namespace ActsExamples
