// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/AlignmentMillePede/ActsSolverFromMille.hpp"

#include "Acts/Definitions/Alignment.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Material/MaterialInteraction.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsAlignment/Kernel/Alignment.hpp"
#include "ActsAlignment/Kernel/detail/AlignmentEngine.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsPlugins/Mille/ActsToMille.hpp"

#include <memory>

#include <Mille/MilleFactory.h>

namespace ActsExamples {

ActsSolverFromMille::ActsSolverFromMille(
    Config cfg, std::unique_ptr<const Acts::Logger> logger)
    : IAlgorithm("ActsSolverFromMille", std::move(logger)),
      m_cfg(std::move(cfg)) {
  // retrieve tracking geo
  m_trackingGeometry = m_cfg.trackingGeometry;

  // instantiate the alignment tool instance
  Acts::Navigator::Config navcfg{m_cfg.trackingGeometry};
  navcfg.resolvePassive = false;
  navcfg.resolveMaterial = true;
  navcfg.resolveSensitive = true;
  Acts::Navigator navigator(navcfg,
                            this->logger().cloneWithSuffix("Navigator"));
  Stepper stepper{m_cfg.magneticField};
  m_align = std::make_shared<Alignment>(
      Fitter(Propagator(stepper, Acts::Navigator(navcfg))));
}

ProcessCode ActsSolverFromMille::execute(
    const AlgorithmContext& /*ClangShutUp*/) const {
  // this algorithm will not do anything at event-time
  return ProcessCode::SUCCESS;
}

ProcessCode ActsSolverFromMille::finalize() {
  auto milleReader = Mille::spawnMilleReader(m_cfg.milleInput);
  auto openStat = milleReader->open(m_cfg.milleInput);
  if (!openStat) {
    ACTS_FATAL("Failed to read the mille binary " << m_cfg.milleInput);
    return ProcessCode::ABORT;
  }

  // Assign indices to the alignable surfaces

  // We wish to have a relation between alignment parameter indices and real
  // geometry. The unordered_map does not give us this - so perform a manual
  // sorting.
  std::vector<std::pair<Acts::GeometryIdentifier, const Acts::Surface*>>
      sortedGeo;
  ActsAlignment::AlignmentResult alignResult;

  sortedGeo.insert(sortedGeo.end(),
                   m_trackingGeometry->geoIdSurfaceMap().begin(),
                   m_trackingGeometry->geoIdSurfaceMap().end());
  std::sort(
      sortedGeo.begin(), sortedGeo.end(),
      [&](const std::pair<Acts::GeometryIdentifier, const Acts::Surface*>& lhs,
          const std::pair<Acts::GeometryIdentifier, const Acts::Surface*>&
              rhs) { return (lhs.first.layer() < rhs.first.layer()); });

  const Acts::Surface* firstSurf = nullptr;
  unsigned int iSurface = 0;
  for (auto& [geoID, surface] : sortedGeo) {
    // only consider sensitive surfaces
    if (!surface->isSensitive()) {
      continue;
    }
    // use the first sensitive surface as trajectory reference in the kalman
    if (firstSurf == nullptr) {
      firstSurf = surface;
    }
    if (!m_cfg.fixModules.contains(geoID)) {
      alignResult.idxedAlignSurfaces.emplace(surface, iSurface);
      iSurface++;
    }
  }

  ACTS_INFO("Performing alignment fit on collected Mille records");
  std::vector<ActsAlignment::detail::TrackAlignmentState> alignmentStates;
  ActsAlignment::detail::TrackAlignmentState state;
  std::size_t iRec = 0;
  while (ActsPlugins::ActsToMille::unpackMilleRecord(
             *milleReader, state, alignResult.idxedAlignSurfaces) ==
         Mille::MilleDecoder::ReadResult::OK) {
    if (++iRec % 10000 == 0) {
      ACTS_INFO("     Reading input record " << iRec);
    }
    alignmentStates.push_back(state);
  }

  /// TODO: Should try a local iteration without track state info.
  /// Can use the linearised info in the Track Alignment State
  /// to calculate approximate track parameter & residual updates
  /// and then repeat the solution. As in Millepede, probably
  /// safe to keep the "big matrix" and only update the right hand side.
  m_align->calculateAlignmentParameters(alignmentStates, alignResult);

  /// in a real experiment, the results would be written out
  /// and stored e.g. in a DB file for further use / validation.
  /// For this initial demo, we just print them out.

  ACTS_INFO("Performed internal alignment. ");
  ACTS_INFO(std::setw(16) << "  Tracks used: " << alignmentStates.size());
  ACTS_INFO(std::setw(16) << "  avg Chi2/NDF = "
                          << alignResult.averageChi2ONdf);
  ACTS_INFO(std::setw(16) << "  Chi2   = " << alignResult.chi2);
  ACTS_INFO(std::setw(16) << "  delta Chi2   = " << alignResult.deltaChi2);
  ACTS_INFO(std::setw(16) << "  Alignment parameter updates: ");
  std::vector<std::string> parLabels{"dx", "dy", "dz", "rx", "ry", "rz"};
  for (auto [surface, index] : alignResult.idxedAlignSurfaces) {
    ACTS_INFO(std::setw(20)
              << " Surface with geo ID " << surface->geometryId() << ": ");
    for (std::size_t i = 0; i < Acts::eAlignmentSize; ++i) {
      std::size_t row = Acts::eAlignmentSize * index + i;
      ACTS_INFO(std::setw(20)
                << parLabels[i] << " = " << std::setw(10)
                << alignResult.deltaAlignmentParameters(row) << std::setw(6)
                << " +/- " << std::setw(10)
                << std::sqrt(alignResult.alignmentCovariance(row, row)));
    }
  }

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
