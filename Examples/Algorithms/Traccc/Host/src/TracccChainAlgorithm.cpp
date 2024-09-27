// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Traccc/Host/TracccChainAlgorithm.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Plugins/Traccc/BarcodeMap.hpp"
#include "Acts/Plugins/Traccc/CellConversion.hpp"
#include "Acts/Plugins/Traccc/SurfaceMap.hpp"
#include "Acts/Plugins/Traccc/TrackConversion.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/Traccc/Conversion/CellMapConversion.hpp"
#include "ActsExamples/Traccc/Conversion/DigitizationConversion.hpp"
#include "ActsExamples/Traccc/Conversion/MeasurementConversion.hpp"
#include "ActsExamples/Traccc/Conversion/MeasurementMap.hpp"
#include "ActsExamples/Traccc/Conversion/SeedConversion.hpp"
#include "ActsExamples/Traccc/Conversion/SpacePointConversion.hpp"
#include "ActsExamples/Traccc/Conversion/TrackConversion.hpp"

#include <algorithm>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <vector>

#include "traccc/edm/measurement.hpp"

ActsExamples::Traccc::Host::TracccChainAlgorithm::TracccChainAlgorithm(
    const Common::TracccChainConfig& cfg, Acts::Logging::Level lvl)
    : IAlgorithm("TracccHostChain", lvl),
      m_cfg(cfg),
      m_field(Acts::CovfiePlugin::covfieField(
          m_cfg.field != nullptr ? *m_cfg.field
                                 : throw std::invalid_argument(
                                       "Missing magnetic field, without which "
                                       "the traccc algorithm cannot run."))),
      m_surfaceTransforms{Acts::TracccPlugin::createSurfaceMap(
          m_cfg.detrayStore != nullptr
              ? m_cfg.detrayStore->detector
              : throw std::invalid_argument(
                    "Missing detray detector, without which the traccc "
                    "algorithm cannot run."))},
      m_barcodeMap{Acts::TracccPlugin::createBarcodeMap(
          m_cfg.detrayStore != nullptr
              ? m_cfg.detrayStore->detector
              : throw std::invalid_argument(
                    "Missing detray detector, without which the traccc "
                    "algorithm cannot run."))},
      m_digitizationConfig{
          Common::Conversion::tracccConfig(m_cfg.digitizationConfigs)},
      clusterizationAlgorithm(m_hostMemoryResource),
      spacepointFormationAlgorithm(m_hostMemoryResource),
      seedingAlgorithm(m_cfg.seedfinderConfig, m_cfg.spacepointGridConfig,
                       m_cfg.seedfilterConfig, m_hostMemoryResource),
      trackParametersEstimationAlgorithm(m_hostMemoryResource),
      findingAlgorithm(m_cfg.findingConfig),
      fittingAlgorithm(m_cfg.fittingConfig),
      ambiguityResolutionAlgorithm(m_cfg.ambiguityResolutionConfig) {
  ACTS_INFO("traccc chain algorithm has a barcode map with "
            << m_barcodeMap.size() << " elements");

  if (m_cfg.digitizationConfigs.empty()) {
    throw std::invalid_argument(
        "Digitization configuration is empty; cannot run any meaningful "
        "reconstruction.");
  }

  if (m_cfg.trackingGeometry == nullptr) {
    throw std::invalid_argument(
        "Missing ACTS tracking geometry, without which the traccc algorithm "
        "cannot run.");
  }

  if (m_cfg.inputCells.empty()) {
    throw std::invalid_argument(
        "Cell input was not given, without which the traccc algorithm cannot "
        "run.");
  }
  m_inputCells.initialize(m_cfg.inputCells);

  if (m_cfg.inputMeasurements.empty()) {
    throw std::invalid_argument(
        "Measurement input was not given, without which the traccc algorithm "
        "cannot run (yet).");
  }
  m_inputMeasurements.initialize(m_cfg.inputMeasurements);

  if (!m_cfg.inputSpacePoints.empty()) {
    if (m_cfg.inputSeeds.empty()) {
      throw std::invalid_argument(
          "Spacepoint input was given, but seed input was not; this is "
          "invalid.");
    }
    m_inputSpacePoints.initialize(m_cfg.inputSpacePoints);
    m_inputSeeds.initialize(m_cfg.inputSeeds);
  }

  if (!m_cfg.outputSpacePoints.empty()) {
    m_outputSpacePoints.initialize(m_cfg.outputSpacePoints);
  }

  if (!m_cfg.outputSeeds.empty()) {
    m_outputSeeds.initialize(m_cfg.outputSeeds);
  }

  if (!m_cfg.outputTracks.empty()) {
    m_outputTracks.initialize(m_cfg.outputTracks);
  }
}

ActsExamples::ProcessCode
ActsExamples::Traccc::Host::TracccChainAlgorithm::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  vecmem::host_memory_resource mr;
  // Read the cells
  const auto cellsMap = m_inputCells(ctx);

  ACTS_INFO("Read " << std::accumulate(
                           cellsMap.begin(), cellsMap.end(), std::size_t{0},
                           [](const std::size_t& acc, const auto& i) {
                             return acc + i.second.size();
                           })
                    << " cells from " << cellsMap.size() << " modules");

  // Convert the cells
  auto tcm = Common::Conversion::tracccCellsMap(cellsMap);
  auto [cells, modules] = Acts::TracccPlugin::createCellsAndModules(
      mr, tcm, m_surfaceTransforms, m_digitizationConfig, m_barcodeMap,
      logger().cloneWithSuffix("CreateCellsAndModules"));

  ACTS_INFO("Beginning reconstruction algorithm; have cell map of size "
            << tcm.size() << ", cell vector of size " << cells.size()
            << ", and module vector of size " << modules.size());

  auto measurements = clusterizationAlgorithm(vecmem::get_data(cells),
                                              vecmem::get_data(modules));
  ACTS_INFO("Ran the clustering algorithm; produced " << measurements.size()
                                                      << " measurements");

  auto spacepoints = spacepointFormationAlgorithm(
      vecmem::get_data(measurements), vecmem::get_data(modules));
  ACTS_INFO("Ran the spacepoint formation algorithm; produced "
            << spacepoints.size() << " spacepoints");

  auto seeds = seedingAlgorithm(spacepoints);
  ACTS_INFO("Ran the seeding algorithm; produced " << seeds.size() << " seeds");

  // Now we have both traccc measurements and acts measurements
  // We have run the traccc digitization and we expect that Acts
  // digitization has also been run, since otherwise we cannot do compare
  // and do truth matching.

  // Read the acts measurements
  auto& actsMeasurements = m_inputMeasurements(ctx);
  // Determine which traccc measurements are equivalent to which Acts
  // measurements. This is needed since we cannot guarantee that the
  // measurements have the same ordering. We run the following to obtain a
  // mapping between the two measurement collections. Note: if the number of
  // measurements don't match during the perhaps mergeCommonCorner or
  // doMerge is false in the digitization algorithm configuration.
  ACTS_INFO("ACTS measurement container contains " << actsMeasurements.size()
                                                   << " measurements");

  std::map<std::size_t, std::size_t> tracccToActsMeasurementIndexMap =
      Common::Conversion::makeTracccToActsMeasurementIndexMap(
          measurements, actsMeasurements, m_cfg.detrayStore->detector);

  ACTS_INFO("traccc to ACTS measurement index map has "
            << tracccToActsMeasurementIndexMap.size() << " keys");

  if (m_inputSpacePoints.isInitialized() && m_inputSeeds.isInitialized()) {
    ACTS_INFO("Using ACTS spacepoints and seeds; discarding traccc results");

    std::map<std::size_t, std::size_t> actsToTracccMeasurementIndexMap =
        Common::Conversion::makeActsToTracccMeasurementIndexMap(
            actsMeasurements, measurements, m_cfg.detrayStore->detector);

    // Read the exteranlly generated spacepoints and convert them.
    ACTS_INFO("Reading and converting ACTS spacepoints to traccc spacepoints");
    auto& actsSpacePoints = m_inputSpacePoints(ctx);
    spacepoints.clear();

    auto [newSpacepoints, actsToTracccSpacepointIndexMap] = ActsExamples::
        Traccc::Common::Conversion::convertActsToTracccSpacePoints(
            actsSpacePoints, actsToTracccMeasurementIndexMap, measurements);

    spacepoints = std::move(newSpacepoints);

    ACTS_INFO("Converted " << actsSpacePoints.size() << " ACTS spacepoints to "
                           << spacepoints.size() << " traccc spacepoints");

    // Read the exteranlly generated seeds and convert them.
    ACTS_INFO("Reading and converting ACTS seeds to traccc seeds");
    auto& actsSeeds = m_inputSeeds(ctx);

    std::vector<traccc::seed, std::pmr::polymorphic_allocator<traccc::seed>>
        newSeeds =
            ActsExamples::Traccc::Common::Conversion::convertActsToTracccSeeds(
                actsSeeds, actsSpacePoints, actsToTracccSpacepointIndexMap);

    seeds = std::move(newSeeds);

    ACTS_INFO("Converted " << actsSeeds.size() << " ACTS seeds to "
                           << seeds.size() << " traccc seeds");
  } else {
    ACTS_INFO("Continuing with traccc data (not using converted ACTS results)");

    if (m_outputSpacePoints.isInitialized() || m_outputSeeds.isInitialized()) {
      ACTS_INFO(
          "Performing conversions to output traccc spacepoint and seed data");
      MeasurementContainer convertedMeasurements =
          Common::Conversion::convertTracccToActsMeasurements(
              m_cfg.detrayStore->detector, measurements.cbegin(),
              measurements.cend());

      ACTS_INFO("Converted traccc measurement container contains "
                << convertedMeasurements.size() << " measurements");

      if (actsMeasurements.size() != convertedMeasurements.size()) {
        ACTS_WARNING(
            "ACTS produced a different number of measurements than traccc!");
      }

      std::map<std::size_t, std::size_t> tmpTracccToActsMeasurementIndexMap =
          Common::Conversion::makeTracccToActsMeasurementIndexMap(
              measurements, convertedMeasurements, m_cfg.detrayStore->detector);

      ACTS_INFO("Converting traccc spacepoints to ACTS spacepoints");

      auto [convertedSpacePoints, spacepointIndexMapTracccToActs] =
          ActsExamples::Traccc::Common::Conversion::
              convertTracccToActsSpacePoints(spacepoints,
                                             tmpTracccToActsMeasurementIndexMap,
                                             convertedMeasurements);

      ACTS_INFO("Converted " << spacepoints.size() << " traccc spacepoints to "
                             << convertedSpacePoints.size()
                             << " ACTS spacepoints");

      if (m_outputSpacePoints.isInitialized()) {
        if (m_outputSeeds.isInitialized()) {
          m_outputSpacePoints(ctx,
                              std::vector<SimSpacePoint>(convertedSpacePoints));
        } else {
          m_outputSpacePoints(ctx, std::move(convertedSpacePoints));
        }
        ACTS_INFO(
            "Wrote ACTS spacepoints (converted from traccc) to whiteboard");
      }

      if (m_outputSeeds.isInitialized()) {
        ACTS_INFO("Converting traccc seeds to ACTS seeds");
        SimSeedContainer convertedSeeds =
            ActsExamples::Traccc::Common::Conversion::convertTracccToActsSeeds(
                seeds, spacepointIndexMapTracccToActs, convertedSpacePoints);
        ACTS_INFO("Converted " << seeds.size() << " traccc seeds to "
                               << convertedSeeds.size() << " ACTS seeds");
        m_outputSeeds(ctx, std::move(convertedSeeds));
        ACTS_INFO("Wrote ACTS seeds (converted from traccc) to whiteboard");
      }
    }
  }

  if (m_outputTracks.isInitialized()) {
    // We run the reconstruction with traccc.
    const field_t::view_t fieldView(m_field);

    // Traccc expects a field vector of a constant field.
    auto bv = fieldView.at(0.f, 0.f, 0.f);

    auto params = trackParametersEstimationAlgorithm(spacepoints, seeds,
                                                     {bv[0], bv[1], bv[2]});
    ACTS_INFO("Ran the track parameter estimation algorithm; estimated "
              << params.size() << " track parameters");

    auto trackCandidates =
        findingAlgorithm(m_cfg.detrayStore->detector, field_t::view_t(m_field),
                         measurements, params);
    ACTS_INFO("Ran the finding algorithm; found " << trackCandidates.size()
                                                  << " track candidates");

    auto tracks = fittingAlgorithm(m_cfg.detrayStore->detector,
                                   field_t::view_t(m_field), trackCandidates);
    ACTS_INFO("Ran the fitting algorithm; fitted " << tracks.size()
                                                   << " tracks");

    if (m_cfg.enableAmbiguityResolution) {
      tracks = ambiguityResolutionAlgorithm(tracks);
      ACTS_INFO("Ran the ambiguity resolution algorithm; kept " << tracks.size()
                                                                << " tracks");
    } else {
      ACTS_INFO("Skipped the ambiguity resolution algorithm");
    }

    // Now we convert the traccc tracks to acts tracks.
    ConstTrackContainer convertedTracks =
        Common::Conversion::convertTracccToActsTracks(
            tracks, tracccToActsMeasurementIndexMap, actsMeasurements,
            *m_cfg.trackingGeometry, m_cfg.detrayStore->detector,
            logger().cloneWithSuffix("TracccToActsTracks"));

    m_outputTracks(ctx, std::move(convertedTracks));
  }

  return ActsExamples::ProcessCode::SUCCESS;
}
