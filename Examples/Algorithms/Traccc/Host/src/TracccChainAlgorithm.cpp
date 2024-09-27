// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Acts Examples include(s)
#include "ActsExamples/Traccc/Host/TracccChainAlgorithm.hpp"

#include "Acts/Plugins/Traccc/BarcodeMap.hpp"
#include "Acts/Plugins/Traccc/CellConversion.hpp"
#include "Acts/Plugins/Traccc/SurfaceMap.hpp"
#include "Acts/Plugins/Traccc/TrackConversion.hpp"
#include "ActsExamples/Traccc/Conversion/CellMapConversion.hpp"
#include "ActsExamples/Traccc/Conversion/DigitizationConversion.hpp"
#include "ActsExamples/Traccc/Conversion/MeasurementMatch.hpp"
#include "ActsExamples/Traccc/Conversion/SeedConversion.hpp"
#include "ActsExamples/Traccc/Conversion/SpacePointConversion.hpp"
#include "ActsExamples/Traccc/Conversion/TrackConversion.hpp"
#include "ActsExamples/Traccc/Io/ReadDetector.hpp"
#include "ActsExamples/Traccc/Measurement/Debug.hpp"
#include "ActsExamples/Traccc/Util/IndexMap.hpp"

#include "vecmem/utils/tuple.hpp"

ActsExamples::Traccc::Host::TracccChainAlgorithm::TracccChainAlgorithm(
    const Common::TracccChainConfig& cfg, Acts::Logging::Level lvl)
    : IAlgorithm("TracccHostChain", lvl),
      m_cfg(cfg),
      m_detector(Common::Io::readDetector<detector_t>(
          &m_hostMemoryResource, cfg.detectorFilePath, cfg.materialFilePath,
          cfg.gridFilePath)),
      m_field(Acts::CovfiePlugin::covfieField(*m_cfg.field)),
      m_surfaceTransforms{Acts::TracccPlugin::createSurfaceMap(m_detector)},
      m_barcodeMap{Acts::TracccPlugin::createBarcodeMap(m_detector)},
      m_digitizationConfig{
          Common::Conversion::tracccConfig(m_cfg.digitizationConfigs)},
      clusterizationAlgorithm(m_hostMemoryResource),
      spacepointFormationAlgorithm(m_hostMemoryResource),
      seedingAlgorithm(m_cfg.seedfinderConfig, m_cfg.spacepointGridConfig,
                       m_cfg.seedfilterConfig, m_hostMemoryResource),
      trackParametersEstimationAlgorithm(m_hostMemoryResource),
      findingAlgorithm(m_cfg.findingConfig),
      fittingAlgorithm(m_cfg.fittingConfig),
      ambiguityResolutionAlgorithm(m_cfg.ambiguityResolutionConfig) {}

std::tuple<vecmem::vector<traccc::measurement>,
           vecmem::vector<traccc::spacepoint>, vecmem::vector<traccc::seed>>
ActsExamples::Traccc::Host::TracccChainAlgorithm::runClusterization(
    const vecmem::vector<traccc::cell>& cells,
    const Types::SiliconDetectorDescriptionType::const_view& detector_desc,
    vecmem::host_memory_resource& mr) const {
  typename Types::ClusterizationAlgorithmType::output_type measurements{&mr};
  typename Types::SpacepointFormationAlgorithmType::output_type spacepoints{
      &mr};
  typename Types::SeedingAlgorithmType::output_type seeds{&mr};

  measurements =
      clusterizationAlgorithm(vecmem::get_data(cells), detector_desc);
  ACTS_INFO("Ran the clusterization algorithm");

  spacepoints = spacepointFormationAlgorithm(vecmem::get_data(measurements),
                                             detector_desc);
  ACTS_INFO("Ran the spacepoint formation algorithm");

  seeds = seedingAlgorithm(spacepoints);
  ACTS_INFO("Ran the seeding algorithm");

  return std::make_tuple(std::move(measurements), std::move(spacepoints),
                         std::move(seeds));
}

traccc::host_container<traccc::fitting_result<traccc::default_algebra>,
                       traccc::track_state<traccc::default_algebra>>
ActsExamples::Traccc::Host::TracccChainAlgorithm::runReconstruction(
    const vecmem::vector<traccc::measurement>& measurements,
    const vecmem::vector<traccc::spacepoint>& spacepoints,
    const vecmem::vector<traccc::seed>& seeds,
    vecmem::host_memory_resource& mr) const {
  typename Types::TrackParametersEstimationAlgorithmType::output_type params{
      &mr};
  typename Types::FindingAlgorithmType::output_type trackCandidates{&mr};
  typename Types::FittingAlgorithmType::output_type tracks{&mr};

  const field_t::view_t fieldView(m_field);

  // Traccc expects a field vector of a constant field.
  auto bv = fieldView.at(0.f, 0.f, 0.f);

  params = trackParametersEstimationAlgorithm(spacepoints, seeds,
                                              {bv[0], bv[1], bv[2]});
  ACTS_INFO("Ran the parameters estimation algorithm");

  trackCandidates = findingAlgorithm(m_detector, field_t::view_t(m_field),
                                     measurements, params);
  ACTS_INFO("Ran the finding algorithm");

  tracks =
      fittingAlgorithm(m_detector, field_t::view_t(m_field), trackCandidates);
  ACTS_INFO("Ran the fitting algorithm");

  if (m_cfg.enableAmbiguityResolution) {
    tracks = ambiguityResolutionAlgorithm(tracks);
    ACTS_INFO("Ran the ambiguity resolution algorithm");
  } else {
    ACTS_INFO("Skipped the ambiguity resolution algorithm");
  }

  return tracks;
}

ActsExamples::ProcessCode
ActsExamples::Traccc::Host::TracccChainAlgorithm::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  vecmem::host_memory_resource mr;

  // Read the cells
  const auto cellsMap = m_inputCells(ctx);

  // Convert the cells
  auto tcm = Common::Conversion::tracccCellsMap(cellsMap);
  auto [cells, modules] = Acts::TracccPlugin::createCellsAndModules(
      mr, tcm, m_surfaceTransforms, m_digitizationConfig, m_barcodeMap);

  // Run the traccc digitization
  auto [measurements, spacepoints, seeds] =
      runClusterization(cells, vecmem::get_data(modules), mr);

  // Now we have both traccc measurements and acts measurements
  // We have run the traccc digitization and we expect that Acts digitization
  // has also been run, since otherwise we cannot do compare and do truth
  // matching.

  // Read the acts measurements
  auto& actsMeasurements = m_inputMeasurements(ctx);
  // Determine which traccc measurements are equivalent to which Acts
  // measurements. This is needed since we cannot guarantee that the
  // measurements have the same ordering. We run the following to obtain a
  // mapping between the two measurement collections. Note: if the number of
  // measurements don't match during the perhaps mergeCommonCorner or doMerge
  // is false in the digitization algorithm configuration.
  MeasurementContainer convertedMeasurements;
  Common::Conversion::convertMeasurements(m_detector, measurements.cbegin(),
                                          measurements.cend(),
                                          convertedMeasurements);

  std::map<std::size_t, std::size_t> indexMap = Common::Util::matchMap(
      convertedMeasurements.cbegin(), convertedMeasurements.cend(),
      actsMeasurements.cbegin(), actsMeasurements.cend(),
      Common::Conversion::MeasurementGeoIDHash{},
      Common::Conversion::MeasurementAproxEquals{});

  auto measurementConv = Common::Util::makeConversionFromIndexMap(
      measurements.cbegin(), measurements.cend(), actsMeasurements.cbegin(),
      actsMeasurements.cend(), indexMap,
      Common::Conversion::TracccMeasurementHash{},
      std::equal_to<traccc::measurement>{});

  ACTS_DEBUG(std::string("Traccc (1) and Acts (2) measurement index pairing "
                         "information:\n") +
             Common::Measurement::pairingStatistics(
                 measurements, actsMeasurements, m_detector));

  // Check if we want to fetch the measurements, spacepoints, and seeds
  // instead of use the ones created by traccc.
  if (!m_cfg.reconstructionOnly) {
    // Convert the traccc spacepoints to traccc space points.
    // Create an empty container to hold the converted space points.
    SimSpacePointContainer convertedSpacePoints;
    auto spacePointConv =
        ActsExamples::Traccc::Common::Conversion::convertSpacePoints(
            spacepoints, measurementConv, convertedSpacePoints);

    // Repeat the process for the traccc seeds to obtain the seeds converted
    // to the Acts edm.
    SimSeedContainer convertedSeeds;
    ActsExamples::Traccc::Common::Conversion::convertSeeds(
        seeds, spacePointConv, convertedSeeds);

    // We have now obtained the traccc seeds as Acts seeds which is useful for
    // comparison. The converted seeds will be outputted along with the
    // converted tracks.

    // We now want to obtain the converted tracks.
    // We run the reconstruction with traccc.
    auto tracks = runReconstruction(measurements, spacepoints, seeds, mr);

    // Now we convert the traccc tracks to acts tracks.
    auto convertedTracks = Common::Conversion::convertTracks(
        tracks, measurementConv, *m_cfg.trackingGeometry, m_detector);

    // Write results.
    m_outputSpacePoints(ctx, std::move(convertedSpacePoints));
    m_outputSeeds(ctx, std::move(convertedSeeds));
    m_outputTracks(ctx, std::move(convertedTracks));

  } else {
    // Use externally generated measurements, spacepoints, and seeds instead
    // of the ones generated by traccc.
    ACTS_INFO(
        "Flag 'reconstruction only' set to true - discarding traccc "
        "digitization data and using external measurements, spacepoints, and "
        "seeds");

    auto invMeasurementConv =
        measurementConv.invert(Common::Conversion::ActsMeasurementHash{},
                               Common::Conversion::ActsMeasurementEquals{});

    // We have previously ensured that the traccc measurement and the
    // externally generated acts measurements are the same and obtained their
    // conversion data. Thus we do not need to convert the acts measurements
    // to traccc measurements.

    // Read the exteranlly generated spacepoints and convert them.
    auto& actsSpacePoints = m_inputSpacePoints(ctx);
    spacepoints.clear();
    auto spacePointConv =
        ActsExamples::Traccc::Common::Conversion::convertSpacePoints(
            actsSpacePoints, invMeasurementConv, spacepoints);

    // Read the exteranlly generated seeds and convert them.
    auto& actsSeeds = m_inputSeeds(ctx);
    seeds.clear();
    ActsExamples::Traccc::Common::Conversion::convertSeeds(
        actsSeeds, spacePointConv, seeds);

    // We run the reconstruction with traccc.
    auto tracks = runReconstruction(measurements, spacepoints, seeds, mr);

    // Now we convert the traccc tracks to acts tracks.
    auto convertedTracks = Common::Conversion::convertTracks(
        tracks, measurementConv, *m_cfg.trackingGeometry, m_detector);

    // Write results.
    m_outputTracks(ctx, std::move(convertedTracks));
  }

  return ActsExamples::ProcessCode::SUCCESS;
}
