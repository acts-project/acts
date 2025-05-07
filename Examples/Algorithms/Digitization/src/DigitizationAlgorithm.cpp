// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Digitization/DigitizationAlgorithm.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "ActsExamples/Digitization/ModuleClusters.hpp"
#include "ActsExamples/EventData/GeometryContainers.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"

#include <algorithm>
#include <array>
#include <limits>
#include <ostream>
#include <stdexcept>
#include <string>
#include <utility>

namespace ActsExamples {

DigitizationAlgorithm::DigitizationAlgorithm(Config config,
                                             Acts::Logging::Level level)
    : IAlgorithm("DigitizationAlgorithm", level), m_cfg(std::move(config)) {
  if (m_cfg.inputSimHits.empty()) {
    throw std::invalid_argument("Missing simulated hits input collection");
  }
  if (m_cfg.surfaceByIdentifier.empty()) {
    throw std::invalid_argument("Missing Surface-GeometryID association map");
  }
  if (!m_cfg.randomNumbers) {
    throw std::invalid_argument("Missing random numbers tool");
  }
  if (m_cfg.digitizationConfigs.empty()) {
    throw std::invalid_argument("Missing digitization configuration");
  }

  if (m_cfg.doClusterization) {
    if (m_cfg.outputMeasurements.empty()) {
      throw std::invalid_argument("Missing measurements output collection");
    }
    if (m_cfg.outputClusters.empty()) {
      throw std::invalid_argument("Missing cluster output collection");
    }
    if (m_cfg.outputMeasurementParticlesMap.empty()) {
      throw std::invalid_argument(
          "Missing hit-to-particles map output collection");
    }
    if (m_cfg.outputMeasurementSimHitsMap.empty()) {
      throw std::invalid_argument(
          "Missing hit-to-simulated-hits map output collection");
    }
    if (m_cfg.outputParticleMeasurementsMap.empty()) {
      throw std::invalid_argument(
          "Missing particle-to-measurements map output collection");
    }
    if (m_cfg.outputSimHitMeasurementsMap.empty()) {
      throw std::invalid_argument(
          "Missing particle-to-simulated-hits map output collection");
    }

    m_outputMeasurements.initialize(m_cfg.outputMeasurements);
    m_outputClusters.initialize(m_cfg.outputClusters);
    m_outputMeasurementParticlesMap.initialize(
        m_cfg.outputMeasurementParticlesMap);
    m_outputMeasurementSimHitsMap.initialize(m_cfg.outputMeasurementSimHitsMap);
    m_outputParticleMeasurementsMap.initialize(
        m_cfg.outputParticleMeasurementsMap);
    m_outputSimHitMeasurementsMap.initialize(m_cfg.outputSimHitMeasurementsMap);
  }

  if (m_cfg.doOutputCells) {
    if (m_cfg.outputCells.empty()) {
      throw std::invalid_argument("Missing cell output collection");
    }

    m_outputCells.initialize(m_cfg.outputCells);
  }

  m_inputHits.initialize(m_cfg.inputSimHits);

  // Create the digitizers from the configuration
  std::vector<std::pair<Acts::GeometryIdentifier, Digitizer>> digitizerInput;

  for (std::size_t i = 0; i < m_cfg.digitizationConfigs.size(); ++i) {
    GeometricConfig geoCfg;
    Acts::GeometryIdentifier geoId = m_cfg.digitizationConfigs.idAt(i);

    const auto& digiCfg = m_cfg.digitizationConfigs.valueAt(i);
    geoCfg = digiCfg.geometricDigiConfig;
    // Copy so we can sort in-place
    SmearingConfig smCfg = digiCfg.smearingDigiConfig;

    std::vector<Acts::BoundIndices> indices;
    for (auto& gcf : smCfg.params) {
      indices.push_back(gcf.index);
    }
    indices.insert(indices.begin(), geoCfg.indices.begin(),
                   geoCfg.indices.end());

    // Make sure the configured input parameter indices are sorted and unique
    std::ranges::sort(indices);

    auto dup = std::adjacent_find(indices.begin(), indices.end());
    if (dup != indices.end()) {
      std::invalid_argument(
          "Digitization configuration contains duplicate parameter indices");
    }

    switch (smCfg.params.size()) {
      case 0u:
        digitizerInput.emplace_back(geoId, makeDigitizer<0u>(digiCfg));
        break;
      case 1u:
        digitizerInput.emplace_back(geoId, makeDigitizer<1u>(digiCfg));
        break;
      case 2u:
        digitizerInput.emplace_back(geoId, makeDigitizer<2u>(digiCfg));
        break;
      case 3u:
        digitizerInput.emplace_back(geoId, makeDigitizer<3u>(digiCfg));
        break;
      case 4u:
        digitizerInput.emplace_back(geoId, makeDigitizer<4u>(digiCfg));
        break;
      default:
        throw std::invalid_argument("Unsupported smearer size");
    }
  }

  m_digitizers = Acts::GeometryHierarchyMap<Digitizer>(digitizerInput);
}

ProcessCode DigitizationAlgorithm::execute(const AlgorithmContext& ctx) const {
  // Retrieve input
  const auto& simHits = m_inputHits(ctx);
  ACTS_DEBUG("Loaded " << simHits.size() << " sim hits");

  // Prepare output containers
  // need list here for stable addresses
  MeasurementContainer measurements;
  ClusterContainer clusters;

  IndexMultimap<SimBarcode> measurementParticlesMap;
  IndexMultimap<Index> measurementSimHitsMap;
  measurements.reserve(simHits.size());
  measurementParticlesMap.reserve(simHits.size());
  measurementSimHitsMap.reserve(simHits.size());

  // Setup random number generator
  auto rng = m_cfg.randomNumbers->spawnGenerator(ctx);

  // Some statistics
  std::size_t skippedHits = 0;

  // Some algorithms do the clusterization themselves such as the traccc chain.
  // Thus we need to store the cell data from the simulation.
  CellsMap cellsMap;

  ACTS_DEBUG("Starting loop over modules ...");
  for (const auto& simHitsGroup : groupByModule(simHits)) {
    // Manual pair unpacking instead of using
    //   auto [moduleGeoId, moduleSimHits] : ...
    // otherwise clang on macos complains that it is unable to capture the local
    // binding in the lambda used for visiting the smearer below.
    Acts::GeometryIdentifier moduleGeoId = simHitsGroup.first;
    const auto& moduleSimHits = simHitsGroup.second;

    auto surfaceItr = m_cfg.surfaceByIdentifier.find(moduleGeoId);

    if (surfaceItr == m_cfg.surfaceByIdentifier.end()) {
      // this is either an invalid geometry id or a misconfigured smearer
      // setup; both cases can not be handled and should be fatal.
      ACTS_ERROR("Could not find surface " << moduleGeoId
                                           << " for configured smearer");
      return ProcessCode::ABORT;
    }

    const Acts::Surface* surfacePtr = surfaceItr->second;

    auto digitizerItr = m_digitizers.find(moduleGeoId);
    if (digitizerItr == m_digitizers.end()) {
      ACTS_VERBOSE("No digitizer present for module " << moduleGeoId);
      continue;
    } else {
      ACTS_VERBOSE("Digitizer found for module " << moduleGeoId);
    }

    // Run the digitizer. Iterate over the hits for this surface inside the
    // visitor so we do not need to lookup the variant object per-hit.
    std::visit(
        [&](const auto& digitizer) {
          ModuleClusters moduleClusters(
              digitizer.geometric.segmentation, digitizer.geometric.indices,
              m_cfg.doMerge, m_cfg.mergeNsigma, m_cfg.mergeCommonCorner);

          for (auto h = moduleSimHits.begin(); h != moduleSimHits.end(); ++h) {
            const auto& simHit = *h;
            const auto simHitIdx = simHits.index_of(h);

            DigitizedParameters dParameters;

            if (simHit.depositedEnergy() < m_cfg.minEnergyDeposit) {
              ACTS_VERBOSE("Skip hit because energy deposit to small");
              continue;
            }

            // Geometric part - 0, 1, 2 local parameters are possible
            if (!digitizer.geometric.indices.empty()) {
              ACTS_VERBOSE("Configured to geometric digitize "
                           << digitizer.geometric.indices.size()
                           << " parameters.");
              const auto& cfg = digitizer.geometric;
              Acts::Vector3 driftDir = cfg.drift(simHit.position(), rng);
              auto channelsRes = m_channelizer.channelize(
                  simHit, *surfacePtr, ctx.geoContext, driftDir,
                  cfg.segmentation, cfg.thickness);
              if (!channelsRes.ok() || channelsRes->empty()) {
                ACTS_DEBUG(
                    "Geometric channelization did not work, skipping this "
                    "hit.");
                continue;
              }
              ACTS_VERBOSE("Activated " << channelsRes->size()
                                        << " channels for this hit.");
              dParameters =
                  localParameters(digitizer.geometric, *channelsRes, rng);
            }

            // Smearing part - (optionally) rest
            if (!digitizer.smearing.indices.empty()) {
              ACTS_VERBOSE("Configured to smear "
                           << digitizer.smearing.indices.size()
                           << " parameters.");
              auto res =
                  digitizer.smearing(rng, simHit, *surfacePtr, ctx.geoContext);
              if (!res.ok()) {
                ++skippedHits;
                ACTS_DEBUG("Problem in hit smearing, skip hit ("
                           << res.error().message() << ")");
                continue;
              }
              const auto& [par, cov] = res.value();
              for (Eigen::Index ip = 0; ip < par.rows(); ++ip) {
                dParameters.indices.push_back(digitizer.smearing.indices[ip]);
                dParameters.values.push_back(par[ip]);
                dParameters.variances.push_back(cov(ip, ip));
              }
            }

            // Check on success - threshold could have eliminated all channels
            if (dParameters.values.empty()) {
              ACTS_VERBOSE(
                  "Parameter digitization did not yield a measurement.");
              continue;
            }

            moduleClusters.add(std::move(dParameters), simHitIdx);
          }

          auto digitizeParametersResult = moduleClusters.digitizedParameters();

          // Store the cell data into a map.
          if (m_cfg.doOutputCells) {
            std::vector<Cluster::Cell> cells;
            for (const auto& [dParameters, simHitsIdxs] :
                 digitizeParametersResult) {
              for (const auto& cell : dParameters.cluster.channels) {
                cells.push_back(cell);
              }
            }
            cellsMap.insert({moduleGeoId, std::move(cells)});
          }

          if (m_cfg.doClusterization) {
            for (auto& [dParameters, simHitsIdxs] : digitizeParametersResult) {
              auto measurement =
                  createMeasurement(measurements, moduleGeoId, dParameters);
              clusters.emplace_back(std::move(dParameters.cluster));

              for (auto [i, simHitIdx] : Acts::enumerate(simHitsIdxs)) {
                measurementParticlesMap.emplace_hint(
                    measurementParticlesMap.end(), measurement.index(),
                    simHits.nth(simHitIdx)->particleId());
                measurementSimHitsMap.emplace_hint(measurementSimHitsMap.end(),
                                                   measurement.index(),
                                                   simHitIdx);
              }
            }
          }
        },
        *digitizerItr);
  }

  if (skippedHits > 0) {
    ACTS_WARNING(
        skippedHits
        << " skipped in Digitization. Enable DEBUG mode to see more details.");
  }

  if (m_cfg.doClusterization) {
    m_outputMeasurements(ctx, std::move(measurements));
    m_outputClusters(ctx, std::move(clusters));

    // invert them before they are moved
    m_outputParticleMeasurementsMap(
        ctx, invertIndexMultimap(measurementParticlesMap));
    m_outputSimHitMeasurementsMap(ctx,
                                  invertIndexMultimap(measurementSimHitsMap));

    m_outputMeasurementParticlesMap(ctx, std::move(measurementParticlesMap));
    m_outputMeasurementSimHitsMap(ctx, std::move(measurementSimHitsMap));
  }

  if (m_cfg.doOutputCells) {
    m_outputCells(ctx, std::move(cellsMap));
  }

  return ProcessCode::SUCCESS;
}

DigitizedParameters DigitizationAlgorithm::localParameters(
    const GeometricConfig& geoCfg,
    const std::vector<ActsFatras::Segmentizer::ChannelSegment>& channels,
    RandomEngine& rng) const {
  DigitizedParameters dParameters;

  const auto& binningData = geoCfg.segmentation.binningData();

  // For digital readout, the weight needs to be split in x and y
  std::array<double, 2u> pos = {0., 0.};
  std::array<double, 2u> totalWeight = {0., 0.};
  std::array<std::size_t, 2u> bmin = {std::numeric_limits<std::size_t>::max(),
                                      std::numeric_limits<std::size_t>::max()};
  std::array<std::size_t, 2u> bmax = {0, 0};

  // The component digital store
  std::array<std::set<std::size_t>, 2u> componentChannels;

  // Combine the channels
  for (const auto& ch : channels) {
    auto bin = ch.bin;
    double charge = geoCfg.charge(ch.activation, rng);
    // Loop and check
    if (charge > geoCfg.threshold) {
      double weight = geoCfg.digital ? 1. : charge;
      for (std::size_t ib = 0; ib < 2; ++ib) {
        if (geoCfg.digital && geoCfg.componentDigital) {
          // only fill component of this row/column if not yet filled
          if (!componentChannels[ib].contains(bin[ib])) {
            totalWeight[ib] += weight;
            pos[ib] += weight * binningData[ib].center(bin[ib]);
            componentChannels[ib].insert(bin[ib]);
          }
        } else {
          totalWeight[ib] += weight;
          pos[ib] += weight * binningData[ib].center(bin[ib]);
        }
        // min max channels
        bmin[ib] = std::min(bmin[ib], static_cast<std::size_t>(bin[ib]));
        bmax[ib] = std::max(bmax[ib], static_cast<std::size_t>(bin[ib]));
      }
      // Create a copy of the channel, as activation may change
      auto chdig = ch;
      chdig.bin = ch.bin;
      chdig.activation = charge;
      dParameters.cluster.channels.push_back(chdig);
    }
  }
  if (totalWeight[0] > 0. && totalWeight[1] > 0.) {
    pos[0] /= totalWeight[0];
    pos[1] /= totalWeight[1];
    dParameters.indices = geoCfg.indices;
    for (auto idx : dParameters.indices) {
      dParameters.values.push_back(pos[idx]);
    }
    std::size_t size0 = (bmax[0] - bmin[0] + 1);
    std::size_t size1 = (bmax[1] - bmin[1] + 1);

    dParameters.variances = geoCfg.variances({size0, size1}, bmin);
    dParameters.cluster.sizeLoc0 = size0;
    dParameters.cluster.sizeLoc1 = size1;
  }

  return dParameters;
}

}  // namespace ActsExamples
