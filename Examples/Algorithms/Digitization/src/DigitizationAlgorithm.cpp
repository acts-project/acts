// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Digitization/DigitizationAlgorithm.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Plugins/Digitization/Clusterization.hpp"
#include "ActsExamples/EventData/GeometryContainers.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsFatras/Digitization/UncorrelatedHitSmearer.hpp"

#include <algorithm>
#include <set>
#include <stdexcept>
#include <string>
#include <type_traits>

// Anonymous namespace that contains types and code for hit merging
namespace {
using hit_t =
    ActsExamples::SimHitContainer::size_type;  // Only keep index into event
                                               // simhit container
using cell_t = ActsFatras::Channelizer::ChannelSegment;

struct SingleCell {
  hit_t simHit;
  cell_t value;

  SingleCell(size_t hitIdx, cell_t cell)
      : simHit(hitIdx), value(std::move(cell)){};

  double depositedEnergy() { return value.activation; }
};

struct TruthCluster {
  std::set<hit_t> simHits;
  std::vector<cell_t> cells;

  TruthCluster(size_t hitIdx, std::vector<cell_t> cells_ = {})
      : simHits{hitIdx}, cells(std::move(cells_)){};

  TruthCluster(SingleCell cell) : TruthCluster(cell.simHit, {cell.value}){};

  TruthCluster() : simHits(), cells(){};

  void add(SingleCell& cell) {
    simHits.insert(cell.simHit);
    cells.push_back(cell.value);
  }

  void add(hit_t hit) { simHits.insert(hit); }

  void add(TruthCluster clus) {
    for (hit_t hit : clus.simHits) {
      simHits.insert(hit);
    }
    for (cell_t& cell : clus.cells) {
      cells.push_back(std::move(cell));
    }
    clus.simHits.clear();
    clus.cells.clear();
  }
};

struct DigitizedCluster {
  std::set<hit_t> simHits;
  ActsExamples::DigitizedParameters params;

  DigitizedCluster(std::set<hit_t> hits) : simHits(std::move(hits)), params(){};

  template <typename T>
  void addSmearedParams(
      const std::unordered_map<hit_t, ActsExamples::DigitizedParameters>&
          measurementMap,
      const T& smearIdx) {
    std::vector<Acts::ActsScalar> values(smearIdx.size(), 0);
    std::vector<Acts::ActsScalar> variances(smearIdx.size(), 0);
    Acts::ActsScalar norm = 1.0 / simHits.size();

    for (hit_t hit : simHits) {
      auto parItr = measurementMap.find(hit);
      if (parItr != measurementMap.end()) {
        for (size_t i = 0; i < smearIdx.size(); i++) {
          values.at(i) += parItr->second.values.at(i) * norm;
          variances.at(i) += parItr->second.variances.at(i) * norm * norm;
        }
      }
    }
    for (size_t i = 0; i < smearIdx.size(); i++) {
      params.indices.push_back(smearIdx.at(i));
      params.values.push_back(values.at(i));
      params.variances.push_back(values.at(i));
    }
  }
};

auto findInMap(std::unordered_map<hit_t, ActsExamples::DigitizedParameters>&
                   measurementMap,
               SingleCell& cell) {
  return measurementMap.find(cell.simHit);
}

// Return first match
auto findInMap(std::unordered_map<hit_t, ActsExamples::DigitizedParameters>&
                   measurementMap,
               TruthCluster& clus) {
  auto end = measurementMap.end();
  for (hit_t hit : clus.simHits) {
    auto it = measurementMap.find(hit);
    if (it != end)
      return it;
  }
  return end;
}

template <typename T>
std::vector<TruthCluster> mergeMeasurements(
    std::vector<T>& cells,
    std::unordered_map<hit_t, ActsExamples::DigitizedParameters>&
        measurementMap,
    double nsigma) {
  std::vector<TruthCluster> clusters;
  std::vector<bool> used(cells.size(), false);

  for (size_t i = 0; i < cells.size(); i++) {
    if (used.at(i))
      continue;

    // Cell has not yet been claimed, so claim it
    clusters.push_back(TruthCluster(cells.at(i)));
    used.at(i) = true;

    // Cells previously visited by index `i' have already been added
    // to a cluster or used to seed a new cluster, so start at the
    // next unseen one
    for (size_t j = i + 1; j < cells.size(); j++) {
      // Still may have already been used, so check it
      if (used.at(j))
        continue;

      // Now, iterate through hits matched to current cluster and see
      // if any of them match the hit associated to the current cell
      auto& [h1, dparams] = *findInMap(measurementMap, cells.at(j)); // TODO check
      bool matched = false;
      for (auto hitIdx : clusters.back().simHits) {
        auto& [h2, dparams_2] = *measurementMap.find(hitIdx); // TODO check

        // Consider the cell matched to the currently considered
        // simhit until we find evidence to the contrary. This way,
        // merging still works when digitization is done by geometry
        // only.
        matched = true;

        // Obvious match when simhits are the same so can skip some
        // computations
        if (h1 == h2) {
          break;
        }

	// Loop over smeared directions Important note: we assume that
	// at this point the digitized parameters do not containt
	// dimensions that will be computed from the cluster
        for (size_t k = 0; k < dparams.values.size() and matched; k++) {
          auto maxdist = nsigma * std::sqrt(dparams.variances.at(k));
          if (std::abs(dparams.values.at(k) - dparams_2.values.at(k)) >
              maxdist) {
            matched = false;
          }
        }
        if (matched) {
          // Cell matched at least one hit in the cluster, no need to
          // keep checking
          break;
        }
      }  // loop over `k' (dimensions)
      if (matched) {
        // Claim cell `j'
        used.at(j) = true;
        clusters.back().add(std::move(cells.at(j)));
      }
    }  // loop over 'j'
  }    // loop over `i'
  return clusters;
}

void mergeClusters(
    std::vector<TruthCluster>& clusters,
    const ActsExamples::GeometricConfig& geoCfg,
    std::unordered_map<hit_t, ActsExamples::DigitizedParameters> measurementMap,
    double nsigma = 1.0) {

  // Start by merging actual clusters when there are any. To do so,
  // first the cell map used by the clusterization code has to be
  // filled
  std::unordered_map<size_t, std::pair<SingleCell, bool>> cellMap;
  for (TruthCluster& clus : clusters) {
    // TODO validate that at this stage only one simhit per cluster?
    for (cell_t& cell : clus.cells) {
      size_t index = cell.bin[0] + geoCfg.segmentation.bins(0) * cell.bin[1];
      // FIXME handle the cases where many hits pass through the same cell
      cellMap.insert(
          {index, {SingleCell(*clus.simHits.begin(), std::move(cell)), false}});
    }
  }
  if (not cellMap.empty()) {
    // Case where we actually have geometric clusters; use the
    // clusterization code from the digitization plugin
    std::vector<std::vector<SingleCell>> merged =
        Acts::createClusters(cellMap, geoCfg.segmentation.bins(0));
    for (std::vector<SingleCell>& cellv : merged) {
      // At this stage, the cellv vector contains cells that form a
      // consistent cluster based on a connected component analysis
      // only. Still have to check if they match based on the smeared
      // indices (a good example of this would a for a timing
      // detector). At this stage, the measurement map contains
      // parameters from smeared dimensions only, so re-merge based on
      // that.
      clusters = mergeMeasurements(cellv, measurementMap, nsigma);
    }
  } else {
    // Smeared measurements only -- merge based on that. There are no
    // cells, so merge the TruthClusters directly.
    clusters = mergeMeasurements(clusters, measurementMap, nsigma);
  }
}
}  // namespace

ActsExamples::DigitizationAlgorithm::DigitizationAlgorithm(
    DigitizationConfig cfg, Acts::Logging::Level lvl)
    : ActsExamples::BareAlgorithm("DigitizationAlgorithm", lvl),
      m_cfg(std::move(cfg)) {
  if (m_cfg.inputSimHits.empty()) {
    throw std::invalid_argument("Missing simulated hits input collection");
  }
  if (m_cfg.outputMeasurements.empty()) {
    throw std::invalid_argument("Missing measurements output collection");
  }
  if (m_cfg.outputSourceLinks.empty()) {
    throw std::invalid_argument("Missing source links output collection");
  }
  if (m_cfg.outputMeasurementParticlesMap.empty()) {
    throw std::invalid_argument(
        "Missing hit-to-particles map output collection");
  }
  if (m_cfg.outputMeasurementSimHitsMap.empty()) {
    throw std::invalid_argument(
        "Missing hit-to-simulated-hits map output collection");
  }
  if (not m_cfg.trackingGeometry) {
    throw std::invalid_argument("Missing tracking geometry");
  }
  if (not m_cfg.randomNumbers) {
    throw std::invalid_argument("Missing random numbers tool");
  }

  if (m_cfg.digitizationConfigs.empty()) {
    throw std::invalid_argument("Missing digitization configuration");
  }

  // create the smearers from the configuration
  std::vector<std::pair<Acts::GeometryIdentifier, Digitizer>> digitizerInput;

  for (size_t i = 0; i < m_cfg.digitizationConfigs.size(); ++i) {
    GeometricConfig geoCfg;
    Acts::GeometryIdentifier geoId = m_cfg.digitizationConfigs.idAt(i);

    const auto& digiCfg = m_cfg.digitizationConfigs.valueAt(i);
    geoCfg = digiCfg.geometricDigiConfig;
    // Copy so we can sort in-place
    SmearingConfig smCfg = digiCfg.smearingDigiConfig;

    std::vector<Acts::BoundIndices> indices;
    for (auto& gcf : smCfg) {
      indices.push_back(gcf.index);
    }
    indices.insert(indices.begin(), geoCfg.indices.begin(),
                   geoCfg.indices.end());

    // Make sure the configured input parameter indices are sorted and unique
    std::sort(indices.begin(), indices.end());

    auto dup = std::adjacent_find(indices.begin(), indices.end());
    if (dup != indices.end()) {
      std::invalid_argument(
          "Digitization configuration contains duplicate parameter indices");
    }

    switch (smCfg.size()) {
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

ActsExamples::ProcessCode ActsExamples::DigitizationAlgorithm::execute(
    const AlgorithmContext& ctx) const {
  // Retrieve input
  const auto& simHits = ctx.eventStore.get<SimHitContainer>(m_cfg.inputSimHits);

  // Prepare output containers
  IndexSourceLinkContainer sourceLinks;
  MeasurementContainer measurements;
  ClusterContainer clusters;
  IndexMultimap<ActsFatras::Barcode> measurementParticlesMap;
  IndexMultimap<Index> measurementSimHitsMap;
  sourceLinks.reserve(simHits.size());
  measurements.reserve(simHits.size());
  measurementParticlesMap.reserve(simHits.size());
  measurementSimHitsMap.reserve(simHits.size());

  // Setup random number generator
  auto rng = m_cfg.randomNumbers->spawnGenerator(ctx);

  ACTS_DEBUG("Starting loop over modules ...");
  for (auto simHitsGroup : groupByModule(simHits)) {
    // Manual pair unpacking instead of using
    //   auto [moduleGeoId, moduleSimHits] : ...
    // otherwise clang on macos complains that it is unable to capture the local
    // binding in the lambda used for visiting the smearer below.
    Acts::GeometryIdentifier moduleGeoId = simHitsGroup.first;
    const auto& moduleSimHits = simHitsGroup.second;

    const Acts::Surface* surfacePtr =
        m_cfg.trackingGeometry->findSurface(moduleGeoId);

    if (not surfacePtr) {
      // this is either an invalid geometry id or a misconfigured smearer
      // setup; both cases can not be handled and should be fatal.
      ACTS_ERROR("Could not find surface " << moduleGeoId
                                           << " for configured smearer");
      return ProcessCode::ABORT;
    }

    auto digitizerItr = m_digitizers.find(moduleGeoId);
    if (digitizerItr == m_digitizers.end()) {
      ACTS_DEBUG("No digitizer present for module " << moduleGeoId);
      continue;
    } else {
      ACTS_DEBUG("Digitizer found for module " << moduleGeoId);
    }

    // Run the digitizer. Iterate over the hits for this surface inside the
    // visitor so we do not need to lookup the variant object per-hit.
    std::visit(
        [&](const auto& digitizer) {
	  // Use TruthCluster objects to keep simhits
	  // association with clusters and/or cells since they can be
	  // merged downstream
          std::vector<::TruthCluster> moduleClusters;

	  // This map will hold initial smeared measurments associated
	  // to all hits on the module. When smearing is not
	  // requested, the map will hold an empty measurement
          std::unordered_map<hit_t, DigitizedParameters> measurementMap;

          for (auto h = moduleSimHits.begin(); h != moduleSimHits.end(); ++h) {
            const auto& simHit = *h;
            const auto simHitIdx = simHits.index_of(h);

            // Geometric part - 0, 1, 2 local parameters are possible
            if (not digitizer.geometric.indices.empty()) {
              ACTS_VERBOSE("Configured to geometric digitize "
                           << digitizer.geometric.indices.size()
                           << " parameters.");
              auto channels = channelizing(digitizer.geometric, simHit,
                                           *surfacePtr, ctx.geoContext, rng);
              if (channels.empty()) {
                ACTS_DEBUG(
                    "Geometric channelization did not work, skipping this hit.")
                continue;
              }
              ACTS_VERBOSE("Activated " << channels.size()
                                        << " channels for this hit.");

              moduleClusters.emplace_back(simHitIdx, channels);
            } else {
	      // In this case, geometric digitization is not performed
	      // so create a pseudo-cluster, i.e. associate the hit
	      // with an empty cell list
              moduleClusters.emplace_back(simHitIdx);
            }

            // Smearing part - (optionally) rest
            DigitizedParameters meas;
            if (not digitizer.smearing.indices.empty()) {
              ACTS_VERBOSE("Configured to smear "
                           << digitizer.smearing.indices.size()
                           << " parameters.");

              auto res =
                  digitizer.smearing(rng, simHit, *surfacePtr, ctx.geoContext);
              if (not res.ok()) {
                ACTS_DEBUG("Problem in hit smearing, skipping this hit.")
                continue;
              }
              const auto& [par, cov] = res.value();
              for (Eigen::Index ip = 0; ip < par.rows(); ++ip) {
                meas.indices.push_back(digitizer.smearing.indices[ip]);
                meas.values.push_back(par[ip]);
                meas.variances.push_back(cov(ip, ip));
              }
            }
	    // Unconditionaly associate a parameter set to each hit,
	    // even if empty
            measurementMap.insert({simHitIdx, meas});
          }

          if (m_cfg.mergeClusters)
            mergeClusters(moduleClusters, digitizer.geometric, measurementMap,
                          m_cfg.nSigmaMerge);

          for (::TruthCluster& cluster : moduleClusters) {
            DigitizedCluster dClus(std::move(cluster.simHits));
	    // N.B.: mergeClusters() assumes it runs before we create
	    // measurements from the cluster with localParameters()
            dClus.params =
                localParameters(digitizer.geometric, cluster.cells, rng);
            dClus.addSmearedParams(measurementMap, digitizer.smearing.indices);

            // Check on success - threshold could have eliminated all channels
            if (dClus.params.values.empty()) {
              ACTS_VERBOSE(
                  "Parameter digitization did not yield a measurement.")
              continue;
            }

            // The measurement container is unordered and the index under which
            // the measurement will be stored is known before adding it.
            Index measurementIdx = measurements.size();
            IndexSourceLink sourceLink(moduleGeoId, measurementIdx);

            // Add to output containers:
            // index map and source link container are geometry-ordered.
            // since the input is also geometry-ordered, new items can
            // be added at the end.
            sourceLinks.emplace_hint(sourceLinks.end(), std::move(sourceLink));
            measurements.emplace_back(
                createMeasurement(dClus.params, sourceLink));
            clusters.emplace_back(std::move(dClus.params.cluster));
            // this digitization does hit merging so there may be more than one
            // mapping entry for each digitized hit.
            for (auto idx : dClus.simHits) {
              measurementParticlesMap.emplace_hint(
                  measurementParticlesMap.end(), measurementIdx,
                  simHits.nth(idx)->particleId());
              measurementSimHitsMap.emplace_hint(measurementSimHitsMap.end(),
                                                 measurementIdx, idx);
            }
          }
        },
        *digitizerItr);
  }

  ctx.eventStore.add(m_cfg.outputSourceLinks, std::move(sourceLinks));
  ctx.eventStore.add(m_cfg.outputMeasurements, std::move(measurements));
  ctx.eventStore.add(m_cfg.outputClusters, std::move(clusters));
  ctx.eventStore.add(m_cfg.outputMeasurementParticlesMap,
                     std::move(measurementParticlesMap));
  ctx.eventStore.add(m_cfg.outputMeasurementSimHitsMap,
                     std::move(measurementSimHitsMap));
  return ProcessCode::SUCCESS;
}

std::vector<ActsFatras::Channelizer::ChannelSegment>
ActsExamples::DigitizationAlgorithm::channelizing(
    const GeometricConfig& geoCfg, const SimHit& hit,
    const Acts::Surface& surface, const Acts::GeometryContext& gctx,
    RandomEngine& rng) const {
  Acts::Vector3 driftDir = geoCfg.drift(hit.position(), rng);

  auto driftedSegment =
      m_surfaceDrift.toReadout(gctx, surface, geoCfg.thickness, hit.position(),
                               hit.unitDirection(), driftDir);
  auto maskedSegmentRes = m_surfaceMask.apply(surface, driftedSegment);
  if (maskedSegmentRes.ok()) {
    auto maskedSegment = maskedSegmentRes.value();
    // Now Channelize
    return m_channelizer.segments(gctx, surface, geoCfg.segmentation,
                                  maskedSegment);
  }
  return {};
}

ActsExamples::DigitizedParameters
ActsExamples::DigitizationAlgorithm::localParameters(
    const GeometricConfig& geoCfg,
    const std::vector<ActsFatras::Channelizer::ChannelSegment>& channels,
    RandomEngine& rng) const {
  DigitizedParameters dParameters;

  const auto& binningData = geoCfg.segmentation.binningData();

  Acts::ActsScalar totalWeight = 0.;
  Acts::Vector2 m(0., 0.);
  size_t b0min = SIZE_MAX;
  size_t b0max = 0;
  size_t b1min = SIZE_MAX;
  size_t b1max = 0;
  // Combine the channels
  for (const auto& ch : channels) {
    auto bin = ch.bin;
    Acts::ActsScalar charge =
        geoCfg.digital ? 1. : geoCfg.charge(ch.activation, ch.activation, rng);
    if (geoCfg.digital or charge > geoCfg.threshold) {
      totalWeight += charge;
      size_t b0 = bin[0];
      size_t b1 = bin[1];
      m += Acts::Vector2(charge * binningData[0].center(b0),
                         charge * binningData[1].center(b1));
      b0min = std::min(b0min, b0);
      b0max = std::max(b0max, b0);
      b1min = std::min(b1min, b1);
      b1max = std::max(b1max, b1);
      // Create a copy of the channel, as activation may change
      auto chdig = ch;
      chdig.bin = ch.bin;
      chdig.activation = charge;
      dParameters.cluster.channels.push_back(chdig);
    }
  }
  if (totalWeight > 0.) {
    m *= 1. / totalWeight;
    dParameters.indices = geoCfg.indices;
    for (auto idx : dParameters.indices) {
      dParameters.values.push_back(m[idx]);
    }
    size_t size0 = static_cast<size_t>(b0max - b0min + 1);
    size_t size1 = static_cast<size_t>(b1max - b1min + 1);
    auto variances = geoCfg.variances(size0, size1, rng);
    if (variances.size() == dParameters.indices.size()) {
      dParameters.variances = variances;
    } else {
      dParameters.variances =
          std::vector<Acts::ActsScalar>(dParameters.indices.size(), -1.);
    }

    dParameters.cluster.sizeLoc0 = size0;
    dParameters.cluster.sizeLoc1 = size1;
  }

  return dParameters;
}
