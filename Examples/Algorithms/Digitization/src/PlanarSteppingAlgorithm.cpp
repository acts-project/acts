// This file is part of the Acts project.
//
// Copyright (C) 2017-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Digitization/PlanarSteppingAlgorithm.hpp"

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Digitization/DigitizationModule.hpp"
#include "Acts/Digitization/PlanarModuleCluster.hpp"
#include "Acts/Digitization/PlanarModuleStepper.hpp"
#include "Acts/Digitization/Segmentation.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/DetectorElementBase.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Plugins/Identification/IdentifiedDetectorElement.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "ActsExamples/EventData/GeometryContainers.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <list>
#include <stdexcept>

ActsExamples::PlanarSteppingAlgorithm::PlanarSteppingAlgorithm(
    ActsExamples::PlanarSteppingAlgorithm::Config config,
    Acts::Logging::Level level)
    : ActsExamples::BareAlgorithm("PlanarSteppingAlgorithm", level),
      m_cfg(std::move(config)) {
  if (m_cfg.inputSimHits.empty()) {
    throw std::invalid_argument("Missing input hits collection");
  }
  if (m_cfg.outputClusters.empty()) {
    throw std::invalid_argument("Missing output clusters collection");
  }
  if (m_cfg.outputSourceLinks.empty()) {
    throw std::invalid_argument("Missing source links output collection");
  }
  if (m_cfg.outputDigiSourceLinks.empty()) {
    throw std::invalid_argument(
        "Missing digitization source links output collection");
  }
  if (m_cfg.outputMeasurements.empty()) {
    throw std::invalid_argument("Missing measurements output collection");
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
  if (!m_cfg.planarModuleStepper) {
    throw std::invalid_argument("Missing planar module stepper");
  }
  if (!m_cfg.randomNumbers) {
    throw std::invalid_argument("Missing random numbers tool");
  }
  // fill the digitizables map to allow lookup by geometry id only
  m_cfg.trackingGeometry->visitSurfaces([this](const Acts::Surface* surface) {
    Digitizable dg;
    // require a valid surface
    dg.surface = surface;
    if (dg.surface == nullptr) {
      return;
    }
    // require an associated detector element
    dg.detectorElement = dynamic_cast<const Acts::IdentifiedDetectorElement*>(
        dg.surface->associatedDetectorElement());
    if (dg.detectorElement == nullptr) {
      return;
    }
    // require an associated digitization module
    dg.digitizer = dg.detectorElement->digitizationModule().get();
    if (dg.digitizer == nullptr) {
      return;
    }
    // record all valid surfaces
    this->m_digitizables.insert_or_assign(surface->geometryId(), dg);
  });
}

ActsExamples::ProcessCode ActsExamples::PlanarSteppingAlgorithm::execute(
    const AlgorithmContext& ctx) const {
  using ClusterContainer =
      ActsExamples::GeometryIdMultimap<Acts::PlanarModuleCluster>;

  // retrieve input
  const auto& simHits = ctx.eventStore.get<SimHitContainer>(m_cfg.inputSimHits);

  // prepare output containers
  ClusterContainer clusters;

  // These are lists because we need stable addresses
  std::list<IndexSourceLink> sourceLinkStorage;
  std::list<Acts::DigitizationSourceLink> digiSourceLinks;

  GeometryIdMultiset<std::reference_wrapper<IndexSourceLink>> sourceLinks;
  MeasurementContainer measurements;
  IndexMultimap<ActsFatras::Barcode> hitParticlesMap;
  IndexMultimap<Index> hitSimHitsMap;
  clusters.reserve(simHits.size());
  measurements.reserve(simHits.size());
  hitParticlesMap.reserve(simHits.size());
  hitSimHitsMap.reserve(simHits.size());

  for (auto&& [moduleGeoId, moduleSimHits] : groupByModule(simHits)) {
    // can only digitize hits on digitizable surfaces
    const auto it = m_digitizables.find(moduleGeoId);
    if (it == m_digitizables.end()) {
      continue;
    }

    const auto& dg = it->second;
    // local intersection / direction
    const auto invTransfrom = dg.surface->transform(ctx.geoContext).inverse();

    // use iterators manually so we can retrieve the hit index in the container
    for (auto ih = moduleSimHits.begin(); ih != moduleSimHits.end(); ++ih) {
      const auto& simHit = *ih;
      const auto simHitIdx = simHits.index_of(ih);

      Acts::Vector2 localIntersect =
          (invTransfrom * simHit.position()).head<2>();
      Acts::Vector3 localDirection =
          invTransfrom.linear() * simHit.unitDirection();

      // compute digitization steps
      const auto thickness = dg.detectorElement->thickness();
      const auto lorentzAngle = dg.digitizer->lorentzAngle();
      auto lorentzShift = thickness * std::tan(lorentzAngle);
      lorentzShift *= -(dg.digitizer->readoutDirection());
      // now calculate the steps through the silicon
      std::vector<Acts::DigitizationStep> dSteps =
          m_cfg.planarModuleStepper->cellSteps(ctx.geoContext, *dg.digitizer,
                                               localIntersect, localDirection);
      // everything under threshold or edge effects
      if (dSteps.empty()) {
        ACTS_VERBOSE("No steps returned from stepper.");
        continue;
      }

      // lets create a cluster - centroid method
      double localX = 0.;
      double localY = 0.;
      double totalPath = 0.;
      // the cells to be used
      std::vector<Acts::DigitizationCell> usedCells;
      usedCells.reserve(dSteps.size());
      // loop over the steps
      for (auto dStep : dSteps) {
        // @todo implement smearing
        localX += dStep.stepLength * dStep.stepCellCenter.x();
        localY += dStep.stepLength * dStep.stepCellCenter.y();
        totalPath += dStep.stepLength;
        usedCells.push_back(Acts::DigitizationCell(dStep.stepCell.channel0,
                                                   dStep.stepCell.channel1,
                                                   dStep.stepLength));
      }
      // divide by the total path
      localX /= totalPath;
      localX += lorentzShift;
      localY /= totalPath;

      // get the segmentation & find the corresponding cell id
      const Acts::Segmentation& segmentation = dg.digitizer->segmentation();
      auto binUtility = segmentation.binUtility();
      Acts::Vector2 localPosition(localX, localY);
      // @todo remove unneccesary conversion
      // size_t bin0 = binUtility.bin(localPosition, 0);
      // size_t bin1 = binUtility.bin(localPosition, 1);
      // size_t binSerialized = binUtility.serialize({{bin0, bin1, 0}});

      // the covariance is currently set to some arbitrary value.
      Acts::SymMatrix3 cov;
      cov << 0.05, 0., 0., 0., 0.05, 0., 0., 0.,
          900. * Acts::UnitConstants::ps * Acts::UnitConstants::ps;
      Acts::Vector3 par(localX, localY, simHit.time());

      // create the planar cluster
      digiSourceLinks.emplace_back(moduleGeoId,
                                   std::vector<std::size_t>{simHitIdx});
      Acts::DigitizationSourceLink& digiSourceLink = digiSourceLinks.back();

      Acts::PlanarModuleCluster cluster(
          dg.surface->getSharedPtr(), digiSourceLink, std::move(cov), localX,
          localY, simHit.time(), std::move(usedCells));

      // the measurement container is unordered and the index under which
      // the measurement will be stored is known before adding it.
      Index hitIdx = measurements.size();
      sourceLinkStorage.emplace_back(moduleGeoId, hitIdx);
      IndexSourceLink& sourceLink = sourceLinkStorage.back();

      sourceLinks.insert(sourceLinks.end(), sourceLink);

      auto meas = Acts::makeMeasurement(sourceLink, par, cov, Acts::eBoundLoc0,
                                        Acts::eBoundLoc1, Acts::eBoundTime);

      // add to output containers. since the input is already geometry-order,
      // new elements in geometry containers can just be appended at the end.
      clusters.emplace_hint(clusters.end(), moduleGeoId, std::move(cluster));
      measurements.emplace_back(std::move(meas));
      // no hit merging -> only one mapping per digitized hit.
      hitParticlesMap.emplace_hint(hitParticlesMap.end(), hitIdx,
                                   simHit.particleId());
      hitSimHitsMap.emplace_hint(hitSimHitsMap.end(), hitIdx, simHitIdx);
    }
  }

  ACTS_DEBUG("digitized " << simHits.size() << " hits into " << clusters.size()
                          << " clusters");

  ctx.eventStore.add(m_cfg.outputClusters, std::move(clusters));
  ctx.eventStore.add(m_cfg.outputSourceLinks, std::move(sourceLinks));
  ctx.eventStore.add(m_cfg.outputSourceLinks + "__storage",
                     std::move(sourceLinkStorage));
  ctx.eventStore.add(m_cfg.outputDigiSourceLinks, std::move(digiSourceLinks));
  ctx.eventStore.add(m_cfg.outputMeasurements, std::move(measurements));
  ctx.eventStore.add(m_cfg.outputMeasurementParticlesMap,
                     std::move(hitParticlesMap));
  ctx.eventStore.add(m_cfg.outputMeasurementSimHitsMap,
                     std::move(hitSimHitsMap));
  return ActsExamples::ProcessCode::SUCCESS;
}
