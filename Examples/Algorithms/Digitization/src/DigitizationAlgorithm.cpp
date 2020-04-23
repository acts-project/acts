// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Digitization/DigitizationAlgorithm.hpp"

#include <iostream>
#include <stdexcept>

#include "ACTFW/EventData/GeometryContainers.hpp"
#include "ACTFW/EventData/SimHit.hpp"
#include "ACTFW/EventData/SimParticle.hpp"
#include "ACTFW/EventData/SimVertex.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/DetectorElementBase.hpp"
#include "Acts/Geometry/GeometryID.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Plugins/Digitization/DigitizationModule.hpp"
#include "Acts/Plugins/Digitization/PlanarModuleCluster.hpp"
#include "Acts/Plugins/Digitization/PlanarModuleStepper.hpp"
#include "Acts/Plugins/Digitization/Segmentation.hpp"
#include "Acts/Plugins/Identification/IdentifiedDetectorElement.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"
#include "Acts/Utilities/Units.hpp"

FW::DigitizationAlgorithm::DigitizationAlgorithm(
    FW::DigitizationAlgorithm::Config cfg, Acts::Logging::Level lvl)
    : FW::BareAlgorithm("DigitizationAlgorithm", lvl), m_cfg(std::move(cfg)) {
  if (m_cfg.inputSimulatedHits.empty()) {
    throw std::invalid_argument("Missing input hits collection");
  }
  if (m_cfg.outputClusters.empty()) {
    throw std::invalid_argument("Missing output clusters collection");
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
    if (not dg.surface) {
      return;
    }
    // require an associated detector element
    dg.detectorElement = dynamic_cast<const Acts::IdentifiedDetectorElement*>(
        dg.surface->associatedDetectorElement());
    if (not dg.detectorElement) {
      return;
    }
    // require an associated digitization module
    dg.digitizer = dg.detectorElement->digitizationModule().get();
    if (not dg.digitizer) {
      return;
    }
    // record all valid surfaces
    this->m_digitizables.insert_or_assign(surface->geoID(), dg);
  });
}

FW::ProcessCode FW::DigitizationAlgorithm::execute(
    const AlgorithmContext& ctx) const {
  // Prepare the input and output collections
  const auto& hits =
      ctx.eventStore.get<SimHitContainer>(m_cfg.inputSimulatedHits);
  FW::GeometryIdMultimap<Acts::PlanarModuleCluster> clusters;

  for (auto&& [moduleGeoId, moduleHits] : groupByModule(hits)) {
    // can only digitize hits on digitizable surfaces
    const auto it = m_digitizables.find(moduleGeoId);
    if (it == m_digitizables.end()) {
      continue;
    }

    const auto& dg = it->second;
    // local intersection / direction
    const auto invTransfrom = dg.surface->transform(ctx.geoContext).inverse();

    // use iterators manually so we can retrieve the hit index in the container
    for (auto ih = moduleHits.begin(); ih != moduleHits.end(); ++ih) {
      const auto& hit = *ih;
      const auto idx = hits.index_of(ih);

      Acts::Vector2D localIntersect = (invTransfrom * hit.position()).head<2>();
      Acts::Vector3D localDirection =
          invTransfrom.linear() * hit.unitDirection();

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
      if (!dSteps.size()) {
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
      Acts::Vector2D localPosition(localX, localY);
      // @todo remove unneccesary conversion
      // size_t bin0 = binUtility.bin(localPosition, 0);
      // size_t bin1 = binUtility.bin(localPosition, 1);
      // size_t binSerialized = binUtility.serialize({{bin0, bin1, 0}});

      // the covariance is currently set to 0.
      Acts::ActsSymMatrixD<3> cov;
      cov << 0.05, 0., 0., 0., 0.05, 0., 0., 0.,
          900. * Acts::UnitConstants::ps * Acts::UnitConstants::ps;

      // create the planar cluster
      Acts::PlanarModuleCluster pCluster(
          dg.surface->getSharedPtr(), Identifier(identifier_type(idx), {idx}),
          std::move(cov), localX, localY, hit.time(), std::move(usedCells));

      // insert into the cluster container. since the input data is already
      // sorted by geoId, we should always be able to add at the end.
      clusters.emplace_hint(clusters.end(), hit.geometryId(),
                            std::move(pCluster));
    }
  }

  ACTS_DEBUG("digitized " << hits.size() << " hits into " << clusters.size()
                          << " clusters");

  // write the clusters to the EventStore
  ctx.eventStore.add(m_cfg.outputClusters, std::move(clusters));
  return FW::ProcessCode::SUCCESS;
}
