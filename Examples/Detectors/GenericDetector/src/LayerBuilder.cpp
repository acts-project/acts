// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/GenericDetector/LayerBuilder.hpp"

#include <span>

namespace ActsExamples::Generic {

namespace {
std::vector<std::shared_ptr<const Acts::Surface>> makeConstSurfaces(
    std::span<const std::shared_ptr<Acts::Surface>> surfaces) {
  std::vector<std::shared_ptr<const Acts::Surface>> constSurfaces;
  constSurfaces.reserve(surfaces.size());
  for (const auto& surface : surfaces) {
    constSurfaces.push_back(surface);
  }
  return constSurfaces;
}
}  // namespace

const Acts::LayerVector LayerBuilder::centralLayers(
    const Acts::GeometryContext& gctx) const {
  // create the vector
  Acts::LayerVector cLayers;
  cLayers.reserve(m_cfg.centralProtoLayers.size());
  // the layer counter
  std::size_t icl = 0;
  for (auto& cpl : m_cfg.centralProtoLayers) {
    // create the layer actually
    Acts::ProtoLayer constProtoLayer{cpl.protoLayer};
    Acts::MutableLayerPtr cLayer = m_cfg.layerCreator->cylinderLayer(
        gctx, makeConstSurfaces(cpl.surfaces), cpl.bins0, cpl.bins1,
        constProtoLayer);

    // the layer is built let's see if it needs material
    if (!m_cfg.centralLayerMaterial.empty()) {
      std::shared_ptr<const Acts::ISurfaceMaterial> layerMaterialPtr =
          m_cfg.centralLayerMaterial.at(icl);
      // central material
      if (m_cfg.centralLayerMaterialConcentration.at(icl) == 0.) {
        // the layer surface is the material surface
        cLayer->surfaceRepresentation().assignSurfaceMaterial(layerMaterialPtr);
        ACTS_VERBOSE("- and material at central layer surface.");
      } else {
        // approach surface material
        // get the approach descriptor - at this stage we know that the
        // approachDescriptor exists
        auto approachSurfaces =
            cLayer->approachDescriptor()->containedSurfaces();
        if (m_cfg.centralLayerMaterialConcentration.at(icl) > 0) {
          auto mutableOuterSurface =
              const_cast<Acts::Surface*>(approachSurfaces.at(1));
          mutableOuterSurface->assignSurfaceMaterial(layerMaterialPtr);
          ACTS_VERBOSE("- and material at outer approach surface");
        } else {
          auto mutableInnerSurface =
              const_cast<Acts::Surface*>(approachSurfaces.at(0));
          mutableInnerSurface->assignSurfaceMaterial(layerMaterialPtr);
          ACTS_VERBOSE("- and material at inner approach surface");
        }
      }
    }
    // push it into the layer vector
    cLayers.push_back(cLayer);
    ++icl;
  }
  return cLayers;
}

const Acts::LayerVector LayerBuilder::negativeLayers(
    const Acts::GeometryContext& gctx) const {
  return constructEndcapLayers(gctx, -1);
}

const Acts::LayerVector LayerBuilder::positiveLayers(
    const Acts::GeometryContext& gctx) const {
  return constructEndcapLayers(gctx, 1);
}

LayerBuilder::LayerBuilder(const Config& cfg,
                           std::unique_ptr<const Acts::Logger> log)
    : Acts::ILayerBuilder(), m_cfg(cfg), m_logger(std::move(log)) {}

const Acts::LayerVector LayerBuilder::constructEndcapLayers(
    const Acts::GeometryContext& gctx, int side) const {
  // The from negative or positive proto layers
  const auto& protoLayers =
      (side < 0) ? m_cfg.negativeProtoLayers : m_cfg.positiveProtoLayers;

  // create the vector
  Acts::LayerVector eLayers;
  eLayers.reserve(protoLayers.size());

  // the layer counter
  std::size_t ipnl = 0;
  // loop over the proto layers and create the actual layers
  for (auto& ple : protoLayers) {
    /// the layer is created
    Acts::ProtoLayer constProtoLayer{ple.protoLayer};
    Acts::MutableLayerPtr eLayer =
        m_cfg.layerCreator->discLayer(gctx, makeConstSurfaces(ple.surfaces),
                                      ple.bins0, ple.bins1, constProtoLayer);

    // the layer is built let's see if it needs material
    if (!m_cfg.posnegLayerMaterial.empty()) {
      std::shared_ptr<const Acts::ISurfaceMaterial> layerMaterialPtr =
          m_cfg.posnegLayerMaterial[ipnl];
      // central material
      if (m_cfg.posnegLayerMaterialConcentration.at(ipnl) == 0.) {
        // assign the surface material - the layer surface is the material
        // surface
        eLayer->surfaceRepresentation().assignSurfaceMaterial(layerMaterialPtr);
        ACTS_VERBOSE("- and material at central layer surface.");
      } else {
        // approach surface material
        // get the approach descriptor - at this stage we know that the
        // approachDescriptor exists
        auto approachSurfaces =
            eLayer->approachDescriptor()->containedSurfaces();
        if (m_cfg.posnegLayerMaterialConcentration.at(ipnl) > 0.) {
          int sf = side < 0 ? 0 : 1;
          auto mutableInnerSurface =
              const_cast<Acts::Surface*>(approachSurfaces.at(sf));
          mutableInnerSurface->assignSurfaceMaterial(layerMaterialPtr);
          ACTS_VERBOSE("- and material at outer approach surfaces.");
        } else {
          int sf = side < 0 ? 1 : 0;
          auto mutableOuterSurface =
              const_cast<Acts::Surface*>(approachSurfaces.at(sf));
          mutableOuterSurface->assignSurfaceMaterial(layerMaterialPtr);
          ACTS_VERBOSE("- and material at inner approach surfaces.");
        }
      }
    }
    // push it into the layer vector
    eLayers.push_back(eLayer);
    ++ipnl;
  }
  return eLayers;
}

}  // namespace ActsExamples::Generic
