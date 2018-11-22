// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// PassiveLayerBuilder.cpp, Acts project
///////////////////////////////////////////////////////////////////

#include "Acts/Tools/PassiveLayerBuilder.hpp"
#include "Acts/Layers/CylinderLayer.hpp"
#include "Acts/Layers/DiscLayer.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/MaterialProperties.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Utilities/Definitions.hpp"

Acts::PassiveLayerBuilder::PassiveLayerBuilder(
    const PassiveLayerBuilder::Config& plConfig,
    std::unique_ptr<const Logger>      logger)
  : m_cfg(), m_logger(std::move(logger)), m_nLayers(), m_cLayers(), m_pLayers()
{
  setConfiguration(plConfig);
}

void
Acts::PassiveLayerBuilder::setConfiguration(
    const PassiveLayerBuilder::Config& plConfig)
{
  //!< @todo add configuration check
  m_cfg = plConfig;
  constructLayers();
}

void
Acts::PassiveLayerBuilder::setLogger(std::unique_ptr<const Logger> newLogger)
{
  m_logger = std::move(newLogger);
}

void
Acts::PassiveLayerBuilder::constructLayers()
{
  // DEBUG: Flush layers in case the class was already initialized before
  m_nLayers.clear();
  m_cLayers.clear();
  m_pLayers.clear();

  // the central layers
  size_t numcLayers = m_cfg.centralLayerRadii.size();
  if (numcLayers != 0u) {
    ACTS_DEBUG("Configured to build " << numcLayers
                                      << " passive central layers.");
    m_cLayers.reserve(numcLayers);
    // loop through
    for (size_t icl = 0; icl < numcLayers; ++icl) {
      // some screen output
      ACTS_VERBOSE("- build layer " << icl << " with radius = "
                                    << m_cfg.centralLayerRadii.at(icl)
                                    << " and halfZ = "
                                    << m_cfg.centralLayerHalflengthZ.at(icl));
      // create the layer and push it back
      auto cBounds = std::make_shared<const CylinderBounds>(
          m_cfg.centralLayerRadii[icl], m_cfg.centralLayerHalflengthZ.at(icl));
      // create the layer
      MutableLayerPtr cLayer = CylinderLayer::create(
          nullptr, cBounds, nullptr, m_cfg.centralLayerThickness.at(icl));
      // assign the material to the layer surface
      std::shared_ptr<const SurfaceMaterial> material = nullptr;
      // create the material from jobOptions
      if (!m_cfg.centralLayerMaterial.empty()) {
        // create homogeneous material
        material = m_cfg.centralLayerMaterial.at(icl);
        // sign it to the surface
        cLayer->surfaceRepresentation().setAssociatedMaterial(material);
      }
      // push it into the layer vector
      m_cLayers.push_back(cLayer);
    }
  }

  // pos/neg layers
  size_t numpnLayers = m_cfg.posnegLayerPositionZ.size();
  if (numpnLayers != 0u) {
    ACTS_DEBUG("Configured to build 2 * "
               << numpnLayers
               << " passive positive/negative side layers.");
    m_pLayers.reserve(numpnLayers);
    m_nLayers.reserve(numpnLayers);
    // loop through
    for (size_t ipnl = 0; ipnl < numpnLayers; ++ipnl) {
      // some screen output
      ACTS_VERBOSE("- build layers " << (2 * ipnl) << " and " << (2 * ipnl) + 1
                                     << " at +/- z = "
                                     << m_cfg.posnegLayerPositionZ.at(ipnl)
                                     << " and rMin/rMax = "
                                     << m_cfg.posnegLayerRmin.at(ipnl)
                                     << " / "
                                     << m_cfg.posnegLayerRmax.at(ipnl));
      // create the share disc bounds
      std::shared_ptr<const DiscBounds> dBounds
          = std::make_shared<const RadialBounds>(
              m_cfg.posnegLayerRmin.at(ipnl), m_cfg.posnegLayerRmax.at(ipnl));
      // create the layer transforms
      const Transform3D* nTransform = new Transform3D(
          Translation3D(0., 0., -m_cfg.posnegLayerPositionZ.at(ipnl)));
      const Transform3D* pTransform = new Transform3D(
          Translation3D(0., 0., m_cfg.posnegLayerPositionZ.at(ipnl)));
      // create the layers
      MutableLayerPtr nLayer
          = DiscLayer::create(std::shared_ptr<const Transform3D>(nTransform),
                              dBounds,
                              nullptr,
                              m_cfg.posnegLayerThickness.at(ipnl));
      MutableLayerPtr pLayer
          = DiscLayer::create(std::shared_ptr<const Transform3D>(pTransform),
                              dBounds,
                              nullptr,
                              m_cfg.posnegLayerThickness.at(ipnl));
      // assign the material to the layer surface
      std::shared_ptr<const SurfaceMaterial> material = nullptr;
      // create the material from jobOptions
      if (!m_cfg.posnegLayerMaterial.empty()) {
        // create homogeneous material
        material = m_cfg.posnegLayerMaterial.at(ipnl);
        // sign it to the surface
        nLayer->surfaceRepresentation().setAssociatedMaterial(material);
        pLayer->surfaceRepresentation().setAssociatedMaterial(material);
      }
      // push it into the layer vector
      m_nLayers.push_back(nLayer);
      m_pLayers.push_back(pLayer);
    }
  }
}
