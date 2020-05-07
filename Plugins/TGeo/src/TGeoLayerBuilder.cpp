// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/TGeo/TGeoLayerBuilder.hpp"
#include <stdio.h>
#include "Acts/Geometry/ProtoLayer.hpp"
#include "Acts/Plugins/TGeo/TGeoDetectorElement.hpp"
#include "TGeoManager.h"
#include "TGeoMatrix.h"

Acts::TGeoLayerBuilder::TGeoLayerBuilder(
    const Acts::TGeoLayerBuilder::Config& config,
    std::unique_ptr<const Logger> logger)
    : m_cfg(), m_logger(std::move(logger)) {
  setConfiguration(config);
}

Acts::TGeoLayerBuilder::~TGeoLayerBuilder() = default;

void Acts::TGeoLayerBuilder::setConfiguration(
    const Acts::TGeoLayerBuilder::Config& config) {
  m_cfg = config;
}

void Acts::TGeoLayerBuilder::setLogger(
    std::unique_ptr<const Logger> newLogger) {
  m_logger = std::move(newLogger);
}

const Acts::LayerVector Acts::TGeoLayerBuilder::negativeLayers(
    const GeometryContext& gctx) const {
  // @todo Remove this hack once the m_elementStore mess is sorted out
  auto mutableThis = const_cast<TGeoLayerBuilder*>(this);
  LayerVector nVector;
  mutableThis->buildLayers(gctx, nVector, -1);
  return nVector;
}

const Acts::LayerVector Acts::TGeoLayerBuilder::centralLayers(
    const GeometryContext& gctx) const {
  // @todo Remove this hack once the m_elementStore mess is sorted out
  auto mutableThis = const_cast<TGeoLayerBuilder*>(this);
  LayerVector cVector;
  mutableThis->buildLayers(gctx, cVector, 0);
  return cVector;
}

const Acts::LayerVector Acts::TGeoLayerBuilder::positiveLayers(
    const GeometryContext& gctx) const {
  // @todo Remove this hack once the m_elementStore mess is sorted out
  auto mutableThis = const_cast<TGeoLayerBuilder*>(this);
  LayerVector pVector;
  mutableThis->buildLayers(gctx, pVector, 1);
  return pVector;
}

void Acts::TGeoLayerBuilder::buildLayers(const GeometryContext& gctx,
                                         LayerVector& layers, int type) {
  // Bail out if you have no gGeoManager
  if (gGeoManager == nullptr) {
    ACTS_WARNING("No gGeoManager found - bailing out.");
    return;
  }

  // Prepare which ones to build
  std::vector<LayerConfig> layerConfigs = m_cfg.layerConfigurations[type + 1];
  std::string layerType = m_layerTypes[type + 1];

  // Appropriate screen output
  std::string addonOutput = m_cfg.layerSplitToleranceR[type + 1] > 0.
                                ? std::string(", splitting in r")
                                : std::string("");
  addonOutput += m_cfg.layerSplitToleranceZ[type + 1] > 0.
                     ? std::string(", splitting in z")
                     : std::string("");
  addonOutput += std::string(".");

  // Screen output of the configuration
  ACTS_DEBUG(layerType << " layers : found " << layerConfigs.size()
                       << " configuration(s)" + addonOutput);
  for (auto layerCfg : layerConfigs) {
    // Prepare the layer surfaces
    using LayerSurfaceVector = std::vector<std::shared_ptr<const Surface>>;
    LayerSurfaceVector layerSurfaces;

    ACTS_DEBUG("- layer configuration found for layer "
               << layerCfg.layerName << " with sensor " << layerCfg.sensorName);
    ACTS_DEBUG("- layers radially bound to rmin/rmax = "
               << layerCfg.parseRangeR.first << "/"
               << layerCfg.parseRangeR.second);
    // Step down from the top volume each time to collect the logical tree
    TGeoVolume* tvolume = gGeoManager->GetTopVolume();
    if (tvolume != nullptr) {
      // Recursively step down
      resolveSensitive(gctx, layerSurfaces, tvolume, nullptr, TGeoIdentity(),
                       layerCfg, type);
      // screen output
      ACTS_DEBUG(
          "- number of senstive sensors found : " << layerSurfaces.size());

      // Helper function to fill the layer
      auto fillLayer = [&](const LayerSurfaceVector lSurfaces,
                           const LayerConfig& lCfg) -> void {
        // Create the layer  - either way as cylinder or disk
        if (type == 0) {
          ProtoLayer pl(gctx, lSurfaces);
          pl.envelope[Acts::binR] = {lCfg.envelope.first, lCfg.envelope.second};
          pl.envelope[Acts::binZ] = {lCfg.envelope.second,
                                     lCfg.envelope.second};
          layers.push_back(m_cfg.layerCreator->cylinderLayer(
              gctx, lSurfaces, lCfg.binsLoc0, lCfg.binsLoc1, pl));
        } else {
          ProtoLayer pl(gctx, lSurfaces);
          pl.envelope[Acts::binR] = {lCfg.envelope.first, lCfg.envelope.second};
          pl.envelope[Acts::binZ] = {lCfg.envelope.second,
                                     lCfg.envelope.second};
          layers.push_back(m_cfg.layerCreator->discLayer(
              gctx, lSurfaces, lCfg.binsLoc0, lCfg.binsLoc1, pl));
        }
      };

      // There is no split to be attempted
      if (layerCfg.splitParametersR.empty() and
          layerCfg.splitParametersZ.empty()) {
        // No splitting to be done, fill and return
        fillLayer(layerSurfaces, layerCfg);
        return;
      }

      std::vector<LayerSurfaceVector> splitLayerSurfaces = {layerSurfaces};

      // Helper method : perform the actual split
      auto splitSurfaces =
          [&](std::string splitValue, BinningValue bValue,
              const std::vector<LayerSurfaceVector>& preSplitSurfaces,
              double splitTolerance,
              std::pair<double, double> splitRange = {0., 0.},
              std::vector<double> splitParameters = {})
          -> std::vector<LayerSurfaceVector> {
        ACTS_DEBUG("- split attempt in " << splitValue);
        ACTS_DEBUG("- split layers seperated by more than " << splitTolerance);
        // Re-evaluate
        bool reevaluate = splitParameters.empty();
        if (reevaluate) {
          ACTS_DEBUG("- split parameters to be re-evaluated");
        }

        // The vector of surfaces after splitting
        std::vector<LayerSurfaceVector> postSplitSurfaces;
        // Go through and split them accordingly
        for (const auto& surfaceSet : preSplitSurfaces) {
          ACTS_DEBUG("- split surface set with " << surfaceSet.size()
                                                 << " surfaces.");
          // The Split parameters are empty, parse again
          if (reevaluate) {
            // Loop over sub set for new splitting range and parameters
            for (const auto& surface : surfaceSet) {
              // Get the surface parameter
              double surfacePar = surface->binningPositionValue(gctx, bValue);
              registerSplit(splitParameters, surfacePar, splitTolerance,
                            splitRange);
            }
          }
          // Output the split range
          ACTS_DEBUG("- split range is = " << splitRange.first << ", "
                                           << splitRange.second);
          ACTS_DEBUG("- number of proposed splits is "
                     << splitParameters.size());
          // Allocate expected sub vector
          std::vector<LayerSurfaceVector> setSplitSurfaces{
              splitParameters.size(), LayerSurfaceVector{}};
          // Filling loop (2nd loop if split in case of re-evaluation)
          for (const auto& surface : surfaceSet) {
            // Get the surface parameter
            double surfacePar = surface->binningPositionValue(gctx, bValue);
            unsigned isplit = 0;
            for (const auto& splitPar : splitParameters) {
              if (std::abs(splitPar - surfacePar) < splitTolerance) {
                setSplitSurfaces[isplit].push_back(surface);
              }
              ++isplit;
            }
          }
          // Insert the split set into the post split set
          postSplitSurfaces.insert(postSplitSurfaces.end(),
                                   setSplitSurfaces.begin(),
                                   setSplitSurfaces.end());
          // reset the split parameters
          if (reevaluate) {
            splitParameters.clear();
          }
        }
        unsigned int iss = 0;
        ACTS_VERBOSE("  - splitting yielded " << postSplitSurfaces.size()
                                              << " surface sets:");
        for (const auto& sset : postSplitSurfaces) {
          ACTS_VERBOSE("  - set " << iss++ << " has " << sset.size()
                                  << " surfaces.");
        }

        // Return them to the callers
        return postSplitSurfaces;
      };

      // Split in R first if split tolerance is set
      if (m_cfg.layerSplitToleranceR[type + 1] > 0.) {
        //  Split the surfaces in R
        splitLayerSurfaces = splitSurfaces(
            "r", binR, splitLayerSurfaces, m_cfg.layerSplitToleranceR[type + 1],
            layerCfg.splitRangeR, layerCfg.splitParametersR);
        // This invalidates the Z parameters and range
        layerCfg.splitParametersZ.clear();
        layerCfg.splitRangeZ = {std::numeric_limits<double>::max(),
                                -std::numeric_limits<double>::max()};
      }

      // Split in Z then if configured to do so
      if (m_cfg.layerSplitToleranceZ[type + 1] > 0.) {
        //  Split the surfaces in Z
        splitLayerSurfaces = splitSurfaces(
            "z", binZ, splitLayerSurfaces, m_cfg.layerSplitToleranceZ[type + 1],
            layerCfg.splitRangeZ, layerCfg.splitParametersZ);
      }

      // Now go through and fill, @todo adapt layer configurations
      unsigned int il = 0;
      for (const auto& slSurfaces : splitLayerSurfaces) {
        ACTS_VERBOSE("  - layer " << il++ << " has " << slSurfaces.size()
                                  << " surfaces.");
        fillLayer(slSurfaces, layerCfg);
      }
    }
  }
}

void Acts::TGeoLayerBuilder::resolveSensitive(
    const GeometryContext& gctx,
    std::vector<std::shared_ptr<const Acts::Surface>>& layerSurfaces,
    TGeoVolume* tgVolume, TGeoNode* tgNode, const TGeoMatrix& tgTransform,
    LayerConfig& layerConfig, int type, bool correctBranch,
    const std::string& offset) {
  if (tgVolume != nullptr) {
    std::string volumeName = tgVolume->GetName();
    /// Some screen output indicating that the volume was found
    if (m_cfg.nodeSearchDebug) {
      ACTS_VERBOSE(offset << "[o] Volume : " << volumeName
                          << " - checking for volume name "
                          << layerConfig.layerName);
    }
    // Once in the current branch stepping down means staying inside the branch
    bool correctVolume = correctBranch;
    if (!correctVolume &&
        (volumeName.find(layerConfig.layerName) != std::string::npos ||
         match(layerConfig.layerName.c_str(), volumeName.c_str()))) {
      correctVolume = true;
      if (m_cfg.nodeSearchDebug) {
        ACTS_VERBOSE(offset << "    triggered current branch!");
      }
    }
    // Loop over the daughters and collect them
    auto daugthers = tgVolume->GetNodes();
    // Screen output
    if (m_cfg.nodeSearchDebug) {
      ACTS_VERBOSE(offset << "has " << tgVolume->GetNdaughters()
                          << " daughters.");
    }
    // A daughter iterator
    TIter iObj(daugthers);
    // While loop over the objects for the recursive parsing
    while (TObject* obj = iObj()) {
      // dynamic_cast to a node
      TGeoNode* node = dynamic_cast<TGeoNode*>(obj);
      if (node != nullptr) {
        resolveSensitive(gctx, layerSurfaces, nullptr, node, tgTransform,
                         layerConfig, type, correctVolume, offset + "  ");
      }
    }
  } else if (m_cfg.nodeSearchDebug) {
    ACTS_VERBOSE("No volume present.");
  }

  /// If you have a node, get the volume and step down further
  if (tgNode != nullptr) {
    // Get the matrix of the current node for positioning
    const TGeoMatrix* tgMatrix = tgNode->GetMatrix();

    // Build the matrix
    TGeoHMatrix parseTransform =
        TGeoCombiTrans(tgTransform) * TGeoCombiTrans(*tgMatrix);

    // The translation of the node for parsing
    const Double_t* translation = parseTransform.GetTranslation();

    double x = m_cfg.unit * translation[0];
    double y = m_cfg.unit * translation[1];
    double z = m_cfg.unit * translation[2];
    double r = std::sqrt(x * x + y * y);

    // The name of the nodefor cross checking
    std::string tNodeName = tgNode->GetName();
    if (m_cfg.nodeSearchDebug) {
      ACTS_VERBOSE(offset << "[>] Node : " << tNodeName
                          << " - checking for sensor name "
                          << layerConfig.sensorName);
    }
    // Find out the branch hit, ingle layer depth supported by sensor==layer
    bool branchHit =
        correctBranch || (layerConfig.sensorName == layerConfig.layerName);
    if (branchHit &&
        (tNodeName.find(layerConfig.sensorName) != std::string::npos ||
         match(layerConfig.sensorName.c_str(), tNodeName.c_str()))) {
      if (m_cfg.nodeSearchDebug) {
        ACTS_VERBOSE(offset << "Sensor name '" << layerConfig.sensorName
                            << "' found in branch '" << layerConfig.layerName
                            << "'.");
      }

      // Create the detector element
      // - check on the type for the side
      // - check for the parsing volume
      bool insideParseRange = r >= layerConfig.parseRangeR.first and
                              r <= layerConfig.parseRangeR.second;

      if (insideParseRange and ((type == 0) || type * z > 0.)) {
        //  Senstive volume found, collect it
        if (m_cfg.nodeSearchDebug) {
          ACTS_VERBOSE(offset << "[>>] accepted.");
        }
        // Create the element
        auto identifier =
            m_cfg.identifierProvider != nullptr
                ? m_cfg.identifierProvider->identify(gctx, *tgNode)
                : Identifier();
        auto tgElement = std::make_shared<const Acts::TGeoDetectorElement>(
            identifier, tgNode, &tgTransform, layerConfig.localAxes,
            m_cfg.unit);
        // Record the element @todo solve with provided cache
        m_elementStore.push_back(tgElement);
        // Register the shared pointer to the surface for layer building
        layerSurfaces.push_back(tgElement->surface().getSharedPtr());

        // Record split range for eventual splitting
        double surfaceR = tgElement->surface().binningPositionValue(gctx, binR);
        double surfaceZ = tgElement->surface().binningPositionValue(gctx, binZ);

        // Split in R if configured to do so
        if (m_cfg.layerSplitToleranceR[type + 1] > 0.) {
          registerSplit(layerConfig.splitParametersR, surfaceR,
                        m_cfg.layerSplitToleranceR[type + 1],
                        layerConfig.splitRangeR);
        }
        // Split in Z if configured to do so
        if (m_cfg.layerSplitToleranceZ[type + 1] > 0.) {
          registerSplit(layerConfig.splitParametersZ, surfaceZ,
                        m_cfg.layerSplitToleranceZ[type + 1],
                        layerConfig.splitRangeZ);
        }
      } else if (type * z < 0 and m_cfg.nodeSearchDebug) {
        ACTS_VERBOSE("[xx] cancelled by side check.");
      } else if (not insideParseRange and m_cfg.nodeSearchDebug) {
        ACTS_VERBOSE("[xx] cancelled by parse range on side " << type);
        ACTS_VERBOSE("     r = " << r << " in ("
                                 << layerConfig.parseRangeR.first << ", "
                                 << layerConfig.parseRangeR.second << "] ?");
      }
    } else {
      // This is not yet the senstive one
      if (m_cfg.nodeSearchDebug) {
        ACTS_VERBOSE(offset << "[<<] not accepted, stepping down.");
      }
      // Build the matrix
      TGeoHMatrix nTransform =
          TGeoCombiTrans(tgTransform) * TGeoCombiTrans(*tgMatrix);
      std::string suffix = "_transform";
      nTransform.SetName((tNodeName + suffix).c_str());
      // If it's not accepted, get the associated volume
      TGeoVolume* nodeVolume = tgNode->GetVolume();
      // Now step down one further
      resolveSensitive(gctx, layerSurfaces, nodeVolume, nullptr, nTransform,
                       layerConfig, type, correctBranch, offset + "  ");
    }
  } else if (m_cfg.nodeSearchDebug) {
    ACTS_VERBOSE("No node present.");
  }
}
