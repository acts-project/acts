// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/GenericDetector/ProtoLayerCreator.hpp"

using Acts::VectorHelpers::phi;

namespace ActsExamples::Generic {

std::vector<ProtoLayerSurfaces> ProtoLayerCreator::centralProtoLayers(
    const Acts::GeometryContext& gctx) const {
  // create the vector
  std::vector<ProtoLayerSurfaces> cpLayers;

  // ----------------------- central layers -------------------------
  // the central layers
  std::size_t numcLayers = m_cfg.centralLayerRadii.size();
  if (numcLayers != 0u) {
    ACTS_DEBUG("Configured to build " << numcLayers
                                      << " active central layers.");
    cpLayers.reserve(numcLayers);
    // loop through
    for (std::size_t icl = 0; icl < numcLayers; ++icl) {
      // layer R/Z
      double layerR = m_cfg.centralLayerRadii.at(icl);
      // some screen output
      ACTS_DEBUG("Build layer " << icl << " with target radius = " << layerR);

      // prepare the Surface vector
      std::vector<std::shared_ptr<Acts::Surface>> sVector;
      // assign the current envelope
      double layerEnvelopeCoverZ =
          !m_cfg.centralLayerEnvelopes.empty()
              ? m_cfg.centralLayerEnvelopes.at(icl).second
              : 0.;
      // module size & tilt
      double modulePhiTilt = m_cfg.centralModuleTiltPhi.at(icl);
      double moduleHalfX = m_cfg.centralModuleHalfX.at(icl);
      double moduleHalfY = m_cfg.centralModuleHalfY.at(icl);
      double moduleThickness = m_cfg.centralModuleThickness.at(icl);
      // create the shared module
      auto moduleBounds =
          std::make_shared<Acts::RectangleBounds>(moduleHalfX, moduleHalfY);
      std::size_t nCentralModules =
          m_cfg.centralModuleBinningSchema.at(icl).first *
          m_cfg.centralModuleBinningSchema.at(icl).second;

      ACTS_DEBUG("- number of modules "
                 << nCentralModules << " ( from "
                 << m_cfg.centralModuleBinningSchema.at(icl).first << " x "
                 << m_cfg.centralModuleBinningSchema.at(icl).second << " )");

      sVector.reserve(nCentralModules);

      // prepartation :
      // create the Module material from input
      std::shared_ptr<const Acts::ISurfaceMaterial> moduleMaterialPtr = nullptr;
      if (!m_cfg.centralModuleMaterial.empty()) {
        // get the sensor material from configuration
        moduleMaterialPtr = m_cfg.centralModuleMaterial.at(icl);
      }

      // confirm
      if (m_cfg.centralModulePositions.at(icl).size() != nCentralModules) {
        ACTS_WARNING("Mismatching module numbers, configuration error!");
        ACTS_WARNING("- Binning schema suggests : " << nCentralModules);
        ACTS_WARNING("- Positions provided are  : "
                     << m_cfg.centralModulePositions.at(icl).size());
      }
      // loop over the position, create the modules
      for (auto& moduleCenter : m_cfg.centralModulePositions.at(icl)) {
        // create the association transform
        double modulePhi = phi(moduleCenter);
        // the local z axis is the normal vector
        Acts::Vector3 moduleLocalZ(cos(modulePhi + modulePhiTilt),
                                   sin(modulePhi + modulePhiTilt), 0.);
        // the local y axis is the global z axis
        Acts::Vector3 moduleLocalY(0., 0., 1);
        // the local x axis the normal to local y,z
        Acts::Vector3 moduleLocalX(-sin(modulePhi + modulePhiTilt),
                                   cos(modulePhi + modulePhiTilt), 0.);
        // create the RotationMatrix
        Acts::RotationMatrix3 moduleRotation;
        moduleRotation.col(0) = moduleLocalX;
        moduleRotation.col(1) = moduleLocalY;
        moduleRotation.col(2) = moduleLocalZ;
        // get the moduleTransform
        std::shared_ptr<Acts::Transform3> mutableModuleTransform =
            std::make_shared<Acts::Transform3>(
                Acts::Translation3(moduleCenter) * moduleRotation);
        // stereo angle if necessary
        if (!m_cfg.centralModuleFrontsideStereo.empty() &&
            m_cfg.centralModuleFrontsideStereo.at(icl) != 0.) {
          // twist by the stereo angle
          double stereo = m_cfg.centralModuleFrontsideStereo.at(icl);
          (*mutableModuleTransform) *=
              Acts::AngleAxis3(-stereo, Acts::Vector3::UnitZ());
        }

        // Finalize the transform
        auto moduleTransform = std::const_pointer_cast<const Acts::Transform3>(
            mutableModuleTransform);
        // create the module
        auto moduleElement = m_cfg.detectorElementFactory(
            moduleTransform, moduleBounds, moduleThickness, moduleMaterialPtr);

        // register the surface
        sVector.push_back(moduleElement->surface().getSharedPtr());
        // IF double modules exist
        // and the backside one (if configured to do so)
        if (!m_cfg.centralModuleBacksideGap.empty()) {
          // create the module identifier

          Acts::Vector3 bsModuleCenter =
              moduleCenter +
              m_cfg.centralModuleBacksideGap.at(icl) * moduleLocalZ;
          mutableModuleTransform = std::make_shared<Acts::Transform3>(
              Acts::Translation3(bsModuleCenter) * moduleRotation);
          // apply the stereo
          if (!m_cfg.centralModuleBacksideStereo.empty()) {
            // twist by the stereo angle
            double stereoBackSide = m_cfg.centralModuleBacksideStereo.at(icl);
            (*mutableModuleTransform) *=
                Acts::AngleAxis3(-stereoBackSide, Acts::Vector3::UnitZ());
          }
          // Finalize the transform
          moduleTransform = std::const_pointer_cast<const Acts::Transform3>(
              mutableModuleTransform);
          // create the backseide moulde
          auto bsModuleElement =
              m_cfg.detectorElementFactory(moduleTransform, moduleBounds,
                                           moduleThickness, moduleMaterialPtr);
        }
      }

      std::size_t phiBins = m_cfg.centralModuleBinningSchema.at(icl).first;
      phiBins *= m_cfg.centralLayerBinMultipliers.first;
      std::size_t zBins = m_cfg.centralModuleBinningSchema.at(icl).second;
      zBins *= m_cfg.centralLayerBinMultipliers.second;
      // create the surface array - it will also fill the accessible binmember
      // cache if available
      Acts::MutableProtoLayer pl(gctx, sVector);
      pl.envelope[Acts::AxisDirection::AxisR] = {m_cfg.approachSurfaceEnvelope,
                                                 m_cfg.approachSurfaceEnvelope};
      pl.envelope[Acts::AxisDirection::AxisZ] = {layerEnvelopeCoverZ,
                                                 layerEnvelopeCoverZ};

      // Record the proto layer and the surfaces for the later layer building
      ProtoLayerSurfaces pls{std::move(pl), sVector, phiBins, zBins};
      cpLayers.push_back(std::move(pls));
    }
  }
  return cpLayers;
}

std::vector<ProtoLayerSurfaces> ProtoLayerCreator::negativeProtoLayers(
    const Acts::GeometryContext& gctx) const {
  return createProtoLayers(gctx, -1);
}

std::vector<ProtoLayerSurfaces> ProtoLayerCreator::positiveProtoLayers(
    const Acts::GeometryContext& gctx) const {
  return createProtoLayers(gctx, 1);
}

ProtoLayerCreator::ProtoLayerCreator(const ProtoLayerCreator::Config& cfg,
                                     std::unique_ptr<const Acts::Logger> log)
    : m_cfg(cfg), m_logger(std::move(log)) {
  if (!m_cfg.detectorElementFactory) {
    throw std::invalid_argument("Detector element factory is not set");
  }
}

std::vector<ProtoLayerSurfaces> ProtoLayerCreator::createProtoLayers(
    const Acts::GeometryContext& gctx, int side) const {
  // Count the current detector modules identifiers
  // the return layers
  std::vector<ProtoLayerSurfaces> epLayers;
  // -------------------------------- endcap type layers
  // pos/neg layers
  std::size_t numpnLayers = m_cfg.posnegLayerPositionsZ.size();
  if (numpnLayers != 0u) {
    ACTS_DEBUG("Configured to build 2 * "
               << numpnLayers << " passive positive/negative side layers.");
    epLayers.reserve(numpnLayers);

    /// this is the loop over the layer positions
    for (std::size_t ipnl = 0; ipnl < numpnLayers; ++ipnl) {
      // some screen output
      ACTS_VERBOSE("- building layer "
                   << ipnl << " and " << numpnLayers + ipnl << " at z = "
                   << side * m_cfg.posnegLayerPositionsZ.at(ipnl));
      /// some preparation work
      // define the layer envelope
      double layerEnvelopeR = m_cfg.posnegLayerEnvelopeR.at(ipnl);
      // prepare for the r binning
      std::vector<std::shared_ptr<Acts::Surface>> esVector;
      // now fill the vectors
      std::size_t ipnR = 0;
      for (auto& discModulePositions : m_cfg.posnegModulePositions.at(ipnl)) {
        ACTS_VERBOSE("- building ring " << ipnR << " for this layer.");
        // now prepare all the shared stuff
        // (0) module specifications
        double moduleThickness = m_cfg.posnegModuleThickness.at(ipnl).at(ipnR);
        double moduleMinHalfX = m_cfg.posnegModuleMinHalfX.at(ipnl).at(ipnR);
        double moduleMaxHalfX = 0.;
        if (m_cfg.posnegModuleMaxHalfX.size() > ipnl &&
            m_cfg.posnegModuleMaxHalfX.at(ipnl).size() > ipnR) {
          moduleMaxHalfX = m_cfg.posnegModuleMaxHalfX.at(ipnl).at(ipnR);
        }
        double moduleHalfY = m_cfg.posnegModuleHalfY.at(ipnl).at(ipnR);
        // (1) module bounds
        // create the bounds
        std::shared_ptr<const Acts::PlanarBounds> moduleBounds;
        if (moduleMaxHalfX != 0. && moduleMinHalfX != moduleMaxHalfX) {
          moduleBounds = std::make_shared<Acts::TrapezoidBounds>(
              moduleMinHalfX, moduleMaxHalfX, moduleHalfY);
        } else {
          moduleBounds = std::make_shared<Acts::RectangleBounds>(moduleMinHalfX,
                                                                 moduleHalfY);
        }
        // (2)) module material
        // create the Module material from input
        std::shared_ptr<const Acts::ISurfaceMaterial> moduleMaterialPtr =
            nullptr;
        if (!m_cfg.posnegModuleMaterial.empty()) {
          // and create the shared pointer
          moduleMaterialPtr = m_cfg.posnegModuleMaterial.at(ipnl).at(ipnR);
        }

        // low loop over the phi positions and build the stuff
        for (auto& ringModulePosition : discModulePositions) {
          // the module transform from the position
          double modulePhi = phi(ringModulePosition);
          // the center position of the modules
          Acts::Vector3 moduleCenter(ringModulePosition);
          moduleCenter.z() *= side;
          // the rotation matrix of the module
          Acts::Vector3 moduleLocalY(cos(modulePhi), sin(modulePhi), 0.);
          // take different axis to have the same readout direction
          Acts::Vector3 moduleLocalZ(0., 0., side * 1.);
          Acts::Vector3 moduleLocalX = moduleLocalY.cross(moduleLocalZ);
          // local rotation matrices
          // create the RotationMatrix - negative side
          Acts::RotationMatrix3 moduleRotation;
          moduleRotation.col(0) = moduleLocalX;
          moduleRotation.col(1) = moduleLocalY;
          moduleRotation.col(2) = moduleLocalZ;
          // the transforms for the two modules
          std::shared_ptr<const Acts::Transform3> moduleTransform =
              std::make_shared<const Acts::Transform3>(
                  Acts::Translation3(moduleCenter) * moduleRotation);

          // create the module
          auto moduleElement =
              m_cfg.detectorElementFactory(moduleTransform, moduleBounds,
                                           moduleThickness, moduleMaterialPtr);

          // now deal with the potential backside
          if (!m_cfg.posnegModuleBacksideGap.empty()) {
            // the new centers
            moduleCenter =
                moduleCenter +
                m_cfg.posnegModuleBacksideGap.at(ipnl).at(ipnR) * moduleLocalZ;
            // the new transforms
            auto mutableModuleTransform = std::make_shared<Acts::Transform3>(
                Acts::Translation3(moduleCenter) * moduleRotation);
            // apply the stereo
            if (!m_cfg.posnegModuleBacksideStereo.empty()) {
              // twist by the stereo angle
              double stereoBackSide =
                  m_cfg.posnegModuleBacksideStereo.at(ipnl).at(ipnR);
              (*mutableModuleTransform) *=
                  Acts::AngleAxis3(-stereoBackSide, Acts::Vector3::UnitZ());
            }
            // Finalize the transform
            moduleTransform = std::const_pointer_cast<const Acts::Transform3>(
                mutableModuleTransform);
            // everything is set for the next module
            auto bsModuleElement = m_cfg.detectorElementFactory(
                moduleTransform, moduleBounds, moduleThickness,
                moduleMaterialPtr);
          }
          // create the surface
          esVector.push_back(moduleElement->surface().getSharedPtr());
        }
        // counter of rings
        ++ipnR;
      }
      // the binning
      std::size_t layerBinsR = m_cfg.posnegModulePhiBins.at(ipnl).size();
      // never multiply 1 single r-bin, does not make sense
      if (layerBinsR > 1) {
        // multiply with the given bin multiplier
        layerBinsR *= m_cfg.posnegLayerBinMultipliers.first;
      }
      std::size_t layerBinsPhi = 0;
      // take the minimum phi bins in that layer
      for (unsigned int phiBins : m_cfg.posnegModulePhiBins.at(ipnl)) {
        layerBinsPhi = phiBins < layerBinsPhi ? phiBins : layerBinsPhi;
        layerBinsPhi *= m_cfg.posnegLayerBinMultipliers.second;
      }
      // create the layers with the surface arrays
      Acts::MutableProtoLayer ple(gctx, esVector);
      ple.envelope[Acts::AxisDirection::AxisR] = {layerEnvelopeR,
                                                  layerEnvelopeR};
      ple.envelope[Acts::AxisDirection::AxisZ] = {
          m_cfg.approachSurfaceEnvelope, m_cfg.approachSurfaceEnvelope};

      // push it into the layer vector
      ProtoLayerSurfaces ples{std::move(ple), esVector, layerBinsR,
                              layerBinsPhi};
      epLayers.push_back(std::move(ples));
    }
  }
  return epLayers;
}

}  // namespace ActsExamples::Generic
