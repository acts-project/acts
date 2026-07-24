// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Digitization/DigitizationConfigurator.hpp"

#include "Acts/Surfaces/AnnulusBounds.hpp"
#include "Acts/Surfaces/DiscTrapezoidBounds.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"
#include "Acts/Utilities/IAxis.hpp"
#include "Acts/Utilities/IMultiAxis.hpp"
#include "Acts/Utilities/Zip.hpp"
#include "ActsExamples/Digitization/SmearingConfig.hpp"
#include "ActsFatras/Digitization/UncorrelatedHitSmearer.hpp"

#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <vector>

namespace ActsExamples {

namespace {

/// @note This does not really compare if the configs are equal, therefore
/// it is no operator==. The contained std::function types cannot really
/// be checked for equality.
bool segmentationMaybeEqual(const std::shared_ptr<const Acts::IMultiAxis> &a,
                            const std::shared_ptr<const Acts::IMultiAxis> &b) {
  if (a == b) {
    return true;
  }
  if (a == nullptr || b == nullptr) {
    return false;
  }
  // The axis equality does not cover the axis directions
  return *a == *b &&
         std::ranges::equal(*a, *b, [](const auto &aa, const auto &ba) {
           return aa.getDirection() == ba.getDirection();
         });
}

bool digiConfigMaybeEqual(DigiComponentsConfig &a, DigiComponentsConfig &b) {
  // Check smearing config
  for (const auto &[as, bs] :
       Acts::zip(a.smearingDigiConfig.params, b.smearingDigiConfig.params)) {
    if (as.index != bs.index) {
      return false;
    }
  }
  if (a.smearingDigiConfig.maxRetries != b.smearingDigiConfig.maxRetries) {
    return false;
  }
  // Check geometric config
  const auto &ag = a.geometricDigiConfig;
  const auto &bg = b.geometricDigiConfig;
  return (ag.indices == bg.indices &&
          segmentationMaybeEqual(ag.segmentation, bg.segmentation) &&
          ag.thickness == bg.thickness && ag.threshold == bg.threshold &&
          ag.digital == bg.digital);
}

}  // namespace

void DigitizationConfigurator::operator()(const Acts::Surface *surface) {
  if (!surface->isSensitive()) {
    return;
  }

  Acts::GeometryIdentifier geoId = surface->geometryId();
  const auto dInputConfig = inputDigiComponents.find(geoId);
  if (dInputConfig == inputDigiComponents.end()) {
    return;
  }

  // The output config, copy over the smearing part
  DigiComponentsConfig dOutputConfig;
  dOutputConfig.smearingDigiConfig = dInputConfig->smearingDigiConfig;

  if (!dInputConfig->geometricDigiConfig.indices.empty()) {
    // Copy over what can be done
    dOutputConfig.geometricDigiConfig.indices =
        dInputConfig->geometricDigiConfig.indices;
    dOutputConfig.geometricDigiConfig.thickness =
        dInputConfig->geometricDigiConfig.thickness;
    dOutputConfig.geometricDigiConfig.chargeSmearer =
        dInputConfig->geometricDigiConfig.chargeSmearer;
    dOutputConfig.geometricDigiConfig.digital =
        dInputConfig->geometricDigiConfig.digital;

    const Acts::SurfaceBounds &sBounds = surface->bounds();
    const std::vector<double> boundValues = sBounds.values();

    const auto &inputSegmentation =
        dInputConfig->geometricDigiConfig.segmentation;
    if (inputSegmentation == nullptr) {
      throw std::invalid_argument(
          "DigitizationConfigurator: geometric digitization configured "
          "without a segmentation");
    }
    std::vector<std::unique_ptr<const Acts::IAxis>> outputAxes;

    switch (sBounds.type()) {
      // The module is a rectangle module
      case Acts::SurfaceBounds::eRectangle: {
        if (inputSegmentation->getAxis(0).getDirection() ==
            Acts::AxisDirection::AxisX) {
          const double minX = boundValues[Acts::RectangleBounds::eMinX];
          const double maxX = boundValues[Acts::RectangleBounds::eMaxX];

          const unsigned int nBins = inputSegmentation->getAxis(0).getNBins();

          outputAxes.push_back(Acts::IAxis::createEquidistant(
              Acts::AxisBoundaryType::Bound, minX, maxX, nBins,
              Acts::AxisDirection::AxisX));
        }
        if (inputSegmentation->getAxis(0).getDirection() ==
                Acts::AxisDirection::AxisY ||
            inputSegmentation->getNAxes() == 2) {
          const unsigned int accessBin =
              inputSegmentation->getNAxes() == 2 ? 1 : 0;

          const double minY = boundValues[Acts::RectangleBounds::eMinY];
          const double maxY = boundValues[Acts::RectangleBounds::eMaxY];

          const unsigned int nBins =
              inputSegmentation->getAxis(accessBin).getNBins();

          outputAxes.push_back(Acts::IAxis::createEquidistant(
              Acts::AxisBoundaryType::Bound, minY, maxY, nBins,
              Acts::AxisDirection::AxisY));
        }
      } break;

      // The module is a trapezoid module
      case Acts::SurfaceBounds::eTrapezoid: {
        if (inputSegmentation->getAxis(0).getDirection() ==
            Acts::AxisDirection::AxisX) {
          const double maxX =
              std::max(boundValues[Acts::TrapezoidBounds::eHalfLengthXnegY],
                       boundValues[Acts::TrapezoidBounds::eHalfLengthXposY]);

          const unsigned int nBins = inputSegmentation->getAxis(0).getNBins();

          outputAxes.push_back(Acts::IAxis::createEquidistant(
              Acts::AxisBoundaryType::Bound, -maxX, maxX, nBins,
              Acts::AxisDirection::AxisX));
        }
        if (inputSegmentation->getAxis(0).getDirection() ==
                Acts::AxisDirection::AxisY ||
            inputSegmentation->getNAxes() == 2) {
          const unsigned int accessBin =
              inputSegmentation->getNAxes() == 2 ? 1 : 0;

          const double maxY = boundValues[Acts::TrapezoidBounds::eHalfLengthY];

          const unsigned int nBins =
              inputSegmentation->getAxis(accessBin).getNBins();

          outputAxes.push_back(Acts::IAxis::createEquidistant(
              Acts::AxisBoundaryType::Bound, -maxY, maxY, nBins,
              Acts::AxisDirection::AxisY));
        }
      } break;

      // The module is an annulus module
      case Acts::SurfaceBounds::eAnnulus: {
        if (inputSegmentation->getAxis(0).getDirection() ==
            Acts::AxisDirection::AxisR) {
          const double minR = boundValues[Acts::AnnulusBounds::eMinR];
          const double maxR = boundValues[Acts::AnnulusBounds::eMaxR];

          const unsigned int nBins = inputSegmentation->getAxis(0).getNBins();

          outputAxes.push_back(Acts::IAxis::createEquidistant(
              Acts::AxisBoundaryType::Bound, minR, maxR, nBins,
              Acts::AxisDirection::AxisR));
        }
        if (inputSegmentation->getAxis(0).getDirection() ==
                Acts::AxisDirection::AxisPhi ||
            inputSegmentation->getNAxes() == 2) {
          const double averagePhi =
              boundValues[Acts::AnnulusBounds::eAveragePhi];
          const double minPhi =
              averagePhi - boundValues[Acts::AnnulusBounds::eMinPhiRel];
          const double maxPhi =
              averagePhi + boundValues[Acts::AnnulusBounds::eMaxPhiRel];

          const unsigned int nBins = inputSegmentation->getAxis(0).getNBins();

          outputAxes.push_back(Acts::IAxis::createEquidistant(
              Acts::AxisBoundaryType::Bound, minPhi, maxPhi, nBins,
              Acts::AxisDirection::AxisPhi));
        }
      } break;

      // The module is a Disc Trapezoid
      case Acts::SurfaceBounds::eDiscTrapezoid: {
        const double minR = boundValues[Acts::DiscTrapezoidBounds::eMinR];
        const double maxR = boundValues[Acts::DiscTrapezoidBounds::eMaxR];

        if (inputSegmentation->getAxis(0).getDirection() ==
            Acts::AxisDirection::AxisR) {
          const unsigned int nBins = inputSegmentation->getAxis(0).getNBins();

          outputAxes.push_back(Acts::IAxis::createEquidistant(
              Acts::AxisBoundaryType::Bound, minR, maxR, nBins,
              Acts::AxisDirection::AxisR));
        }
        if (inputSegmentation->getAxis(0).getDirection() ==
                Acts::AxisDirection::AxisPhi ||
            inputSegmentation->getNAxes() == 2) {
          const unsigned int accessBin =
              inputSegmentation->getNAxes() == 2 ? 1 : 0;

          const double hxMinR =
              boundValues[Acts::DiscTrapezoidBounds::eHalfLengthXminR];
          const double hxMaxR =
              boundValues[Acts::DiscTrapezoidBounds::eHalfLengthXmaxR];

          const double averagePhi =
              boundValues[Acts::DiscTrapezoidBounds::eAveragePhi];
          const double alphaMinR = std::atan2(minR, hxMinR);
          const double alphaMaxR = std::atan2(maxR, hxMaxR);
          const double alpha = std::max(alphaMinR, alphaMaxR);

          const unsigned int nBins =
              inputSegmentation->getAxis(accessBin).getNBins();

          outputAxes.push_back(Acts::IAxis::createEquidistant(
              Acts::AxisBoundaryType::Bound, averagePhi - alpha,
              averagePhi + alpha, nBins, Acts::AxisDirection::AxisPhi));
        }
      } break;

      case Acts::SurfaceBounds::eDisc: {
        if (inputSegmentation->getAxis(0).getDirection() ==
            Acts::AxisDirection::AxisR) {
          const double minR = boundValues[Acts::RadialBounds::eMinR];
          const double maxR = boundValues[Acts::RadialBounds::eMaxR];

          const unsigned int nBins = inputSegmentation->getAxis(0).getNBins();

          outputAxes.push_back(Acts::IAxis::createEquidistant(
              Acts::AxisBoundaryType::Bound, minR, maxR, nBins,
              Acts::AxisDirection::AxisR));
        }
        if (inputSegmentation->getAxis(0).getDirection() ==
                Acts::AxisDirection::AxisPhi ||
            inputSegmentation->getNAxes() == 2) {
          const unsigned int accessBin =
              inputSegmentation->getNAxes() == 2 ? 1 : 0;

          const double averagePhi =
              boundValues[Acts::RadialBounds::eAveragePhi];
          const double halfPhiSector =
              boundValues[Acts::RadialBounds::eHalfPhiSector];
          const double minPhi = averagePhi - halfPhiSector;
          const double maxPhi = averagePhi + halfPhiSector;

          const unsigned int nBins =
              inputSegmentation->getAxis(accessBin).getNBins();

          outputAxes.push_back(Acts::IAxis::createEquidistant(
              Acts::AxisBoundaryType::Bound, minPhi, maxPhi, nBins,
              Acts::AxisDirection::AxisPhi));
        }
      } break;

      default:
        break;
    }

    // Set the adapted segmentation class
    if (outputAxes.size() == 1) {
      dOutputConfig.geometricDigiConfig.segmentation =
          Acts::IMultiAxis::create(*outputAxes[0]);
    } else if (outputAxes.size() == 2) {
      dOutputConfig.geometricDigiConfig.segmentation =
          Acts::IMultiAxis::create(*outputAxes[0], *outputAxes[1]);
    }
  }

  // Compactify the output map where possible
  if (compactify) {
    // Check for a representing volume configuration, insert if not
    // present
    Acts::GeometryIdentifier volGeoId =
        Acts::GeometryIdentifier().withVolume(geoId.volume());

    auto volRep = volumeLayerComponents.find(volGeoId);
    if (volRep != volumeLayerComponents.end() &&
        digiConfigMaybeEqual(dOutputConfig, volRep->second)) {
      // return if the volume representation already covers this one
      return;
    }
    volumeLayerComponents[volGeoId] = dOutputConfig;
    outputDigiComponents.push_back({volGeoId, dOutputConfig});

    // Check for a representing layer configuration, insert if not present
    Acts::GeometryIdentifier volLayGeoId =
        Acts::GeometryIdentifier(volGeoId).withLayer(geoId.layer());
    const auto volLayRep = volumeLayerComponents.find(volLayGeoId);
    if (volLayRep != volumeLayerComponents.end() &&
        digiConfigMaybeEqual(dOutputConfig, volLayRep->second)) {
      return;
    }
    volumeLayerComponents[volLayGeoId] = dOutputConfig;
    outputDigiComponents.push_back({volLayGeoId, dOutputConfig});
  }

  // Insert into the output list
  outputDigiComponents.push_back({geoId, dOutputConfig});
}

}  // namespace ActsExamples
