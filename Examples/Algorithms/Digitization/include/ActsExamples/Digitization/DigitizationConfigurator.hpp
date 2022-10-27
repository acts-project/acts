// This file is part of the Acts project.
//
// Copyright (C) 2020-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryHierarchyMap.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Surfaces/AnnulusBounds.hpp"
#include "Acts/Surfaces/DiscTrapezoidBounds.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/BinningData.hpp"
#include "ActsExamples/Digitization/DigitizationAlgorithm.hpp"
#include "ActsExamples/Digitization/DigitizationConfig.hpp"

#include <map>
#include <memory>

namespace ActsExamples {

/// Helper configurator that takes a simplified (per volume, per layer)
/// Input digitization file and creates a full fletched per module
/// digitization configuration file.
///
/// It acts as a visitor and then builds a fully developed digitization file
/// for the geometric digitization, filling in the correct dimensions and
/// numer of bins.
///
/// The simplified file is assumed to have just one bin for the geometric
/// digitization, which is then used to calculate the number of bins with
/// respect to the bounds range.
struct DigitizationConfigurator {
  /// Simplified input components for digitization
  Acts::GeometryHierarchyMap<DigiComponentsConfig> inputDigiComponents;

  /// Final collection of output components
  std::vector<std::pair<Acts::GeometryIdentifier, DigiComponentsConfig>>
      outputDigiComponents;

  /// This tries to compactify the output map
  bool compactify = false;

  /// High level reference configurations for compactification
  std::map<Acts::GeometryIdentifier, DigiComponentsConfig>
      volumeLayerComponents;

  /// The visitor call for the geometry
  ///
  /// @param surface is the surfaces that is visited
  ///
  /// it adds an appropriate entry into the digitisation configuration
  void operator()(const Acts::Surface* surface) {
    if (surface->associatedDetectorElement() != nullptr) {
      Acts::GeometryIdentifier geoId = surface->geometryId();
      auto dInputConfig = inputDigiComponents.find((geoId));
      if (dInputConfig != inputDigiComponents.end()) {
        // The output config, copy over the smearing part
        DigiComponentsConfig dOutputConfig;
        dOutputConfig.smearingDigiConfig = dInputConfig->smearingDigiConfig;

        if (not dInputConfig->geometricDigiConfig.indices.empty()) {
          // Copy over what can be done
          dOutputConfig.geometricDigiConfig.indices =
              dInputConfig->geometricDigiConfig.indices;
          dOutputConfig.geometricDigiConfig.drift =
              dInputConfig->geometricDigiConfig.drift;
          dOutputConfig.geometricDigiConfig.thickness =
              dInputConfig->geometricDigiConfig.thickness;
          dOutputConfig.geometricDigiConfig.charge =
              dInputConfig->geometricDigiConfig.charge;
          dOutputConfig.geometricDigiConfig.digital =
              dInputConfig->geometricDigiConfig.digital;
          dOutputConfig.geometricDigiConfig.variances =
              dInputConfig->geometricDigiConfig.variances;

          const Acts::SurfaceBounds& sBounds = surface->bounds();
          auto boundValues = sBounds.values();

          const auto& inputSegmentation =
              dInputConfig->geometricDigiConfig.segmentation;
          Acts::BinUtility outputSegmentation;

          switch (sBounds.type()) {
            // The module is a rectangle module
            case Acts::SurfaceBounds::eRectangle: {
              if (inputSegmentation.binningData()[0].binvalue == Acts::binX) {
                Acts::ActsScalar minX =
                    boundValues[Acts::RectangleBounds::eMinX];
                Acts::ActsScalar maxX =
                    boundValues[Acts::RectangleBounds::eMaxX];
                unsigned int nBins = std::round(
                    (maxX - minX) / inputSegmentation.binningData()[0].step);
                outputSegmentation +=
                    Acts::BinUtility(nBins, minX, maxX, Acts::open, Acts::binX);
              }
              if (inputSegmentation.binningData()[0].binvalue == Acts::binY or
                  inputSegmentation.dimensions() == 2) {
                unsigned int accessBin =
                    inputSegmentation.dimensions() == 2 ? 1 : 0;
                Acts::ActsScalar minY =
                    boundValues[Acts::RectangleBounds::eMinY];
                Acts::ActsScalar maxY =
                    boundValues[Acts::RectangleBounds::eMaxY];
                unsigned int nBins =
                    std::round((maxY - minY) /
                               inputSegmentation.binningData()[accessBin].step);
                outputSegmentation +=
                    Acts::BinUtility(nBins, minY, maxY, Acts::open, Acts::binY);
              }
            } break;

            // The module is a trapezoid module
            case Acts::SurfaceBounds::eTrapezoid: {
              if (inputSegmentation.binningData()[0].binvalue == Acts::binX) {
                Acts::ActsScalar maxX = std::max(
                    boundValues[Acts::TrapezoidBounds::eHalfLengthXnegY],
                    boundValues[Acts::TrapezoidBounds::eHalfLengthXposY]);
                unsigned int nBins = std::round(
                    2 * maxX / inputSegmentation.binningData()[0].step);
                outputSegmentation += Acts::BinUtility(nBins, -maxX, maxX,
                                                       Acts::open, Acts::binX);
              }
              if (inputSegmentation.binningData()[0].binvalue == Acts::binY or
                  inputSegmentation.dimensions() == 2) {
                unsigned int accessBin =
                    inputSegmentation.dimensions() == 2 ? 1 : 0;
                Acts::ActsScalar maxY =
                    boundValues[Acts::TrapezoidBounds::eHalfLengthY];
                unsigned int nBins =
                    std::round((2 * maxY) /
                               inputSegmentation.binningData()[accessBin].step);
                outputSegmentation += Acts::BinUtility(nBins, -maxY, maxY,
                                                       Acts::open, Acts::binY);
              }
            } break;

            // The module is an annulus module
            case Acts::SurfaceBounds::eAnnulus: {
              if (inputSegmentation.binningData()[0].binvalue == Acts::binR) {
                Acts::ActsScalar minR = boundValues[Acts::AnnulusBounds::eMinR];
                Acts::ActsScalar maxR = boundValues[Acts::AnnulusBounds::eMaxR];
                unsigned int nBins = std::round(
                    (maxR - minR) / inputSegmentation.binningData()[0].step);
                outputSegmentation +=
                    Acts::BinUtility(nBins, minR, maxR, Acts::open, Acts::binR);
              }
              if (inputSegmentation.binningData()[0].binvalue == Acts::binPhi or
                  inputSegmentation.dimensions() == 2) {
                unsigned int accessBin =
                    inputSegmentation.dimensions() == 2 ? 1 : 0;
                Acts::ActsScalar averagePhi =
                    boundValues[Acts::AnnulusBounds::eAveragePhi];
                Acts::ActsScalar minPhi =
                    averagePhi - boundValues[Acts::AnnulusBounds::eMinPhiRel];
                Acts::ActsScalar maxPhi =
                    averagePhi + boundValues[Acts::AnnulusBounds::eMaxPhiRel];
                unsigned int nBins =
                    std::round((maxPhi - minPhi) /
                               inputSegmentation.binningData()[accessBin].step);
                outputSegmentation += Acts::BinUtility(
                    nBins, minPhi, maxPhi, Acts::open, Acts::binPhi);
              }

            } break;

            // The module is a Disc Trapezoid
            case Acts::SurfaceBounds::eDiscTrapezoid: {
              Acts::ActsScalar minR =
                  boundValues[Acts::DiscTrapezoidBounds::eMinR];
              Acts::ActsScalar maxR =
                  boundValues[Acts::DiscTrapezoidBounds::eMaxR];

              if (inputSegmentation.binningData()[0].binvalue == Acts::binR) {
                unsigned int nBins = std::round(
                    (maxR - minR) / inputSegmentation.binningData()[0].step);
                outputSegmentation +=
                    Acts::BinUtility(nBins, minR, maxR, Acts::open, Acts::binR);
              }
              if (inputSegmentation.binningData()[0].binvalue == Acts::binPhi or
                  inputSegmentation.dimensions() == 2) {
                unsigned int accessBin =
                    inputSegmentation.dimensions() == 2 ? 1 : 0;
                Acts::ActsScalar hxMinR =
                    boundValues[Acts::DiscTrapezoidBounds::eHalfLengthXminR];
                Acts::ActsScalar hxMaxR =
                    boundValues[Acts::DiscTrapezoidBounds::eHalfLengthXmaxR];

                Acts::ActsScalar averagePhi =
                    boundValues[Acts::DiscTrapezoidBounds::eAveragePhi];
                Acts::ActsScalar alphaMinR = std::atan2(minR, hxMinR);
                Acts::ActsScalar alphaMaxR = std::atan2(maxR, hxMaxR);
                Acts::ActsScalar alpha = std::max(alphaMinR, alphaMaxR);
                unsigned int nBins =
                    std::round(2 * alpha /
                               inputSegmentation.binningData()[accessBin].step);
                outputSegmentation += Acts::BinUtility(
                    nBins, averagePhi - alpha, averagePhi + alpha, Acts::open,
                    Acts::binPhi);
              }

            } break;

            case Acts::SurfaceBounds::eDisc: {
              if (inputSegmentation.binningData()[0].binvalue == Acts::binR) {
                Acts::ActsScalar minR = boundValues[Acts::RadialBounds::eMinR];
                Acts::ActsScalar maxR = boundValues[Acts::RadialBounds::eMaxR];
                unsigned int nBins = std::round(
                    (maxR - minR) / inputSegmentation.binningData()[0].step);
                outputSegmentation +=
                    Acts::BinUtility(nBins, minR, maxR, Acts::open, Acts::binR);
              }
              if (inputSegmentation.binningData()[0].binvalue == Acts::binPhi or
                  inputSegmentation.dimensions() == 2) {
                unsigned int accessBin =
                    inputSegmentation.dimensions() == 2 ? 1 : 0;

                Acts::ActsScalar averagePhi =
                    boundValues[Acts::RadialBounds::eAveragePhi];
                Acts::ActsScalar halfPhiSector =
                    boundValues[Acts::RadialBounds::eHalfPhiSector];
                Acts::ActsScalar minPhi = averagePhi - halfPhiSector;
                Acts::ActsScalar maxPhi = averagePhi + halfPhiSector;

                unsigned int nBins =
                    std::round((maxPhi - minPhi) /
                               inputSegmentation.binningData()[accessBin].step);
                outputSegmentation += Acts::BinUtility(
                    nBins, minPhi, maxPhi, Acts::open, Acts::binPhi);
              }

            } break;

            default:
              break;
          }
          // Set the adapted segmentation class
          dOutputConfig.geometricDigiConfig.segmentation = outputSegmentation;
        }

        // Compactify the output map where possible
        if (compactify) {
          // Check for a representing volume configuration, insert if not
          // present
          Acts::GeometryIdentifier volGeoId =
              Acts::GeometryIdentifier().setVolume(geoId.volume());

          auto volRep = volumeLayerComponents.find(volGeoId);
          if (volRep != volumeLayerComponents.end() and
              dOutputConfig == volRep->second) {
            // return if the volume representation already covers this one
            return;
          } else {
            volumeLayerComponents[volGeoId] = dOutputConfig;
            outputDigiComponents.push_back({volGeoId, dOutputConfig});
          }

          // Check for a representing layer configuration, insert if not present
          Acts::GeometryIdentifier volLayGeoId =
              Acts::GeometryIdentifier(volGeoId).setLayer(geoId.layer());
          auto volLayRep = volumeLayerComponents.find(volLayGeoId);

          if (volLayRep != volumeLayerComponents.end() and
              dOutputConfig == volLayRep->second) {
            return;
          } else {
            volumeLayerComponents[volLayGeoId] = dOutputConfig;
            outputDigiComponents.push_back({volLayGeoId, dOutputConfig});
          }
        }

        // Insert into the output list
        outputDigiComponents.push_back({geoId, dOutputConfig});
      }
    }
  }
};
}  // namespace ActsExamples
