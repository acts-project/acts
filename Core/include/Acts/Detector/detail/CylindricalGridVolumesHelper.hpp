// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Detector/Portal.hpp"
#include "Acts/Detector/PortalGenerators.hpp"
#include "Acts/Detector/detail/CylindricalDetectorHelper.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Navigation/DetectorVolumeFinders.hpp"
#include "Acts/Navigation/SurfaceCandidatesUpdators.hpp"

#include <limits>
#include <map>
#include <memory>
#include <vector>

namespace Acts {

namespace Experimental {

namespace detail {
namespace CylindricalGridVolumesHelper {

struct Options {
  /// The volume naming schema
  std::string volumeBaseName = "GridVolume";
  /// Boolean to indicate if common surfaces should be generated
  bool generateCommonPortals = false;
  /// If this is set to some reasonable value, a polygon approximation
  /// will be used for the detector volume, and the full description will
  /// be planar surfaces
  int polygonApproximation = -1;
  /// An additional transform
  Transform3 transform = Transform3::Identity();
};

/// @brief Construct the volumes for this grid
///
/// @tparam grid_type the type of the grid used, it encapsulates the
/// need actual bound, closed, variable, equidistant, etc. types
///
/// @note no checking is done if the grid is actually a meaningful z,r,phi grid
/// @note no checking is done on the dimensionality of the grid, has to be done
/// upstream
///
/// @param gctx the geometry context
/// @param cylindricalGrid the cylindrical grid
/// @param rootVolumesGrid [pin, out] the root volumes grid
/// @param detectorSurfaces the detector surfaces to be filled
///
/// @return a vector of detector voluemes
template <typename grid_type, typename root_volumes_grid_type>
std::vector<std::shared_ptr<DetectorVolume>> buildVolumes(
    const GeometryContext& gctx, const grid_type& cylindricalGrid,
    root_volumes_grid_type& rootVolumesGrid,
    const std::vector<std::shared_ptr<Surface>>& detectorSurfaces = {},
    const Options& options = Options{}) {
  // Shorthand for the grid point
  using GridPoint = typename root_volumes_grid_type::point_t;

  // Get the axis and
  auto axes = cylindricalGrid.axes();
  // Generate the axis bin edges
  auto edgesZ = axes[0]->getBinEdges();
  auto edgesR = axes[1]->getBinEdges();
  auto edgesPhi = axes[2]->getBinEdges();
  // The detector volumes to be returned
  std::vector<std::shared_ptr<DetectorVolume>> detectorVolumes;
  detectorVolumes.reserve((edgesZ.size() - 1) * (edgesR.size() - 1) *
                          (edgesPhi.size() - 1));

  /// Default portal generation code
  auto portalGenerator = defaultPortalGenerator();

  // Remeber the last r-phi disc
  std::vector<std::shared_ptr<DetectorVolume>> lastRPhiDisc = {};

  for (auto [iz, z] : enumerate(edgesZ)) {
    if (iz > 0) {
      // Collect the full disc
      std::vector<std::shared_ptr<DetectorVolume>> rPhiDisc = {};
      rPhiDisc.reserve((edgesR.size() - 1) * (edgesPhi.size() - 1));

      // The z position of the volume
      ActsScalar zPos = 0.5 * (z + edgesZ[iz - 1]);
      ActsScalar zHalfRange = 0.5 * (z - edgesZ[iz - 1]);
      // Create the bin transform
      Transform3 binTransform = Transform3::Identity();
      binTransform.pretranslate(Vector3(0., 0., zPos));
      // remember the last phi ring
      std::vector<std::shared_ptr<DetectorVolume>> lastPhiRing = {};

      for (auto [ir, r] : enumerate(edgesR)) {
        if (ir > 0) {
          // The r range of the volume
          ActsScalar rMin = edgesR[ir - 1];
          ActsScalar rMax = r;
          // Portal connection: phi ring for common portals option
          std::vector<std::shared_ptr<DetectorVolume>> phiRing = {};
          phiRing.reserve(edgesPhi.size() - 1);

          // Portal connection: last/first for individual portal option
          std::shared_ptr<DetectorVolume> lastVolume = nullptr;
          std::shared_ptr<DetectorVolume> firstVolume = nullptr;
          for (auto [iphi, phi] : enumerate(edgesPhi)) {
            if (iphi > 0) {
              // Construct the volume base name
              std::string volumeName = options.volumeBaseName + "_z" +
                                       std::to_string(iz - 1) + "_r" +
                                       std::to_string(ir - 1) + "_phi" +
                                       std::to_string(iphi - 1);
              // The phi range of the volume
              ActsScalar phiPos = 0.5 * (phi + edgesPhi[iphi - 1]);
              ActsScalar phiHalfRange = 0.5 * (phi - edgesPhi[iphi - 1]);
              // Construct the bounds
              auto binBounds = std::make_unique<CylinderVolumeBounds>(
                  rMin, rMax, zHalfRange, phiHalfRange, phiPos);
              // The bin volume to be created, either with surfaces or empty
              std::shared_ptr<DetectorVolume> binVolume = nullptr;
              // Pull the surfaces from the grid
              GridPoint p = {zPos, 0.5 * (rMin + rMax), phiPos};
              auto binSurfaceIndices = cylindricalGrid.atPosition(p);
              // Empty surfaces
              if (binSurfaceIndices.empty()) {
                // Contruct the bin volume
                binVolume = DetectorVolumeFactory::construct(
                    portalGenerator, gctx, volumeName, binTransform,
                    std::move(binBounds), tryAllPortals());
              } else {
                // The surfaces to be filled
                std::vector<std::shared_ptr<Surface>> binSurfaces;
                binSurfaces.reserve(binSurfaceIndices.size());
                for (auto s : binSurfaceIndices) {
                  binSurfaces.push_back(detectorSurfaces.at(s));
                }
                // Contruct the bin volume
                binVolume = DetectorVolumeFactory::construct(
                    portalGenerator, gctx, volumeName, binTransform,
                    std::move(binBounds), binSurfaces, {}, tryNoVolumes(),
                    tryAllPortalsAndSurfaces());
              }
              // Register the detector volume in the grid
              rootVolumesGrid.atPosition(p) = detectorVolumes.size();
              // Add the volume to the phi ring for the r fusing
              phiRing.push_back(binVolume);
              // And to the full disc for the z fusing
              rPhiDisc.push_back(binVolume);
              if (not options.generateCommonPortals) {
                if (lastVolume != nullptr) {
                  CylindricalDetectorHelper::fuseInPhi(*lastVolume, *binVolume);
                } else {
                  firstVolume = binVolume;
                }
                // Register the last volume
                lastVolume = binVolume;
              }
              // and store it
              detectorVolumes.push_back(binVolume);
            }
          }
          // Connect the phi ring - generate common portals
          if (options.generateCommonPortals) {
            CylindricalDetectorHelper::connectInPhi(gctx, phiRing);
          } else {
            CylindricalDetectorHelper::fuseInPhi(*lastVolume, *firstVolume);
          }
          if (not lastPhiRing.empty()) {
            // Connect the phi rings
            for (auto [iphi, v] : enumerate(lastPhiRing)) {
              auto keepCoverVolume = v;
              auto wasteCoverVolume = phiRing.at(iphi);
              CylindricalDetectorHelper::fuseInR(*keepCoverVolume,
                                                 *wasteCoverVolume);
            }
          }
          lastPhiRing = phiRing;
        }
      }  // end of r loop
      if (not lastRPhiDisc.empty()) {
        // Connect the r-phi discs
        for (auto [ir, v] : enumerate(lastRPhiDisc)) {
          auto keepEndplateVolume = v;
          auto wastEndplateVolume = rPhiDisc.at(ir);
          CylindricalDetectorHelper::fuseInZ(*keepEndplateVolume,
                                             *wastEndplateVolume);
        }
      }
      lastRPhiDisc = rPhiDisc;
    }
  }
  return detectorVolumes;
}

}  // namespace CylindricalGridVolumesHelper
}  // namespace detail
}  // namespace Experimental
}  // namespace Acts
