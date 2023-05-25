// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Detector/detail/SupportHelper.hpp"

#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"

#include <stdexcept>

std::vector<std::shared_ptr<Acts::Surface>>
Acts::Experimental::detail::SupportHelper::cylindricalSupport(
    const Transform3& transform, const std::array<ActsScalar, 6u>& bounds,
    unsigned int splits) {
  // Return vector preparation
  std::vector<std::shared_ptr<Acts::Surface>> cSupport;
  if (splits == 1u) {
    // No splitting is done in this case
    cSupport.push_back(Surface::makeShared<CylinderSurface>(
        transform, std::make_shared<CylinderBounds>(bounds)));
  } else {
    // Split into n(splits) planar surfaces, prep work:
    ActsScalar r = bounds[0u];
    ActsScalar halfZ = bounds[1u];
    ActsScalar minPhi = bounds[3u] - bounds[2u];
    ActsScalar maxPhi = bounds[3u] + bounds[2u];
    ActsScalar dHalfPhi = (maxPhi - minPhi) / (2 * splits);
    ActsScalar cosPhiHalf = std::cos(dHalfPhi);
    ActsScalar sinPhiHalf = std::sin(dHalfPhi);
    ActsScalar planeR = r * cosPhiHalf;
    ActsScalar planeHalfX = r * sinPhiHalf;
    ActsScalar planeZ = transform.translation().z();

    auto sRectangle =
        std::make_shared<Acts::RectangleBounds>(planeHalfX, halfZ);
    // Now create the Trapezoids
    for (unsigned int iphi = 0; iphi < splits; ++iphi) {
      // Get the moduleTransform
      ActsScalar phi = -M_PI + (iphi + 0.5) * 2 * dHalfPhi;
      ActsScalar cosPhi = std::cos(phi);
      ActsScalar sinPhi = std::sin(phi);
      ActsScalar planeX = planeR * cosPhi;
      ActsScalar planeY = planeR * sinPhi;

      Acts::Vector3 planeCenter(planeX, planeY, planeZ);
      Acts::Vector3 planeAxisZ(cosPhi, sinPhi, 0.);
      Acts::Vector3 planeAxisY(0., 0., 1.);
      Acts::Vector3 planeAxisX = planeAxisY.cross(planeAxisZ);

      RotationMatrix3 planeRotation;
      planeRotation.col(0) = planeAxisX;
      planeRotation.col(1) = planeAxisY;
      planeRotation.col(2) = planeAxisZ;

      Transform3 sTransform{planeRotation};
      sTransform.pretranslate(planeCenter);
      // Place it
      cSupport.push_back(
          Surface::makeShared<PlaneSurface>(sTransform, sRectangle));
    }
  }

  return cSupport;
}

std::vector<std::shared_ptr<Acts::Surface>>
Acts::Experimental::detail::SupportHelper::discSupport(
    const Transform3& transform, const std::array<ActsScalar, 4u>& bounds,
    unsigned int splits) {
  // Return vector
  std::vector<std::shared_ptr<Acts::Surface>> dSupport;
  if (splits == 1u) {
    // No splitting is done in this case
    dSupport.push_back(Surface::makeShared<DiscSurface>(
        transform, std::make_shared<RadialBounds>(bounds)));
  } else {
    // Split into n(splits) planar surfaces in phi, prep work:
    ActsScalar minR = bounds[0u];
    ActsScalar maxR = bounds[1u];
    ActsScalar minPhi = bounds[3u] - bounds[2u];
    ActsScalar maxPhi = bounds[3u] + bounds[2u];
    ActsScalar dHalfPhi = (maxPhi - minPhi) / (2 * splits);
    ActsScalar cosPhiHalf = std::cos(dHalfPhi);
    ActsScalar sinPhiHalf = std::sin(dHalfPhi);
    ActsScalar maxLocY = maxR * cosPhiHalf;
    ActsScalar minLocY = minR * cosPhiHalf;
    ActsScalar hR = 0.5 * (maxLocY + minLocY);
    ActsScalar hY = 0.5 * (maxLocY - minLocY);
    ActsScalar hXminY = minR * sinPhiHalf;
    ActsScalar hXmaxY = maxR * sinPhiHalf;
    // Split trapezoid
    auto sTrapezoid =
        std::make_shared<Acts::TrapezoidBounds>(hXminY, hXmaxY, hY);
    Vector3 zAxis = transform.rotation().col(2);
    ActsScalar zPosition = transform.translation().z();
    // Now create the Trapezoids
    for (unsigned int iphi = 0; iphi < splits; ++iphi) {
      // Create the split module transform
      ActsScalar phi = -M_PI + (iphi + 0.5) * 2 * dHalfPhi;
      auto sTransform = Transform3(
          Translation3(hR * std::cos(phi), hR * std::sin(phi), zPosition) *
          AngleAxis3(phi - 0.5 * M_PI, zAxis));
      // Place it
      dSupport.push_back(
          Surface::makeShared<PlaneSurface>(sTransform, sTrapezoid));
    }
  }
  return dSupport;
}

void Acts::Experimental::detail::SupportHelper::addSupport(
    std::vector<std::shared_ptr<Surface>>& layerSurfaces,
    std::vector<size_t>& assignToAll, const Extent& layerExtent,
    Surface::SurfaceType layerRepresentation,
    const std::array<ActsScalar, 5u>& layerSupportValues,
    std::optional<Transform3> layerTransform, unsigned int supportSplits) {
  // Cylinder and Disc section
  if (layerRepresentation == Surface::SurfaceType::Cylinder or
      layerRepresentation == Surface::SurfaceType::Disc) {
    // Bail out if you have no measure of R, Z
    if (not layerExtent.constrains(binZ) or not layerExtent.constrains(binR)) {
      throw std::runtime_error(
          "SupportHelper::addSupport(...) - z or phi are not constrained.");
    }
    /// Overall parameters of support surfaces
    ActsScalar minZ = layerExtent.min(binZ);
    ActsScalar maxZ = layerExtent.max(binZ);
    ActsScalar minR = layerExtent.min(binR);
    ActsScalar maxR = layerExtent.max(binR);
    ActsScalar minPhi = -M_PI;
    ActsScalar maxPhi = M_PI;
    bool sectoral = false;
    bool concentric = false;
    // Check if concentric
    if (layerTransform.has_value() and
        layerTransform.value().isApprox(Transform3::Identity())) {
      concentric = true;
    }
    // Check if we are dealing with a sectoral setup
    if (layerExtent.constrains(binPhi)) {
      minPhi = layerExtent.min(binPhi);
      maxPhi = layerExtent.max(binPhi);
      sectoral = true;
    }

    // Get the main support parameters:
    // - doff .. offset (in r.z)
    // - demin,d emax .. envelop min, max (in z,r)
    // - dphimin, dphimin .. envelop min, max (in phi)
    auto [doff, demin, demax, dphimin, dphimax] = layerSupportValues;
    // phi treatment is common between the cylinders and discs
    if (sectoral) {
      minPhi -= std::abs(demin);
      maxPhi += std::abs(demax);
    }
    // Average phi and half phi
    ActsScalar avgPhi = 0.5 * (maxPhi + minPhi);
    ActsScalar halfPhi = 0.5 * (maxPhi - minPhi);
    // Now specify into Cylinder or disc
    if (layerRepresentation == Surface::SurfaceType::Cylinder) {
      ActsScalar layerR = doff < 0 ? minR + doff : maxR + doff;
      minZ -= std::abs(demin);
      maxZ += std::abs(demax);
      ActsScalar midZ = 0.5 * (minZ + maxZ);
      ActsScalar halfZ = 0.5 * (maxZ - minZ);
      // midZ / halfZ are overwritten if the cylinder
      // is chosen to be concentric
      Transform3 sTransform = Transform3::Identity();
      if (concentric) {
        midZ = 0.;
        halfZ = std::max(std::abs(minZ), std::abs(maxZ));
      } else {
        sTransform.pretranslate(Vector3(0., 0., midZ));
      }
      auto cSupport = SupportHelper::cylindricalSupport(
          sTransform, {layerR, halfZ, halfPhi, avgPhi, 0., 0.}, supportSplits);
      // Remember the surfaces to be assigned to all bins, once the
      // support surfaces are split they enter the standard bin assignment
      if (supportSplits == 1u and cSupport.size() == 1u) {
        assignToAll.push_back(layerSurfaces.size());
      }
      // Add those to the layer surfaces
      layerSurfaces.insert(layerSurfaces.end(), cSupport.begin(),
                           cSupport.end());

    } else {
      // Disc section
      ActsScalar layerZ = doff < 0 ? minZ + doff : maxZ + doff;
      minR -= std::abs(demin);
      maxR += std::abs(demax);
      Transform3 sTransform = Transform3::Identity();
      sTransform.pretranslate(Vector3(0., 0., layerZ));
      auto dSupport = SupportHelper::discSupport(
          sTransform, {minR, maxR, halfPhi, avgPhi}, supportSplits);
      // Remember the surfaces to be assigned to all bins, once the
      // support surfaces are split they enter the standard bin assignment
      if (supportSplits == 1u and dSupport.size() == 1u) {
        assignToAll.push_back(layerSurfaces.size());
      }
      // Add those to the layer surfaces
      layerSurfaces.insert(layerSurfaces.end(), dSupport.begin(),
                           dSupport.end());
    }
  } else {
    throw std::invalid_argument(
        "SupportHelper: currently only cylindrical/disc support building "
        "possible.");
  }
}
