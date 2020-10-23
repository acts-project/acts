// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsFatras/Digitization/PlanarSurfaceMask.hpp"

#include "ActsFatras/Digitization/DigitizationError.hpp"
#include <Acts/Surfaces/AnnulusBounds.hpp>
#include <Acts/Surfaces/DiscTrapezoidBounds.hpp>
#include <Acts/Surfaces/PlanarBounds.hpp>
#include <Acts/Surfaces/RadialBounds.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/Surfaces/detail/IntersectionHelper2D.hpp>
#include <Acts/Utilities/Helpers.hpp>

Acts::Result<std::array<Acts::Vector2D, 2>>
ActsFatras::PlanarSurfaceMask::apply(const Acts::Surface& surface,
                                     const Segment2D& segment) const {
  auto surfaceType = surface.type();
  Segment2D clipped(segment);

  // Plane surface section -------------------
  if (surfaceType == Acts::Surface::Plane or
      surface.bounds().type() == Acts::SurfaceBounds::eDiscTrapezoid) {
    Acts::Vector2D localStart =
        (surfaceType == Acts::Surface::Plane)
            ? segment[0]
            : Acts::Vector2D(Acts::VectorHelpers::perp(segment[0]),
                             Acts::VectorHelpers::phi(segment[0]));

    Acts::Vector2D localEnd =
        (surfaceType == Acts::Surface::Plane)
            ? segment[1]
            : Acts::Vector2D(Acts::VectorHelpers::perp(segment[1]),
                             Acts::VectorHelpers::phi(segment[1]));

    bool startInside = surface.bounds().inside(localStart, true);
    bool endInside = surface.bounds().inside(localEnd, true);

    if (startInside and endInside) {
      return segment;
    }

    // It's either planar or disc trapezoid bounds
    const Acts::PlanarBounds* planarBounds = nullptr;
    const Acts::DiscTrapezoidBounds* dtbBounds = nullptr;
    if (surfaceType == Acts::Surface::Plane) {
      planarBounds =
          static_cast<const Acts::PlanarBounds*>(&(surface.bounds()));
      if (planarBounds->type() == Acts::SurfaceBounds::eEllipse) {
        return DigitizationError::UndefinedSurface;
      }
    } else {
      dtbBounds =
          static_cast<const Acts::DiscTrapezoidBounds*>(&(surface.bounds()));
    }
    auto vertices = planarBounds != nullptr ? planarBounds->vertices(1)
                                            : dtbBounds->vertices(1);
    if (not startInside) {
      auto maskedR = polygonMask(segment[0], segment[1], vertices);
      if (maskedR.ok()) {
        clipped[0] = maskedR.value();
      } else {
        return DigitizationError::MaskingError;
      }
    }
    if (not endInside) {
      auto maskedR = polygonMask(segment[1], segment[0], vertices);
      if (maskedR.ok()) {
        clipped[1] = maskedR.value();
      } else {
        return DigitizationError::MaskingError;
      }
    }
    return clipped;

    // Disc surface section --------------------
  } else if (surfaceType == Acts::Surface::Disc) {
    // Polar coordinates
    Acts::Vector2D sPolar(Acts::VectorHelpers::perp(segment[0]),
                          Acts::VectorHelpers::phi(segment[0]));
    Acts::Vector2D ePolar(Acts::VectorHelpers::perp(segment[1]),
                          Acts::VectorHelpers::phi(segment[1]));

    bool startInside = surface.bounds().inside(sPolar, true);
    bool endInside = surface.bounds().inside(ePolar, true);

    if (startInside and endInside) {
      return segment;
    }

    // DiscTrapezoidalBounds are already excluded
    std::array<double, 2> radialE;
    std::array<std::pair<double, double>, 2> polarE;
    std::array<Segment2D, 2> sectorEdges;
    if (surface.bounds().type() == Acts::SurfaceBounds::eDisc) {
      auto rBounds =
          static_cast<const Acts::RadialBounds*>(&(surface.bounds()));
      double rMin = rBounds->get(Acts::RadialBounds::eMinR);
      double rMax = rBounds->get(Acts::RadialBounds::eMaxR);
      radialE = {rMin, rMax};
      double hPhi = rBounds->get(Acts::RadialBounds::eHalfPhiSector);
      double aPhi = rBounds->get(Acts::RadialBounds::eAveragePhi);
      std::pair<double, double> phiMinMax = {aPhi - hPhi, aPhi + hPhi};
      polarE = {phiMinMax, phiMinMax};
      double cphiMin = std::cos(aPhi - hPhi);
      double sphiMin = std::sin(aPhi - hPhi);
      double cphiMax = std::cos(aPhi + hPhi);
      double sphiMax = std::sin(aPhi + hPhi);
      Segment2D minEdge = {Acts::Vector2D(rMin * cphiMin, rMin * sphiMin),
                           Acts::Vector2D(rMax * cphiMin, rMax * sphiMin)};
      Segment2D maxEdge = {Acts::Vector2D(rMin * cphiMax, rMin * sphiMax),
                           Acts::Vector2D(rMax * cphiMax, rMax * sphiMax)};
      sectorEdges = {minEdge, maxEdge};

    } else if (surface.bounds().type() == Acts::SurfaceBounds::eAnnulus) {
      auto aBounds =
          static_cast<const Acts::AnnulusBounds*>(&(surface.bounds()));
      double rMin = aBounds->rMin();
      double rMax = aBounds->rMax();
      radialE = {rMin, rMax};

      auto polarEdges = aBounds->corners();
      auto polarCenter = aBounds->moduleOrigin();

      std::vector<Acts::Vector2D> cartesianEdges;
      for (const auto& pe : polarEdges) {
        double R = pe[Acts::eBoundLoc0];
        double phi = pe[Acts::eBoundLoc1];
        Acts::Vector2D ce(R * std::cos(phi) - polarCenter.x(),
                          R * std::sin(phi) - polarCenter.y());
        cartesianEdges.push_back(ce);
      }
      std::sort(cartesianEdges.begin(), cartesianEdges.end(),
                [](const Acts::Vector2D& a, const Acts::Vector2D& b) {
                  return (Acts::VectorHelpers::perp(a) <
                          Acts::VectorHelpers::perp(b));
                });
      std::sort(
          cartesianEdges.begin(), cartesianEdges.begin() + 1,
          [](const Acts::Vector2D& a, const Acts::Vector2D& b) {
            return (Acts::VectorHelpers::phi(a) < Acts::VectorHelpers::phi(b));
          });
      std::sort(
          cartesianEdges.begin() + 2, cartesianEdges.begin() + 3,
          [](const Acts::Vector2D& a, const Acts::Vector2D& b) {
            return (Acts::VectorHelpers::phi(a) < Acts::VectorHelpers::phi(b));
          });
      Segment2D minEdge = {cartesianEdges[0], cartesianEdges[3]};
      Segment2D maxEdge = {cartesianEdges[1], cartesianEdges[2]};
      sectorEdges = {minEdge, maxEdge};
      std::pair<double, double> phisMinR = {
          Acts::VectorHelpers::phi(cartesianEdges[0]),
          Acts::VectorHelpers::phi(cartesianEdges[2])};
      std::pair<double, double> phisMaxR = {
          Acts::VectorHelpers::phi(cartesianEdges[1]),
          Acts::VectorHelpers::phi(cartesianEdges[3])};
      polarE = {phisMinR, phisMaxR};
    }

    Acts::detail::IntersectionHelper2D intersector;

    std::array<PolarSegment, 2> segEdges = {
        PolarSegment(segment[0], segment[0], segment[1], startInside),
        PolarSegment(segment[1], segment[1], segment[0], endInside)};

    // Clip to radial bound
    for (auto& se : segEdges) {
      Acts::Intersection2D solution;
      if (not se.inside) {
        Acts::Vector2D sedir = (se.end - se.start).normalized();
        for (const auto& sector : sectorEdges) {
          auto tSolution = intersector.intersectSegment(sector[0], sector[1],
                                                        se.start, sedir, true);
          if (tSolution and tSolution.pathLength > 0. and
              tSolution.pathLength < solution.pathLength) {
            se.clipped = tSolution.position;
            solution = tSolution;
          }
        }
        size_t ir = 0;
        // Check radial first
        for (auto r : radialE) {
          auto tSolution = intersector.intersectCircleSegment(
              r, polarE[ir].first, polarE[ir].second, se.start, sedir);
          if (tSolution and tSolution.pathLength > 0. and
              tSolution.pathLength < solution.pathLength) {
            se.clipped = tSolution.position;
            solution = tSolution;
          }
          ++ir;
        }
        // No clipping happened
        if (not solution) {
          return DigitizationError::MaskingError;
        }
      }
    }
    clipped = {segEdges[0].clipped, segEdges[1].clipped};
    return clipped;
  }
  return DigitizationError::UndefinedSurface;
}

Acts::Result<Acts::Vector2D> ActsFatras::PlanarSurfaceMask::polygonMask(
    const Acts::Vector2D& outside, const Acts::Vector2D& inside,
    const std::vector<Acts::Vector2D>& vertices) const {
  Acts::Intersection2D solution;
  Acts::detail::IntersectionHelper2D intersector;

  for (size_t iv = 0; iv < vertices.size(); ++iv) {
    const Acts::Vector2D& s0 = vertices[iv];
    const Acts::Vector2D& s1 =
        (iv + 1) < vertices.size() ? vertices[iv + 1] : vertices[0];

    Acts::Vector2D lineSegment2D = (inside - outside).normalized();
    auto intersection =
        intersector.intersectSegment(s0, s1, outside, lineSegment2D, true);
    if (intersection and intersection.pathLength < solution.pathLength) {
      solution = intersection;
    }
  }

  if (solution.status == Acts::Intersection2D::Status::reachable) {
    return Acts::Result<Acts::Vector2D>::success(solution.position);
  }
  return Acts::Result<Acts::Vector2D>::failure(DigitizationError::MaskingError);
}
