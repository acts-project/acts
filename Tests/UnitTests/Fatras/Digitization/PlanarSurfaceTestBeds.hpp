// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Surfaces/AnnulusBounds.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/DiscTrapezoidBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/MathHelpers.hpp"
#include "Acts/Utilities/ProtoAxis.hpp"

#include <numbers>
#include <tuple>
#include <vector>

#include "BoundRandomValues.hpp"

namespace ActsFatras {

using Randomizer = std::function<Acts::Vector2(double, double)>;

using PlanarTestBed =
    std::tuple<std::string, std::shared_ptr<const Acts::Surface>,
               std::vector<Acts::DirectedProtoAxis>, Randomizer>;

/// Helper struct to create a testbed for Digitization steps
struct PlanarSurfaceTestBeds {
  /// Call operator for creating a testbed of different surfaces.
  /// It returns a testbed for all planar surface types and the
  /// corresponding cartesian/polar segmentation.
  ///
  /// @param rScale is a parameter how far the random numbers
  /// should be generated (1 -> inside to boundary, 1.1 -> 10% outside)
  std::vector<PlanarTestBed> operator()(double rScale) const {
    double irScale = (2. - rScale);

    // Pixel test in Rectangle
    double xhalf = 3.;
    double yhalf = 6.5;
    auto rectangle = std::make_shared<Acts::RectangleBounds>(xhalf, yhalf);
    auto rSurface = Acts::Surface::makeShared<Acts::PlaneSurface>(
        Acts::Transform3::Identity(), rectangle);
    std::vector<Acts::DirectedProtoAxis> pixelated;
    pixelated.emplace_back(Acts::AxisDirection::AxisX,
                           Acts::AxisBoundaryType::Open, -xhalf, xhalf, 15);
    pixelated.emplace_back(Acts::AxisDirection::AxisY,
                           Acts::AxisBoundaryType::Open, -yhalf, yhalf, 26);
    RectangleRandom rRandom(xhalf * rScale, yhalf * rScale);

    // Cartesian strip test in Trapezoid
    double xhalfminy = 2.;
    double xhalfmaxy = 3.5;
    yhalf = 4.;
    auto trapezoid =
        std::make_shared<Acts::TrapezoidBounds>(xhalfminy, xhalfmaxy, yhalf);
    auto tSurface = Acts::Surface::makeShared<Acts::PlaneSurface>(
        Acts::Transform3::Identity(), trapezoid);
    std::vector<Acts::DirectedProtoAxis> stripsX;
    stripsX.emplace_back(Acts::AxisDirection::AxisX,
                         Acts::AxisBoundaryType::Open, -xhalfmaxy, xhalfmaxy,
                         35);
    stripsX.emplace_back(Acts::AxisDirection::AxisY,
                         Acts::AxisBoundaryType::Open, -yhalf, yhalf, 1);
    TrapezoidRandom tRandom(xhalfminy * rScale, xhalfmaxy * rScale,
                            yhalf * rScale);

    // Phi strip test in DiscTrapezoid
    double rmin = 2.;
    double rmax = 7.5;
    double xmin = 2.;
    double xmax = 3.5;
    double ymax = Acts::fastCathetus(rmax, xmax);
    double alpha = std::max(std::atan2(xmin, rmin), std::atan2(xmax, ymax));

    auto discTrapezoid =
        std::make_shared<Acts::DiscTrapezoidBounds>(xmin, xmax, rmin, rmax);
    auto dtSurface = Acts::Surface::makeShared<Acts::DiscSurface>(
        Acts::Transform3::Identity(), discTrapezoid);
    std::vector<Acts::DirectedProtoAxis> stripsPhi;
    stripsPhi.emplace_back(Acts::AxisDirection::AxisR,
                           Acts::AxisBoundaryType::Open, rmin, rmax, 1);
    stripsPhi.emplace_back(
        Acts::AxisDirection::AxisPhi, Acts::AxisBoundaryType::Open,
        std::numbers::pi / 2. - alpha, std::numbers::pi / 2. + alpha, 25);
    TrapezoidRandom dtRandom(xmin * rScale, xmax * rScale, rmin * irScale,
                             ymax * rScale);

    // Raidal disc test
    auto discRadial = std::make_shared<Acts::RadialBounds>(
        rmin, rmax, std::numbers::pi / 4., std::numbers::pi / 2.);
    auto dSurface = Acts::Surface::makeShared<Acts::DiscSurface>(
        Acts::Transform3::Identity(), discRadial);
    std::vector<Acts::DirectedProtoAxis> rphiseg;
    rphiseg.emplace_back(Acts::AxisDirection::AxisR,
                         Acts::AxisBoundaryType::Open, rmin, rmax, 10);
    rphiseg.emplace_back(Acts::AxisDirection::AxisPhi,
                         Acts::AxisBoundaryType::Open,
                         (std::numbers::pi / 2. - std::numbers::pi / 4.),
                         (std::numbers::pi / 2. + std::numbers::pi / 4.), 20);

    DiscRandom dRandom(
        rmin * irScale, rmax * rScale,
        (std::numbers::pi / 2. - std::numbers::pi / 4.) * irScale,
        (std::numbers::pi / 2. + std::numbers::pi / 4.) * rScale);

    // Annulus disc test
    rmin = 2.5;
    rmax = 5.5;
    Acts::Vector2 aorigin(0.1, -0.3);
    double phimin = -0.25;
    double phimax = 0.38;
    auto annulus = std::make_shared<Acts::AnnulusBounds>(rmin, rmax, phimin,
                                                         phimax, aorigin);
    auto aSurface = Acts::Surface::makeShared<Acts::DiscSurface>(
        Acts::Transform3::Identity() *
            Acts::Translation3(-aorigin.x(), -aorigin.y(), 0.),
        annulus);

    auto vertices = annulus->vertices(72);
    std::ranges::for_each(vertices, [&](Acts::Vector2& v) {
      double r = Acts::VectorHelpers::perp(v);
      rmin = std::min(rmin, r);
      rmax = std::max(rmax, r);
    });

    std::vector<Acts::DirectedProtoAxis> stripsPhiA;
    stripsPhiA.emplace_back(Acts::AxisDirection::AxisR,
                            Acts::AxisBoundaryType::Open, rmin, rmax, 1);
    stripsPhiA.emplace_back(Acts::AxisDirection::AxisPhi,
                            Acts::AxisBoundaryType::Open, phimin, phimax, 12);
    AnnulusRandom aRandom(rmin * irScale, rmax * rScale, phimin * rScale,
                          phimax * rScale, aorigin.x(), aorigin.y());

    return {{"Rectangle", std::move(rSurface), pixelated, rRandom},
            {"Trapezoid", std::move(tSurface), stripsX, tRandom},
            {"DiscTrapezoid", std::move(dtSurface), stripsPhi, dtRandom},
            {"DiscRadial", std::move(dSurface), rphiseg, dRandom},
            {"Annulus", std::move(aSurface), stripsPhiA, aRandom}};
  }
};

}  // namespace ActsFatras
