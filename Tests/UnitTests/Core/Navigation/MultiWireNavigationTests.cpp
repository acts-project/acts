// This file is part of the Acts project.
//
// Copyright (C) 2022-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Detector/Detector.hpp"
#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Detector/MultiWireStructureBuilder.hpp"
#include "Acts/Navigation/DetectorVolumeFinders.hpp"
#include "Acts/Navigation/NavigationState.hpp"
#include "Acts/Navigation/NavigationStateFillers.hpp"
#include "Acts/Navigation/NavigationStateUpdators.hpp"
#include "Acts/Navigation/SurfaceCandidatesUpdators.hpp"
#include "Acts/Plugins/ActSVG/GridSvgConverter.hpp"
#include "Acts/Plugins/ActSVG/SurfaceSvgConverter.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/StrawSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"
#include "Acts/Utilities/detail/Grid.hpp"
#include "actsvg/display/grids.hpp"

#include <fstream>
#include <memory>
#include <string>
#include <vector>

using namespace Acts;
using namespace Acts::Experimental;
using namespace Acts::detail;

GeometryContext tContext;

BOOST_AUTO_TEST_SUITE(Experimental)

BOOST_AUTO_TEST_CASE(Navigation_in_Indexed_Surfaces) {
  using GlobalBin = size_t;
  using LocalBin = std::array<size_t, 2u>;

  std::vector<std::shared_ptr<Acts::Surface>> strawSurfaces = {};

  // svg file to save visualized objects
  actsvg::svg::file gFile;

  // Set the number of surfaces along each dimension of the multi wire structure
  // aligned along z axis
  std::size_t nSurfacesY = 4;
  std::size_t nSurfacesX = 15;

  double radius = 15.;
  double halfZ = 250.;

  // The transform of the 1st surface
  Vector3 ipos = {-0.5 * nSurfacesX * 2 * radius + radius,
                  -0.5 * nSurfacesY * 2 * radius + radius, 0.};

  Vector3 pos = ipos;

  Svg::SurfaceConverter::Options sOptions;

  // Generate the surfaces
  for (std::size_t i = 0; i < nSurfacesY; i++) {
    for (std::size_t j = 0; j < nSurfacesX; j++) {
      auto surface = Surface::makeShared<StrawSurface>(
          Transform3(Translation3(pos)), radius, halfZ);
      strawSurfaces.push_back(surface);
      pos.x() = ipos.x() + 2 * j * radius;
      pos.y() = ipos.y() + 2 * i * radius;

      auto svStraw =
          Svg::SurfaceConverter::convert(tContext, *surface, sOptions);
      auto xyObject = Svg::View::xy(
          svStraw, "svg_surface_" + std::to_string(i) + std::to_string(j));
      gFile.add_object(xyObject);
    }
  }

  std::vector<ActsScalar> vBounds = {0.5 * nSurfacesX * 2 * radius,
                                     0.5 * nSurfacesY * 2 * radius, halfZ};

  MultiWireStructureBuilder::Config mlCfg;
  mlCfg.name = "Multi_Layer_With_Wires";
  mlCfg.mlSurfaces = strawSurfaces;

  mlCfg.mlBinning = {
      ProtoBinning(Acts::binX, Acts::detail::AxisBoundaryType::Bound,
                   -vBounds[0], vBounds[0], nSurfacesX, 1u),
      ProtoBinning(Acts::binY, Acts::detail::AxisBoundaryType::Bound,
                   -vBounds[1], vBounds[1], nSurfacesY, 0u)};
  mlCfg.mlBounds = vBounds;

  MultiWireStructureBuilder mlBuilder(mlCfg);
  auto [volumes, portals, roots] = mlBuilder.construct(tContext);

  Acts::Experimental::NavigationState nState;
  nState.position = Acts::Vector3(0., -60., 0.);
  nState.direction = Acts::Vector3(-1., 1., 0.);

  nState.currentVolume = volumes.front().get();
  nState.currentVolume->updateNavigationState(tContext, nState);

  // check the surface candidates after update (12 surfaces + 6 portals)
  BOOST_CHECK(nState.surfaceCandidates.size() == 18u);

  // Construct the grid for visualization with the track and the highlighted
  // bins
  Axis<AxisType::Equidistant, AxisBoundaryType::Bound> axisX(
      -vBounds[0], vBounds[0], nSurfacesX);
  Axis<AxisType::Equidistant, AxisBoundaryType::Bound> axisY(
      -vBounds[1], vBounds[1], nSurfacesY);
  Grid<std::tuple<GlobalBin, LocalBin>, decltype(axisX), decltype(axisY)>
      gridXY({axisX, axisY});

  Svg::GridConverter::Options cOptions;
  auto pGrid = Svg::GridConverter::convert(gridXY, {binX, binY}, cOptions);
  auto oGrid = actsvg::display::grid("MultiWireGrid", pGrid);
  gFile.add_object(oGrid);

  actsvg::style::fill fill_style{{{200, 150, 50}}};
  actsvg::style::stroke stroke_style{{{200, 150, 50}}};
  actsvg::style::stroke stroke_black = actsvg::style::stroke();

  Acts::Svg::Style candidates_style;
  candidates_style.fillColor = {200, 100, 100};

  sOptions.style = candidates_style;

  for (const auto sfc : nState.surfaceCandidates) {
    if (!sfc.portal) {
      auto activatedBin = actsvg::draw::rectangle(
          "bin",
          {static_cast<float>(sfc.surface->center(tContext).x()),
           static_cast<float>(sfc.surface->center(tContext).y())},
          0.5 * axisX.getBinWidth(1), 0.5 * axisY.getBinWidth(1), fill_style,
          stroke_black);
      auto pcSurf =
          Svg::SurfaceConverter::convert(tContext, *sfc.surface, sOptions);

      gFile.add_object(activatedBin);
      gFile.add_object(Svg::View::xy(pcSurf, "svg_surface_candidate"));
    }
  }

  auto dy = axisY.getMax() - nState.position.y();
  auto tanphi = tan(Acts::VectorHelpers::phi(nState.direction));
  auto dx = dy / tanphi;

  actsvg::point2 start = {static_cast<float>(nState.position.x()),
                          static_cast<float>(nState.position.y())};
  actsvg::point2 end = {static_cast<float>(nState.position.x() + dx),
                        static_cast<float>(nState.position.y() + dy)};

  auto track_line = actsvg::draw::line("track_line", start, end);
  gFile.add_object(track_line);

  std::ofstream gFile_ofstream;
  gFile_ofstream.open("grid.svg");
  gFile_ofstream << gFile;
  gFile_ofstream.close();
}

BOOST_AUTO_TEST_SUITE_END()