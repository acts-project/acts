// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Detector/detail/IndexedGridFiller.hpp"
#include "Acts/Detector/detail/IndexedSurfacesGenerator.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Plugins/ActSVG/IndexedSurfacesSvgConverter.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Tests/CommonHelpers/CylindricalTrackingGeometry.hpp"
#include "Acts/Utilities/Grid.hpp"
#include "Acts/Utilities/GridAxisGenerators.hpp"

#include <tuple>

using namespace Acts;
using namespace Acts::Svg;
using namespace Acts::Test;
using namespace Acts::Experimental;
using namespace Acts::Experimental::detail;

GeometryContext tContext;
CylindricalTrackingGeometry cGeometry = CylindricalTrackingGeometry(tContext);

IndexedSurfacesConverter::Options generateDrawOptions() {
  // The converter options
  IndexedSurfacesConverter::Options isOptions;
  // Sensitive surface style
  Style sensitiveStyle;
  sensitiveStyle.fillColor = {51, 153, 255};
  sensitiveStyle.fillOpacity = 0.9;
  sensitiveStyle.highlightColor = {255, 153, 51};
  sensitiveStyle.highlights = {"onmouseover", "onmouseout"};
  sensitiveStyle.strokeWidth = 0.5;
  sensitiveStyle.strokeColor = {0, 0, 0};
  sensitiveStyle.nSegments = 72u;
  std::pair<GeometryIdentifier, Style> allSensitives = {GeometryIdentifier(0u),
                                                        sensitiveStyle};

  // Hierarchy map of styles
  GeometryHierarchyMap<Style> surfaceStyles({allSensitives});
  isOptions.surfaceStyles = surfaceStyles;

  // The grid style
  GridConverter::Options gridOptions;
  Style gridStyle;
  gridStyle.fillOpacity = 0.;
  gridStyle.strokeColor = {0, 0, 255};
  gridStyle.strokeWidth = 1.;
  gridStyle.highlightStrokeWidth = 3;
  gridStyle.highlightStrokeColor = {255, 0, 0};
  gridOptions.style = gridStyle;

  isOptions.gridOptions = gridOptions;
  return isOptions;
}

auto drawOptions = generateDrawOptions();

BOOST_AUTO_TEST_SUITE(ActSvg)

BOOST_AUTO_TEST_CASE(RingDisc1D) {
  // A single ring
  CylindricalTrackingGeometry::DetectorStore dStore;
  auto rSurfaces = cGeometry.surfacesRing(dStore, 6.4, 12.4, 36., 0.125, 0.,
                                          55., 0., 2., 22u);

  IndexedSurfacesGenerator<decltype(rSurfaces), IndexedSurfacesImpl> irSurfaces{
      rSurfaces, {}, {binPhi}};

  GridAxisGenerators::EqClosed aGenerator{{-M_PI, M_PI}, 44u};
  PolyhedronReferenceGenerator<1u, true> rGenerator;

  auto indexedRing = irSurfaces(tContext, aGenerator, rGenerator);
  // The displaying
  auto pIndexedRing = IndexedSurfacesConverter::convert(
      tContext, rSurfaces, indexedRing, drawOptions);
  auto pIndexRingView = View::xy(pIndexedRing, "RingDisc1D");
  toFile({pIndexRingView}, pIndexRingView._id + ".svg");
}

BOOST_AUTO_TEST_CASE(RingDisc1DWithSupport) {
  // A single ring
  CylindricalTrackingGeometry::DetectorStore dStore;
  auto rSurfaces = cGeometry.surfacesRing(dStore, 6.4, 12.4, 36., 0.125, 0.,
                                          55., 0., 2., 22u);

  auto rBounds = std::make_shared<RadialBounds>(20., 20.);
  auto dSurface = Surface::makeShared<DiscSurface>(Transform3::Identity(),
                                                   std::move(rBounds));
  rSurfaces.push_back(dSurface.get());

  IndexedSurfacesGenerator<decltype(rSurfaces), IndexedSurfacesImpl> irSurfaces{
      rSurfaces, {rSurfaces.size() - 1u}, {binPhi}};

  GridAxisGenerators::EqClosed aGenerator{{-M_PI, M_PI}, 44u};
  PolyhedronReferenceGenerator<1u, true> rGenerator;

  auto indexedRing = irSurfaces(tContext, aGenerator, rGenerator);
  // The displaying
  auto pIndexedRing = IndexedSurfacesConverter::convert(
      tContext, rSurfaces, indexedRing, drawOptions);
  auto pIndexRingView = View::xy(pIndexedRing, "RingDisc1DWithSupport");
  toFile({pIndexRingView}, pIndexRingView._id + ".svg");
}

BOOST_AUTO_TEST_CASE(RingDisc2D) {
  // Two rings to make a disc
  CylindricalTrackingGeometry::DetectorStore dStore;
  auto rSurfacesR0 = cGeometry.surfacesRing(dStore, 6.4, 12.4, 18., 0.125, 0.,
                                            42., 0., 2., 22u);

  auto rSurfacesR1 = cGeometry.surfacesRing(dStore, 12.4, 20.4, 30., 0.125, 0.,
                                            80., 0., 2., 22u);

  decltype(rSurfacesR0) rSurfaces = rSurfacesR0;
  rSurfaces.insert(rSurfaces.end(), rSurfacesR1.begin(), rSurfacesR1.end());

  IndexedSurfacesGenerator<decltype(rSurfaces), IndexedSurfacesImpl> irSurfaces{
      rSurfaces, {}, {binR, binPhi}};

  GridAxisGenerators::VarBoundEqClosed aGenerator{
      {24., 74., 110.}, {-M_PI, M_PI}, 44u};
  PolyhedronReferenceGenerator<1u, true> rGenerator;

  auto indexedRing = irSurfaces(tContext, aGenerator, rGenerator);
  // The displaying
  auto pIndexedRing = IndexedSurfacesConverter::convert(
      tContext, rSurfaces, indexedRing, drawOptions);
  auto pIndexRingView = View::xy(pIndexedRing, "RingDisc2D");

  toFile({pIndexRingView}, pIndexRingView._id + ".svg");
}

BOOST_AUTO_TEST_CASE(RingDisc2DFine) {
  // Three rings to make a disc
  CylindricalTrackingGeometry::DetectorStore dStore;
  auto rSurfacesR0 = cGeometry.surfacesRing(dStore, 6.4, 12.4, 18., 0.125, 0.,
                                            42., 0., 2., 22u);

  auto rSurfacesR1 = cGeometry.surfacesRing(dStore, 12.4, 20.4, 30., 0.125, 0.,
                                            80., 0., 2., 22u);

  auto rSurfacesR2 = cGeometry.surfacesRing(dStore, 18.4, 28.4, 30., 0.125, 0.,
                                            122., 0., 2., 36u);

  decltype(rSurfacesR0) rSurfaces = rSurfacesR0;
  rSurfaces.insert(rSurfaces.end(), rSurfacesR1.begin(), rSurfacesR1.end());
  rSurfaces.insert(rSurfaces.end(), rSurfacesR2.begin(), rSurfacesR2.end());

  IndexedSurfacesGenerator<decltype(rSurfaces), IndexedSurfacesImpl> irSurfaces{
      rSurfaces, {}, {binR, binPhi}};

  GridAxisGenerators::EqBoundEqClosed aGenerator{
      {24., 152}, 8u, {-M_PI, M_PI}, 88u};

  PolyhedronReferenceGenerator<1u, true> rGenerator;
  auto indexedRing = irSurfaces(tContext, aGenerator, rGenerator);
  // The displaying
  auto pIndexedRing = IndexedSurfacesConverter::convert(
      tContext, rSurfaces, indexedRing, drawOptions);
  auto pIndexRingView = View::xy(pIndexedRing, "RingDisc2DFine");

  toFile({pIndexRingView}, pIndexRingView._id + ".svg");
}

BOOST_AUTO_TEST_CASE(RingDisc2DFineExpanded) {
  // Three rings to make a disc
  CylindricalTrackingGeometry::DetectorStore dStore;
  auto rSurfacesR0 = cGeometry.surfacesRing(dStore, 6.4, 12.4, 18., 0.125, 0.,
                                            42., 0., 2., 22u);

  auto rSurfacesR1 = cGeometry.surfacesRing(dStore, 12.4, 20.4, 30., 0.125, 0.,
                                            80., 0., 2., 22u);

  auto rSurfacesR2 = cGeometry.surfacesRing(dStore, 18.4, 28.4, 30., 0.125, 0.,
                                            122., 0., 2., 36u);

  decltype(rSurfacesR0) rSurfaces = rSurfacesR0;
  rSurfaces.insert(rSurfaces.end(), rSurfacesR1.begin(), rSurfacesR1.end());
  rSurfaces.insert(rSurfaces.end(), rSurfacesR2.begin(), rSurfacesR2.end());

  IndexedSurfacesGenerator<decltype(rSurfaces), IndexedSurfacesImpl> irSurfaces{
      rSurfaces, {}, {binR, binPhi}, {2u, 4u}};

  GridAxisGenerators::EqBoundEqClosed aGenerator{
      {24., 152}, 8u, {-M_PI, M_PI}, 88u};
  PolyhedronReferenceGenerator<1u, true> rGenerator;

  auto indexedRing = irSurfaces(tContext, aGenerator, rGenerator);
  // The displaying
  auto pIndexedRing = IndexedSurfacesConverter::convert(
      tContext, rSurfaces, indexedRing, drawOptions);
  auto pIndexRingView = View::xy(pIndexedRing, "RingDisc2DFineExpanded");

  toFile({pIndexRingView}, pIndexRingView._id + ".svg");
}

BOOST_AUTO_TEST_CASE(Cylinder2D) {
  CylindricalTrackingGeometry::DetectorStore dStore;
  auto surfaces = cGeometry.surfacesCylinder(dStore, 8.4, 36., 0.15, 0.145,
                                             116., 3., 2., {52, 14});

  IndexedSurfacesGenerator<decltype(surfaces), IndexedSurfacesImpl> icSurfaces{
      surfaces, {}, {binZ, binPhi}, {1u, 1u}};

  GridAxisGenerators::EqBoundEqClosed aGenerator{
      {-500., 500}, 28, {-M_PI, M_PI}, 52u};
  PolyhedronReferenceGenerator<1u, true> rGenerator;

  auto indexedCylinder = icSurfaces(tContext, aGenerator, rGenerator);
  // The displaying
  auto pIndexeCylinder = IndexedSurfacesConverter::convert(
      tContext, surfaces, indexedCylinder, drawOptions);
  auto pIndexCylinderView = View::zphi(pIndexeCylinder, "Cylinder");

  toFile({pIndexCylinderView}, pIndexCylinderView._id + ".svg");
}

BOOST_AUTO_TEST_SUITE_END()
