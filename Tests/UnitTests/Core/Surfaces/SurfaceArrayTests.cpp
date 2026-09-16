// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"
#include "Acts/Utilities/Axis.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/Diagnostics.hpp"
#include "Acts/Utilities/Helpers.hpp"

#include <algorithm>
#include <array>
#include <cfenv>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numbers>
#include <ranges>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <boost/format.hpp>

using Acts::VectorHelpers::phi;

using namespace Acts;

namespace ActsTests {

// Create a test context
GeometryContext tgContext = GeometryContext::dangerouslyDefaultConstruct();

using SrfVec = std::vector<std::shared_ptr<const Surface>>;
struct SurfaceArrayFixture {
  std::vector<std::shared_ptr<const Surface>> m_surfaces;

  SurfaceArrayFixture() { BOOST_TEST_MESSAGE("setup fixture"); }
  ~SurfaceArrayFixture() { BOOST_TEST_MESSAGE("teardown fixture"); }

  SrfVec fullPhiTestSurfacesEC(std::size_t n = 10, double shift = 0,
                               double zbase = 0, double r = 10) {
    SrfVec res;

    double phiStep = 2 * std::numbers::pi / n;
    for (std::size_t i = 0; i < n; ++i) {
      double z = zbase + ((i % 2 == 0) ? 1 : -1) * 0.2;

      Transform3 trans;
      trans.setIdentity();
      trans.rotate(Eigen::AngleAxisd(i * phiStep + shift, Vector3(0, 0, 1)));
      trans.translate(Vector3(r, 0, z));

      auto bounds = std::make_shared<const RectangleBounds>(2, 1);
      std::shared_ptr<const Surface> srf =
          Surface::makeShared<PlaneSurface>(trans, bounds);

      res.push_back(srf);
      m_surfaces.push_back(
          std::move(srf));  // keep shared, will get destroyed at the end
    }

    return res;
  }

  SrfVec fullPhiTestSurfacesBRL(int n = 10, double shift = 0, double zbase = 0,
                                double incl = std::numbers::pi / 9.,
                                double w = 2, double h = 1.5) {
    SrfVec res;

    double phiStep = 2 * std::numbers::pi / n;
    for (int i = 0; i < n; ++i) {
      double z = zbase;

      Transform3 trans;
      trans.setIdentity();
      trans.rotate(Eigen::AngleAxisd(i * phiStep + shift, Vector3(0, 0, 1)));
      trans.translate(Vector3(10, 0, z));
      trans.rotate(Eigen::AngleAxisd(incl, Vector3(0, 0, 1)));
      trans.rotate(Eigen::AngleAxisd(std::numbers::pi / 2., Vector3(0, 1, 0)));

      auto bounds = std::make_shared<const RectangleBounds>(w, h);
      std::shared_ptr<const Surface> srf =
          Surface::makeShared<PlaneSurface>(trans, bounds);

      res.push_back(srf);
      m_surfaces.push_back(
          std::move(srf));  // keep shared, will get destroyed at the end
    }

    return res;
  }

  SrfVec straightLineSurfaces(
      std::size_t n = 10., double step = 3, const Vector3& origin = {0, 0, 1.5},
      const Transform3& pretrans = Transform3::Identity(),
      const Vector3& dir = {0, 0, 1}) {
    SrfVec res;
    for (std::size_t i = 0; i < n; ++i) {
      Transform3 trans;
      trans.setIdentity();
      trans.translate(origin + dir * step * i);
      // trans.rotate(AngleAxis3(std::numbers::pi/9., Vector3(0, 0, 1)));
      trans.rotate(AngleAxis3(std::numbers::pi / 2., Vector3(1, 0, 0)));
      trans = trans * pretrans;

      auto bounds = std::make_shared<const RectangleBounds>(2, 1.5);

      std::shared_ptr<const Surface> srf =
          Surface::makeShared<PlaneSurface>(trans, bounds);

      res.push_back(srf);
      m_surfaces.push_back(
          std::move(srf));  // keep shared, will get destroyed at the end
    }

    return res;
  }

  SrfVec makeBarrel(int nPhi, int nZ, double w, double h) {
    double z0 = -(nZ - 1) * w;
    SrfVec res;

    for (int i = 0; i < nZ; i++) {
      double z = i * w * 2 + z0;
      SrfVec ring =
          fullPhiTestSurfacesBRL(nPhi, 0, z, std::numbers::pi / 9., w, h);
      res.insert(res.end(), ring.begin(), ring.end());
    }

    return res;
  }

  void draw_surfaces(const SrfVec& surfaces, const std::string& fname) {
    std::ofstream os;
    os.open(fname);

    os << std::fixed << std::setprecision(4);

    std::size_t nVtx = 0;
    for (const auto& srfx : surfaces) {
      std::shared_ptr<const PlaneSurface> srf =
          std::dynamic_pointer_cast<const PlaneSurface>(srfx);
      const PlanarBounds* bounds =
          dynamic_cast<const PlanarBounds*>(&srf->bounds());

      for (const auto& vtxloc : bounds->vertices()) {
        Vector3 vtx = srf->localToGlobalTransform(tgContext) *
                      Vector3(vtxloc.x(), vtxloc.y(), 0);
        os << "v " << vtx.x() << " " << vtx.y() << " " << vtx.z() << "\n";
      }

      // connect them
      os << "f";
      for (std::size_t i = 1; i <= bounds->vertices().size(); ++i) {
        os << " " << nVtx + i;
      }
      os << "\n";

      nVtx += bounds->vertices().size();
    }

    os.close();
  }
};

BOOST_AUTO_TEST_SUITE(SurfacesSuite)

BOOST_FIXTURE_TEST_CASE(SurfaceArray_create, SurfaceArrayFixture) {
  GeometryContext tgContext = GeometryContext::dangerouslyDefaultConstruct();

  SrfVec brl = makeBarrel(30, 7, 2, 1);
  std::vector<const Surface*> brlRaw = unpackSmartPointers(brl);
  draw_surfaces(brl, "SurfaceArray_create_BRL_1.obj");

  Axis<AxisType::Equidistant, AxisBoundaryType::Closed> phiAxis(
      -std::numbers::pi, std::numbers::pi, 30u);
  Axis<AxisType::Equidistant, AxisBoundaryType::Bound> zAxis(-14, 14, 7u);

  double R = 10;
  auto itransform = [R](const Vector2& loc) {
    return Vector3(R * std::cos(loc[0]), R * std::sin(loc[0]), loc[1]);
  };

  auto cylinder =
      Surface::makeShared<CylinderSurface>(Transform3::Identity(), R, 10);
  SurfaceArray sa(tgContext, brl, cylinder, 1., std::tuple{phiAxis, zAxis});

  // let's see if we can access all surfaces
  sa.toStream(tgContext, std::cout);

  // the bin a surface's own centre falls into has to hold it
  for (const auto& srf : brl) {
    const Vector3 ctr = srf->referencePosition(tgContext, AxisDirection::AxisR);
    const Vector3 normal = srf->normal(tgContext, ctr, Vector3::UnitZ());
    const auto binContent = sa.at(tgContext, ctr, normal);

    BOOST_CHECK(std::ranges::find(binContent, srf.get()) != binContent.end());
  }

  const Vector3 crossing = itransform(Vector2(0, 0));
  const auto axes = sa.getAxes();
  const std::array<std::size_t, 2> crossingBins{axes.at(0)->getBin(0.),
                                                axes.at(1)->getBin(0.)};

  // at normal incidence the track does not slide, so the lookup is the bin
  BOOST_CHECK(std::ranges::equal(
      sa.neighbors(tgContext, crossing, crossing.normalized()),
      sa.neighbors(crossingBins, {0, 0})));

  // an inclined track slides along z only, so only z widens. the z bins are 4
  // wide against a tolerance of 1: a slope of 3 reaches one bin, 7 two.
  for (const auto& [slope, distance] :
       std::vector<std::pair<double, std::uint8_t>>{{3., 1}, {7., 2}}) {
    const Vector3 direction = Vector3(1, 0, slope).normalized();
    BOOST_CHECK(std::ranges::equal(sa.neighbors(tgContext, crossing, direction),
                                   sa.neighbors(crossingBins, {0, distance})));
    BOOST_CHECK(
        !std::ranges::equal(sa.neighbors(tgContext, crossing, direction),
                            sa.neighbors(crossingBins, {distance, distance})));
  }

  // the floor is served regardless of the crossing angle
  SurfaceArray floored(tgContext, brl, cylinder, 1., std::tuple{phiAxis, zAxis},
                       {{1, 1}, {1, 2}});
  BOOST_CHECK(std::ranges::equal(
      floored.neighbors(tgContext, crossing, crossing.normalized()),
      floored.neighbors(crossingBins, {1, 1})));
  BOOST_CHECK_THROW(SurfaceArray(tgContext, brl, cylinder, 1.,
                                 std::tuple{phiAxis, zAxis}, {{2, 2}, {1, 2}}),
                    std::invalid_argument);

  // nothing is cached past the bound, so asking for it is an error
  BOOST_CHECK_THROW(floored.neighbors(crossingBins, {2, 0}), std::out_of_range);
  BOOST_CHECK_THROW(floored.neighbors(crossingBins, {0, 3}), std::out_of_range);
  BOOST_CHECK_THROW(floored.neighbors({10000, 0}, {0, 0}), std::out_of_range);

  // a scalar bound is the isotropic window it used to describe
  ACTS_PUSH_IGNORE_DEPRECATED()
  const SurfaceArray scalarBound(tgContext, brl, cylinder, 1.,
                                 std::tuple{phiAxis, zAxis}, std::uint8_t{1});
  BOOST_CHECK_EQUAL(scalarBound.maxNeighborDistance(), 1u);
  ACTS_POP_IGNORE_DEPRECATED()
  const SurfaceArray::NeighborWindow scalarWindow =
      scalarBound.neighborWindow();
  BOOST_CHECK((scalarWindow.min == std::array<std::uint8_t, 2>{1, 1}));
  BOOST_CHECK((scalarWindow.max == std::array<std::uint8_t, 2>{1, 1}));
  BOOST_CHECK(std::ranges::equal(
      scalarBound.neighbors(tgContext, crossing, crossing.normalized()),
      scalarBound.neighbors(crossingBins, {1, 1})));
}

BOOST_AUTO_TEST_CASE(SurfaceArray_overfill) {
  const auto representative = Surface::makeShared<PlaneSurface>(
      Transform3::Identity(), std::make_shared<RectangleBounds>(3.5, 3.5));
  // The footprint covers three cells, so expansion has overlapping
  // contributions.
  const auto module = Surface::makeShared<PlaneSurface>(
      Transform3::Identity(), std::make_shared<RectangleBounds>(0.6, 0.2));
  const auto check = [&]<AxisBoundaryType boundary>() {
    const Axis<AxisType::Equidistant, boundary> axis(-3.5, 3.5, 7);
    for (const std::uint8_t radius : {0, 1, 2, 3, 255}) {
      const SurfaceArray array(tgContext, {module}, representative, 0.,
                               {axis, axis}, {{0, 0}, {0, 0}}, radius);
      for (int x = -3; x <= 3; ++x) {
        for (int y = -3; y <= 3; ++y) {
          const auto content =
              array.at(tgContext, Vector3(x, y, 0), Vector3::UnitZ());
          const bool expected =
              std::abs(x) <= radius + 1 && std::abs(y) <= radius;
          BOOST_CHECK_EQUAL(content.size(), expected ? 1u : 0u);
          if (expected) {
            BOOST_CHECK_EQUAL(content.front(), module.get());
          }
        }
      }
      for (std::size_t bin = 0; bin < array.size(); ++bin) {
        if (!array.isValidBin(bin)) {
          BOOST_CHECK(array.at(bin).empty());
        }
      }
      BOOST_CHECK(!array.isValidBin(array.size() + 10));
    }
  };
  check.template operator()<AxisBoundaryType::Bound>();
  check.template operator()<AxisBoundaryType::Open>();
}

BOOST_FIXTURE_TEST_CASE(SurfaceArray_overfillPhiSeam, SurfaceArrayFixture) {
  const auto modules = makeBarrel(12, 3, 2, 0.2);
  const auto cylinder =
      Surface::makeShared<CylinderSurface>(Transform3::Identity(), 10., 10.);
  const Axis<AxisType::Variable, AxisBoundaryType::Closed> phiAxis(
      {-std::numbers::pi, -3., -2., -1., 0., 1., 2., 3., std::numbers::pi});
  const Axis<AxisType::Equidistant, AxisBoundaryType::Bound> zAxis(-6., 6., 3);
  const SurfaceArray nominal(tgContext, modules, cylinder, 1., {phiAxis, zAxis},
                             {{0, 0}, {3, 3}});
  for (const std::uint8_t radius : {0, 1, 2, 3}) {
    const SurfaceArray expanded(tgContext, modules, cylinder, 1.,
                                {phiAxis, zAxis}, {{0, 0}, {0, 0}}, radius);
    for (std::size_t phiBin = 1; phiBin <= phiAxis.getNBins(); ++phiBin) {
      for (std::size_t zBin = 1; zBin <= zAxis.getNBins(); ++zBin) {
        BOOST_CHECK(
            std::ranges::equal(expanded.neighbors({phiBin, zBin}, 0),
                               nominal.neighbors({phiBin, zBin}, radius)));
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(SurfaceArray_maximumWindow) {
  const auto plane = Surface::makeShared<PlaneSurface>(
      Transform3::Identity(), std::make_shared<RectangleBounds>(1., 1.));
  const Axis<AxisType::Equidistant, AxisBoundaryType::Bound> axis(-1., 1., 1);
  // Both counters must terminate at the largest representable distance.
  const SurfaceArray array(tgContext, {plane}, plane, 0., {axis, axis},
                           {{0, 0}, {255, 255}});
  BOOST_CHECK_EQUAL(array.neighbors({1, 1}, {255, 255}).size(), 1u);
}

BOOST_FIXTURE_TEST_CASE(SurfaceArray_periodicWindow, SurfaceArrayFixture) {
  const auto modules = makeBarrel(30, 1, 2, 0.2);
  const auto cylinder =
      Surface::makeShared<CylinderSurface>(Transform3::Identity(), 10., 10.);
  const Axis<AxisType::Equidistant, AxisBoundaryType::Closed> phiAxis(
      -std::numbers::pi, std::numbers::pi, 30);
  const Axis<AxisType::Equidistant, AxisBoundaryType::Bound> zAxis(-3., 3., 1);
  const SurfaceArray array(tgContext, modules, cylinder, 1., {phiAxis, zAxis},
                           {{0, 0}, {2, 0}});
  const double phi = 0.1;
  const Vector3 crossing(10. * std::cos(phi), 10. * std::sin(phi), 0.);
  for (const double slide : {0., 2. * std::numbers::pi - 0.01,
                             2. * std::numbers::pi, 4. * std::numbers::pi}) {
    const Vector3 direction =
        Vector3(std::cos(phi) - 10. * slide * std::sin(phi),
                std::sin(phi) + 10. * slide * std::cos(phi), 0.)
            .normalized();
    const auto expected =
        array.neighbors({phiAxis.getBin(phi), 1},
                        {static_cast<std::uint8_t>(slide == 0. ? 0 : 2), 0});
    BOOST_CHECK(std::ranges::equal(
        array.neighbors(tgContext, crossing, direction), expected));
    BOOST_CHECK(std::ranges::equal(
        array.neighbors(tgContext, crossing, -direction), expected));
  }
}

BOOST_FIXTURE_TEST_CASE(SurfaceArray_variablePeriodicWindow,
                        SurfaceArrayFixture) {
  const auto modules = fullPhiTestSurfacesEC(30);
  const auto disc =
      Surface::makeShared<DiscSurface>(Transform3::Identity(), 0., 20.);
  const Axis<AxisType::Equidistant, AxisBoundaryType::Bound> rAxis(0., 20., 1);
  const Axis<AxisType::Variable, AxisBoundaryType::Closed> phiAxis(
      {-std::numbers::pi, -3., -2.9, -2.8, -2.7, -1., 0., 1., 2., 3.,
       std::numbers::pi});
  const SurfaceArray array(tgContext, modules, disc, 1., {rAxis, phiAxis},
                           {{0, 0}, {0, 4}});
  const double phi = 3.05;
  const Vector3 crossing(10. * std::cos(phi), 10. * std::sin(phi), 0.);
  const Vector3 direction =
      Vector3(-8. * std::sin(phi), 8. * std::cos(phi), 1.).normalized();
  BOOST_CHECK(
      std::ranges::equal(array.neighbors(tgContext, crossing, direction),
                         array.neighbors({1, phiAxis.getBin(phi)}, {0, 4})));
  // At the polar singularity, an undefined derivative must not enter getBin.
  std::feclearexcept(FE_INVALID | FE_DIVBYZERO);
  BOOST_CHECK(std::ranges::equal(
      array.neighbors(tgContext, Vector3::Zero(), Vector3::UnitZ()),
      array.neighbors({1, phiAxis.getBin(0.)}, {0, 4})));
  BOOST_CHECK(
      std::ranges::equal(array.neighbors(tgContext, Vector3::Zero(), direction),
                         array.neighbors({1, phiAxis.getBin(0.)}, {0, 4})));
  BOOST_CHECK_EQUAL(std::fetestexcept(FE_INVALID | FE_DIVBYZERO), 0);
}

BOOST_AUTO_TEST_CASE(SurfaceArray_singleElement) {
  const double w = 3;
  const double h = 4;
  const auto bounds = std::make_shared<const RectangleBounds>(w, h);
  auto srf = Surface::makeShared<PlaneSurface>(Transform3::Identity(), bounds);

  SurfaceArray sa(srf);

  const auto binContent =
      sa.at(tgContext, Vector3(42, 42, 42), Vector3::UnitX());
  BOOST_CHECK_EQUAL(binContent.size(), 1u);
  BOOST_CHECK_EQUAL(binContent[0], srf.get());
  BOOST_CHECK_EQUAL(sa.surfaces().size(), 1u);
  BOOST_CHECK_EQUAL(sa.surfaces().at(0), srf.get());
}

BOOST_AUTO_TEST_CASE(SurfaceArrayToStreamPreservesStreamState) {
  const auto bounds = std::make_shared<const RectangleBounds>(3., 4.);
  auto surface =
      Surface::makeShared<PlaneSurface>(Transform3::Identity(), bounds);
  SurfaceArray surfaceArray(surface);

  std::ostringstream stream;
  stream << std::scientific << std::showpos << std::setfill('#')
         << std::setprecision(3);
  stream.width(17);

  const auto flags = stream.flags();
  const auto precision = stream.precision();
  const auto width = stream.width();
  const auto fill = stream.fill();

  surfaceArray.toStream(tgContext, stream);

  BOOST_CHECK(stream.flags() == flags);
  BOOST_CHECK_EQUAL(stream.precision(), precision);
  BOOST_CHECK_EQUAL(stream.width(), width);
  BOOST_CHECK_EQUAL(stream.fill(), fill);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
