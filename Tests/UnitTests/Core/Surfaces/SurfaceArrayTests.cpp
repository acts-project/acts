// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/detail/Axis.hpp"
#include "Acts/Utilities/detail/AxisFwd.hpp"
#include "Acts/Utilities/detail/Grid.hpp"
#include "Acts/Utilities/detail/grid_helper.hpp"

#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <boost/format.hpp>

using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;

namespace bdata = boost::unit_test::data;
namespace tt = boost::test_tools;

namespace Acts {

namespace Test {

// Create a test context
GeometryContext tgContext = GeometryContext();

using SrfVec = std::vector<std::shared_ptr<const Surface>>;
struct SurfaceArrayFixture {
  std::vector<std::shared_ptr<const Surface>> m_surfaces;

  SurfaceArrayFixture() { BOOST_TEST_MESSAGE("setup fixture"); }
  ~SurfaceArrayFixture() { BOOST_TEST_MESSAGE("teardown fixture"); }

  SrfVec fullPhiTestSurfacesEC(size_t n = 10, double shift = 0,
                               double zbase = 0, double r = 10) {
    SrfVec res;

    double phiStep = 2 * M_PI / n;
    for (size_t i = 0; i < n; ++i) {
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
                                double incl = M_PI / 9., double w = 2,
                                double h = 1.5) {
    SrfVec res;

    double phiStep = 2 * M_PI / n;
    for (int i = 0; i < n; ++i) {
      double z = zbase;

      Transform3 trans;
      trans.setIdentity();
      trans.rotate(Eigen::AngleAxisd(i * phiStep + shift, Vector3(0, 0, 1)));
      trans.translate(Vector3(10, 0, z));
      trans.rotate(Eigen::AngleAxisd(incl, Vector3(0, 0, 1)));
      trans.rotate(Eigen::AngleAxisd(M_PI / 2., Vector3(0, 1, 0)));

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
      size_t n = 10., double step = 3, const Vector3& origin = {0, 0, 1.5},
      const Transform3& pretrans = Transform3::Identity(),
      const Vector3& dir = {0, 0, 1}) {
    SrfVec res;
    for (size_t i = 0; i < n; ++i) {
      Transform3 trans;
      trans.setIdentity();
      trans.translate(origin + dir * step * i);
      // trans.rotate(AngleAxis3(M_PI/9., Vector3(0, 0, 1)));
      trans.rotate(AngleAxis3(M_PI / 2., Vector3(1, 0, 0)));
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
      // std::cout << "z=" << z << std::endl;
      SrfVec ring = fullPhiTestSurfacesBRL(nPhi, 0, z, M_PI / 9., w, h);
      res.insert(res.end(), ring.begin(), ring.end());
    }

    return res;
  }

  void draw_surfaces(const SrfVec& surfaces, const std::string& fname) {
    std::ofstream os;
    os.open(fname);

    os << std::fixed << std::setprecision(4);

    size_t nVtx = 0;
    for (const auto& srfx : surfaces) {
      std::shared_ptr<const PlaneSurface> srf =
          std::dynamic_pointer_cast<const PlaneSurface>(srfx);
      const PlanarBounds* bounds =
          dynamic_cast<const PlanarBounds*>(&srf->bounds());

      for (const auto& vtxloc : bounds->vertices()) {
        Vector3 vtx =
            srf->transform(tgContext) * Vector3(vtxloc.x(), vtxloc.y(), 0);
        os << "v " << vtx.x() << " " << vtx.y() << " " << vtx.z() << "\n";
      }

      // connect them
      os << "f";
      for (size_t i = 1; i <= bounds->vertices().size(); ++i) {
        os << " " << nVtx + i;
      }
      os << "\n";

      nVtx += bounds->vertices().size();
    }

    os.close();
  }
};

BOOST_AUTO_TEST_SUITE(Surfaces)

BOOST_FIXTURE_TEST_CASE(SurfaceArray_create, SurfaceArrayFixture) {
  GeometryContext tgContext = GeometryContext();

  SrfVec brl = makeBarrel(30, 7, 2, 1);
  std::vector<const Surface*> brlRaw = unpack_shared_vector(brl);
  draw_surfaces(brl, "SurfaceArray_create_BRL_1.obj");

  detail::Axis<detail::AxisType::Equidistant, detail::AxisBoundaryType::Closed>
      phiAxis(-M_PI, M_PI, 30u);
  detail::Axis<detail::AxisType::Equidistant, detail::AxisBoundaryType::Bound>
      zAxis(-14, 14, 7u);

  double angleShift = 2 * M_PI / 30. / 2.;
  auto transform = [angleShift](const Vector3& pos) {
    return Vector2(phi(pos) + angleShift, pos.z());
  };
  double R = 10;
  auto itransform = [angleShift, R](const Vector2& loc) {
    return Vector3(R * std::cos(loc[0] - angleShift),
                   R * std::sin(loc[0] - angleShift), loc[1]);
  };

  auto sl = std::make_unique<
      SurfaceArray::SurfaceGridLookup<decltype(phiAxis), decltype(zAxis)>>(
      transform, itransform,
      std::make_tuple(std::move(phiAxis), std::move(zAxis)));
  sl->fill(tgContext, brlRaw);
  SurfaceArray sa(std::move(sl), brl);

  // let's see if we can access all surfaces
  sa.toStream(tgContext, std::cout);

  for (const auto& srf : brl) {
    Vector3 ctr = srf->binningPosition(tgContext, binR);
    std::vector<const Surface*> binContent = sa.at(ctr);

    BOOST_CHECK_EQUAL(binContent.size(), 1u);
    BOOST_CHECK_EQUAL(srf.get(), binContent.at(0));
  }

  std::vector<const Surface*> neighbors =
      sa.neighbors(itransform(Vector2(0, 0)));
  BOOST_CHECK_EQUAL(neighbors.size(), 9u);

  auto sl2 = std::make_unique<
      SurfaceArray::SurfaceGridLookup<decltype(phiAxis), decltype(zAxis)>>(
      transform, itransform,
      std::make_tuple(std::move(phiAxis), std::move(zAxis)));
  // do NOT fill, only completebinning
  sl2->completeBinning(tgContext, brlRaw);
  SurfaceArray sa2(std::move(sl2), brl);
  sa.toStream(tgContext, std::cout);
  for (const auto& srf : brl) {
    Vector3 ctr = srf->binningPosition(tgContext, binR);
    std::vector<const Surface*> binContent = sa2.at(ctr);

    BOOST_CHECK_EQUAL(binContent.size(), 1u);
    BOOST_CHECK_EQUAL(srf.get(), binContent.at(0));
  }
}

BOOST_AUTO_TEST_CASE(SurfaceArray_singleElement) {
  double w = 3, h = 4;
  auto bounds = std::make_shared<const RectangleBounds>(w, h);
  auto srf = Surface::makeShared<PlaneSurface>(Transform3::Identity(), bounds);

  SurfaceArray sa(srf);

  auto binContent = sa.at(Vector3(42, 42, 42));
  BOOST_CHECK_EQUAL(binContent.size(), 1u);
  BOOST_CHECK_EQUAL(binContent.at(0), srf.get());
  BOOST_CHECK_EQUAL(sa.surfaces().size(), 1u);
  BOOST_CHECK_EQUAL(sa.surfaces().at(0), srf.get());
}

BOOST_AUTO_TEST_CASE(SurfaceArray_manyElementsSingleLookup) {
  double w = 3, h = 4;
  auto bounds = std::make_shared<const RectangleBounds>(w, h);
  auto srf0 = Surface::makeShared<PlaneSurface>(Transform3::Identity(), bounds);
  auto srf1 = Surface::makeShared<PlaneSurface>(Transform3::Identity(), bounds);

  std::vector<const Surface*> sfPointers = {srf0.get(), srf1.get()};
  std::vector<std::shared_ptr<const Surface>> surfaces = {srf0, srf1};

  auto singleLookUp =
      std::make_unique<Acts::SurfaceArray::SingleElementLookup>(sfPointers);

  SurfaceArray sa(std::move(singleLookUp), surfaces);

  auto binContent = sa.at(Vector3(42, 42, 42));
  BOOST_CHECK_EQUAL(binContent.size(), 2u);
  BOOST_CHECK_EQUAL(sa.surfaces().size(), 2u);
}

BOOST_AUTO_TEST_SUITE_END()
}  // namespace Test

}  // namespace Acts
