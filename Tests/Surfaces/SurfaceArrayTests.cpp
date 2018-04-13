// This file is part of the ACTS project.
//
// Copyright (C) 2017-2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE SurfaceArray
#include <boost/test/included/unit_test.hpp>

#include <boost/format.hpp>
#include <boost/test/data/test_case.hpp>

#include "ACTS/Surfaces/CylinderBounds.hpp"
#include "ACTS/Surfaces/PlaneSurface.hpp"
#include "ACTS/Surfaces/RadialBounds.hpp"
#include "ACTS/Surfaces/RectangleBounds.hpp"
#include "ACTS/Surfaces/SurfaceArray.hpp"
#include "ACTS/Tools/LayerCreator.hpp"
#include "ACTS/Tools/SurfaceArrayCreator.hpp"
#include "ACTS/Utilities/BinningType.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/detail/Grid.hpp"

namespace bdata = boost::unit_test::data;
namespace tt    = boost::test_tools;

namespace Acts {

namespace Test {

  using SrfVec = std::vector<const Surface*>;
  struct SurfaceArrayFixture
  {
    std::vector<std::unique_ptr<const Surface>> m_surfaces;

    SurfaceArrayFixture() { BOOST_TEST_MESSAGE("setup fixture"); }
    ~SurfaceArrayFixture() { BOOST_TEST_MESSAGE("teardown fixture"); }

    SrfVec
    fullPhiTestSurfacesEC(size_t n     = 10,
                          double shift = 0,
                          double zbase = 0,
                          double r     = 10)
    {

      SrfVec res;

      double phiStep = 2 * M_PI / n;
      for (size_t i = 0; i < n; ++i) {

        double z = zbase + ((i % 2 == 0) ? 1 : -1) * 0.2;

        Transform3D trans;
        trans.setIdentity();
        trans.rotate(Eigen::AngleAxisd(i * phiStep + shift, Vector3D(0, 0, 1)));
        trans.translate(Vector3D(r, 0, z));

        auto bounds = std::make_shared<const RectangleBounds>(2, 1);

        auto transptr = std::make_shared<const Transform3D>(trans);
        auto srf      = std::make_unique<const PlaneSurface>(transptr, bounds);

        res.push_back(srf.get());  // use raw pointer
        m_surfaces.push_back(
            std::move(srf));  // keep unique, will get destroyed at the end
      }

      return res;
    }

    SrfVec
    fullPhiTestSurfacesBRL(int    n     = 10,
                           double shift = 0,
                           double zbase = 0,
                           double incl  = M_PI / 9.,
                           double w     = 2,
                           double h     = 1.5)
    {

      SrfVec res;

      double phiStep = 2 * M_PI / n;
      for (int i = 0; i < n; ++i) {

        double z = zbase;

        Transform3D trans;
        trans.setIdentity();
        trans.rotate(Eigen::AngleAxisd(i * phiStep + shift, Vector3D(0, 0, 1)));
        trans.translate(Vector3D(10, 0, z));
        trans.rotate(Eigen::AngleAxisd(incl, Vector3D(0, 0, 1)));
        trans.rotate(Eigen::AngleAxisd(M_PI / 2., Vector3D(0, 1, 0)));

        auto bounds = std::make_shared<const RectangleBounds>(w, h);

        auto transptr = std::make_shared<const Transform3D>(trans);
        auto srf      = std::make_unique<const PlaneSurface>(transptr, bounds);

        res.push_back(srf.get());  // use raw pointer
        m_surfaces.push_back(
            std::move(srf));  // keep unique, will get destroyed at the end
      }

      return res;
    }

    SrfVec
    makeBarrel(int nPhi, int nZ, double w, double h)
    {
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

    void
    draw_surfaces(SrfVec surfaces, std::string fname)
    {

      std::ofstream os;
      os.open(fname);

      os << std::fixed << std::setprecision(4);

      size_t nVtx = 0;
      for (const auto& srfx : surfaces) {
        const PlaneSurface* srf = dynamic_cast<const PlaneSurface*>(srfx);
        const PlanarBounds* bounds
            = dynamic_cast<const PlanarBounds*>(&srf->bounds());

        for (const auto& vtxloc : bounds->vertices()) {
          Vector3D vtx = srf->transform() * Vector3D(vtxloc.x(), vtxloc.y(), 0);
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

  BOOST_FIXTURE_TEST_CASE(SurfaceArray_create, SurfaceArrayFixture)
  {
    SrfVec brl = makeBarrel(30, 7, 2, 1);
    draw_surfaces(brl, "SurfaceArray_create_BRL_1.obj");

    detail::Axis<detail::AxisType::Equidistant,
                 detail::AxisBoundaryType::Closed>
        phiAxis(-M_PI, M_PI, 30u);
    detail::Axis<detail::AxisType::Equidistant, detail::AxisBoundaryType::Bound>
        zAxis(-14, 14, 7u);

    double angleShift = 2 * M_PI / 30. / 2.;
    auto transform    = [angleShift](const Vector3D& pos) {
      return Vector2D(pos.phi() + angleShift, pos.z());
    };
    double R        = 10;
    auto itransform = [angleShift, R](const Vector2D& loc) {
      return Vector3D(R * std::cos(loc[0] - angleShift),
                      R * std::sin(loc[0] - angleShift),
                      loc[1]);
    };

    auto sl
        = std::make_unique<SurfaceArray::SurfaceGridLookup<decltype(phiAxis),
                                                           decltype(zAxis)>>(
            transform,
            itransform,
            std::make_tuple(std::move(phiAxis), std::move(zAxis)));
    sl->fill(brl);
    SurfaceArray sa(std::move(sl), brl);

    // let's see if we can access all surfaces
    sa.dump(std::cout);

    for (const auto& srf : brl) {
      Vector3D ctr        = srf->binningPosition(binR);
      SrfVec   binContent = sa.at(ctr);

      BOOST_TEST(binContent.size() == 1);
      BOOST_TEST(srf == binContent.at(0));
    }

    std::vector<const Surface*> neighbors
        = sa.neighbors(itransform(Vector2D(0, 0)));
    BOOST_TEST(neighbors.size() == 9);

    auto sl2
        = std::make_unique<SurfaceArray::SurfaceGridLookup<decltype(phiAxis),
                                                           decltype(zAxis)>>(
            transform,
            itransform,
            std::make_tuple(std::move(phiAxis), std::move(zAxis)));
    // do NOT fill, only completebinning
    sl2->completeBinning(brl);
    SurfaceArray sa2(std::move(sl2), brl);
    sa.dump(std::cout);
    for (const auto& srf : brl) {
      Vector3D ctr        = srf->binningPosition(binR);
      SrfVec   binContent = sa2.at(ctr);

      BOOST_TEST(binContent.size() == 1);
      BOOST_TEST(srf == binContent.at(0));
    }
  }

  BOOST_AUTO_TEST_CASE(SurfaceArray_singleElement)
  {
    double w = 3, h = 4;
    auto   bounds = std::make_shared<const RectangleBounds>(w, h);
    auto   transptr
        = std::make_shared<const Transform3D>(Transform3D::Identity());
    auto srf = std::make_unique<const PlaneSurface>(transptr, bounds);

    SurfaceArray sa(srf.get());

    auto binContent = sa.at(Vector3D(42, 42, 42));
    BOOST_TEST(binContent.size() == 1);
    BOOST_TEST(binContent.at(0) == srf.get());
    BOOST_TEST(sa.surfaces().size() == 1);
    BOOST_TEST(sa.surfaces().at(0) == srf.get());
  }

  BOOST_AUTO_TEST_SUITE_END()
}  // end of namespace Test

}  // end of namespace Acts
