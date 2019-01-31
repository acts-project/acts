// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE SurfaceArray
#include <boost/test/included/unit_test.hpp>

#include <boost/format.hpp>
#include <boost/test/data/test_case.hpp>

#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"
#include "Acts/Tools/LayerCreator.hpp"
#include "Acts/Tools/SurfaceArrayCreator.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/VariantData.hpp"
#include "Acts/Utilities/detail/Grid.hpp"

using Acts::VectorHelpers::phi;
using Acts::VectorHelpers::perp;

namespace bdata = boost::unit_test::data;
namespace tt    = boost::test_tools;

namespace Acts {

namespace Test {

  using SrfVec = std::vector<std::shared_ptr<const Surface>>;
  struct SurfaceArrayFixture
  {
    std::vector<std::shared_ptr<const Surface>> m_surfaces;

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
        std::shared_ptr<const Surface> srf
            = Surface::makeShared<PlaneSurface>(transptr, bounds);

        res.push_back(srf);
        m_surfaces.push_back(
            std::move(srf));  // keep shared, will get destroyed at the end
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
        std::shared_ptr<const Surface> srf
            = Surface::makeShared<PlaneSurface>(transptr, bounds);

        res.push_back(srf);
        m_surfaces.push_back(
            std::move(srf));  // keep shared, will get destroyed at the end
      }

      return res;
    }

    SrfVec
    straightLineSurfaces(size_t             n        = 10.,
                         double             step     = 3,
                         const Vector3D&    origin   = {0, 0, 1.5},
                         const Transform3D& pretrans = Transform3D::Identity(),
                         const Vector3D&    dir      = {0, 0, 1})
    {
      SrfVec res;
      for (size_t i = 0; i < n; ++i) {
        Transform3D trans;
        trans.setIdentity();
        trans.translate(origin + dir * step * i);
        // trans.rotate(AngleAxis3D(M_PI/9., Vector3D(0, 0, 1)));
        trans.rotate(AngleAxis3D(M_PI / 2., Vector3D(1, 0, 0)));
        trans = trans * pretrans;

        auto bounds = std::make_shared<const RectangleBounds>(2, 1.5);

        auto transptr = std::make_shared<const Transform3D>(trans);
        std::shared_ptr<const Surface> srf
            = Surface::makeShared<PlaneSurface>(transptr, bounds);

        res.push_back(srf);
        m_surfaces.push_back(
            std::move(srf));  // keep shared, will get destroyed at the end
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
    draw_surfaces(const SrfVec& surfaces, const std::string& fname)
    {

      std::ofstream os;
      os.open(fname);

      os << std::fixed << std::setprecision(4);

      size_t nVtx = 0;
      for (const auto& srfx : surfaces) {
        std::shared_ptr<const PlaneSurface> srf
            = std::dynamic_pointer_cast<const PlaneSurface>(srfx);
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
    SrfVec                      brl    = makeBarrel(30, 7, 2, 1);
    std::vector<const Surface*> brlRaw = unpack_shared_vector(brl);
    draw_surfaces(brl, "SurfaceArray_create_BRL_1.obj");

    detail::Axis<detail::AxisType::Equidistant,
                 detail::AxisBoundaryType::Closed>
        phiAxis(-M_PI, M_PI, 30u);
    detail::Axis<detail::AxisType::Equidistant, detail::AxisBoundaryType::Bound>
        zAxis(-14, 14, 7u);

    double angleShift = 2 * M_PI / 30. / 2.;
    auto transform    = [angleShift](const Vector3D& pos) {
      return Vector2D(phi(pos) + angleShift, pos.z());
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
    sl->fill(brlRaw);
    SurfaceArray sa(std::move(sl), brl);

    // let's see if we can access all surfaces
    sa.dump(std::cout);

    for (const auto& srf : brl) {
      Vector3D                    ctr        = srf->binningPosition(binR);
      std::vector<const Surface*> binContent = sa.at(ctr);

      BOOST_CHECK_EQUAL(binContent.size(), 1);
      BOOST_CHECK_EQUAL(srf.get(), binContent.at(0));
    }

    std::vector<const Surface*> neighbors
        = sa.neighbors(itransform(Vector2D(0, 0)));
    BOOST_CHECK_EQUAL(neighbors.size(), 9);

    auto sl2
        = std::make_unique<SurfaceArray::SurfaceGridLookup<decltype(phiAxis),
                                                           decltype(zAxis)>>(
            transform,
            itransform,
            std::make_tuple(std::move(phiAxis), std::move(zAxis)));
    // do NOT fill, only completebinning
    sl2->completeBinning(brlRaw);
    SurfaceArray sa2(std::move(sl2), brl);
    sa.dump(std::cout);
    for (const auto& srf : brl) {
      Vector3D                    ctr        = srf->binningPosition(binR);
      std::vector<const Surface*> binContent = sa2.at(ctr);

      BOOST_CHECK_EQUAL(binContent.size(), 1);
      BOOST_CHECK_EQUAL(srf.get(), binContent.at(0));
    }
  }

  BOOST_AUTO_TEST_CASE(SurfaceArray_singleElement)
  {
    double w = 3, h = 4;
    auto   bounds = std::make_shared<const RectangleBounds>(w, h);
    auto   transptr
        = std::make_shared<const Transform3D>(Transform3D::Identity());
    auto srf = Surface::makeShared<PlaneSurface>(transptr, bounds);

    SurfaceArray sa(srf);

    auto binContent = sa.at(Vector3D(42, 42, 42));
    BOOST_CHECK_EQUAL(binContent.size(), 1);
    BOOST_CHECK_EQUAL(binContent.at(0), srf.get());
    BOOST_CHECK_EQUAL(sa.surfaces().size(), 1);
    BOOST_CHECK_EQUAL(sa.surfaces().at(0), srf.get());
  }

  BOOST_FIXTURE_TEST_CASE(SurfaceArray_toVariantData, SurfaceArrayFixture)
  {
    SrfVec                      brl    = makeBarrel(30, 7, 2, 1);
    std::vector<const Surface*> brlRaw = unpack_shared_vector(brl);

    detail::Axis<detail::AxisType::Equidistant,
                 detail::AxisBoundaryType::Closed>
                        phiAxis(-M_PI, M_PI, 30u);
    std::vector<double> zAxis_bin_edges_exp = {-14, -10, 3, 5, 8, 14};
    detail::Axis<detail::AxisType::Variable, detail::AxisBoundaryType::Bound>
        zAxis(zAxis_bin_edges_exp);

    double angleShift = 2 * M_PI / 30. / 2.;
    auto transform    = [angleShift](const Vector3D& pos) {
      return Vector2D(phi(pos) + angleShift, pos.z());
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
    sl->fill(brlRaw);
    SurfaceArray sa(std::move(sl), brl);
    sa.dump(std::cout);

    variant_data data = sa.toVariantData();
    // std::cout << data << std::endl;

    const variant_map& var_map = boost::get<variant_map>(data);
    BOOST_CHECK_EQUAL(var_map.get<std::string>("type"), "SurfaceArray");
    const variant_map& sa_var_pl = var_map.get<variant_map>("payload");
    BOOST_CHECK_EQUAL(sa_var_pl.count("surfacegridlookup"), 1);
    const variant_map& sgl_var_pl
        = sa_var_pl.get<variant_map>("surfacegridlookup")
              .get<variant_map>("payload");
    BOOST_CHECK_EQUAL(sgl_var_pl.get<int>("dimensions"), 2);
    const variant_vector& axes = sgl_var_pl.get<variant_vector>("axes");
    BOOST_CHECK_EQUAL(axes.size(), 2);

    const variant_map& phiAxis_pl
        = axes.get<variant_map>(0).get<variant_map>("payload");
    BOOST_CHECK_EQUAL(phiAxis_pl.get<std::string>("axisboundarytype"),
                      "closed");
    BOOST_CHECK_EQUAL(phiAxis_pl.get<std::string>("axistype"), "equidistant");
    BOOST_CHECK_EQUAL(phiAxis_pl.get<double>("min"), -M_PI);
    BOOST_CHECK_EQUAL(phiAxis_pl.get<double>("max"), M_PI);
    BOOST_CHECK_EQUAL(phiAxis_pl.get<int>("nbins"), 30);

    const variant_map& zAxis_pl
        = axes.get<variant_map>(1).get<variant_map>("payload");
    BOOST_CHECK_EQUAL(zAxis_pl.get<std::string>("axisboundarytype"), "bound");
    BOOST_CHECK_EQUAL(zAxis_pl.get<std::string>("axistype"), "variable");
    const variant_vector& zAxis_bin_edges
        = zAxis_pl.get<variant_vector>("bin_edges");
    BOOST_CHECK_EQUAL(zAxis_bin_edges.size(), 6);
    for (size_t i = 0; i < zAxis_bin_edges.size(); i++) {
      BOOST_CHECK_EQUAL(zAxis_bin_edges.get<double>(i),
                        zAxis_bin_edges_exp.at(i));
    }

    SurfaceArray sa2(data, transform, itransform);
    sa2.dump(std::cout);

    std::ostringstream dumpExp_os;
    sa.dump(dumpExp_os);
    std::string                           dumpExp = dumpExp_os.str();
    boost::test_tools::output_test_stream dumpAct;
    sa2.dump(dumpAct);
    BOOST_CHECK(dumpAct.is_equal(dumpExp));
  }

  BOOST_FIXTURE_TEST_CASE(SurfaceArray_toVariantData_1D, SurfaceArrayFixture)
  {
    detail::Axis<detail::AxisType::Equidistant, detail::AxisBoundaryType::Bound>
         zAxis(0, 30, 10);
    auto transform = [](const Vector3D& pos) {
      return std::array<double, 1>({{pos.z()}});
    };
    auto itransform = [](const std::array<double, 1>& loc) {
      return Vector3D(0, 0, loc[0]);
    };
    auto sl
        = std::make_unique<SurfaceArray::SurfaceGridLookup<decltype(zAxis)>>(
            transform, itransform, std::make_tuple(zAxis));

    // same thing in 1D
    SrfVec                      line    = straightLineSurfaces();
    std::vector<const Surface*> lineRaw = unpack_shared_vector(line);
    sl->fill(lineRaw);
    SurfaceArray sa(std::move(sl), line);

    sa.dump(std::cout);

    variant_data data = sa.toVariantData();
    // std::cout << data << std::endl;

    const variant_map& var_map = boost::get<variant_map>(data);
    BOOST_CHECK_EQUAL(var_map.get<std::string>("type"), "SurfaceArray");
    const variant_map& sa_var_pl = var_map.get<variant_map>("payload");
    BOOST_CHECK_EQUAL(sa_var_pl.count("surfacegridlookup"), 1);
    const variant_map& sgl_var_pl
        = sa_var_pl.get<variant_map>("surfacegridlookup")
              .get<variant_map>("payload");
    BOOST_CHECK_EQUAL(sgl_var_pl.get<int>("dimensions"), 1);
    const variant_vector& axes = sgl_var_pl.get<variant_vector>("axes");
    BOOST_CHECK_EQUAL(axes.size(), 1);

    const variant_map& zAxis_pl
        = axes.get<variant_map>(0).get<variant_map>("payload");
    BOOST_CHECK_EQUAL(zAxis_pl.get<std::string>("axisboundarytype"), "bound");
    BOOST_CHECK_EQUAL(zAxis_pl.get<std::string>("axistype"), "equidistant");
    BOOST_CHECK_EQUAL(zAxis_pl.get<double>("min"), 0);
    BOOST_CHECK_EQUAL(zAxis_pl.get<double>("max"), 30);
    BOOST_CHECK_EQUAL(zAxis_pl.get<int>("nbins"), 10);

    SurfaceArray sa2(data, transform, itransform);
    sa2.dump(std::cout);

    std::ostringstream dumpExp_os;
    sa.dump(dumpExp_os);
    std::string                           dumpExp = dumpExp_os.str();
    boost::test_tools::output_test_stream dumpAct;
    sa2.dump(dumpAct);
    BOOST_CHECK(dumpAct.is_equal(dumpExp));
  }

  BOOST_AUTO_TEST_SUITE_END()
}  // namespace Test

}  // namespace Acts
