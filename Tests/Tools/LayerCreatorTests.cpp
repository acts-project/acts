// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE SurfaceArrayCreator
#include <boost/test/included/unit_test.hpp>

#include <boost/format.hpp>
#include <boost/test/data/test_case.hpp>

#include "ACTS/Layers/CylinderLayer.hpp"
#include "ACTS/Layers/DiscLayer.hpp"
#include "ACTS/Surfaces/CylinderBounds.hpp"
#include "ACTS/Surfaces/PlaneSurface.hpp"
#include "ACTS/Surfaces/RadialBounds.hpp"
#include "ACTS/Surfaces/RectangleBounds.hpp"
#include "ACTS/Tools/LayerCreator.hpp"
#include "ACTS/Tools/SurfaceArrayCreator.hpp"
#include "ACTS/Utilities/BinningType.hpp"

#include "ACTS/Utilities/Definitions.hpp"

#include <fstream>
#include <random>

namespace bdata = boost::unit_test::data;
namespace tt    = boost::test_tools;

namespace Acts {

namespace Test {

#define CHECK_ROTATION_ANGLE(t, a, tolerance)                                  \
  {                                                                            \
    Vector3D v = (*t) * Vector3D(1, 0, 0);                                     \
    BOOST_CHECK_SMALL(v.phi() - (a), tolerance);                               \
  }

  using SrfVec = std::vector<const Surface*>;

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

  struct LayerCreatorFixture
  {
    std::shared_ptr<const SurfaceArrayCreator>  p_SAC;
    std::shared_ptr<LayerCreator>               p_LC;
    std::vector<std::unique_ptr<const Surface>> m_surfaces;

    LayerCreatorFixture()
    {
      p_SAC = std::make_shared<const SurfaceArrayCreator>(
          SurfaceArrayCreator::Config(),
          Acts::getDefaultLogger("SurfaceArrayCreator",
                                 Acts::Logging::VERBOSE));
      LayerCreator::Config cfg;
      cfg.surfaceArrayCreator = p_SAC;
      p_LC                    = std::make_shared<LayerCreator>(
          cfg, Acts::getDefaultLogger("LayerCreator", Acts::Logging::VERBOSE));
    }

    template <typename... Args>
    bool
    checkBinning(Args&&... args)
    {
      return p_LC->checkBinning(std::forward<Args>(args)...);
    }

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
      for (size_t i = 0; i < n; ++i) {

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

      for (size_t i = 0; i < nZ; i++) {
        double z = i * w * 2 + z0;
        std::cout << "z=" << z << std::endl;
        SrfVec ring = fullPhiTestSurfacesBRL(nPhi, 0, z, M_PI / 9., w, h);
        res.insert(res.end(), ring.begin(), ring.end());
      }

      return res;
    }
  };

  BOOST_AUTO_TEST_SUITE(Tools);

  BOOST_FIXTURE_TEST_CASE(LayerCreator_createCylinderLayer, LayerCreatorFixture)
  {
    std::vector<const Surface*> srf;

    srf = makeBarrel(30, 7, 2, 1.5);
    draw_surfaces(srf, "LayerCreator_createCylinderLayer_BRL_1.obj");

    double                         envR = 0.1, envZ = 0.5;
    std::shared_ptr<CylinderLayer> layer
        = std::dynamic_pointer_cast<CylinderLayer>(
            p_LC->cylinderLayer(srf, envR, envZ, equidistant, equidistant));

    double rMax = 10.6071, rMin = 9.59111;  // empirical
    BOOST_TEST(layer->thickness() == (rMax - rMin) + 2 * envR,
               tt::tolerance(1e-3));

    const CylinderBounds* bounds = &layer->bounds();
    BOOST_TEST(bounds->r() == (rMax + rMin) / 2., tt::tolerance(1e-3));
    BOOST_TEST(bounds->halflengthZ() == 14 + envZ, tt::tolerance(1e-3));
    BOOST_TEST(checkBinning(layer->surfaceArray()->objectGrid(), srf));
    auto bu = layer->surfaceArray()->binUtility();
    BOOST_TEST(bu->bins(0) == 30);
    BOOST_TEST(bu->bins(1) == 7);
    BOOST_TEST(bu->binningData().at(0).min == -M_PI, tt::tolerance(1e-3));
    BOOST_TEST(bu->binningData().at(0).max == M_PI, tt::tolerance(1e-3));
    BOOST_TEST(bu->binningData().at(1).min == -14, tt::tolerance(1e-3));
    BOOST_TEST(bu->binningData().at(1).max == 14, tt::tolerance(1e-3));

    layer = std::dynamic_pointer_cast<CylinderLayer>(
        p_LC->cylinderLayer(srf, envR, envZ, 30, 7));
    BOOST_TEST(layer->thickness() == (rMax - rMin) + 2 * envR,
               tt::tolerance(1e-3));
    bounds = &layer->bounds();
    BOOST_TEST(bounds->r() == (rMax + rMin) / 2., tt::tolerance(1e-3));
    BOOST_TEST(bounds->halflengthZ() == 14 + envZ, tt::tolerance(1e-3));
    BOOST_TEST(checkBinning(layer->surfaceArray()->objectGrid(), srf));
    bu = layer->surfaceArray()->binUtility();
    BOOST_TEST(bu->bins(0) == 30);
    BOOST_TEST(bu->bins(1) == 7);

    BOOST_TEST(bu->binningData().at(0).min == -M_PI, tt::tolerance(1e-3));
    BOOST_TEST(bu->binningData().at(0).max == M_PI, tt::tolerance(1e-3));

    BOOST_TEST(bu->binningData().at(1).min == -14, tt::tolerance(1e-3));
    BOOST_TEST(bu->binningData().at(1).max == 14, tt::tolerance(1e-3));

    layer = std::dynamic_pointer_cast<CylinderLayer>(
        p_LC->cylinderLayer(srf, 1, 20, 25, equidistant, equidistant));
    BOOST_TEST(layer->thickness() == 19, tt::tolerance(1e-3));
    bounds = &layer->bounds();
    BOOST_TEST(bounds->r() == 10.5, tt::tolerance(1e-3));
    BOOST_TEST(bounds->halflengthZ() == 25, tt::tolerance(1e-3));

    // this should fail, b/c it's a completely inconvenient binning
    BOOST_TEST(!checkBinning(layer->surfaceArray()->objectGrid(), srf));

    bu = layer->surfaceArray()->binUtility();

    BOOST_TEST(bu->bins(0) == 30);
    BOOST_TEST(bu->bins(1) == 7);

    BOOST_TEST(bu->binningData().at(0).min == -M_PI, tt::tolerance(1e-3));
    BOOST_TEST(bu->binningData().at(0).max == M_PI, tt::tolerance(1e-3));

    BOOST_TEST(bu->binningData().at(1).min == -25, tt::tolerance(1e-3));
    BOOST_TEST(bu->binningData().at(1).max == 25, tt::tolerance(1e-3));
  }

  BOOST_FIXTURE_TEST_CASE(LayerCreator_createDiscLayer, LayerCreatorFixture)
  {
    std::vector<const Surface*> surfaces;
    auto                        ringa = fullPhiTestSurfacesEC(30, 0, 0, 10);
    surfaces.insert(surfaces.end(), ringa.begin(), ringa.end());
    auto ringb = fullPhiTestSurfacesEC(30, 0, 0, 15);
    surfaces.insert(surfaces.end(), ringb.begin(), ringb.end());
    auto ringc = fullPhiTestSurfacesEC(30, 0, 0, 20);
    surfaces.insert(surfaces.end(), ringc.begin(), ringc.end());
    draw_surfaces(surfaces, "LayerCreator_createDiscLayer_EC_1.obj");

    std::shared_ptr<DiscLayer> layer = std::dynamic_pointer_cast<DiscLayer>(
        p_LC->discLayer(surfaces, -10, 10, 5, 25, equidistant, equidistant));
    BOOST_TEST(layer->thickness() == 20, tt::tolerance(1e-3));
    const RadialBounds* bounds
        = dynamic_cast<const RadialBounds*>(&layer->bounds());
    BOOST_TEST(bounds->rMin() == 5, tt::tolerance(1e-3));
    BOOST_TEST(bounds->rMax() == 25, tt::tolerance(1e-3));
    BOOST_TEST(checkBinning(layer->surfaceArray()->objectGrid(), surfaces));
    auto bu = layer->surfaceArray()->binUtility();
    BOOST_TEST(bu->bins(0) == 3);
    BOOST_TEST(bu->bins(1) == 30);
    BOOST_TEST(bu->binningData().at(0).min == 5, tt::tolerance(1e-3));
    BOOST_TEST(bu->binningData().at(0).max == 25, tt::tolerance(1e-3));
    BOOST_TEST(bu->binningData().at(1).min == -M_PI, tt::tolerance(1e-3));
    BOOST_TEST(bu->binningData().at(1).max == M_PI, tt::tolerance(1e-3));
    // check that it's applying a rotation transform to improve phi binning
    BOOST_TEST(bu->transform() != nullptr);
    double actAngle = ((*bu->transform()) * Vector3D(1, 0, 0)).phi();
    double expAngle = -2 * M_PI / 30 / 2.;
    BOOST_TEST(actAngle == expAngle, tt::tolerance(1e-3));

    double envMinR = 1, envMaxR = 2, envZ = 5;
    size_t nBinsR = 3, nBinsPhi = 30;
    layer = std::dynamic_pointer_cast<DiscLayer>(
        p_LC->discLayer(surfaces, envMinR, envMaxR, envZ, nBinsR, nBinsPhi));

    double rMin = 8, rMax = 22.0227;
    BOOST_TEST(layer->thickness() == 0.4 + 2 * envZ, tt::tolerance(1e-3));
    bounds = dynamic_cast<const RadialBounds*>(&layer->bounds());
    BOOST_TEST(bounds->rMin() == rMin - envMinR, tt::tolerance(1e-3));
    BOOST_TEST(bounds->rMax() == rMax + envMaxR, tt::tolerance(1e-3));
    BOOST_TEST(checkBinning(layer->surfaceArray()->objectGrid(), surfaces));
    bu = layer->surfaceArray()->binUtility();
    BOOST_TEST(bu->bins(0) == nBinsR);
    BOOST_TEST(bu->bins(1) == nBinsPhi);
    BOOST_TEST(bu->binningData().at(0).min == rMin, tt::tolerance(1e-3));
    BOOST_TEST(bu->binningData().at(0).max == rMax, tt::tolerance(1e-3));
    BOOST_TEST(bu->binningData().at(1).min == -M_PI, tt::tolerance(1e-3));
    BOOST_TEST(bu->binningData().at(1).max == M_PI, tt::tolerance(1e-3));
    // check that it's applying a rotation transform to improve phi binning
    BOOST_TEST(bu->transform() != nullptr);
    actAngle = ((*bu->transform()) * Vector3D(1, 0, 0)).phi();
    expAngle = -2 * M_PI / 30 / 2.;
    BOOST_TEST(actAngle == expAngle, tt::tolerance(1e-3));

    layer = std::dynamic_pointer_cast<DiscLayer>(p_LC->discLayer(
        surfaces, envMinR, envMaxR, envZ, equidistant, equidistant));
    BOOST_TEST(layer->thickness() == 0.4 + 2 * envZ, tt::tolerance(1e-3));
    bounds = dynamic_cast<const RadialBounds*>(&layer->bounds());
    BOOST_TEST(bounds->rMin() == rMin - envMinR, tt::tolerance(1e-3));
    BOOST_TEST(bounds->rMax() == rMax + envMaxR, tt::tolerance(1e-3));
    BOOST_TEST(checkBinning(layer->surfaceArray()->objectGrid(), surfaces));
    bu = layer->surfaceArray()->binUtility();
    BOOST_TEST(bu->bins(0) == nBinsR);
    BOOST_TEST(bu->bins(1) == nBinsPhi);
    BOOST_TEST(bu->binningData().at(0).min == rMin, tt::tolerance(1e-3));
    BOOST_TEST(bu->binningData().at(0).max == rMax, tt::tolerance(1e-3));
    BOOST_TEST(bu->binningData().at(1).min == -M_PI, tt::tolerance(1e-3));
    BOOST_TEST(bu->binningData().at(1).max == M_PI, tt::tolerance(1e-3));
    // check that it's applying a rotation transform to improve phi binning
    BOOST_TEST(bu->transform() != nullptr);
    actAngle = ((*bu->transform()) * Vector3D(1, 0, 0)).phi();
    expAngle = -2 * M_PI / 30 / 2.;
    BOOST_TEST(actAngle == expAngle, tt::tolerance(1e-3));
  }

  BOOST_AUTO_TEST_SUITE_END();
}  // end of namespace Test

}  // end of namespace Acts
