// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// clang-format off
#define BOOST_TEST_MODULE LayerCreator
#include <boost/test/included/unit_test.hpp>
#include <boost/format.hpp>
#include <boost/test/data/test_case.hpp>
// clang-format on

#include <fstream>
#include <random>

#include "Acts/Layers/CylinderLayer.hpp"
#include "Acts/Layers/DiscLayer.hpp"
#include "Acts/Layers/ProtoLayer.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Tools/LayerCreator.hpp"
#include "Acts/Tools/SurfaceArrayCreator.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/VariantData.hpp"

namespace bdata = boost::unit_test::data;
namespace tt    = boost::test_tools;

namespace Acts {

namespace Test {

#define CHECK_ROTATION_ANGLE(t, a, tolerance)                                  \
  {                                                                            \
    Vector3D v = (*t) * Vector3D(1, 0, 0);                                     \
    CHECK_CLOSE_ABS(VectorHelpers::phi(v), (a), tolerance);                    \
  }

  using SrfVec = std::vector<std::shared_ptr<const Surface>>;

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

  struct LayerCreatorFixture
  {
    std::shared_ptr<const SurfaceArrayCreator>  p_SAC;
    std::shared_ptr<LayerCreator>               p_LC;
    std::vector<std::shared_ptr<const Surface>> m_surfaces;

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

    bool
    checkBinContentSize(const SurfaceArray* sArray, size_t n)
    {
      size_t nBins  = sArray->size();
      bool   result = true;
      for (size_t i = 0; i < nBins; ++i) {
        if (!sArray->isValidBin(i)) {
          continue;
        }
        std::vector<const Surface*> binContent = sArray->at(i);
        BOOST_TEST_INFO("Bin: " << i);
        BOOST_CHECK_EQUAL(binContent.size(), n);
        result = result && binContent.size() == n;
      }

      return result;
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
        std::shared_ptr<PlaneSurface> srf
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
        std::shared_ptr<PlaneSurface> srf
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
        std::cout << "z=" << z << std::endl;
        SrfVec ring = fullPhiTestSurfacesBRL(nPhi, 0, z, M_PI / 9., w, h);
        res.insert(res.end(), ring.begin(), ring.end());
      }

      return res;
    }

    std::pair<SrfVec, std::vector<std::pair<const Surface*, const Surface*>>>
    makeBarrelStagger(int    nPhi,
                      int    nZ,
                      double shift = 0,
                      double incl  = M_PI / 9.,
                      double w     = 2,
                      double h     = 1.5)
    {
      double z0 = -(nZ - 1) * w;
      SrfVec res;

      std::vector<std::pair<const Surface*, const Surface*>> pairs;

      for (int i = 0; i < nZ; i++) {
        double z = i * w * 2 + z0;

        double phiStep = 2 * M_PI / nPhi;
        for (int j = 0; j < nPhi; ++j) {

          Transform3D trans;
          trans.setIdentity();
          trans.rotate(
              Eigen::AngleAxisd(j * phiStep + shift, Vector3D(0, 0, 1)));
          trans.translate(Vector3D(10, 0, z));
          trans.rotate(Eigen::AngleAxisd(incl, Vector3D(0, 0, 1)));
          trans.rotate(Eigen::AngleAxisd(M_PI / 2., Vector3D(0, 1, 0)));

          auto bounds = std::make_shared<const RectangleBounds>(w, h);

          auto transAptr = std::make_shared<const Transform3D>(trans);

          std::shared_ptr<PlaneSurface> srfA
              = Surface::makeShared<PlaneSurface>(transAptr, bounds);

          Vector3D    nrm    = srfA->normal();
          Transform3D transB = trans;
          transB.pretranslate(nrm * 0.1);
          auto transBptr = std::make_shared<const Transform3D>(transB);
          std::shared_ptr<PlaneSurface> srfB
              = Surface::makeShared<PlaneSurface>(transBptr, bounds);

          pairs.push_back(std::make_pair(srfA.get(), srfB.get()));

          res.push_back(srfA);
          res.push_back(srfB);
          m_surfaces.push_back(std::move(srfA));
          m_surfaces.push_back(std::move(srfB));
        }
      }

      return std::make_pair(res, pairs);
    }
  };

  BOOST_AUTO_TEST_SUITE(Tools)

  BOOST_FIXTURE_TEST_CASE(LayerCreator_createCylinderLayer, LayerCreatorFixture)
  {
    std::vector<std::shared_ptr<const Surface>> srf;

    srf = makeBarrel(30, 7, 2, 1.5);
    draw_surfaces(srf, "LayerCreator_createCylinderLayer_BRL_1.obj");

    // CASE I
    double     envR = 0.1, envZ = 0.5;
    ProtoLayer pl(srf);
    pl.envR = {envR, envR};
    pl.envZ = {envZ, envZ};
    std::shared_ptr<CylinderLayer> layer
        = std::dynamic_pointer_cast<CylinderLayer>(
            p_LC->cylinderLayer(srf, equidistant, equidistant, pl));

    double rMax = 10.6071, rMin = 9.59111;  // empirical
    CHECK_CLOSE_REL(layer->thickness(), (rMax - rMin) + 2 * envR, 1e-3);

    const CylinderBounds* bounds = &layer->bounds();
    CHECK_CLOSE_REL(bounds->r(), (rMax + rMin) / 2., 1e-3);
    CHECK_CLOSE_REL(bounds->halflengthZ(), 14 + envZ, 1e-3);
    BOOST_CHECK(checkBinning(*layer->surfaceArray()));
    auto axes = layer->surfaceArray()->getAxes();
    BOOST_CHECK_EQUAL(axes.at(0)->getNBins(), 30);
    BOOST_CHECK_EQUAL(axes.at(1)->getNBins(), 7);
    CHECK_CLOSE_REL(axes.at(0)->getMin(), -M_PI, 1e-3);
    CHECK_CLOSE_REL(axes.at(0)->getMax(), M_PI, 1e-3);
    CHECK_CLOSE_REL(axes.at(1)->getMin(), -14, 1e-3);
    CHECK_CLOSE_REL(axes.at(1)->getMax(), 14, 1e-3);

    // CASE II

    ProtoLayer pl2(srf);
    pl2.envR = {envR, envR};
    pl2.envZ = {envZ, envZ};
    layer    = std::dynamic_pointer_cast<CylinderLayer>(
        p_LC->cylinderLayer(srf, 30, 7, pl2));
    CHECK_CLOSE_REL(layer->thickness(), (rMax - rMin) + 2 * envR, 1e-3);
    bounds = &layer->bounds();
    CHECK_CLOSE_REL(bounds->r(), (rMax + rMin) / 2., 1e-3);
    CHECK_CLOSE_REL(bounds->halflengthZ(), 14 + envZ, 1e-3);
    BOOST_CHECK(checkBinning(*layer->surfaceArray()));
    axes = layer->surfaceArray()->getAxes();
    BOOST_CHECK_EQUAL(axes.at(0)->getNBins(), 30);
    BOOST_CHECK_EQUAL(axes.at(1)->getNBins(), 7);
    CHECK_CLOSE_REL(axes.at(0)->getMin(), -M_PI, 1e-3);
    CHECK_CLOSE_REL(axes.at(0)->getMax(), M_PI, 1e-3);
    CHECK_CLOSE_REL(axes.at(1)->getMin(), -14, 1e-3);
    CHECK_CLOSE_REL(axes.at(1)->getMax(), 14, 1e-3);

    layer = std::dynamic_pointer_cast<CylinderLayer>(
        p_LC->cylinderLayer(srf, 13, 3, pl2));
    CHECK_CLOSE_REL(layer->thickness(), (rMax - rMin) + 2 * envR, 1e-3);
    bounds = &layer->bounds();
    CHECK_CLOSE_REL(bounds->r(), (rMax + rMin) / 2., 1e-3);
    CHECK_CLOSE_REL(bounds->halflengthZ(), 14 + envZ, 1e-3);
    // this succeeds despite sub-optimal binning
    // since we now have multientry bins
    BOOST_CHECK(checkBinning(*layer->surfaceArray()));
    axes = layer->surfaceArray()->getAxes();
    BOOST_CHECK_EQUAL(axes.at(0)->getNBins(), 13);
    BOOST_CHECK_EQUAL(axes.at(1)->getNBins(), 3);
    CHECK_CLOSE_REL(axes.at(0)->getMin(), -M_PI, 1e-3);
    CHECK_CLOSE_REL(axes.at(0)->getMax(), M_PI, 1e-3);
    CHECK_CLOSE_REL(axes.at(1)->getMin(), -14, 1e-3);
    CHECK_CLOSE_REL(axes.at(1)->getMax(), 14, 1e-3);

    // CASE III
    ProtoLayer pl3;
    pl3.minR = 1;
    pl3.maxR = 20;
    pl3.minZ = -25;
    pl3.maxZ = 25;
    layer    = std::dynamic_pointer_cast<CylinderLayer>(
        p_LC->cylinderLayer(srf, equidistant, equidistant, pl3));
    CHECK_CLOSE_REL(layer->thickness(), 19, 1e-3);
    bounds = &layer->bounds();
    CHECK_CLOSE_REL(bounds->r(), 10.5, 1e-3);
    CHECK_CLOSE_REL(bounds->halflengthZ(), 25, 1e-3);

    // this should fail, b/c it's a completely inconvenient binning
    // but it succeeds despite sub-optimal binning
    // since we now have multientry bins
    BOOST_CHECK(checkBinning(*layer->surfaceArray()));

    axes = layer->surfaceArray()->getAxes();
    BOOST_CHECK_EQUAL(axes.at(0)->getNBins(), 30);
    BOOST_CHECK_EQUAL(axes.at(1)->getNBins(), 7);
    CHECK_CLOSE_REL(axes.at(0)->getMin(), -M_PI, 1e-3);
    CHECK_CLOSE_REL(axes.at(0)->getMax(), M_PI, 1e-3);
    CHECK_CLOSE_REL(axes.at(1)->getMin(), -25, 1e-3);
    CHECK_CLOSE_REL(axes.at(1)->getMax(), 25, 1e-3);
  }

  BOOST_FIXTURE_TEST_CASE(LayerCreator_createDiscLayer, LayerCreatorFixture)
  {
    std::vector<std::shared_ptr<const Surface>> surfaces;
    auto ringa = fullPhiTestSurfacesEC(30, 0, 0, 10);
    surfaces.insert(surfaces.end(), ringa.begin(), ringa.end());
    auto ringb = fullPhiTestSurfacesEC(30, 0, 0, 15);
    surfaces.insert(surfaces.end(), ringb.begin(), ringb.end());
    auto ringc = fullPhiTestSurfacesEC(30, 0, 0, 20);
    surfaces.insert(surfaces.end(), ringc.begin(), ringc.end());
    draw_surfaces(surfaces, "LayerCreator_createDiscLayer_EC_1.obj");

    ProtoLayer pl(surfaces);
    pl.minZ                          = -10;
    pl.maxZ                          = 10;
    pl.minR                          = 5;
    pl.maxR                          = 25;
    std::shared_ptr<DiscLayer> layer = std::dynamic_pointer_cast<DiscLayer>(
        p_LC->discLayer(surfaces, equidistant, equidistant, pl));
    CHECK_CLOSE_REL(layer->thickness(), 20, 1e-3);
    const RadialBounds* bounds
        = dynamic_cast<const RadialBounds*>(&layer->bounds());
    CHECK_CLOSE_REL(bounds->rMin(), 5, 1e-3);
    CHECK_CLOSE_REL(bounds->rMax(), 25, 1e-3);
    BOOST_CHECK(checkBinning(*layer->surfaceArray()));
    auto axes = layer->surfaceArray()->getAxes();
    BOOST_CHECK_EQUAL(axes.at(0)->getNBins(), 3);
    BOOST_CHECK_EQUAL(axes.at(1)->getNBins(), 30);
    CHECK_CLOSE_REL(axes.at(0)->getMin(), 5, 1e-3);
    CHECK_CLOSE_REL(axes.at(0)->getMax(), 25, 1e-3);
    CHECK_CLOSE_REL(axes.at(1)->getMin(), -M_PI, 1e-3);
    CHECK_CLOSE_REL(axes.at(1)->getMax(), M_PI, 1e-3);
    checkBinContentSize(layer->surfaceArray(), 1);

    // check that it's applying a rotation transform to improve phi binning
    // BOOST_CHECK_NE(bu->transform(), nullptr);
    // double actAngle = ((*bu->transform()) * Vector3D(1, 0, 0)).phi();
    // double expAngle = -2 * M_PI / 30 / 2.;
    // CHECK_CLOSE_REL(actAngle, expAngle, 1e-3);

    double     envMinR = 1, envMaxR = 1, envZ = 5;
    size_t     nBinsR = 3, nBinsPhi = 30;
    ProtoLayer pl2(surfaces);
    pl2.envR = {envMinR, envMaxR};
    pl2.envZ = {envZ, envZ};
    layer    = std::dynamic_pointer_cast<DiscLayer>(
        p_LC->discLayer(surfaces, nBinsR, nBinsPhi, pl2));

    double rMin = 8, rMax = 22.0227;
    CHECK_CLOSE_REL(layer->thickness(), 0.4 + 2 * envZ, 1e-3);
    bounds = dynamic_cast<const RadialBounds*>(&layer->bounds());
    CHECK_CLOSE_REL(bounds->rMin(), rMin - envMinR, 1e-3);
    CHECK_CLOSE_REL(bounds->rMax(), rMax + envMaxR, 1e-3);
    BOOST_CHECK(checkBinning(*layer->surfaceArray()));
    axes = layer->surfaceArray()->getAxes();
    BOOST_CHECK_EQUAL(axes.at(0)->getNBins(), nBinsR);
    BOOST_CHECK_EQUAL(axes.at(1)->getNBins(), nBinsPhi);
    CHECK_CLOSE_REL(axes.at(0)->getMin(), rMin, 1e-3);
    CHECK_CLOSE_REL(axes.at(0)->getMax(), rMax, 1e-3);
    CHECK_CLOSE_REL(axes.at(1)->getMin(), -M_PI, 1e-3);
    CHECK_CLOSE_REL(axes.at(1)->getMax(), M_PI, 1e-3);
    checkBinContentSize(layer->surfaceArray(), 1);

    // check that it's applying a rotation transform to improve phi binning
    // BOOST_CHECK_NE(bu->transform(), nullptr);
    // actAngle = ((*bu->transform()) * Vector3D(1, 0, 0)).phi();
    // expAngle = -2 * M_PI / 30 / 2.;
    // CHECK_CLOSE_REL(actAngle, expAngle, 1e-3);

    layer = std::dynamic_pointer_cast<DiscLayer>(
        p_LC->discLayer(surfaces, equidistant, equidistant, pl2));
    CHECK_CLOSE_REL(layer->thickness(), 0.4 + 2 * envZ, 1e-3);
    bounds = dynamic_cast<const RadialBounds*>(&layer->bounds());
    CHECK_CLOSE_REL(bounds->rMin(), rMin - envMinR, 1e-3);
    CHECK_CLOSE_REL(bounds->rMax(), rMax + envMaxR, 1e-3);
    BOOST_CHECK(checkBinning(*layer->surfaceArray()));
    axes = layer->surfaceArray()->getAxes();
    BOOST_CHECK_EQUAL(axes.at(0)->getNBins(), nBinsR);
    BOOST_CHECK_EQUAL(axes.at(1)->getNBins(), nBinsPhi);
    CHECK_CLOSE_REL(axes.at(0)->getMin(), rMin, 1e-3);
    CHECK_CLOSE_REL(axes.at(0)->getMax(), rMax, 1e-3);
    CHECK_CLOSE_REL(axes.at(1)->getMin(), -M_PI, 1e-3);
    CHECK_CLOSE_REL(axes.at(1)->getMax(), M_PI, 1e-3);
    checkBinContentSize(layer->surfaceArray(), 1);

    // check that it's applying a rotation transform to improve phi binning
    // BOOST_CHECK_NE(bu->transform(), nullptr);
    // actAngle = ((*bu->transform()) * Vector3D(1, 0, 0)).phi();
    // expAngle = -2 * M_PI / 30 / 2.;
    // CHECK_CLOSE_REL(actAngle, expAngle, 1e-3);
  }

  BOOST_FIXTURE_TEST_CASE(LayerCreator_barrelStagger, LayerCreatorFixture)
  {

    auto barrel = makeBarrelStagger(30, 7, 0, M_PI / 9.);
    auto brl    = barrel.first;
    draw_surfaces(brl, "LayerCreator_barrelStagger.obj");

    double     envR = 0, envZ = 0;
    ProtoLayer pl(brl);
    pl.envR = {envR, envR};
    pl.envZ = {envZ, envZ};
    std::shared_ptr<CylinderLayer> layer
        = std::dynamic_pointer_cast<CylinderLayer>(
            p_LC->cylinderLayer(brl, equidistant, equidistant, pl));

    std::cout << (*layer->surfaceArray()) << std::endl;
    auto axes = layer->surfaceArray()->getAxes();
    BOOST_CHECK_EQUAL(axes.at(0)->getNBins(), 30);
    BOOST_CHECK_EQUAL(axes.at(1)->getNBins(), 7);

    // check if binning is good!
    for (const auto& pr : barrel.second) {
      auto A = pr.first;
      auto B = pr.second;

      // std::cout << A->center().phi() << " ";
      // std::cout << B->center().phi() << std::endl;
      // std::cout << "dPHi = " << A->center().phi() - B->center().phi() <<
      // std::endl;

      Vector3D ctr        = A->binningPosition(binR);
      auto     binContent = layer->surfaceArray()->at(ctr);
      BOOST_CHECK_EQUAL(binContent.size(), 2);
      std::set<const Surface*> act;
      act.insert(binContent[0]);
      act.insert(binContent[1]);

      std::set<const Surface*> exp;
      exp.insert(A);
      exp.insert(B);
      BOOST_CHECK(exp == act);
    }

    // checkBinning should also report everything is fine
    checkBinning(*layer->surfaceArray());
  }

  BOOST_FIXTURE_TEST_CASE(LayerCreator_Cylinder_toVariantData,
                          LayerCreatorFixture)
  {

    auto barrel = makeBarrelStagger(30, 7, 0, M_PI / 9.);
    auto brl    = barrel.first;
    draw_surfaces(brl, "LayerCreator_barrelStagger.obj");

    ProtoLayer                     pl(brl);
    std::shared_ptr<CylinderLayer> layer
        = std::dynamic_pointer_cast<CylinderLayer>(
            p_LC->cylinderLayer(brl, equidistant, equidistant, pl));

    std::cout << (*layer->surfaceArray()) << std::endl;

    const variant_data var_layer = layer->toVariantData();
    // std::cout << var_layer << std::endl;

    auto layer2 = std::dynamic_pointer_cast<CylinderLayer>(
        CylinderLayer::create(var_layer));
    std::cout << (*layer2->surfaceArray()) << std::endl;

    auto sa  = layer->surfaceArray();
    auto sa2 = layer2->surfaceArray();

    BOOST_CHECK(sa);
    BOOST_CHECK(sa2);

    CHECK_CLOSE_OR_SMALL(sa->transform(), sa2->transform(), 1e-6, 1e-9);

    // let's make sure the binning is really ok
    // we check that lookup at center of input surface centers
    // gives the same number of bin content surfaces
    // which also have compatible transforms.
    // This is as close to "ok" as we can get I think.
    for (const auto& pr : barrel.second) {
      auto A = pr.first;

      Vector3D ctr = A->binningPosition(binR);

      std::vector<const Surface*> bc1 = sa->at(ctr);
      std::vector<const Surface*> bc2 = sa2->at(ctr);

      BOOST_CHECK_EQUAL(bc1.size(), bc2.size());

      for (size_t i = 0; i < bc1.size(); i++) {
        auto srf1 = bc1.at(i);
        auto srf2 = bc2.at(i);

        CHECK_CLOSE_OR_SMALL(srf1->transform(), srf2->transform(), 1e-6, 1e-9);
      }
    }
  }

  BOOST_FIXTURE_TEST_CASE(LayerCreator_Disc_toVariantData, LayerCreatorFixture)
  {
    std::vector<std::shared_ptr<const Surface>> surfaces;
    auto ringa = fullPhiTestSurfacesEC(30, 0, 0, 10);
    surfaces.insert(surfaces.end(), ringa.begin(), ringa.end());
    auto ringb = fullPhiTestSurfacesEC(30, 0, 0, 15);
    surfaces.insert(surfaces.end(), ringb.begin(), ringb.end());
    auto ringc = fullPhiTestSurfacesEC(30, 0, 0, 20);
    surfaces.insert(surfaces.end(), ringc.begin(), ringc.end());

    ProtoLayer                 pl(surfaces);
    std::shared_ptr<DiscLayer> layer = std::dynamic_pointer_cast<DiscLayer>(
        p_LC->discLayer(surfaces, equidistant, equidistant, pl));

    const variant_data var_layer = layer->toVariantData();

    std::cout << (*layer->surfaceArray()) << std::endl;
    auto layer2
        = std::dynamic_pointer_cast<DiscLayer>(DiscLayer::create(var_layer));
    std::cout << (*layer2->surfaceArray()) << std::endl;

    auto sa  = layer->surfaceArray();
    auto sa2 = layer2->surfaceArray();

    BOOST_CHECK(sa);
    BOOST_CHECK(sa2);

    CHECK_CLOSE_OR_SMALL(sa->transform(), sa2->transform(), 1e-6, 1e-9);

    for (const auto& srfRef : surfaces) {

      Vector3D ctr = srfRef->binningPosition(binR);

      std::vector<const Surface*> bc1 = sa->at(ctr);
      std::vector<const Surface*> bc2 = sa2->at(ctr);

      BOOST_CHECK_EQUAL(bc1.size(), bc2.size());

      for (size_t i = 0; i < bc1.size(); i++) {
        auto srf1 = bc1.at(i);
        auto srf2 = bc2.at(i);

        CHECK_CLOSE_OR_SMALL(srf1->transform(), srf2->transform(), 1e-6, 1e-9);
      }
    }
  }

  BOOST_AUTO_TEST_SUITE_END()
}  // namespace Test

}  // namespace Acts
