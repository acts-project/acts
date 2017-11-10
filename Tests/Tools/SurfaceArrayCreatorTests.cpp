// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE SurfaceArrayCreator
#include <boost/test/included/unit_test.hpp>

#include <boost/test/data/test_case.hpp>

#include "ACTS/Tools/SurfaceArrayCreator.hpp"
#include "ACTS/Utilities/BinningType.hpp"
#include "ACTS/Surfaces/PlaneSurface.hpp"
#include "ACTS/Surfaces/RectangleBounds.hpp"

#include "ACTS/Utilities/Definitions.hpp"

#include <fstream>
#include <random>

namespace bdata = boost::unit_test::data;
namespace tt    = boost::test_tools;

namespace Acts {

namespace Test {

  #define CHECK_CLOSE_COLLECTION(aa, bb, tolerance) { \
    using std::distance; \
    using std::begin; \
    using std::end; \
    auto a = begin(aa), ae = end(aa); \
    auto b = begin(bb); \
    BOOST_REQUIRE_EQUAL(distance(a, ae), distance(b, end(bb))); \
    for(; a != ae; ++a, ++b) { \
        BOOST_CHECK_CLOSE(*a, *b, tolerance); \
    } \
  }

  #define CHECK_ROTATION_ANGLE(t, a, tolerance) { \
    Vector3D v = (*t)*Vector3D(1, 0, 0);\
    BOOST_CHECK_SMALL(v.phi()-(a), tolerance);\
  }



  using SrfVec = std::vector<const Surface*>;

  struct SurfaceArrayCreatorFixture {
    SurfaceArrayCreator m_SAC;
    std::vector<std::unique_ptr<const Surface>> m_surfaces;

    SurfaceArrayCreatorFixture() : m_SAC(Acts::getDefaultLogger("SurfaceArrayCreator", Acts::Logging::VERBOSE)) { 
      BOOST_TEST_MESSAGE( "setup fixture" );
    }
    ~SurfaceArrayCreatorFixture() { BOOST_TEST_MESSAGE( "teardown fixture" ); }

    template<typename ...Args>
    BinUtility createBinUtility(Args&&... args) {
      return m_SAC.createBinUtility(std::forward<Args>(args)...);
    }
    
    template<typename ...Args>
    BinUtility createEquidistantBinUtility(Args&&... args) {
      return m_SAC.createEquidistantBinUtility(std::forward<Args>(args)...);
    }

    SrfVec fullPhiTestSurfacesEC(size_t n = 10, double shift = 0, double zbase = 0) {

      SrfVec res;

      double phiStep = 2*M_PI / n;
      for (size_t i=0;i<n;++i) {

        double z = zbase + ((i%2 == 0) ? 1 : -1) * 0.2;

        Transform3D trans;
        trans.setIdentity();
        trans.rotate(Eigen::AngleAxisd(i*phiStep + shift, Vector3D(0, 0, 1)));
        trans.translate(Vector3D(10, 0, z));
        
        auto bounds = std::make_shared<const RectangleBounds>(2, 1);

        auto transptr = std::make_shared<const Transform3D>(trans);
        auto srf = std::make_unique<const PlaneSurface>(transptr, bounds);

        res.push_back(srf.get()); // use raw pointer
        m_surfaces.push_back(std::move(srf)); // keep unique, will get destroyed at the end
      }

      return res;
    }

    SrfVec fullPhiTestSurfacesBRL(size_t n = 10, double shift = 0, double zbase = 0, double incl = M_PI/9.) {

      SrfVec res;

      double phiStep = 2*M_PI / n;
      for (size_t i=0;i<n;++i) {

        double z = zbase;

        Transform3D trans;
        trans.setIdentity();
        trans.rotate(Eigen::AngleAxisd(i*phiStep + shift, Vector3D(0, 0, 1)));
        trans.translate(Vector3D(10, 0, z));
        trans.rotate(Eigen::AngleAxisd(incl, Vector3D(0, 0, 1)));
        trans.rotate(Eigen::AngleAxisd(M_PI/2., Vector3D(0, 1, 0)));
        
        auto bounds = std::make_shared<const RectangleBounds>(2, 1.5);

        auto transptr = std::make_shared<const Transform3D>(trans);
        auto srf = std::make_unique<const PlaneSurface>(transptr, bounds);

        res.push_back(srf.get()); // use raw pointer
        m_surfaces.push_back(std::move(srf)); // keep unique, will get destroyed at the end
      }

      return res;
    }

  };

  void draw_surfaces(SrfVec surfaces, std::string fname)
  {
    
    std::ofstream os;
    os.open(fname);

    os << std::fixed << std::setprecision(4);

    size_t nVtx = 0;
    for (const auto &srfx : surfaces) {
      const PlaneSurface *srf = dynamic_cast<const PlaneSurface*>(srfx);
      const PlanarBounds* bounds = dynamic_cast<const PlanarBounds*>(&srf->bounds());

      for (const auto &vtxloc : bounds->vertices()) {
        Vector3D vtx = srf->transform() * Vector3D(vtxloc.x(), vtxloc.y(), 0);
        os << "v " << vtx.x() << " " << vtx.y() << " " << vtx.z() << "\n";
      }

      // connect them
      os << "f";
      for (size_t i=1;i<=bounds->vertices().size();++i) {
        os << " " << nVtx+i;
      }
      os << "\n";

      nVtx += bounds->vertices().size();
    }

    os.close();
  }
  


  BOOST_AUTO_TEST_SUITE(Tools);

  BOOST_FIXTURE_TEST_CASE( SurfaceArrayCreator_createBinUtility, SurfaceArrayCreatorFixture)
  {
    std::vector<const Surface*> srf;
    auto bu = createBinUtility(srf, BinningValue::binPhi, equidistant, 10, -M_PI, M_PI);
  }
  
  BOOST_FIXTURE_TEST_CASE( SurfaceArrayCreator_createEquidistantBinUtility, SurfaceArrayCreatorFixture)
  {
    // fail on empty srf vector
    BOOST_CHECK_THROW( createEquidistantBinUtility(SrfVec(), BinningValue::binPhi), std::logic_error );

    std::vector<float> bdExp = {
      -3.14159, -2.93215, -2.72271, -2.51327, -2.30383, -2.0944, -1.88496, -1.67552, -1.46608, 
      -1.25664, -1.0472, -0.837758, -0.628319, -0.418879, -0.20944, 0, 0.20944, 0.418879, 
      0.628319, 0.837758, 1.0472, 1.25664, 1.46608, 1.67552, 1.88496, 2.09439, 2.30383, 
      2.51327, 2.72271, 2.93215, 3.14159
    };
    
    double step = 2*M_PI / 30.;

    // endcap style modules

    for (int i=-1;i<=2;i+=2) {
      double z = 10 * i;
      // case 1: one module sits at pi / -pi
      double angleShift = step / 2.;
      auto surfaces = fullPhiTestSurfacesEC(30, angleShift, z);
      auto bu = createEquidistantBinUtility(surfaces, BinningValue::binPhi);
      auto bd = bu.binningData().at(0);
      draw_surfaces(surfaces, "SurfaceArrayCreator_createEquidistantBinUtility_EC_1.obj");
      CHECK_CLOSE_COLLECTION(bdExp, bd.boundaries(), 0.001);
      CHECK_ROTATION_ANGLE(bu.transform(), 0, 1e-3);

      // case 2: two modules sit symmetrically around pi / -pi
      angleShift = 0.;
      surfaces = fullPhiTestSurfacesEC(30, angleShift, z);
      bu = createEquidistantBinUtility(surfaces, BinningValue::binPhi);
      bd = bu.binningData().at(0);
      draw_surfaces(surfaces, "SurfaceArrayCreator_createEquidistantBinUtility_EC_2.obj");
      CHECK_CLOSE_COLLECTION(bdExp, bd.boundaries(), 0.001);
      CHECK_ROTATION_ANGLE(bu.transform(), -0.5*step, 1e-3);
      
      // case 3: two modules sit asymmetrically around pi / -pi shifted up
      angleShift = step / -4.;
      surfaces = fullPhiTestSurfacesEC(30, angleShift, z);
      bu = createEquidistantBinUtility(surfaces, BinningValue::binPhi);
      bd = bu.binningData().at(0);
      draw_surfaces(surfaces, "SurfaceArrayCreator_createEquidistantBinUtility_EC_3.obj");
      CHECK_CLOSE_COLLECTION(bdExp, bd.boundaries(), 0.001);
      CHECK_ROTATION_ANGLE(bu.transform(), step / -4., 1e-3);
      
      // case 4: two modules sit asymmetrically around pi / -pi shifted down
      angleShift = step / 4.;
      surfaces = fullPhiTestSurfacesEC(30, angleShift, z);
      bu = createEquidistantBinUtility(surfaces, BinningValue::binPhi);
      bd = bu.binningData().at(0);
      draw_surfaces(surfaces, "SurfaceArrayCreator_createEquidistantBinUtility_EC_4.obj");
      CHECK_CLOSE_COLLECTION(bdExp, bd.boundaries(), 0.001);
      CHECK_ROTATION_ANGLE(bu.transform(), step / 4., 1e-3);
    }
    

    for (int i=-1;i<=2;i+=2) {
      double z = 10 * i;
      // case 1: one module sits at pi / -pi
      double angleShift = step / 2.;
      auto surfaces = fullPhiTestSurfacesBRL(30, angleShift, z);
      auto bu = createEquidistantBinUtility(surfaces, BinningValue::binPhi);
      auto bd = bu.binningData().at(0);
      draw_surfaces(surfaces, "SurfaceArrayCreator_createEquidistantBinUtility_BRL_1.obj");
      CHECK_CLOSE_COLLECTION(bdExp, bd.boundaries(), 0.001);
      CHECK_ROTATION_ANGLE(bu.transform(), 0, 1e-3);

      // case 2: two modules sit symmetrically around pi / -pi
      angleShift = 0.;
      surfaces = fullPhiTestSurfacesBRL(30, angleShift, z);
      bu = createEquidistantBinUtility(surfaces, BinningValue::binPhi);
      bd = bu.binningData().at(0);
      draw_surfaces(surfaces, "SurfaceArrayCreator_createEquidistantBinUtility_BRL_2.obj");
      CHECK_CLOSE_COLLECTION(bdExp, bd.boundaries(), 0.001);
      CHECK_ROTATION_ANGLE(bu.transform(), -0.5*step, 1e-3);
      
      // case 3: two modules sit asymmetrically around pi / -pi shifted up
      angleShift = step / -4.;
      surfaces = fullPhiTestSurfacesBRL(30, angleShift, z);
      bu = createEquidistantBinUtility(surfaces, BinningValue::binPhi);
      bd = bu.binningData().at(0);
      draw_surfaces(surfaces, "SurfaceArrayCreator_createEquidistantBinUtility_BRL_3.obj");
      CHECK_CLOSE_COLLECTION(bdExp, bd.boundaries(), 0.001);
      CHECK_ROTATION_ANGLE(bu.transform(), step / -4., 1e-3);
      
      // case 4: two modules sit asymmetrically around pi / -pi shifted down
      angleShift = step / 4.;
      surfaces = fullPhiTestSurfacesBRL(30, angleShift, z);
      bu = createEquidistantBinUtility(surfaces, BinningValue::binPhi);
      bd = bu.binningData().at(0);
      draw_surfaces(surfaces, "SurfaceArrayCreator_createEquidistantBinUtility_BRL_4.obj");
      CHECK_CLOSE_COLLECTION(bdExp, bd.boundaries(), 0.001);
      CHECK_ROTATION_ANGLE(bu.transform(), step / 4., 1e-3);
    }

    /*
    std::mt19937 gen(42);
    std::uniform_real_distribution<> dis(-M_PI, M_PI);

    double lo = bd.boundaries().front();
    double up = bd.boundaries().back();

    auto check = [&bu, &up, &lo](double phi) {
      Vector3D pos( 10*std::cos(phi), 10*std::sin(phi), 0);
      auto bt = bu.binTriple(pos);

      return bt;
    };
    
    size_t n = 1e7;
    size_t n5 = n / 20;
    for(size_t i=0;i<n;i++) {
      double phi = dis(gen);

      check(phi);
      
      if(i%n5 == 0) {
        std::cout << ( i / double(n) * 100. ) << "%" << std::endl;
      }

    }
    */

  }

  BOOST_AUTO_TEST_SUITE_END();
}  // end of namespace Test

}  // end of namespace Acts
