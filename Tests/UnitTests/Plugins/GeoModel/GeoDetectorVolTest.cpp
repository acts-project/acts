#include <boost/test/unit_test.hpp>
#include "Acts/Plugins/GeoModel/GeoModelDetectorObjectFactory.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Plugins/GeoModel/GeoModelReader.hpp"

#include <GeoModelKernel/GeoBox.h>
#include <GeoModelKernel/GeoFullPhysVol.h>
#include <GeoModelKernel/GeoLogVol.h>
#include <GeoModelKernel/GeoMaterial.h>

#include "Acts/Utilities/Logger.hpp"

BOOST_AUTO_TEST_SUITE(GeoModelDetObj)

BOOST_AUTO_TEST_CASE(GeoModelDetectorObjectFactory) {
  /*
  auto material = new GeoMaterial("Material", 1.0);
  // Let's create a GeoFullPhysVol object

  // (BOX object)
  auto boxXY = new GeoBox(100, 200, 2);
  auto logXY = new GeoLogVol("LogVolumeXY", boxXY, material);
  auto fphysXY = new GeoFullPhysVol(logXY);
  auto rBounds = std::make_shared<Acts::RectangleBounds>(100, 200);
  */

  //create pars for constructor
  //std::unique_ptr<const Actst::Logger> log;
  Acts::GeoModelDetectorObjectFactory::Config gmConfig;
  gmConfig.materialList = {"Aluminium"};
  //gmConfig.nameList = {};
  //gmConfig.convertFpv = {};
  //gmConfig.convertSubVolumnes = {};

  //create factory instance
  Acts::GeoModelDetectorObjectFactory factory = Acts::GeoModelDetectorObjectFactory(gmConfig);

  std::cout <<"!!!!!!!! test" << std::endl;
  //convert shape
  Acts::GeoModelDetectorObjectFactory::Options gmOptions;
  gmOptions.queries = {"Muon"};
  Acts::GeoModelDetectorObjectFactory::Cache gmCache;
  Acts::GeometryContext gContext;
  std::cout <<"!!!!!!!! test cache" << std::endl;
  Acts::GeoModelTree gmTree = Acts::GeoModelReader::readFromDb("ATLAS-R3-MUONTEST_v3.db");
  std::cout <<"!!!!!!!! test Db" << std::endl;
  //factory.construct(gmCache, gContext, gmTree, gmOptions);
  std::cout <<"!!!!!!!! test2" << std::endl;


  //BOOST_CHECK(gmCache.sensitiveSurfaces.size() == 100);
  std::cout <<"!!!!!!!! test3" << std::endl;

  /*
  PVConstLink physXY{fphysXY};
  auto elementXY =
      Acts::GeoModelDetectorElement::createDetectorElement<Acts::PlaneSurface>(
          physXY, rBounds, Acts::Transform3::Identity(), 2.0);

  const Acts::Surface& surface = elementXY->surface();
  BOOST_CHECK(surface.type() == Acts::Surface::SurfaceType::Plane);
  */

}

BOOST_AUTO_TEST_SUITE_END()
