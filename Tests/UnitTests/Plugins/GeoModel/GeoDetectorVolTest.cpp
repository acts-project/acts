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
  auto material = new GeoMaterial("Material", 1.0);
  // Let's create a GeoFullPhysVol object

  // (BOX object)
  auto boxXY = new GeoBox(100, 200, 2);
  auto logXY = new GeoLogVol("LogVolumeXY", boxXY, material);
  auto fphysXY = new GeoFullPhysVol(logXY);
  PVConstLink physVol{fphysXY};
  auto rBounds = std::make_shared<Acts::RectangleBounds>(100, 200);
  //int shapeId = shape->typeID();
  //std::cout << "!!!!!" <<physVol->typeID() << std::endl;

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
  //TODO find the right data type for the surface converter
  //factory.convertFpv("LogVolumeXY", fphysXY, gmCache, gContext);
}

BOOST_AUTO_TEST_SUITE_END()
