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
#include <typeinfo>

#include "Acts/Utilities/Logger.hpp"

BOOST_AUTO_TEST_SUITE(GeoModelDetObj)

BOOST_AUTO_TEST_CASE(GeoModelDetectorObjectFactory) {

  // Let's create a GeoFullPhysVol object
  auto material = new GeoMaterial("Material", 1.0);
  double gmhlx = 100;
  double gmhly = 200;
  double gmhlz = 2;
  auto boxXY = new GeoBox(gmhlx, gmhly, gmhlz);
  auto logXY = new GeoLogVol("LogVolumeXY", boxXY, material);
  auto fphysXY = new GeoFullPhysVol(logXY);
  PVConstLink physVol{fphysXY};
  auto rBounds = std::make_shared<Acts::RectangleBounds>(100, 200);

  //create pars for conversion
  Acts::GeoModelDetectorObjectFactory::Config gmConfig;
  Acts::GeometryContext gContext;
  Acts::GeoModelDetectorObjectFactory::Cache gmCache;

  //create factory instance
  Acts::GeoModelDetectorObjectFactory factory = Acts::GeoModelDetectorObjectFactory(gmConfig);

  //convert GeoFullPhysVol
  factory.convertFpv("LogVolumeXY", fphysXY, gmCache, gContext);


  //perform checks
  //TODO check achlz
  for (auto surface : gmCache.sensitiveSurfaces){
    auto ss = std::get<1>(surface);
    BOOST_CHECK(ss->type() == Acts::Surface::SurfaceType::Plane);
    const Acts::SurfaceBounds& bounds = ss->bounds();
    const Acts::RectangleBounds* rectBounds = dynamic_cast<const Acts::RectangleBounds*>(&bounds);

    double achlx = rectBounds->halfLengthX();
    double achly = rectBounds->halfLengthY();
    //double achlz = rectBounds->halfLengthZ();
    BOOST_CHECK(gmhlx == achlx);
    BOOST_CHECK(gmhly == achly);
    //BOOST_CHECK(gmhlz == achlz);
  }
}

BOOST_AUTO_TEST_SUITE_END()
