#include "Acts/Plugins/GeoModel/GeoModelDetectorVolumeFactory.hpp"

#include "Acts/Detector/GeometryIdGenerator.hpp"
#include "Acts/Detector/PortalGenerators.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/CutoutCylinderVolumeBounds.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/TrapezoidVolumeBounds.hpp"
#include "Acts/Navigation/InternalNavigation.hpp"

#include "Acts/Plugins/GeoModel/IGeoShapeConverter.hpp"
#include <GeoModelKernel/GeoBox.h>
#include <GeoModelKernel/GeoPcon.h>
#include <GeoModelKernel/GeoShapeShift.h>
#include <GeoModelKernel/GeoShapeSubtraction.h>
#include <GeoModelKernel/GeoShapeUnion.h>
#include <GeoModelKernel/GeoSimplePolygonBrep.h>
#include <GeoModelKernel/GeoTrd.h>
#include <GeoModelKernel/GeoTube.h>
#include <GeoModelKernel/GeoTubs.h>
#include <iostream>


namespace Acts {
//TODO adapt constructor, build construct function, 
Acts::GeoModelDetectorVolumeFactory::GeoModelDetectorVolumeFactory(const Config& cfg, std::unique_ptr<const Logger> mlogger){}
  

void Acts::GeoModelDetectorVolumeFactory::construct(Cache& cache, const GeometryContext& gctx, const GeoModelTree& geoModelTree, const Options& options){
  for (const auto &q : options.queries) {
    std::cout << q << std::endl;
  }
}


/*
void Acts::GeoModelDetectorVolumeFactory::print(){
  std::cout << "sdfdsf " << name << std::endl;
}
*/

}
