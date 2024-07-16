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
#include <typeinfo>


namespace Acts {
//TODO build construct function, 
Acts::GeoModelDetectorVolumeFactory::GeoModelDetectorVolumeFactory(const Config& cfg, std::unique_ptr<const Logger> mlogger) : m_cfg(cfg), m_logger(std::move(mlogger)) {}
  

void Acts::GeoModelDetectorVolumeFactory::construct(Cache& cache, const GeometryContext& gctx, const GeoModelTree& geoModelTree, const Options& options){
  if (geoModelTree.geoReader == nullptr) {
    throw std::invalid_argument("GeoModelTree has no GeoModelReader");
  }
  for (const auto &q : options.queries) {
    //ACTS_VERBOSE("Constructing detector elements for query " << q);
    //load data from database according to querie (Muon)
    auto qFPV = geoModelTree.geoReader->getPublishedNodes<std::string, GeoFullPhysVol *>(q);
    //TODO apply some conditions

    auto matches = [&](const std::string &name, PVConstLink physvol) {

      if (m_cfg.nameList.empty() && m_cfg.materialList.empty()) {
        return true;
      }

      auto matchName = std::any_of(
          m_cfg.nameList.begin(), m_cfg.nameList.end(),
          [&](const auto &n) { return name.find(n) != std::string::npos; });

      std::string matStr = physvol->getLogVol()->getMaterial()->getName();

      auto matchMaterial = std::any_of(
          m_cfg.materialList.begin(), m_cfg.materialList.end(),
          [&](const auto &m) { return matStr.find(m) != std::string::npos; });

      bool match = matchMaterial && matchName;
      GeoIntrusivePtr<const GeoVFullPhysVol> fullVol =
          dynamic_pointer_cast<const GeoVFullPhysVol>(physvol);

      // for the fullphysvol we only check the name

      if (m_cfg.nameList.empty()) {
        return matchMaterial;
      }

      if (m_cfg.materialList.empty() || !(fullVol == nullptr)) {
        return matchName;
      }

      return match;
    };









    //go through each fpv
    for (auto &[name, fpv] : qFPV) {
      PVConstLink physVol{fpv};
      if (!matches(name, physVol)) {
        continue;
      }
      const GeoLogVol *logVol = physVol->getLogVol();//get logVol for the shape of the volume
      const GeoShape *shape = logVol->getShape();//get shape
      const Acts::Transform3 &transform = fpv->getAbsoluteTransform(nullptr);

      //use GeoModelToDetVol to get the bounding boxes
      std::shared_ptr<Experimental::DetectorVolume> box = Acts::GeoModel::convertVolume(gctx, shape, name, transform);
      //TODO save to bounding boxes in a way that is accessible from py script
      cache.boundingBoxes.push_back(box);
    }
  }
}

}
