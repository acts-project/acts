#include "Acts/Plugins/GeoModel/GeoModelDetectorObjectFactory.hpp"

#include "Acts/Detector/GeometryIdGenerator.hpp"
#include "Acts/Detector/PortalGenerators.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/CutoutCylinderVolumeBounds.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/TrapezoidVolumeBounds.hpp"
#include "Acts/Navigation/InternalNavigation.hpp"

#include <GeoModelKernel/GeoBox.h>
#include <GeoModelKernel/GeoPcon.h>
#include <GeoModelKernel/GeoShapeShift.h>
#include <GeoModelKernel/GeoShapeSubtraction.h>
#include <GeoModelKernel/GeoShapeUnion.h>
#include <GeoModelKernel/GeoSimplePolygonBrep.h>
#include <GeoModelKernel/GeoTrd.h>
#include <GeoModelKernel/GeoTube.h>
#include <GeoModelKernel/GeoTubs.h>
#include "Acts/Plugins/GeoModel/GeoModelConverters.hpp"
#include "Acts/Plugins/GeoModel/IGeoShapeConverter.hpp"
#include <iostream>
#include <typeinfo>

namespace {
std::string recType(const GeoShapeShift &gshift);
std::string recType(const GeoShapeUnion &gunion);
std::string recType(const GeoShape &gshape);

std::string recType(const GeoShapeShift &gshift) {
  return "Shift[" + recType(*gshift.getOp()) + "]";
}
std::string recType(const GeoShapeUnion &gunion) {
  return "Union[" + recType(*gunion.getOpA()) + ", " +
         recType(*gunion.getOpB()) + "]";
}
std::string recType(const GeoShape &gshape) {
  if (auto ptr = dynamic_cast<const GeoShapeUnion *>(&gshape); ptr != nullptr) {
    return recType(*ptr);
  }
  if (auto ptr = dynamic_cast<const GeoShapeShift *>(&gshape); ptr != nullptr) {
    return recType(*ptr);
  }
  return gshape.type();
}

}  // namespace

namespace Acts {
Acts::GeoModelDetectorObjectFactory::GeoModelDetectorObjectFactory(
    const Config& cfg, std::unique_ptr<const Logger> mlogger)
    : m_cfg(cfg), m_logger(std::move(mlogger)) {}
  

void Acts::GeoModelDetectorObjectFactory::construct(Cache& cache, const GeometryContext& gctx, const GeoModelTree& geoModelTree, const Options& options){
  if (geoModelTree.geoReader == nullptr) {
    throw std::invalid_argument("GeoModelTree has no GeoModelReader");
  }
  for (const auto &q : options.queries) {
    ACTS_VERBOSE("Constructing detector elements for query " << q);
    //load data from database according to querie (Muon)
    auto qFPV = geoModelTree.geoReader->getPublishedNodes<std::string, GeoFullPhysVol *>(q);

    //lambda to determine if object fits querry
    auto matches = [&](const std::string &name, PVConstLink physvol) {

      if (m_cfg.nameList.empty() && m_cfg.materialList.empty()) {
        return true;
      }

      auto matchName = std::any_of(
          m_cfg.nameList.begin(), m_cfg.nameList.end(),
          [&](const auto &n) {return name.find(n) != std::string::npos; });

      std::string matStr = physvol->getLogVol()->getMaterial()->getName();

      auto matchMaterial = std::any_of(
          m_cfg.materialList.begin(), m_cfg.materialList.end(),
          [&](const auto &m) {return matStr.find(m) != std::string::npos; });

      bool match = matchMaterial && matchName;
      GeoIntrusivePtr<const GeoVFullPhysVol> fullVol =
          dynamic_pointer_cast<const GeoVFullPhysVol>(physvol);

      // for the fullphysvol we only check the name
      if (m_cfg.nameList.empty()) {
        return matchMaterial;
      }

      //if no material specified or we're looking at fpv judge by name only
      if (m_cfg.materialList.empty() || !(fullVol == nullptr)) {
        return matchName;
      }
      return match;
    };

    //go through each fpv
    for (auto &[name, fpv] : qFPV) {
      PVConstLink physVol{fpv};
      //if the match lambda returns false skip the rest of the loop
      if (!matches(name, physVol)) {
        continue;
      }


      //get children
      std::vector<GeoChildNodeWithTrf> subvolumes = getChildrenWithRef(physVol, false);
      if(subvolumes.size()>0){

        //vector containing all subvolumes to be converted to surfaces
        std::vector<GeoChildNodeWithTrf> surfaces = findAllSubVolumes(physVol, matches);
        std::vector<GeoModelSensitiveSurface> sensitives;
  
        for (auto surface : surfaces){
          const Transform3 &transform = fpv->getAbsoluteTransform() * surface.transform;
          convertSensitive(surface.volume, transform, sensitives);
        }
        cache.sensitiveSurfaces.insert(cache.sensitiveSurfaces.end(), sensitives.begin(), sensitives.end());
        if(convertFpv(name)){
          const GeoLogVol *logVol = physVol->getLogVol();//get logVol for the shape of the volume
          const GeoShape *shape = logVol->getShape();//get shape
          const Acts::Transform3 &fpvtransform = fpv->getAbsoluteTransform(nullptr);

          //convert bounding boxes with surfaces inside
          std::shared_ptr<Experimental::DetectorVolume> box = Acts::GeoModel::convertVolume(gctx, shape, name, fpvtransform, sensitives);
          cache.boundingBoxes.push_back(box);
        }
      }
      else{
        //convert fpvs to surfaces
        const Transform3 &transform = fpv->getAbsoluteTransform();
        convertSensitive(fpv, transform, cache.sensitiveSurfaces);
      }
    }
  }
}

void Acts::GeoModelDetectorObjectFactory::convertSensitive(PVConstLink geoPV, const Acts::Transform3 &transform, std::vector<GeoModelSensitiveSurface> &sensitives) {
  const GeoLogVol *logVol = geoPV->getLogVol();
  const GeoShape *shape = logVol->getShape();
  int shapeId = shape->typeID();
  std::string name = logVol->getName();
  std::shared_ptr<const Acts::IGeoShapeConverter> converter = Acts::GeoShapesConverters(shapeId);
  if (converter == nullptr) {
    throw std::runtime_error("The converter for " + recType(*shape) +
                             " is nullptr");

    return;
  }
  auto converted = converter->toSensitiveSurface(geoPV, transform);
  if (converted.ok()) {
    sensitives.push_back(converted.value());
    const auto &[el, sf] = converted.value();

    ACTS_VERBOSE("(successfully converted: " << name << " / " << recType(*shape) << " / " << logVol->getMaterial()->getName() << ")");

    if (!el || !sf) {
      throw std::runtime_error("The Detector Element or the Surface is nllptr");
    }
    return;
  }
  ACTS_ERROR(name << " / " << recType(*shape)<< ") could not be converted by any converter");
}

std::vector<GeoChildNodeWithTrf> Acts::GeoModelDetectorObjectFactory::findAllSubVolumes(PVConstLink vol, std::function<bool(std::string, PVConstLink)> matchFunc){
  std::vector<GeoChildNodeWithTrf> subvolumes = getChildrenWithRef(vol, false);
  std::vector<GeoChildNodeWithTrf> sensitives;
  for (auto subvolume : subvolumes){
    if (matchFunc(subvolume.nodeName, subvolume.volume)) {
      sensitives.push_back(subvolume);
    }
    std::vector<GeoChildNodeWithTrf> senssubsubvolumes = findAllSubVolumes(subvolume.volume, matchFunc);
  std::transform(std::make_move_iterator(senssubsubvolumes.begin()),
                 std::make_move_iterator(senssubsubvolumes.end()),
                 std::back_inserter(sensitives),
                 [&subvolume](GeoChildNodeWithTrf &&vol) {
                   vol.transform = subvolume.transform * vol.transform;
                   return vol;
                 });
  }
  return sensitives;
}
bool Acts::GeoModelDetectorObjectFactory::convertBox(std::string name){
  return (name.find(m_cfg.convertBox) != std::string::npos);
}
void Acts::GeoModelDetectorObjectFactory::convertFpv(std::string name, PVConstLink fpv, Cache& cache, const GeometryContext& gctx){
  std::cout << "pars work" << std::endl;

      //get children
      std::vector<GeoChildNodeWithTrf> subvolumes = getChildrenWithRef(fpv, false);
      if(subvolumes.size()>0){
        std::cout << "subvols" << std::endl;

        //vector containing all subvolumes to be converted to surfaces
        std::vector<GeoChildNodeWithTrf> surfaces = findAllSubVolumes(fpv);
        std::vector<GeoModelSensitiveSurface> sensitives;
  
        for (auto surface : surfaces){
          const Transform3 &transform = fpv->getAbsoluteTransform() * surface.transform;
          convertSensitive(surface.volume, transform, sensitives);
        }
        cache.sensitiveSurfaces.insert(cache.sensitiveSurfaces.end(), sensitives.begin(), sensitives.end());
        if(convertBox(name)){
          const GeoLogVol *logVol = fpv->getLogVol();//get logVol for the shape of the volume
          const GeoShape *shape = logVol->getShape();//get shape
          const Acts::Transform3 &fpvtransform = fpv->getAbsoluteTransform(nullptr);

          //convert bounding boxes with surfaces inside
          std::shared_ptr<Experimental::DetectorVolume> box = Acts::GeoModel::convertVolume(gctx, shape, name, fpvtransform, sensitives);
          cache.boundingBoxes.push_back(box);
        }
      }
      else{
        std::cout << "no subvols" << std::endl;
        //convert fpvs to surfaces
        const Transform3 &transform = fpv->getAbsoluteTransform();
        std::cout << "transforms work" << std::endl;
        convertSensitive(fpv, transform, cache.sensitiveSurfaces);
        std::cout << "converts work" << std::endl;
      }
}
}
