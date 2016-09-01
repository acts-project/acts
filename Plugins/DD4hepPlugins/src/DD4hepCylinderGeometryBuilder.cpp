// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTS/Plugins/DD4hepPlugins/DD4hepCylinderGeometryBuilder.hpp"
#include "ACTS/Detector/TrackingGeometry.hpp"
#include "ACTS/Detector/TrackingVolume.hpp"
#include "ACTS/Layers/CylinderLayer.hpp"
#include "ACTS/Layers/DiscLayer.hpp"
#include "ACTS/Layers/GenericApproachDescriptor.hpp"
#include "ACTS/Material/Material.hpp"
#include "ACTS/Material/SurfaceMaterialProxy.hpp"
#include "ACTS/Plugins/DD4hepPlugins/DD4hepDetElement.hpp"
#include "ACTS/Plugins/DD4hepPlugins/DetExtension.hpp"
#include "ACTS/Plugins/DD4hepPlugins/IDetExtension.hpp"
#include "ACTS/Surfaces/CylinderBounds.hpp"
#include "ACTS/Surfaces/RadialBounds.hpp"
#include "ACTS/Surfaces/Surface.hpp"
#include "ACTS/Surfaces/TrapezoidBounds.hpp"
#include "ACTS/Utilities/BinnedArrayXD.hpp"
#include "ACTS/Utilities/Units.hpp"
#include "ACTS/Volumes/CylinderVolumeBounds.hpp"
#include "ACTS/Volumes/Volume.hpp"
#include "TGeoManager.h"
#include "TGeoMatrix.h"

using namespace Acts::units;

Acts::DD4hepCylinderGeometryBuilder::DD4hepCylinderGeometryBuilder(
    const Config            dgbConfig,
    std::unique_ptr<Logger> logger)
  : m_logger(std::move(logger))
{
  setConfiguration(std::move(dgbConfig));
}

// configuration
void
Acts::DD4hepCylinderGeometryBuilder::setConfiguration(
    const Acts::DD4hepCylinderGeometryBuilder::Config dgbConfig)
{
  // @TODO check consistency
  // copy the configuration
  m_cfg = std::move(dgbConfig);
}

void
Acts::DD4hepCylinderGeometryBuilder::setLogger(
    std::unique_ptr<Logger> newLogger)
{
  m_logger = std::move(newLogger);
}

std::unique_ptr<Acts::TrackingGeometry>
Acts::DD4hepCylinderGeometryBuilder::trackingGeometry() const
{
  ACTS_INFO("Translating DD4hep geometry into ACTS geometry");
  // the return geometry -- and the highest volume
  std::unique_ptr<Acts::TrackingGeometry> trackingGeometry = nullptr;
  Acts::TrackingVolumePtr                 highestVolume    = nullptr;
  Acts::TrackingVolumePtr                 beamPipeVolume   = nullptr;

  // get the sub detectors
  std::vector<DD4hep::Geometry::DetElement>     detElements;
  const DD4hep::Geometry::DetElement::Children& children
      = m_cfg.detWorld.children();
  for (auto& detElement : children) detElements.push_back(detElement.second);
  // sort by id to build detector from bottom to top
  sort(detElements.begin(),
       detElements.end(),
       [](const DD4hep::Geometry::DetElement& a,
          const DD4hep::Geometry::DetElement& b) { return (a.id() < b.id()); });
  // loop over the volumes
  for (auto& detElement : detElements) {
    if (detElement.type() == "beamtube") {
      ACTS_VERBOSE("BeamPipe is being built");
      // extract material
      DD4hep::Geometry::Material mat = detElement.volume().material();
      // create the tracking volume
      beamPipeVolume = Acts::TrackingVolume::create(
          convertTransform(detElement.placement().ptr()->GetMatrix()),
          extractVolumeBounds(detElement),
          std::make_shared<Acts::Material>(mat.radLength(),
                                           mat.intLength(),
                                           mat.A(),
                                           mat.Z(),
                                           mat.density()),
          nullptr,
          {},
          {},
          {},
          "BeamTube");
    } else {
      // assign a new highest volume (and potentially wrap around the given
      // highest volume so far)
      LayerTriple  layerTriple;
      VolumeTriple volumeTriple;
      createSubVolumes(detElement, layerTriple, volumeTriple);
      highestVolume
          = m_cfg.volumeBuilder->trackingVolume(highestVolume,
                                                extractVolumeBounds(detElement),
                                                new LayerTriple(layerTriple),
                                                new VolumeTriple(volumeTriple));
    }
  }
  // if you have a highest volume, stuff it into a TrackingGeometry
  if (highestVolume) {
    // see if the beampipe needs to be wrapped
    if (beamPipeVolume)
      highestVolume = m_cfg.volumeHelper->createContainerTrackingVolume(
          {beamPipeVolume, highestVolume});
    // create the TrackingGeometry
    trackingGeometry = std::make_unique<Acts::TrackingGeometry>(highestVolume);
  }
  // return the geometry to the service
  return trackingGeometry;
}

std::shared_ptr<Acts::Transform3D>
Acts::DD4hepCylinderGeometryBuilder::convertTransform(
    const TGeoMatrix* tGeoTrans) const
{
  // get the placement and orientation in respect to its mother
  const Double_t* rotation    = tGeoTrans->GetRotationMatrix();
  const Double_t* translation = tGeoTrans->GetTranslation();
  auto            transform   = std::make_shared<Acts::Transform3D>(
      Acts::Vector3D(rotation[0], rotation[3], rotation[6]),
      Acts::Vector3D(rotation[1], rotation[4], rotation[7]),
      Acts::Vector3D(rotation[2], rotation[5], rotation[8]),
      Acts::Vector3D(
          translation[0] * _cm, translation[1] * _cm, translation[2] * _cm));
  return (transform);
}

std::shared_ptr<const Acts::VolumeBounds>
Acts::DD4hepCylinderGeometryBuilder::extractVolumeBounds(
    DD4hep::Geometry::DetElement& detElement) const
{
  TGeoShape*   geoShape = detElement.placement().ptr()->GetVolume()->GetShape();
  TGeoConeSeg* tube     = dynamic_cast<TGeoConeSeg*>(geoShape);
  if (!tube) ACTS_ERROR("Volume has wrong shape - needs to be TGeoConeSeg!");
  // get the dimension of TGeo and convert lengths
  double rMin  = tube->GetRmin1() * _cm;
  double rMax  = tube->GetRmax1() * _cm;
  double halfZ = tube->GetDz() * _cm;
  ACTS_DEBUG(
      "[V] Extracting cylindrical volume bounds ( rmin / rmax / halfZ )=  ( "
      << rMin
      << " / "
      << rMax
      << " / "
      << halfZ
      << " )");
  auto cylinderBounds
      = std::make_shared<const Acts::CylinderVolumeBounds>(rMin, rMax, halfZ);
  return cylinderBounds;
}

void
Acts::DD4hepCylinderGeometryBuilder::createSubVolumes(
    DD4hep::Geometry::DetElement& detElement,
    LayerTriple&                  layerTriple,
    VolumeTriple&                 volumeTriple) const
{
  // the negative endcap volume of the current hierarchy
  VolumePtr nEndcapVolume = nullptr;
  // the barrel volume of the current hierarchy
  VolumePtr barrelVolume = nullptr;
  // the positive endcap volume of the current hierarchy
  VolumePtr pEndcapVolume = nullptr;
  // possible layers of the negative end cap
  Acts::LayerVector negativeLayers;
  // possible layers of the central barrel
  Acts::LayerVector centralLayers;
  // possible layers of the positive end cap
  Acts::LayerVector positiveLayers;

  if (detElement.type() == "compound") {
    ACTS_VERBOSE("[V] Volume : "
                 << detElement.name()
                 << " is a compound volume -> resolve now the sub volumes");
    // create tracking volume of compound type
    const DD4hep::Geometry::DetElement::Children& compoundChildren
        = detElement.children();
    for (auto& compoundChild : compoundChildren) {
      DD4hep::Geometry::DetElement compoundDetElement = compoundChild.second;
      // extract the transformation
      TGeoMatrix* transform = compoundDetElement.placement().ptr()->GetMatrix();
      std::shared_ptr<Acts::Transform3D> actsTransform
          = convertTransform(transform);
      // distinguish between TGeoConeSeg used as a cylinder (barrel) and as a
      // disc (end caps)
      Acts::IDetExtension* detExtension
          = compoundDetElement.extension<Acts::IDetExtension>();
      // create disc layers in case of a disc volume, otherwise create
      // cylindrical layers
      if (detExtension->shape() == Acts::ShapeType::Disc) {
        ACTS_VERBOSE("[V] Subvolume : "
                     << compoundDetElement.name()
                     << " is a disc volume -> handling as an endcap");

        if (actsTransform->translation().z() < 0.) {
          ACTS_VERBOSE("[V]       ->is negative endcap");
          nEndcapVolume = std::make_shared<const Volume>(
              actsTransform, extractVolumeBounds(compoundDetElement));
          createDiscLayers(compoundDetElement, negativeLayers, transform);
        } else {
          ACTS_VERBOSE("[V]       ->is positive endcap");
          pEndcapVolume = std::make_shared<const Volume>(
              actsTransform, extractVolumeBounds(compoundDetElement));
          createDiscLayers(compoundDetElement, positiveLayers, transform);
        }
      } else {
        ACTS_VERBOSE("[V] Subvolume : "
                     << compoundDetElement.name()
                     << " is a cylinder volume -> handling as a barrel");
        barrelVolume = std::make_shared<const Volume>(
            actsTransform, extractVolumeBounds(compoundDetElement));
        createCylinderLayers(compoundDetElement, centralLayers, transform);
      }
    }  // compoundchildren
  }    // compoundtype
  else {
    ACTS_VERBOSE("[V] Volume : "
                 << detElement.name()
                 << " is not of compound type -> handling as a barrel");
    // support structure
    // create cylindrical layers
    createCylinderLayers(
        detElement, centralLayers, detElement.placement().ptr()->GetMatrix());
  }
  volumeTriple = VolumeTriple(nEndcapVolume, barrelVolume, pEndcapVolume);
  // set the triples
  layerTriple = LayerTriple(Acts::LayerVector(negativeLayers),
                            Acts::LayerVector(centralLayers),
                            Acts::LayerVector(positiveLayers));
}

void
Acts::DD4hepCylinderGeometryBuilder::createCylinderLayers(
    DD4hep::Geometry::DetElement& motherDetElement,
    Acts::LayerVector&            layers,
    const TGeoMatrix*             motherTransform) const
{
  // get possible layers
  const DD4hep::Geometry::DetElement::Children& children
      = motherDetElement.children();
  // check if volume has layers
  if (!children.empty()) {
    ACTS_VERBOSE("[V] Volume containes layers -> creating cylindrical layers");
    // go through layers
    for (auto& child : children) {
      // get the detector element of the layer
      DD4hep::Geometry::DetElement detElement = child.second;
      // build the matrix
      TGeoMatrix* locTrans = detElement.placement().ptr()->GetMatrix();
      // make the transformation global
      TGeoHMatrix* transform
          = new TGeoHMatrix((*motherTransform) * (*locTrans));
      // get the shape of the layer
      TGeoShape* geoShape
          = detElement.placement().ptr()->GetVolume()->GetShape();
      TGeoConeSeg* tube = dynamic_cast<TGeoConeSeg*>(geoShape);
      if (!tube)
        ACTS_ERROR(
            "[L] Cylinder layer has wrong shape - needs to be TGeoConeSeg!");
      // extract the boundaries
      double rMin           = tube->GetRmin1() * _cm;
      double rMax           = tube->GetRmax1() * _cm;
      double halfZ          = tube->GetDz() * _cm;
      auto   cylinderBounds = std::make_shared<const Acts::CylinderBounds>(
          0.5 * (rMin + rMax), halfZ);
      double thickness = fabs(rMin - rMax);
      ACTS_DEBUG("[L] Creating cylinder layer with dimensions ( rmin / rmax / "
                 "halfZ ) = ( "
                 << rMin
                 << " / "
                 << rMax
                 << " / "
                 << halfZ
                 << " )");
      // if necessary receive the modules contained by the layer and create the
      // layer, otherwise create an empty layer
      Acts::IDetExtension* detExtension
          = detElement.extension<Acts::IDetExtension>();
      std::vector<DD4hep::Geometry::DetElement> modules(
          detExtension->modules());
      // access the axis orienation of the modules
      const std::string axes = detExtension->axes();
      // create the two dimensional BinUtility for the material map of the layer
      Acts::BinUtility*                           materialBinUtil = nullptr;
      std::shared_ptr<const SurfaceMaterialProxy> materialProxy(nullptr);
      // the approachdescriptor telling where the material sits on the layer
      Acts::ApproachDescriptor* approachDescriptor = nullptr;
      // material position on the layer can be inner, outer or center and will
      // be accessed from the Extensions
      Acts::LayerMaterialPos layerPos = LayerMaterialPos::inner;
      // access transformation matrix of the layer
      std::shared_ptr<Acts::Transform3D> actsTransform
          = convertTransform(transform);
      // check if layer should have material
      if (detExtension->hasSupportMaterial()) {
        std::pair<size_t, size_t> materialBins = detExtension->materialBins();
        size_t bins1    = materialBins.first;
        size_t bins2    = materialBins.second;
        materialBinUtil = new Acts::BinUtility(
            bins1, -M_PI, M_PI, Acts::closed, Acts::binPhi);
        (*materialBinUtil) += Acts::BinUtility(
            bins2, -halfZ, halfZ, Acts::open, Acts::binZ, actsTransform);
        // and create material proxy to mark layer for material mapping
        materialProxy
            = std::make_shared<const SurfaceMaterialProxy>(*materialBinUtil);
        // access the material position
        Acts::LayerMaterialPos layerPos = detExtension->layerMaterialPos();
        ACTS_VERBOSE(
            "[L] Layer is marked to carry support material on Surface ( "
            "inner=0 / center=1 / outer=2 ) :   "
            << layerPos
            << "    with binning: ["
            << bins1
            << ", "
            << bins2
            << "]");
        // Create an approachdescriptor for the layer
        // create the new surfaces for the approachdescriptor
        std::vector<const Acts::Surface*> aSurfaces;
        // create the inner boundary surface
        Acts::CylinderSurface* innerBoundary
            = new Acts::CylinderSurface(actsTransform, rMin, halfZ);
        // create outer boundary surface
        Acts::CylinderSurface* outerBoundary
            = new Acts::CylinderSurface(actsTransform, rMax, halfZ);
        // check if the material should be set to the inner or outer boundary
        // and set it in case
        if (layerPos == Acts::LayerMaterialPos::inner)
          innerBoundary->setAssociatedMaterial(materialProxy);
        if (layerPos == Acts::LayerMaterialPos::outer)
          outerBoundary->setAssociatedMaterial(materialProxy);
        // collect the surfaces
        aSurfaces.push_back(innerBoundary);
        aSurfaces.push_back(outerBoundary);
        // create an ApproachDescriptor with standard surfaces surfaces - these
        // will
        // be deleted by the approach descriptor
        approachDescriptor
            = new Acts::GenericApproachDescriptor<const Acts::Surface>(
                aSurfaces);
      }
      if (modules.empty()) {
        auto cylLayer = Acts::CylinderLayer::create(convertTransform(transform),
                                                    cylinderBounds,
                                                    nullptr,
                                                    thickness,
                                                    nullptr,
                                                    approachDescriptor,
                                                    Acts::passive);
        // assign layer to approach surfaces if approach descriptor is given
        if (approachDescriptor) {
          for (auto& sIter : (approachDescriptor->containedSurfaces())) {
            if (sIter) sIter->associateLayer(*cylLayer.get());
          }
        }
        // hand over the possible material if it should be in the center
        if (layerPos == Acts::LayerMaterialPos::central)
          cylLayer->surfaceRepresentation().setAssociatedMaterial(
              materialProxy);
        layers.push_back(cylLayer);

      } else {
        ACTS_VERBOSE(
            "[L] Layer containes modules -> resolving them as surfaces");
        // create surfaces binned in phi and z
        auto surfaceArray
            = createSurfaceArray(modules, binZ, motherTransform, axes);
        auto cylLayer = Acts::CylinderLayer::create(convertTransform(transform),
                                                    cylinderBounds,
                                                    std::move(surfaceArray),
                                                    thickness,
                                                    nullptr,
                                                    approachDescriptor,
                                                    Acts::active);
        // assign layer to approach surfaces if approach descriptor is given
        if (approachDescriptor) {
          for (auto& sIter : (approachDescriptor->containedSurfaces())) {
            if (sIter) sIter->associateLayer(*cylLayer.get());
          }
        }
        // hand over the possible material if it should be in the center
        if (layerPos == Acts::LayerMaterialPos::central)
          cylLayer->surfaceRepresentation().setAssociatedMaterial(
              materialProxy);
        layers.push_back(cylLayer);
      }
    }  // for children
  }    // volume has layers
  ACTS_VERBOSE("[V] Volume has no layers");
}

void
Acts::DD4hepCylinderGeometryBuilder::createDiscLayers(
    DD4hep::Geometry::DetElement& motherDetElement,
    Acts::LayerVector&            layers,
    const TGeoMatrix*             motherTransform) const
{
  // get possible layers
  const DD4hep::Geometry::DetElement::Children& children
      = motherDetElement.children();
  // check if volume has layers
  if (!children.empty()) {
    ACTS_VERBOSE("[V] Volume containes layers -> creating disc layers");
    for (auto& child : children) {
      // get the detector element of the layer
      DD4hep::Geometry::DetElement detElement = child.second;
      // build the matrix
      TGeoMatrix* locTrans = detElement.placement().ptr()->GetMatrix();
      // make the transformation global
      TGeoHMatrix* transform
          = new TGeoHMatrix((*motherTransform) * (*locTrans));
      std::shared_ptr<Acts::Transform3D> actsTransform
          = convertTransform(transform);
      // get the shape of the layer
      TGeoShape* geoShape
          = detElement.placement().ptr()->GetVolume()->GetShape();
      TGeoConeSeg* disc = dynamic_cast<TGeoConeSeg*>(geoShape);
      if (!disc)
        throw "Cylinder layer has wrong shape - needs to be TGeoConeSeg!";
      // extract the boundaries
      double rMin      = disc->GetRmin1() * _cm;
      double rMax      = disc->GetRmax1() * _cm;
      double thickness = 2. * disc->GetDz() * _cm;
      auto discBounds  = std::make_shared<const Acts::RadialBounds>(rMin, rMax);

      ACTS_DEBUG(
          "[L] Creating disc layer with dimensions ( rmin / rmax / halfZ ) = ( "
          << rMin
          << " / "
          << rMax
          << " / "
          << thickness
          << " )");
      // if necessary receive the modules contained by the layer and create the
      // layer, otherwise create empty layer
      Acts::IDetExtension* detExtension
          = detElement.extension<Acts::IDetExtension>();
      std::vector<DD4hep::Geometry::DetElement> modules(
          detExtension->modules());
      // access the axis orienation of the modules
      const std::string axes = detExtension->axes();
      // create the two dimensional BinUtility for the material map of the layer
      Acts::BinUtility*                           materialBinUtil = nullptr;
      std::shared_ptr<const SurfaceMaterialProxy> materialProxy   = nullptr;
      // the approachdescriptor telling where the material sits on the layer
      Acts::ApproachDescriptor* approachDescriptor = nullptr;
      // access the material position
      Acts::LayerMaterialPos layerPos = detExtension->layerMaterialPos();
      // check if layer should have material
      if (detExtension->hasSupportMaterial()) {
        std::pair<size_t, size_t> materialBins = detExtension->materialBins();
        size_t bins1    = materialBins.first;
        size_t bins2    = materialBins.second;
        materialBinUtil = new Acts::BinUtility(
            bins1, -M_PI, M_PI, Acts::closed, Acts::binPhi, actsTransform);
        (*materialBinUtil)
            += Acts::BinUtility(bins2, rMin, rMax, Acts::open, Acts::binR);
        // and create material proxy to mark layer for material mapping
        materialProxy
            = std::make_shared<const SurfaceMaterialProxy>(*materialBinUtil);
        // access the material position
        Acts::LayerMaterialPos layerPos = detExtension->layerMaterialPos();
        ACTS_VERBOSE(
            "[L] Layer is marked to carry support material on Surface ( "
            "inner=0 / center=1 / outer=2 ) :   "
            << layerPos
            << "    with binning: ["
            << bins1
            << ", "
            << bins2
            << "]");
        // Create an approachdescriptor for the layer
        // create the new surfaces - positions first
        double         thickness = fabs(rMax - rMin);
        const Vector3D center    = actsTransform->translation();
        const Vector3D normal    = actsTransform->rotation().col(2);
        Vector3D       aspPosition(center + 0.5 * thickness * normal);
        Vector3D       asnPosition(center - 0.5 * thickness * normal);
        auto           asnTransform
            = std::make_shared<Transform3D>(Translation3D(asnPosition));
        auto aspTransform
            = std::make_shared<Transform3D>(Translation3D(aspPosition));
        // create the vector
        std::vector<const Surface*> aSurfaces;
        // create the inner boundary surface
        std::shared_ptr<const DiscBounds> bounds
            = std::make_shared<const RadialBounds>(rMin, rMax);
        Acts::DiscSurface* innerBoundary
            = new DiscSurface(asnTransform, bounds);
        // create outer boundary surface
        Acts::DiscSurface* outerBoundary
            = new DiscSurface(aspTransform, bounds);
        // check if the material should be set to the inner or outer boundary
        // and set it in case
        if (layerPos == Acts::LayerMaterialPos::inner)
          innerBoundary->setAssociatedMaterial(materialProxy);
        if (layerPos == Acts::LayerMaterialPos::outer)
          outerBoundary->setAssociatedMaterial(materialProxy);
        // collect the surfaces
        aSurfaces.push_back(innerBoundary);
        aSurfaces.push_back(outerBoundary);
        // create an ApproachDescriptor with standard surfaces surfaces - these
        // will
        // be deleted by the approach descriptor
        approachDescriptor
            = new Acts::GenericApproachDescriptor<const Acts::Surface>(
                aSurfaces);
      }

      if (modules.empty()) {
        auto discLayer = Acts::DiscLayer::create(actsTransform,
                                                 discBounds,
                                                 nullptr,
                                                 thickness,
                                                 nullptr,
                                                 approachDescriptor,
                                                 Acts::passive);
        // assign layer to approach surfaces if approach descriptor is given
        if (approachDescriptor) {
          for (auto& sIter : (approachDescriptor->containedSurfaces())) {
            if (sIter) sIter->associateLayer(*discLayer.get());
          }
        }
        // hand over the possible material if it should be in the center
        if (layerPos == Acts::LayerMaterialPos::central)
          discLayer->surfaceRepresentation().setAssociatedMaterial(
              materialProxy);
        layers.push_back(discLayer);
      } else {
        ACTS_VERBOSE(
            "[L] Layer containes modules -> resolving them as surfaces");
        // create surfaces binned in phi and r
        auto surfaceArray
            = createSurfaceArray(modules, binR, motherTransform, axes);
        auto discLayer = Acts::DiscLayer::create(actsTransform,
                                                 discBounds,
                                                 std::move(surfaceArray),
                                                 thickness,
                                                 nullptr,
                                                 approachDescriptor,
                                                 Acts::active);
        // assign layer to approach surfaces if approach descriptor is given
        if (approachDescriptor) {
          for (auto& sIter : (approachDescriptor->containedSurfaces())) {
            if (sIter) sIter->associateLayer(*discLayer.get());
          }
        }
        // hand over the possible material if it should be in the center
        if (layerPos == Acts::LayerMaterialPos::central)
          discLayer->surfaceRepresentation().setAssociatedMaterial(
              materialProxy);
        layers.push_back(discLayer);
      }
    }  // for children
  }    // volume has layers
  ACTS_VERBOSE("[V] Volume has no layers");
}
std::unique_ptr<Acts::SurfaceArray>
Acts::DD4hepCylinderGeometryBuilder::createSurfaceArray(
    std::vector<DD4hep::Geometry::DetElement>& modules,
    Acts::BinningValue                         lValue,
    const TGeoMatrix*                          motherTransform,
    const std::string&                         axes) const
{
  ACTS_VERBOSE("[L] Creating surface array of the layer");
  std::vector<const Acts::Surface*> surfaces;
  for (auto& detElement : modules) {
    // make here the material mapping
    DD4hep::Geometry::Segmentation segmentation;
    // extract segmentation //change later
    if (detElement.volume().isSensitive()) {
      Acts::IDetExtension* detExtension
          = detElement.extension<Acts::IDetExtension>();
      segmentation = detExtension->segmentation();
      if (!segmentation)
        ACTS_ERROR("[S] Detector element is sensitive but Segmentation was not "
                   "handed over in geometry constructor, can not access "
                   "segmentation");
    } else
      ACTS_ERROR("[S] Detector element is not declared sensitive, can not "
                 "access segmentation");
    Acts::DD4hepDetElement* dd4hepDetElement = new Acts::DD4hepDetElement(
        detElement, segmentation, motherTransform, axes, _cm);
    // add surface to surface vector
    surfaces.push_back(&(dd4hepDetElement->surface()));
  }

  return (std::move(binnedSurfaceArray2DPhiL(surfaces, lValue)));
}

// creating a surface array binned in phi and a longitudinal direction which can
// either be z or r
std::unique_ptr<Acts::SurfaceArray>
Acts::DD4hepCylinderGeometryBuilder::binnedSurfaceArray2DPhiL(
    const std::vector<const Acts::Surface*> surfaces,
    Acts::BinningValue                      lValue) const
{
  if (surfaces.empty()) ACTS_ERROR("[L] Active layer has no surfaces");
  // boundaries in r, first value minimum radius and second value maximum radius
  // of the current cylinder layer
  std::vector<std::pair<float, float>> lBoundaries;
  // key values in phi
  std::vector<const Acts::Surface*> keys_dedupP;
  // key values in r
  std::vector<const Acts::Surface*> keys_dedupR;
  // for creating the binned array a vector of surfaces plus their corresponding
  // position is needed
  std::vector<std::pair<const Acts::Surface*, Acts::Vector3D>> posSurfaces;
  for (auto& surface : surfaces) {
    // fill the position and surface vector
    posSurfaces.push_back(std::make_pair(surface, surface->center()));
    // receive the bounds
    const PlanarBounds* planarBounds
        = dynamic_cast<const PlanarBounds*>(&(surface->bounds()));
    if (!planarBounds)
      ACTS_ERROR("[S] Given SurfaceBounds are not planar - not implemented for "
                 "other bounds yet! ");
    // get the boundaries of the longitudinal coordinate
    std::vector<Acts::Vector2D> vertices = planarBounds->vertices();
    if (vertices.empty()) ACTS_ERROR("[S] Vertices of current module empty!");
    // accessing any entry to access the longitudinal coordinate is sufficient
    // make local coordinate global
    // double longCoord = vertices.front().perp();
    double longCoord = fabs(vertices.front().y());
    double lposition = 0;
    // surfaces can be binned in z or r
    if (lValue == Acts::binZ)
      lposition = surface->center().z();
    else
      lposition = surface->center().perp();
    lBoundaries.push_back(std::make_pair<float, float>(lposition - longCoord,
                                                       lposition + longCoord));
  }
  // get the number of keys in phi
  std::unique_copy(begin(surfaces),
                   end(surfaces),
                   back_inserter(keys_dedupP),
                   [](const Acts::Surface* a, const Acts::Surface* b) {
                     return (a->center().phi() == b->center().phi());
                   });
  size_t binsPhi = keys_dedupP.size();
  double step    = 2. * M_PI
      / binsPhi;  // assume that it is 2PI otherwise take from bounds
  double minPhi          = -M_PI;
  double maxPhi          = M_PI;
  double minPhiCorrected = minPhi - 0.5 * step;
  double maxPhiCorrected = maxPhi + 0.5 * step;
  if (minPhiCorrected < -M_PI) {
    minPhiCorrected += step;
    maxPhiCorrected += step;
  }
  // eliminate double values and sort values
  std::vector<float> lValues(createBinValues(lBoundaries));
  // create the 2D bin utility
  auto binUtility
      = std::make_unique<Acts::BinUtility>(lValues, Acts::open, lValue);
  (*binUtility) += Acts::BinUtility(
      binsPhi, minPhiCorrected, maxPhiCorrected, Acts::closed, Acts::binPhi);
  // create the binned array of surfaces
  return (std::move(std::make_unique<Acts::BinnedArrayXD<const Acts::Surface*>>(
      posSurfaces, std::move(binUtility))));
}

std::vector<float>
Acts::DD4hepCylinderGeometryBuilder::createBinValues(
    std::vector<std::pair<float, float>> old) const
{
  sort(old.begin(), old.end(), sortFloatPairs);
  std::vector<float> newlValues;
  std::pair<float, float> current;
  std::pair<float, float> next;
  // eliminate doubles
  std::vector<std::pair<float, float>> oldKeys;
  std::unique_copy(begin(old),
                   end(old),
                   back_inserter(oldKeys),
                   [](std::pair<float, float>& a, std::pair<float, float>& b) {
                     return a == b;
                   });
  for (std::vector<std::pair<float, float>>::iterator it = oldKeys.begin();
       it != (oldKeys.end() - 1);
       ++it) {
    current = *it;
    next    = *(it + 1);
    if (it == oldKeys.begin()) newlValues.push_back(current.first);
    if (next.first <= current.second)
      newlValues.push_back((current.second + next.first) * 0.5);
    else
      newlValues.push_back(current.second);
    if (it == (oldKeys.end() - 2)) newlValues.push_back(next.second);
  }
  return (newlValues);
}

bool
Acts::DD4hepCylinderGeometryBuilder::sortFloatPairs(std::pair<float, float> ap,
                                                    std::pair<float, float> bp)
{
  float a = (ap.second + ap.first) * 0.5;
  float b = (bp.second + bp.first) * 0.5;
  return a < b;
}
