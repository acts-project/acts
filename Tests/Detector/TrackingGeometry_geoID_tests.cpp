///  This file is part of the ACTS project.
/// 
///  Copyright (C) 2016 ACTS project team
/// 
///  This Source Code Form is subject to the terms of the Mozilla Public
///  License, v. 2.0. If a copy of the MPL was not distributed with this
///  file, You can obtain one at http:/// mozilla.org/MPL/2.0/.

///  Boost include(s)
#define BOOST_TEST_MODULE GeometryID Tests
#include <boost/test/included/unit_test.hpp>
#include "ACTS/Utilities/Units.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/BinUtility.hpp"
#include "ACTS/Utilities/BinnedArrayXD.hpp"
#include "ACTS/Surfaces/CylinderSurface.hpp"
#include "ACTS/Surfaces/CylinderBounds.hpp"
#include "ACTS/Layers/CylinderLayer.hpp"
#include "ACTS/Volumes/CylinderVolumeBounds.hpp"
#include "ACTS/Detector/TrackingVolume.hpp"
#include "ACTS/Detector/TrackingGeometry.hpp"

namespace Acts {

///  helper function 
TrackingVolumePtr constructCylinderVolume(double surfaceHalfLengthZ,
                                          double surfaceRadius,
                                          double surfaceRstagger,
                                          double surfaceZoverlap,
                                          double layerEnvelope,
                                          double volumeEnvelope,
                                          double innerVolumeR,
                                          double outerVolumeR,
                                          const std::string& name){
  
  ///  the surface transforms
  auto sfnPosition  = Vector3D(0., 0., -3*surfaceHalfLengthZ-surfaceZoverlap);
  auto sfnTransform = std::make_shared<Transform3D>(Translation3D(sfnPosition));
  auto sfcTransform = nullptr;
  auto sfpPosition  = Vector3D(0., 0., 3*surfaceHalfLengthZ-surfaceZoverlap);
  auto sfpTransform = std::make_shared<Transform3D>(Translation3D(sfpPosition));
  ///  the surfaces
  auto sfn = new CylinderSurface(sfnTransform,surfaceRadius-0.5*surfaceRstagger,surfaceHalfLengthZ);
  auto sfc = new CylinderSurface(sfcTransform,surfaceRadius+0.5*surfaceRstagger,surfaceHalfLengthZ);
  auto sfp = new CylinderSurface(sfpTransform,surfaceRadius-0.5*surfaceRstagger,surfaceHalfLengthZ);
  
  ///  prepare the surfaces
  typedef std::pair<const Surface*, Vector3D> TAP;
  std::vector<TAP> surfaces = { {sfn, sfn->binningPosition(binZ)},
                                {sfc, sfc->binningPosition(binZ)},
                                {sfp, sfp->binningPosition(binZ)} };
  
  ///  make the binned array
  double bUmin  = sfnPosition.z()-surfaceHalfLengthZ;
  double bUmax  = sfpPosition.z()+surfaceHalfLengthZ;                              
  auto bUtility = std::make_unique<BinUtility>(surfaces.size(),bUmin,bUmax,open,binZ);
  std::unique_ptr<SurfaceArray> bArray = std::make_unique<BinnedArrayXD<const Surface*> >(surfaces, std::move(bUtility));
  
  ///  now create the Layer
  auto layer0bounds                       = std::make_shared<CylinderBounds>(surfaceRadius, bUmax);
  auto layer0                             = CylinderLayer::create(nullptr,layer0bounds,std::move(bArray),surfaceRstagger+2*layerEnvelope) ;
  std::unique_ptr<LayerArray> layerArray  = std::make_unique< BinnedArrayXD<LayerPtr> >(layer0); 
    
  ///  create the volume
  auto   volumeBounds      = std::make_shared<CylinderVolumeBounds>(innerVolumeR,outerVolumeR,bUmax+volumeEnvelope);
  TrackingVolumePtr volume = TrackingVolume::create(nullptr,
                                                    volumeBounds,
                                                    nullptr,
                                                    std::move(layerArray),
                                                    {}, {}, {},
                                                    name);
  ///  return the volume
  return volume;
  
}


namespace Test {
  
  ///  create three cylinder surfaces
  ///  the surface radius (will also be the layer radius)
  double iv_surfaceHalfLengthZ     = 50  * Acts::units::_mm;
  double iv_surfaceRadius          = 25. * Acts::units::_mm;
  double iv_surfaceRstagger        = 5.  * Acts::units::_mm;
  double iv_surfaceZoverlap        = 10. * Acts::units::_mm; 
  double iv_layerEnvelope          = 0.5 * Acts::units::_mm; 
  double iv_volumeEnvelope         = 10. * Acts::units::_mm; 
  double iv_volumeRadius           = iv_surfaceRadius+0.5*iv_surfaceRstagger+iv_layerEnvelope+iv_volumeEnvelope;
  
  ///  the surface radius (will also be the layer radius)
  double ov_surfaceHalfLengthZ     = 50.  * Acts::units::_mm;
  double ov_surfaceRadius          = 100. * Acts::units::_mm;
  double ov_surfaceRstagger        = 5.   * Acts::units::_mm;
  double ov_surfaceZoverlap        = 10.  * Acts::units::_mm; 
  double ov_layerEnvelope          = 0.5  * Acts::units::_mm;
  double ov_volumeEnvelope         = 10.  * Acts::units::_mm;
  double ov_volumeRadius           = ov_surfaceRadius+0.5*ov_surfaceRstagger+ov_layerEnvelope+ov_volumeEnvelope;
    
  ///  inner volume 
  auto iVolume = constructCylinderVolume(iv_surfaceHalfLengthZ,
                                         iv_surfaceRadius,     
                                         iv_surfaceRstagger,   
                                         iv_surfaceZoverlap,   
                                         iv_layerEnvelope,
                                         iv_volumeEnvelope,
                                         0.,     
                                         iv_volumeRadius,
                                         "InnerVolume");
  
  BOOST_AUTO_TEST_CASE(GeometryID_innervolume_before_closure_test){
    BOOST_CHECK_EQUAL(0, iVolume->geoID().value());
    // check the boundary surfaces
    for (auto bSf : iVolume->boundarySurfaces()){
      BOOST_CHECK_EQUAL(0, bSf->surfaceRepresentation().geoID().value());
      for (auto lay : iVolume->confinedLayers()->arrayObjects()){
        BOOST_CHECK_EQUAL(0, lay->geoID().value());
        // check the approach surfaces  
        for (auto asf : lay->approachDescriptor()->containedSurfaces() )
            BOOST_CHECK_EQUAL(0, asf->geoID().value());
        // check the layer surface array
        for (auto ssf : lay->surfaceArray()->arrayObjects())
             BOOST_CHECK_EQUAL(0, ssf->geoID().value());
      }      
    }                                       
  }
  
  ///  outer volume 
  auto oVolume = constructCylinderVolume(ov_surfaceHalfLengthZ,
                                         ov_surfaceRadius,     
                                         ov_surfaceRstagger,   
                                         ov_surfaceZoverlap,   
                                         ov_layerEnvelope,     
                                         ov_volumeEnvelope,
                                         iv_volumeRadius,
                                         ov_volumeRadius,
                                         "OuterVolume");
                                         
                                         
  BOOST_AUTO_TEST_CASE(GeometryID_outervolume_before_closure_test){
    BOOST_CHECK_EQUAL(0, oVolume->geoID().value());
    // check the boundary surfaces
    for (auto bSf : iVolume->boundarySurfaces()){
      BOOST_CHECK_EQUAL(0, bSf->surfaceRepresentation().geoID().value());
      for (auto lay : oVolume->confinedLayers()->arrayObjects()){
        BOOST_CHECK_EQUAL(0, lay->geoID().value());
        // check the approach surfaces  
        for (auto asf : lay->approachDescriptor()->containedSurfaces() )
            BOOST_CHECK_EQUAL(0, asf->geoID().value());
        // check the layer surface array
        for (auto ssf : lay->surfaceArray()->arrayObjects())
             BOOST_CHECK_EQUAL(0, ssf->geoID().value());
      }      
    }                                       
  }                                      
                                         
  ///  create the volume array
  typedef std::pair<TrackingVolumePtr, Vector3D> VAP;
  std::vector<VAP> volumes = { {iVolume, iVolume->binningPosition(binR)},
                               {oVolume, oVolume->binningPosition(binR)} };
  ///  the bounds for the container
  double hVolumeHalflength = (4*iv_surfaceHalfLengthZ-iv_surfaceZoverlap+iv_volumeEnvelope);
  auto hVolumeBounds = std::make_shared<CylinderVolumeBounds>(0.,ov_volumeRadius,hVolumeHalflength);
  ///  create the BinUtility & the BinnedArray                                         
  auto vUtility =  std::make_unique<BinUtility>(volumes.size(),0.,ov_volumeRadius,open,binR);                                       
  std::shared_ptr<const TrackingVolumeArray> vArray 
    = std::make_shared<const BinnedArrayXD<TrackingVolumePtr> >(volumes,std::move(vUtility));
  ///  create the container volume
  auto hVolume = TrackingVolume::create(nullptr, hVolumeBounds, vArray, "Container");
  
  ///  pre-check on GeometryID
  BOOST_AUTO_TEST_CASE(GeometryID_before_closure_test){
    ///  let's check that the geometry ID values are all 0
    BOOST_CHECK_EQUAL(0, hVolume->geoID().value());
    for (auto cVol : hVolume->confinedVolumes()->arrayObjects()){
      /// let's check everything is set to 0
      BOOST_CHECK_EQUAL(0, cVol->geoID().value());
      // check the boundary surfaces
      for (auto bSf : cVol->boundarySurfaces()){
        BOOST_CHECK_EQUAL(0, bSf->surfaceRepresentation().geoID().value());
      }  
      for (auto lay : cVol->confinedLayers()->arrayObjects()){
        BOOST_CHECK_EQUAL(0, lay->geoID().value());
        // check the approach surfaces  
        for (auto asf : lay->approachDescriptor()->containedSurfaces() )
            BOOST_CHECK_EQUAL(0, asf->geoID().value());
        // check the layer surface array
        for (auto ssf : lay->surfaceArray()->arrayObjects())
             BOOST_CHECK_EQUAL(0, ssf->geoID().value());
      }      
    }
  }     
  
  /// creating a TrackingGeometry closes the geometry
  TrackingGeometry trackingGeometry(hVolume);
  
  ///  after-check on GeometryID
  BOOST_AUTO_TEST_CASE(GeometryID_after_closure_test){
    ///  let's check that the geometry ID values are all 0
    BOOST_CHECK_EQUAL(0, hVolume->geoID().value());

  }  
 



}  //  end of namespace Test
}  //  end of namespace Acts
  
