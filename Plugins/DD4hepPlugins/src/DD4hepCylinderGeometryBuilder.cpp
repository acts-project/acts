// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTS/Plugins/DD4hepPlugins/DD4hepCylinderGeometryBuilder.hpp"
// Geometry module
#include "ACTS/Detector/TrackingGeometry.hpp"
#include "ACTS/Detector/TrackingVolume.hpp"
#include "ACTS/Material/Material.hpp"
// DD4hepPlugin
#include "ACTS/Plugins/DD4hepPlugins/DD4hepGeometryHelper.hpp"

Acts::DD4hepCylinderGeometryBuilder::DD4hepCylinderGeometryBuilder(
    const Config dgbConfig)
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
  m_config = std::move(dgbConfig);
}

Acts::DD4hepCylinderGeometryBuilder::~DD4hepCylinderGeometryBuilder()
{
}

std::unique_ptr<Acts::TrackingGeometry>
Acts::DD4hepCylinderGeometryBuilder::trackingGeometry() const
{
  // the return geometry -- and the highest volume
  std::unique_ptr<Acts::TrackingGeometry> trackingGeometry = nullptr;
  Acts::TrackingVolumePtr                 highestVolume    = nullptr;
  Acts::TrackingVolumePtr                 beamPipeVolume   = nullptr;

  // get the sub detectors
  std::vector<DD4hep::Geometry::DetElement>     detElements;
  const DD4hep::Geometry::DetElement::Children& children
      = m_config.detWorld.children();
  for (auto& detElement : children) detElements.push_back(detElement.second);
  // sort by id to build detector from bottom to top
  sort(detElements.begin(),
       detElements.end(),
       [](const DD4hep::Geometry::DetElement& a,
          const DD4hep::Geometry::DetElement& b) { return (a.id() < b.id()); });
  // loop over the volumes
  for (auto& detElement : detElements) {
    if (detElement.type() == "beamtube") {
      // MSG_DEBUG("BeamPipe is being built");
      // extract material
      DD4hep::Geometry::Material mat = detElement.volume().material();
      // create the tracking volume
      beamPipeVolume = Acts::TrackingVolume::create(
          Acts::DD4hepGeometryHelper::extractTransform(detElement),
          Acts::DD4hepGeometryHelper::extractVolumeBounds(detElement),
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
      Acts::DD4hepGeometryHelper::createSubVolumes(
          detElement, layerTriple, volumeTriple);
      highestVolume = m_config.volumeBuilder->trackingVolume(
          highestVolume,
          Acts::DD4hepGeometryHelper::extractVolumeBounds(detElement),
          new LayerTriple(layerTriple),
          new VolumeTriple(volumeTriple));
    }
  }
  // if you have a highest volume, stuff it into a TrackingGeometry
  if (highestVolume) {
    // see if the beampipe needs to be wrapped
    if (beamPipeVolume)
      highestVolume = m_config.volumeHelper->createContainerTrackingVolume(
          {beamPipeVolume, highestVolume});
    // create the TrackingGeometry
    trackingGeometry = std::make_unique<Acts::TrackingGeometry>(highestVolume);
  }
  // return the geometry to the service
  return trackingGeometry;
}
