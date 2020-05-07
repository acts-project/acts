// This file is part of the Acts project.
//
// Copyright (C) 2017-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/GeometryID.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"

using json = nlohmann::json;

template <typename object_t>
Acts::JsonHierarchicalObjectConverter<object_t>::
    JsonHierarchicalObjectConverter(
        const Acts::JsonHierarchicalObjectConverter<object_t>::Config& cfg)
    : m_cfg(std::move(cfg)) {
  // Validate the configuration
  if (!m_cfg.logger) {
    throw std::invalid_argument("Missing logger");
  }
}

template <typename object_t>
Acts::HierarchicalGeometryContainer<object_t>
Acts::JsonHierarchicalObjectConverter<object_t>::jsonToHierarchicalContainer(
    const nlohmann::json& map,
    std::function<object_t(const nlohmann::json&)> fromJson) const {
  auto& j = map;
  std::vector<object_t> elements;
  // The return maps
  ACTS_VERBOSE("j2a: Reading maps from json file.");
  ACTS_VERBOSE("j2a: Found entries for " << j.count(m_cfg.volkey)
                                         << " volume(s).");
  // Structured binding
  for (auto& [key, value] : j.items()) {
    // Check if this the volume key
    if (key == m_cfg.volkey) {
      // Get the volume json
      auto volj = value;
      for (auto& [vkey, vvalue] : volj.items()) {
        ACTS_VERBOSE("j2a: -> Found Volume " << vkey);
        // Loop through the information in the volume
        for (auto& [vckey, vcvalue] : vvalue.items()) {
          if (m_cfg.processBoundaries and vckey == m_cfg.boukey and
              not vcvalue.empty()) {
            ACTS_VERBOSE("j2a: --> BoundarySurface(s) to be parsed");
            for (auto& [bkey, bvalue] : vcvalue.items()) {
              ACTS_VERBOSE("j2a: ---> Found boundary surface " << bkey);
              elements.push_back(fromJson(bvalue));
            }
          } else if (vckey == m_cfg.laykey) {
            ACTS_VERBOSE("j2a: --> Layer(s) to be parsed");
            // Loop over layers and repeat
            auto layj = vcvalue;
            for (auto& [lkey, lvalue] : layj.items()) {
              ACTS_VERBOSE("j2a: ---> Found Layer " << lkey);
              // Finally loop over layer components
              for (auto& [lckey, lcvalue] : lvalue.items()) {
                if (m_cfg.processRepresenting and lckey == m_cfg.repkey and
                    not lcvalue.empty()) {
                  ACTS_VERBOSE("j2a: ----> Found representing surface");
                  elements.push_back(fromJson(lcvalue));
                } else if (m_cfg.processApproaches and lckey == m_cfg.appkey and
                           not lcvalue.empty()) {
                  ACTS_VERBOSE("j2a: ----> Found approach surface(s)");
                  // Loop over approach surfaces
                  for (auto& [askey, asvalue] : lcvalue.items()) {
                    ACTS_VERBOSE("j2a: -----> Approach surface " << askey);
                    elements.push_back(fromJson(asvalue));
                  }
                } else if (m_cfg.processSensitives and lckey == m_cfg.senkey and
                           not lcvalue.empty()) {
                  ACTS_VERBOSE("j2a: ----> Found sensitive surface(s)");
                  // Loop over sensitive surfaces
                  for (auto& [sskey, ssvalue] : lcvalue.items()) {
                    ACTS_VERBOSE("j2a: -----> Sensitive surface " << sskey);
                    elements.push_back(fromJson(ssvalue));
                  }
                }
              }
            }
          } else if (m_cfg.processVolumes and vckey == m_cfg.datakey and
                     not vcvalue.empty()) {
            ACTS_VERBOSE("--> Volume Object to be parsed");
            elements.push_back(fromJson(vcvalue));
          }
        }
      }
    }
  }
  HierarchicalGeometryContainer<object_t> container(std::move(elements));
  // Return the filled container
  return container;
}

template <typename object_t>
nlohmann::json
Acts::JsonHierarchicalObjectConverter<object_t>::hierarchicalObjectToJson(
    const HierarchicalGeometryContainer<object_t>& hObject,
    std::function<nlohmann::json(const object_t&)> toJson) const {
  DetectorRep detRep;
  convertToRep(detRep, hObject);
  return detectorRepToJson(detRep, toJson);
}

template <typename object_t>
nlohmann::json
Acts::JsonHierarchicalObjectConverter<object_t>::trackingGeometryToJson(
    const TrackingGeometry& tGeometry,
    std::function<nlohmann::json(const object_t&)> toJson,
    std::function<object_t(const GeometryID&)> initialise) const {
  DetectorRep detRep;
  convertToRep(detRep, *tGeometry.highestTrackingVolume(), initialise);
  return detectorRepToJson(detRep, toJson);
}

template <typename object_t>
void Acts::JsonHierarchicalObjectConverter<object_t>::convertToRep(
    DetectorRep& detRep, const Acts::TrackingVolume& tVolume,
    std::function<object_t(const GeometryID&)> initialise) const {
  // Create a volume representation
  VolumeRep volRep;
  volRep.volumeName = tVolume.volumeName();
  // there are confined volumes
  if (tVolume.confinedVolumes() != nullptr) {
    // get through the volumes
    auto& volumes = tVolume.confinedVolumes()->arrayObjects();
    // loop over the volumes
    for (auto& vol : volumes) {
      // recursive call
      convertToRep(detRep, *vol, initialise);
    }
  }
  // there are dense volumes
  if (!tVolume.denseVolumes().empty()) {
    // loop over the volumes
    for (auto& vol : tVolume.denseVolumes()) {
      // recursive call
      convertToRep(detRep, *vol, initialise);
    }
  }
  // Get the volume Id
  Acts::GeometryID volumeID = tVolume.geoID();
  Acts::GeometryID::Value vid = volumeID.volume();
  volRep.volume = initialise(volumeID);

  // there are confied layers
  if (tVolume.confinedLayers() != nullptr) {
    // get the layers
    auto& layers = tVolume.confinedLayers()->arrayObjects();
    // loop of the volumes
    for (auto& lay : layers) {
      auto layRep = convertToRep(*lay, initialise);
      // it's a valid representation so let's go with it
      Acts::GeometryID layerID = lay->geoID();
      Acts::GeometryID::Value lid = layerID.layer();
      volRep.layers.insert({lid, std::move(layRep)});
    }
  }
  // Let's finally check the boundaries
  for (auto& bsurf : tVolume.boundarySurfaces()) {
    // the surface representation
    auto& bssfRep = bsurf->surfaceRepresentation();
    Acts::GeometryID boundaryID = bssfRep.geoID();
    Acts::GeometryID::Value bid = boundaryID.boundary();
    volRep.boundaries[bid] = initialise(boundaryID);
  }
  volRep.volumeName = tVolume.volumeName();
  volRep.volumeID = volumeID;
  detRep.volumes.insert({vid, std::move(volRep)});
  return;
}

template <typename object_t>
typename Acts::JsonHierarchicalObjectConverter<object_t>::LayerRep
Acts::JsonHierarchicalObjectConverter<object_t>::convertToRep(
    const Acts::Layer& tLayer,
    std::function<object_t(const GeometryID&)> initialise) const {
  // Layer representation
  LayerRep layRep;
  // fill layer ID information
  layRep.layerID = tLayer.geoID();
  if (tLayer.surfaceArray() != nullptr) {
    for (auto& ssf : tLayer.surfaceArray()->surfaces()) {
      if (ssf != nullptr) {
        Acts::GeometryID sensitiveID = ssf->geoID();
        Acts::GeometryID::Value sid = sensitiveID.sensitive();
        layRep.sensitives.insert({sid, initialise(sensitiveID)});
      }
    }
  }
  // the representing
  if (!(tLayer.surfaceRepresentation().geoID() == GeometryID())) {
    layRep.representing = initialise(tLayer.surfaceRepresentation().geoID());
  }
  // the approach
  if (tLayer.approachDescriptor() != nullptr) {
    for (auto& asf : tLayer.approachDescriptor()->containedSurfaces()) {
      // get the surface
      Acts::GeometryID approachID = asf->geoID();
      Acts::GeometryID::Value aid = approachID.approach();
      layRep.approaches.insert({aid, initialise(approachID)});
    }
  }
  // return the layer representation
  return layRep;
}

template <typename object_t>
void Acts::JsonHierarchicalObjectConverter<object_t>::convertToRep(
    DetectorRep& detRep,
    const Acts::HierarchicalGeometryContainer<object_t>& hObject) const {
  size_t index = 0;
  Acts::GeometryID::Value volume = 0;
  Acts::GeometryID::Value layer = 0;
  // Loop through the sorted geometry elements
  for (auto it = hObject.begin(); it != hObject.end(); ++it) {
    index = it - hObject.begin();
    Acts::GeometryID geoId = hObject.idAt(index);
    // The Volumes
    if (geoId.volume() != volume && geoId.boundary() == 0 &&
        geoId.layer() == 0 && geoId.approach() == 0 && geoId.sensitive() == 0) {
      volume = geoId.volume();
      layer = 0;
      VolumeRep volRep;
      volRep.volumeID = geoId;
      volRep.volume = *it;
      detRep.volumes.insert({volume, std::move(volRep)});
    }
    // The Layers
    if (geoId.layer() != layer && geoId.boundary() == 0 &&
        geoId.approach() == 0 && geoId.sensitive() == 0) {
      layer = geoId.layer();
      LayerRep layRep;
      layRep.layerID = geoId;
      layRep.representing = *it;
      detRep.volumes[volume].layers.insert({layer, std::move(layRep)});
    }
    // The approach
    if (geoId.approach() != 0) {
      detRep.volumes[volume].layers[layer].approaches.insert(
          {geoId.approach(), *it});
    }
    // The sensitive element
    if (geoId.sensitive() != 0) {
      detRep.volumes[volume].layers[layer].sensitives.insert(
          {geoId.sensitive(), *it});
    }
    // The approach boundary
    if (geoId.boundary() != 0) {
      detRep.volumes[volume].boundaries.insert({geoId.boundary(), *it});
    }
  }
  return;
}

/// Create Json from a detector represenation
template <typename object_t>
json Acts::JsonHierarchicalObjectConverter<object_t>::detectorRepToJson(
    const DetectorRep& detRep,
    std::function<nlohmann::json(const object_t&)> toJson) const {
  json detectorj;
  ACTS_VERBOSE("a2j: Writing json from detector representation");
  ACTS_VERBOSE("a2j: Found entries for " << detRep.volumes.size()
                                         << " volume(s).");

  json volumesj;
  for (auto& [key, value] : detRep.volumes) {
    json volj;
    ACTS_VERBOSE("a2j: -> Writing Volume: " << key);
    volj[m_cfg.namekey] = value.volumeName;
    if (m_cfg.processVolumes)
      volj[m_cfg.datakey] = toJson(value.volume);
    // Write the boundary surfaces
    if (m_cfg.processBoundaries && not value.boundaries.empty()) {
      ACTS_VERBOSE("a2j: ---> Found " << value.boundaries.size()
                                      << " boundary/ies ");
      json boundariesj;
      for (auto& [bkey, bvalue] : value.boundaries) {
        ACTS_VERBOSE("a2j: ----> Convert boundary " << bkey);
        boundariesj[std::to_string(bkey)] = toJson(bvalue);
      }
      volj[m_cfg.boukey] = boundariesj;
    }
    // Write the layers
    if (not value.layers.empty()) {
      ACTS_VERBOSE("a2j: ---> Found " << value.layers.size() << " layer(s) ");
      json layersj;
      for (auto& [lkey, lvalue] : value.layers) {
        ACTS_VERBOSE("a2j: ----> Convert layer " << lkey);
        json layj;
        // Check for representing
        if (m_cfg.processRepresenting) {
          ACTS_VERBOSE("a2j: ------> Convert representing surface ");
          layj[m_cfg.repkey] = toJson(lvalue.representing);
        }
        // First check for approaches
        if (not lvalue.approaches.empty() and m_cfg.processApproaches) {
          ACTS_VERBOSE("a2j: -----> Found " << lvalue.approaches.size()
                                            << " approach surface(s)");
          json approachesj;
          for (auto& [akey, avalue] : lvalue.approaches) {
            ACTS_VERBOSE("a2j: ------> Convert approach surface " << akey);
            approachesj[std::to_string(akey)] = toJson(avalue);
          }
          // Add to the layer json
          layj[m_cfg.appkey] = approachesj;
        }
        // Then check for sensitive
        if (not lvalue.sensitives.empty() and m_cfg.processSensitives) {
          ACTS_VERBOSE("a2j: -----> Found " << lvalue.sensitives.size()
                                            << " sensitive surface(s)");
          json sensitivesj;
          for (auto& [skey, svalue] : lvalue.sensitives) {
            ACTS_VERBOSE("a2j: ------> Convert sensitive surface " << skey);
            sensitivesj[std::to_string(skey)] = toJson(svalue);
          }
          // Add to the layer json
          layj[m_cfg.senkey] = sensitivesj;
        }
        layersj[std::to_string(lkey)] = layj;
        volj[m_cfg.laykey] = layersj;
      }
    }
    volumesj[std::to_string(key)] = volj;
  }
  // Assign the volume json to the detector json
  detectorj[m_cfg.volkey] = volumesj;
  return detectorj;
}
