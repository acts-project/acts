// This file is part of the Acts project.
//
// Copyright (C) 2017-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

using json = nlohmann::json;

template <typename object_t>
Acts::HierarchicalGeometryContainer<object_t>
Acts::JsonHierarchicalObjectConverter<object_t>::jsonToHierarchicalContainer(
    const nlohmann::json& map,
    std::function<object_t(const Acts::GeometryID&, const nlohmann::json&)>
        fromJson) const {
  ACTS_LOCAL_LOGGER(
      Acts::getDefaultLogger("jsonToHierarchicalContainer", logLevel));
  auto& j = map;
  std::vector<object_t> elements;
  // The return maps
  ACTS_VERBOSE("j2a: Reading maps from json file.");
  ACTS_VERBOSE("j2a: Found entries for " << j.count(m_volkey) << " volume(s).");
  // Structured binding
  for (auto& [key, value] : j.items()) {
    // Check if this the volume key
    if (key == m_volkey) {
      // Get the volume json
      auto volj = value;
      for (auto& [vkey, vvalue] : volj.items()) {
        ACTS_VERBOSE("j2a: -> Found Volume " << vkey);
        // Loop through the information in the volume
        for (auto& [vckey, vcvalue] : vvalue.items()) {
          if (vckey == m_boukey and not vcvalue.empty()) {
            ACTS_VERBOSE("j2a: --> BoundarySurface(s) to be parsed");
            for (auto& [bkey, bvalue] : vcvalue.items()) {
              ACTS_VERBOSE("j2a: ---> Found boundary surface " << bkey);
              for (auto& [okey, ovalue] : bvalue.items()) {
                if (okey == datakey) {
                  Acts::GeometryID id = makeId(vkey, bkey, "0", "0", "0");
                  elements.push_back(fromJson(id, ovalue));
                }
              }
            }
          } else if (vckey == m_laykey) {
            ACTS_VERBOSE("j2a: --> Layer(s) to be parsed");
            // Loop over layers and repeat
            auto layj = vcvalue;
            for (auto& [lkey, lvalue] : layj.items()) {
              ACTS_VERBOSE("j2a: ---> Found Layer " << lkey);
              // Finally loop over layer components
              for (auto& [lckey, lcvalue] : lvalue.items()) {
                if (lckey == m_repkey and not lcvalue.empty()) {
                  ACTS_VERBOSE("j2a: ----> Found representing surface");
                  for (auto& [okey, ovalue] : lcvalue.items()) {
                    if (okey == datakey) {
                      Acts::GeometryID id = makeId(vkey, "0", lkey, "0", "0");
                      elements.push_back(fromJson(id, ovalue));
                    }
                  }
                } else if (lckey == m_appkey and not lcvalue.empty()) {
                  ACTS_VERBOSE("j2a: ----> Found approach surface(s)");
                  // Loop over approach surfaces
                  for (auto& [askey, asvalue] : lcvalue.items()) {
                    ACTS_VERBOSE("j2a: -----> Approach surface " << askey);
                    for (auto& [okey, ovalue] : asvalue.items()) {
                      if (okey == datakey) {
                        Acts::GeometryID id =
                            makeId(vkey, "0", lkey, askey, "0");
                        elements.push_back(fromJson(id, ovalue));
                      }
                    }
                  }
                } else if (lckey == m_senkey and not lcvalue.empty()) {
                  ACTS_VERBOSE("j2a: ----> Found sensitive surface(s)");
                  // Loop over sensitive surfaces
                  for (auto& [sskey, ssvalue] : lcvalue.items()) {
                    ACTS_VERBOSE("j2a: -----> Sensitive surface " << sskey);
                    for (auto& [okey, ovalue] : ssvalue.items()) {
                      if (okey == datakey) {
                        Acts::GeometryID id =
                            makeId(vkey, "0", lkey, "0", sskey);
                        elements.push_back(fromJson(id, ovalue));
                      }
                    }
                  }
                }
              }
            }
          } else if (vckey == datakey and not vcvalue.empty()) {
            ACTS_VERBOSE("j2a: --> Volume Object to be parsed");
            Acts::GeometryID id = makeId(vkey, "0", "0", "0", "0");
            elements.push_back(fromJson(id, vcvalue));
          }
        }
      }
    }
  }
  if (elements.empty()) {
    ACTS_WARNING("No element found in Json file with key " << datakey);
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
  ACTS_LOCAL_LOGGER(
      Acts::getDefaultLogger("hierarchicalObjectToJson", logLevel));
  json detectorj;
  json volumesj;
  json boundaryj;
  json layersj;
  json approachesj;
  json sensitivesj;
  size_t index = 0;
  Acts::GeometryID::Value volume = 0;
  Acts::GeometryID::Value layer = 0;
  // Loop through the container elements backward
  for (auto it = hObject.end(); it-- != hObject.begin();) {
    index = it - hObject.begin();
    Acts::GeometryID geoId = hObject.idAt(index);
    // The Volumes
    if (geoId.volume() != volume && geoId.boundary() == 0 &&
        geoId.layer() == 0 && geoId.approach() == 0 && geoId.sensitive() == 0) {
      json volj;
      volume = geoId.volume();
      ACTS_VERBOSE("c2j: --> Found Volume " << volume);
      layer = 0;
      volj[datakey] = toJson(*it);
      if (boundaryj != json()) {
        ACTS_VERBOSE("c2j: --> Store boundaries in volume");
        volj[m_boukey] = boundaryj;
        boundaryj = json();
      }
      if (layersj != json()) {
        ACTS_VERBOSE("c2j: --> Store layers in volume");
        volj[m_laykey] = layersj;
        layersj = json();
      }
      ACTS_VERBOSE("c2j: -> Store volume " << volume << " in volume list");
      volumesj[std::to_string(volume)] = volj;
    }
    // The Layers
    if (geoId.layer() != layer && geoId.boundary() == 0 &&
        geoId.approach() == 0 && geoId.sensitive() == 0) {
      json layj;
      layer = geoId.layer();
      ACTS_VERBOSE("c2j: ----> Found Layer / representing surface " << layer);
      layj[m_repkey][datakey] = toJson(*it);
      if (approachesj != json()) {
        ACTS_VERBOSE("c2j: ----> Store approaches surfaces in layer");
        layj[m_appkey] = approachesj;
        approachesj = json();
      }
      if (sensitivesj != json()) {
        ACTS_VERBOSE("c2j: ----> Store sensitives surfaces in layer");
        layj[m_senkey] = sensitivesj;
        sensitivesj = json();
      }
      ACTS_VERBOSE("c2j: ---> Store layer " << layer << " in layer list");
      layersj[std::to_string(layer)] = layj;
    }
    // The approach
    if (geoId.approach() != 0) {
      ACTS_VERBOSE("c2j: ------> Found approach surface " << geoId.approach());
      approachesj[std::to_string(geoId.approach())][datakey] = toJson(*it);
    }
    // The sensitive element
    if (geoId.sensitive() != 0) {
      ACTS_VERBOSE("c2j: ------> Found sensitive surface "
                   << geoId.sensitive());
      sensitivesj[std::to_string(geoId.sensitive())][datakey] = toJson(*it);
    }
    // The approach boundary
    if (geoId.boundary() != 0) {
      ACTS_VERBOSE("c2j: ----> Found boundary surface " << geoId.boundary());
      boundaryj[std::to_string(geoId.boundary())][datakey] = toJson(*it);
    }
  }
  detectorj[m_volkey] = volumesj;
  return detectorj;
}
