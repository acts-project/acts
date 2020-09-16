// This file is part of the Acts project.
//
// Copyright (C) 2017-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Json/JsonGeometryConverter.hpp"

#include "Acts/Geometry/ApproachDescriptor.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/CutoutCylinderVolumeBounds.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Material/BinnedSurfaceMaterial.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/HomogeneousVolumeMaterial.hpp"
#include "Acts/Material/InterpolatedMaterialMap.hpp"
#include "Acts/Material/MaterialGridHelper.hpp"
#include "Acts/Material/ProtoSurfaceMaterial.hpp"
#include "Acts/Material/ProtoVolumeMaterial.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include <Acts/Surfaces/AnnulusBounds.hpp>
#include <Acts/Surfaces/CylinderBounds.hpp>
#include <Acts/Surfaces/RadialBounds.hpp>
#include <Acts/Surfaces/SurfaceBounds.hpp>

#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/finder.hpp>
#include <boost/algorithm/string/iter_find.hpp>

namespace {

using json = nlohmann::json;

// helper functions to encode/decode indefinite material
//
// encoded either as `null` for vacuum or to an array of material parameters

json encodeMaterial(const Acts::Material& material) {
  if (!material) {
    return nullptr;
  }
  json encoded = json::array();
  for (unsigned i = 0; i < material.parameters().size(); ++i) {
    encoded.push_back(material.parameters()[i]);
  }
  return encoded;
}

Acts::Material decodeMaterial(const json& encoded) {
  if (encoded.is_null()) {
    return {};
  }
  Acts::Material::ParametersVector params =
      Acts::Material::ParametersVector::Zero();
  for (auto i = params.size(); 0 < i--;) {
    // .at(...) ensures bound checks
    params[i] = encoded.at(i);
  }
  return Acts::Material(params);
}

// helper functions to encode/decode concrete material slabs
//
// encoded as an object w/ two entries: `material` and `thickness`

json encodeMaterialSlab(const Acts::MaterialSlab& slab) {
  return {
      {"material", encodeMaterial(slab.material())},
      {"thickness", slab.thickness()},
  };
}

Acts::MaterialSlab decodeMaterialSlab(const json& encoded) {
  return Acts::MaterialSlab(decodeMaterial(encoded.at("material")),
                            encoded.at("thickness").get<float>());
}

}  // namespace

Acts::JsonGeometryConverter::JsonGeometryConverter(
    const Acts::JsonGeometryConverter::Config& cfg)
    : m_cfg(std::move(cfg)) {
  // Validate the configuration
  if (!m_cfg.logger) {
    throw std::invalid_argument("Missing logger");
  }
}

std::pair<std::map<Acts::GeometryIdentifier,
                   std::shared_ptr<const Acts::ISurfaceMaterial>>,
          std::map<Acts::GeometryIdentifier,
                   std::shared_ptr<const Acts::IVolumeMaterial>>>
Acts::JsonGeometryConverter::jsonToMaterialMaps(const json& materialmaps) {
  auto& j = materialmaps;
  // The return maps
  std::pair<SurfaceMaterialMap, VolumeMaterialMap> maps;
  ACTS_VERBOSE("j2a: Reading material maps from json file.");
  ACTS_VERBOSE("j2a: Found entries for " << j.count(m_cfg.volkey)
                                         << " volume(s).");
  // Structured binding
  for (auto& [key, value] : j.items()) {
    // Check if this the volume key
    if (key == m_cfg.volkey) {
      // Get the volume json
      auto volj = value;
      for (auto& [vkey, vvalue] : volj.items()) {
        // Create the volume id
        int vid = std::stoi(vkey);
        Acts::GeometryIdentifier volumeID;
        volumeID.setVolume(vid);
        ACTS_VERBOSE("j2a: -> Found Volume " << vid);
        // Loop through the information in the volume
        for (auto& [vckey, vcvalue] : vvalue.items()) {
          if (vckey == m_cfg.boukey and m_cfg.processBoundaries and
              not vcvalue.empty()) {
            ACTS_VERBOSE("j2a: --> BoundarySurface(s) to be parsed");
            for (auto& [bkey, bvalue] : vcvalue.items()) {
              // Create the boundary id
              int bid = std::stoi(bkey);
              Acts::GeometryIdentifier boundaryID(volumeID);
              boundaryID.setBoundary(bid);
              ACTS_VERBOSE("j2a: ---> Found boundary surface " << bid);
              if (bvalue[m_cfg.mapkey] == true) {
                auto boumat = jsonToSurfaceMaterial(bvalue);
                maps.first[boundaryID] =
                    std::shared_ptr<const ISurfaceMaterial>(boumat);
              }
            }
          } else if (vckey == m_cfg.laykey) {
            ACTS_VERBOSE("j2a: --> Layer(s) to be parsed");
            // Loop over layers and repeat
            auto layj = vcvalue;
            for (auto& [lkey, lvalue] : layj.items()) {
              // Create the layer id
              int lid = std::stoi(lkey);
              Acts::GeometryIdentifier layerID(volumeID);
              layerID.setLayer(lid);
              ACTS_VERBOSE("j2a: ---> Found Layer " << lid);
              // Finally loop over layer components
              for (auto& [lckey, lcvalue] : lvalue.items()) {
                if (lckey == m_cfg.repkey and m_cfg.processRepresenting and
                    not lcvalue.empty()) {
                  ACTS_VERBOSE("j2a: ----> Found representing surface");
                  if (lcvalue[m_cfg.mapkey] == true) {
                    auto repmat = jsonToSurfaceMaterial(lcvalue);
                    maps.first[layerID] =
                        std::shared_ptr<const Acts::ISurfaceMaterial>(repmat);
                  }
                } else if (lckey == m_cfg.appkey and m_cfg.processApproaches and
                           not lcvalue.empty()) {
                  ACTS_VERBOSE("j2a: ----> Found approach surface(s)");
                  // Loop over approach surfaces
                  for (auto& [askey, asvalue] : lcvalue.items()) {
                    // Create the layer id, todo set to max value
                    int aid = (askey == "*") ? 0 : std::stoi(askey);
                    Acts::GeometryIdentifier approachID(layerID);
                    approachID.setApproach(aid);
                    ACTS_VERBOSE("j2a: -----> Approach surface " << askey);
                    if (asvalue[m_cfg.mapkey] == true) {
                      auto appmat = jsonToSurfaceMaterial(asvalue);
                      maps.first[approachID] =
                          std::shared_ptr<const Acts::ISurfaceMaterial>(appmat);
                    }
                  }
                } else if (lckey == m_cfg.senkey and m_cfg.processSensitives and
                           not lcvalue.empty()) {
                  ACTS_VERBOSE("j2a: ----> Found sensitive surface(s)");
                  // Loop over sensitive surfaces
                  for (auto& [sskey, ssvalue] : lcvalue.items()) {
                    // Create the layer id, todo set to max value
                    int sid = (sskey == "*") ? 0 : std::stoi(sskey);
                    Acts::GeometryIdentifier senisitiveID(layerID);
                    senisitiveID.setSensitive(sid);
                    ACTS_VERBOSE("j2a: -----> Sensitive surface " << sskey);
                    if (ssvalue[m_cfg.mapkey] == true) {
                      auto senmat = jsonToSurfaceMaterial(ssvalue);
                      maps.first[senisitiveID] =
                          std::shared_ptr<const Acts::ISurfaceMaterial>(senmat);
                    }
                  }
                }
              }
            }

          } else if (m_cfg.processVolumes and vckey == m_cfg.matkey and
                     not vcvalue.empty()) {
            ACTS_VERBOSE("--> VolumeMaterial to be parsed");
            if (vcvalue[m_cfg.mapkey] == true) {
              auto intermat = jsonToVolumeMaterial(vcvalue);
              maps.second[volumeID] =
                  std::shared_ptr<const Acts::IVolumeMaterial>(intermat);
            }
          }
        }
      }
    } else if (key == m_cfg.geoversion) {
      ACTS_VERBOSE("Detector version: " << m_cfg.geoversion);
    }
  }

  // Return the filled maps
  return maps;
}

/// Convert method
///
json Acts::JsonGeometryConverter::materialMapsToJson(
    const DetectorMaterialMaps& maps) {
  DetectorRep detRep;
  // Collect all GeometryIdentifiers per VolumeID for the formatted output
  for (auto& [key, value] : maps.first) {
    geo_id_value vid = key.volume();
    auto volRep = detRep.volumes.find(vid);
    if (volRep == detRep.volumes.end()) {
      detRep.volumes.insert({vid, VolumeRep()});
      volRep = detRep.volumes.find(vid);
      volRep->second.volumeID = key;
    }
    geo_id_value lid = key.layer();
    if (lid != 0) {
      // we are on a layer, get the layer rep
      auto layRep = volRep->second.layers.find(lid);
      if (layRep == volRep->second.layers.end()) {
        volRep->second.layers.insert({lid, LayerRep()});
        layRep = volRep->second.layers.find(lid);
        layRep->second.layerID = key;
      }
      // now insert appropriately
      geo_id_value sid = key.sensitive();
      geo_id_value aid = key.approach();
      if (sid != 0) {
        layRep->second.sensitives.insert({sid, value.get()});
      } else if (aid != 0) {
        layRep->second.approaches.insert({aid, value.get()});
      } else {
        layRep->second.representing = value.get();
      }

    } else {
      // not on a layer can only be a boundary surface
      geo_id_value bid = key.boundary();
      volRep->second.boundaries.insert({bid, value.get()});
    }
  }
  for (auto& [key, value] : maps.second) {
    // find the volume representation
    geo_id_value vid = key.volume();
    auto volRep = detRep.volumes.find(vid);
    if (volRep == detRep.volumes.end()) {
      detRep.volumes.insert({vid, VolumeRep()});
      volRep = detRep.volumes.find(vid);
      volRep->second.volumeID = key;
    }
    volRep->second.material = value.get();
  }
  // convert the detector representation to json format
  return detectorRepToJson(detRep);
}

/// Create Json from a detector represenation
json Acts::JsonGeometryConverter::detectorRepToJson(const DetectorRep& detRep) {
  json detectorj;
  ACTS_VERBOSE("a2j: Writing json from detector representation");
  ACTS_VERBOSE("a2j: Found entries for " << detRep.volumes.size()
                                         << " volume(s).");

  json volumesj;
  for (auto& [key, value] : detRep.volumes) {
    json volj;
    ACTS_VERBOSE("a2j: -> Writing Volume: " << key);
    volj[m_cfg.namekey] = value.volumeName;
    std::ostringstream svolumeID;
    svolumeID << value.volumeID;
    volj[m_cfg.geometryidkey] = svolumeID.str();
    if (m_cfg.processVolumes && value.material) {
      volj[m_cfg.matkey] = volumeMaterialToJson(*value.material);
    }
    // Write the layers
    if (not value.layers.empty()) {
      ACTS_VERBOSE("a2j: ---> Found " << value.layers.size() << " layer(s) ");
      json layersj;
      for (auto& [lkey, lvalue] : value.layers) {
        ACTS_VERBOSE("a2j: ----> Convert layer " << lkey);
        json layj;
        std::ostringstream slayerID;
        slayerID << lvalue.layerID;
        layj[m_cfg.geometryidkey] = slayerID.str();
        // First check for approaches
        if (not lvalue.approaches.empty() and m_cfg.processApproaches) {
          ACTS_VERBOSE("a2j: -----> Found " << lvalue.approaches.size()
                                            << " approach surface(s)");
          json approachesj;
          for (auto& [akey, avalue] : lvalue.approaches) {
            ACTS_VERBOSE("a2j: ------> Convert approach surface " << akey);
            approachesj[std::to_string(akey)] = surfaceMaterialToJson(*avalue);
            if (lvalue.approacheSurfaces.find(akey) !=
                lvalue.approacheSurfaces.end())
              addSurfaceToJson(approachesj[std::to_string(akey)],
                               lvalue.approacheSurfaces.at(akey));
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
            sensitivesj[std::to_string(skey)] = surfaceMaterialToJson(*svalue);
            if (lvalue.sensitiveSurfaces.find(skey) !=
                lvalue.sensitiveSurfaces.end())
              addSurfaceToJson(sensitivesj[std::to_string(skey)],
                               lvalue.sensitiveSurfaces.at(skey));
          }
          // Add to the layer json
          layj[m_cfg.senkey] = sensitivesj;
        }
        // Finally check for representing
        if (lvalue.representing != nullptr and m_cfg.processRepresenting) {
          ACTS_VERBOSE("a2j: ------> Convert representing surface ");
          layj[m_cfg.repkey] = surfaceMaterialToJson(*lvalue.representing);
          if (lvalue.representingSurface != nullptr)
            addSurfaceToJson(layj[m_cfg.repkey], lvalue.representingSurface);
        }
        layersj[std::to_string(lkey)] = layj;
      }
      volj[m_cfg.laykey] = layersj;
    }
    // Write the boundary surfaces
    if (not value.boundaries.empty()) {
      ACTS_VERBOSE("a2j: ---> Found " << value.boundaries.size()
                                      << " boundary/ies ");
      json boundariesj;
      for (auto& [bkey, bvalue] : value.boundaries) {
        ACTS_VERBOSE("a2j: ----> Convert boundary " << bkey);
        boundariesj[std::to_string(bkey)] = surfaceMaterialToJson(*bvalue);
        if (value.boundarySurfaces.find(bkey) != value.boundarySurfaces.end())
          addSurfaceToJson(boundariesj[std::to_string(bkey)],
                           value.boundarySurfaces.at(bkey));
      }
      volj[m_cfg.boukey] = boundariesj;
    }

    volumesj[std::to_string(key)] = volj;
  }
  // Assign the volume json to the detector json
  detectorj[m_cfg.volkey] = volumesj;

  return detectorj;
}

/// Create the Surface Material
const Acts::ISurfaceMaterial*
Acts::JsonGeometryConverter::jsonToSurfaceMaterial(const json& material) {
  Acts::ISurfaceMaterial* sMaterial = nullptr;
  // The bin utility for deescribing the data
  Acts::BinUtility bUtility;
  for (auto& [key, value] : material.items()) {
    if (key == m_cfg.transfokeys and not value.empty()) {
      bUtility = Acts::BinUtility(jsonToTransform(value));
      break;
    }
  }
  // Convert the material
  Acts::MaterialSlabMatrix mpMatrix;
  // Structured binding
  for (auto& [key, value] : material.items()) {
    // Check json keys
    if (key == m_cfg.bin0key and not value.empty()) {
      bUtility += jsonToBinUtility(value);
    } else if (key == m_cfg.bin1key and not value.empty()) {
      bUtility += jsonToBinUtility(value);
    }
    if (key == m_cfg.datakey and not value.empty()) {
      mpMatrix = jsonToMaterialMatrix(value);
    }
  }

  // We have protoMaterial
  if (mpMatrix.empty()) {
    sMaterial = new Acts::ProtoSurfaceMaterial(bUtility);
  } else if (bUtility.bins() == 1) {
    sMaterial = new Acts::HomogeneousSurfaceMaterial(mpMatrix[0][0]);
  } else {
    sMaterial = new Acts::BinnedSurfaceMaterial(bUtility, mpMatrix);
  }
  // return what you have
  return sMaterial;
}

/// Create the Volume Material
const Acts::IVolumeMaterial* Acts::JsonGeometryConverter::jsonToVolumeMaterial(
    const json& material) {
  Acts::IVolumeMaterial* vMaterial = nullptr;
  // The bin utility for deescribing the data
  Acts::BinUtility bUtility;
  for (auto& [key, value] : material.items()) {
    if (key == m_cfg.transfokeys and not value.empty()) {
      bUtility = Acts::BinUtility(jsonToTransform(value));
      break;
    }
  }
  // Convert the material
  std::vector<Material> mmat;
  // Structured binding
  for (auto& [key, value] : material.items()) {
    // Check json keys
    if (key == m_cfg.bin0key and not value.empty()) {
      bUtility += jsonToBinUtility(value);
    } else if (key == m_cfg.bin1key and not value.empty()) {
      bUtility += jsonToBinUtility(value);
    } else if (key == m_cfg.bin2key and not value.empty()) {
      bUtility += jsonToBinUtility(value);
    }
    if (key == m_cfg.datakey and not value.empty()) {
      for (const auto& bin : value) {
        mmat.push_back(decodeMaterial(bin));
      }
    }
  }

  // We have protoMaterial
  if (mmat.empty()) {
    vMaterial = new Acts::ProtoVolumeMaterial(bUtility);
  } else if (mmat.size() == 1) {
    vMaterial = new Acts::HomogeneousVolumeMaterial(mmat[0]);
  } else {
    if (bUtility.dimensions() == 2) {
      std::function<Acts::Vector2D(Acts::Vector3D)> transfoGlobalToLocal;
      Acts::Grid2D grid = createGrid2D(bUtility, transfoGlobalToLocal);

      Acts::Grid2D::point_t min = grid.minPosition();
      Acts::Grid2D::point_t max = grid.maxPosition();
      Acts::Grid2D::index_t nBins = grid.numLocalBins();

      Acts::EAxis axis1(min[0], max[0], nBins[0]);
      Acts::EAxis axis2(min[1], max[1], nBins[1]);

      // Build the grid and fill it with data
      MaterialGrid2D mGrid(std::make_tuple(axis1, axis2));

      for (size_t bin = 0; bin < mmat.size(); bin++) {
        mGrid.at(bin) = mmat[bin].parameters();
      }
      MaterialMapper<MaterialGrid2D> matMap(transfoGlobalToLocal, mGrid);
      vMaterial =
          new Acts::InterpolatedMaterialMap<MaterialMapper<MaterialGrid2D>>(
              std::move(matMap), bUtility);
    } else if (bUtility.dimensions() == 3) {
      std::function<Acts::Vector3D(Acts::Vector3D)> transfoGlobalToLocal;
      Acts::Grid3D grid = createGrid3D(bUtility, transfoGlobalToLocal);

      Acts::Grid3D::point_t min = grid.minPosition();
      Acts::Grid3D::point_t max = grid.maxPosition();
      Acts::Grid3D::index_t nBins = grid.numLocalBins();

      Acts::EAxis axis1(min[0], max[0], nBins[0]);
      Acts::EAxis axis2(min[1], max[1], nBins[1]);
      Acts::EAxis axis3(min[2], max[2], nBins[2]);

      // Build the grid and fill it with data
      MaterialGrid3D mGrid(std::make_tuple(axis1, axis2, axis3));

      for (size_t bin = 0; bin < mmat.size(); bin++) {
        mGrid.at(bin) = mmat[bin].parameters();
      }
      MaterialMapper<MaterialGrid3D> matMap(transfoGlobalToLocal, mGrid);
      vMaterial =
          new Acts::InterpolatedMaterialMap<MaterialMapper<MaterialGrid3D>>(
              std::move(matMap), bUtility);
    }
  }
  // return what you have
  return vMaterial;
}

json Acts::JsonGeometryConverter::trackingGeometryToJson(
    const Acts::TrackingGeometry& tGeometry) {
  DetectorRep detRep;
  convertToRep(detRep, *tGeometry.highestTrackingVolume());
  return detectorRepToJson(detRep);
}

void Acts::JsonGeometryConverter::convertToRep(
    DetectorRep& detRep, const Acts::TrackingVolume& tVolume) {
  // The writer reader volume representation
  VolumeRep volRep;
  volRep.volumeName = tVolume.volumeName();
  // there are confined volumes
  if (tVolume.confinedVolumes() != nullptr) {
    // get through the volumes
    auto& volumes = tVolume.confinedVolumes()->arrayObjects();
    // loop over the volumes
    for (auto& vol : volumes) {
      // recursive call
      convertToRep(detRep, *vol);
    }
  }
  // there are dense volumes
  if (m_cfg.processDenseVolumes && !tVolume.denseVolumes().empty()) {
    // loop over the volumes
    for (auto& vol : tVolume.denseVolumes()) {
      // recursive call
      convertToRep(detRep, *vol);
    }
  }
  // Get the volume Id
  Acts::GeometryIdentifier volumeID = tVolume.geometryId();
  geo_id_value vid = volumeID.volume();

  // Write the material if there's one
  if (tVolume.volumeMaterial() != nullptr) {
    volRep.material = tVolume.volumeMaterial();
  } else if (m_cfg.processnonmaterial == true) {
    Acts::BinUtility bUtility = DefaultBin(tVolume);
    Acts::IVolumeMaterial* bMaterial = new Acts::ProtoVolumeMaterial(bUtility);
    volRep.material = bMaterial;
  }
  // there are confied layers
  if (tVolume.confinedLayers() != nullptr) {
    // get the layers
    auto& layers = tVolume.confinedLayers()->arrayObjects();
    // loop of the volumes
    for (auto& lay : layers) {
      auto layRep = convertToRep(*lay);
      if (layRep) {
        // it's a valid representation so let's go with it
        Acts::GeometryIdentifier layerID = lay->geometryId();
        geo_id_value lid = layerID.layer();
        volRep.layers.insert({lid, std::move(layRep)});
      }
    }
  }
  // Let's finally check the boundaries
  for (auto& bsurf : tVolume.boundarySurfaces()) {
    // the surface representation
    auto& bssfRep = bsurf->surfaceRepresentation();
    if (bssfRep.surfaceMaterial() != nullptr) {
      Acts::GeometryIdentifier boundaryID = bssfRep.geometryId();
      geo_id_value bid = boundaryID.boundary();
      // Ignore if the volumeID is not correct (i.e. shared boundary)
      // if (boundaryID.value(Acts::GeometryIdentifier::volume_mask) == vid){
      volRep.boundaries[bid] = bssfRep.surfaceMaterial();
      volRep.boundarySurfaces[bid] = &bssfRep;
      // }
    } else if (m_cfg.processnonmaterial == true) {
      // if no material suface exist add a default one for the mapping
      // configuration
      Acts::GeometryIdentifier boundaryID = bssfRep.geometryId();
      geo_id_value bid = boundaryID.boundary();
      Acts::BinUtility bUtility = DefaultBin(bssfRep);
      Acts::ISurfaceMaterial* bMaterial =
          new Acts::ProtoSurfaceMaterial(bUtility);
      volRep.boundaries[bid] = bMaterial;
      volRep.boundarySurfaces[bid] = &bssfRep;
    }
  }
  // Write if it's good
  if (volRep) {
    volRep.volumeName = tVolume.volumeName();
    volRep.volumeID = volumeID;
    detRep.volumes.insert({vid, std::move(volRep)});
  }
  return;
}

Acts::JsonGeometryConverter::LayerRep Acts::JsonGeometryConverter::convertToRep(
    const Acts::Layer& tLayer) {
  LayerRep layRep;
  // fill layer ID information
  layRep.layerID = tLayer.geometryId();
  if (m_cfg.processSensitives and tLayer.surfaceArray() != nullptr) {
    for (auto& ssf : tLayer.surfaceArray()->surfaces()) {
      if (ssf != nullptr && ssf->surfaceMaterial() != nullptr) {
        Acts::GeometryIdentifier sensitiveID = ssf->geometryId();
        geo_id_value sid = sensitiveID.sensitive();
        layRep.sensitives.insert({sid, ssf->surfaceMaterial()});
        layRep.sensitiveSurfaces.insert({sid, ssf});
      } else if (m_cfg.processnonmaterial == true) {
        // if no material suface exist add a default one for the mapping
        // configuration
        Acts::GeometryIdentifier sensitiveID = ssf->geometryId();
        geo_id_value sid = sensitiveID.sensitive();
        Acts::BinUtility sUtility = DefaultBin(*ssf);
        Acts::ISurfaceMaterial* sMaterial =
            new Acts::ProtoSurfaceMaterial(sUtility);
        layRep.sensitives.insert({sid, sMaterial});
        layRep.sensitiveSurfaces.insert({sid, ssf});
      }
    }
  }
  // the representing
  if (!(tLayer.surfaceRepresentation().geometryId() == GeometryIdentifier())) {
    if (tLayer.surfaceRepresentation().surfaceMaterial() != nullptr) {
      layRep.representing = tLayer.surfaceRepresentation().surfaceMaterial();
      layRep.representingSurface = &tLayer.surfaceRepresentation();
    } else if (m_cfg.processnonmaterial == true) {
      // if no material suface exist add a default one for the mapping
      // configuration
      Acts::BinUtility rUtility = DefaultBin(tLayer.surfaceRepresentation());
      Acts::ISurfaceMaterial* rMaterial =
          new Acts::ProtoSurfaceMaterial(rUtility);
      layRep.representing = rMaterial;
      layRep.representingSurface = &tLayer.surfaceRepresentation();
    }
  }
  // the approach
  if (tLayer.approachDescriptor() != nullptr) {
    for (auto& asf : tLayer.approachDescriptor()->containedSurfaces()) {
      // get the surface and check for material
      if (asf->surfaceMaterial() != nullptr) {
        Acts::GeometryIdentifier approachID = asf->geometryId();
        geo_id_value aid = approachID.approach();
        layRep.approaches.insert({aid, asf->surfaceMaterial()});
        layRep.approacheSurfaces.insert({aid, asf});
      } else if (m_cfg.processnonmaterial == true) {
        // if no material suface exist add a default one for the mapping
        // configuration
        Acts::GeometryIdentifier approachID = asf->geometryId();
        geo_id_value aid = approachID.approach();
        Acts::BinUtility aUtility = DefaultBin(*asf);
        Acts::ISurfaceMaterial* aMaterial =
            new Acts::ProtoSurfaceMaterial(aUtility);
        layRep.approaches.insert({aid, aMaterial});
        layRep.approacheSurfaces.insert({aid, asf});
      }
    }
  }
  // return the layer representation
  return layRep;
}

json Acts::JsonGeometryConverter::surfaceMaterialToJson(
    const Acts::ISurfaceMaterial& sMaterial) {
  json smj;
  // A bin utility needs to be written
  const Acts::BinUtility* bUtility = nullptr;
  // Check if we have a proto material
  auto psMaterial = dynamic_cast<const Acts::ProtoSurfaceMaterial*>(&sMaterial);
  if (psMaterial != nullptr) {
    // Type is proto material
    smj[m_cfg.typekey] = "proto";
    // by default the protoMaterial is not used for mapping
    smj[m_cfg.mapkey] = false;
    bUtility = &(psMaterial->binUtility());
  } else {
    // Now check if we have a homogeneous material
    auto hsMaterial =
        dynamic_cast<const Acts::HomogeneousSurfaceMaterial*>(&sMaterial);
    if (hsMaterial != nullptr) {
      // type is homogeneous
      smj[m_cfg.typekey] = "homogeneous";
      smj[m_cfg.mapkey] = true;
      if (m_cfg.writeData) {
        smj[m_cfg.datakey] = json::array({
            json::array({
                encodeMaterialSlab(hsMaterial->materialSlab(0, 0)),
            }),
        });
      }
    } else {
      // Only option remaining: BinnedSurface material
      auto bsMaterial =
          dynamic_cast<const Acts::BinnedSurfaceMaterial*>(&sMaterial);
      if (bsMaterial != nullptr) {
        // type is binned
        smj[m_cfg.typekey] = "binned";
        smj[m_cfg.mapkey] = true;
        bUtility = &(bsMaterial->binUtility());
        // convert the data
        // get the material matrix
        if (m_cfg.writeData) {
          json mmat = json::array();
          for (const auto& mpVector : bsMaterial->fullMaterial()) {
            json mvec = json::array();
            for (const auto& mp : mpVector) {
              mvec.push_back(encodeMaterialSlab(mp));
            }
            mmat.push_back(std::move(mvec));
          }
          smj[m_cfg.datakey] = std::move(mmat);
        }
      }
    }
  }
  // add the bin utility
  if (bUtility != nullptr && !bUtility->binningData().empty()) {
    std::vector<std::string> binkeys = {m_cfg.bin0key, m_cfg.bin1key};
    // loop over dimensions and write
    auto& binningData = bUtility->binningData();
    // loop over the dimensions
    for (size_t ibin = 0; ibin < binningData.size(); ++ibin) {
      json binj;
      auto cbData = binningData[ibin];
      binj.push_back(Acts::binningValueNames[cbData.binvalue]);
      if (cbData.option == Acts::closed) {
        binj.push_back("closed");
      } else {
        binj.push_back("open");
      }
      binj.push_back(cbData.bins());
      // If protoMaterial has a non uniform binning (non default) then it is
      // used by default in the mapping
      if (smj[m_cfg.typekey] == "proto" && cbData.bins() > 1)
        smj[m_cfg.mapkey] = true;
      // If it's not a proto map, write min / max
      if (smj[m_cfg.typekey] != "proto") {
        std::pair<double, double> minMax = {cbData.min, cbData.max};
        binj.push_back(minMax);
      }
      smj[binkeys[ibin]] = binj;
    }
    std::vector<double> transfo;
    Acts::Transform3D transfo_matrix = bUtility->transform();
    if (not transfo_matrix.isApprox(Acts::Transform3D::Identity())) {
      for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
          transfo.push_back(transfo_matrix(j, i));
        }
      }
      smj[m_cfg.transfokeys] = transfo;
    }
  }
  return smj;
}

json Acts::JsonGeometryConverter::volumeMaterialToJson(
    const Acts::IVolumeMaterial& vMaterial) {
  json smj;
  // A bin utility needs to be written
  const Acts::BinUtility* bUtility = nullptr;
  // Check if we have a proto material
  auto pvMaterial = dynamic_cast<const Acts::ProtoVolumeMaterial*>(&vMaterial);
  if (pvMaterial != nullptr) {
    // Type is proto material
    smj[m_cfg.typekey] = "proto";
    // by default the protoMaterial is not used for mapping
    smj[m_cfg.mapkey] = false;
    bUtility = &(pvMaterial->binUtility());
  } else {
    // Now check if we have a homogeneous material
    auto hvMaterial =
        dynamic_cast<const Acts::HomogeneousVolumeMaterial*>(&vMaterial);
    if (hvMaterial != nullptr) {
      // type is homogeneous
      smj[m_cfg.typekey] = "homogeneous";
      smj[m_cfg.mapkey] = true;
      if (m_cfg.writeData) {
        // array of encoded materials w/ one entry
        smj[m_cfg.datakey] = json::array({
            encodeMaterial(hvMaterial->material({0, 0, 0})),
        });
      }
    } else {
      // Only option remaining: material map
      auto bvMaterial2D = dynamic_cast<
          const Acts::InterpolatedMaterialMap<MaterialMapper<MaterialGrid2D>>*>(
          &vMaterial);
      // Now check if we have a 2D map
      if (bvMaterial2D != nullptr) {
        // type is binned
        smj[m_cfg.typekey] = "interpolated2D";
        smj[m_cfg.mapkey] = true;
        bUtility = &(bvMaterial2D->binUtility());
        // convert the data
        if (m_cfg.writeData) {
          json mmat = json::array();
          MaterialGrid2D grid = bvMaterial2D->getMapper().getGrid();
          for (size_t bin = 0; bin < grid.size(); bin++) {
            mmat.push_back(encodeMaterial(grid.at(bin)));
          }
          smj[m_cfg.datakey] = std::move(mmat);
        }
      } else {
        // Only option remaining: material map
        auto bvMaterial3D = dynamic_cast<const Acts::InterpolatedMaterialMap<
            MaterialMapper<MaterialGrid3D>>*>(&vMaterial);
        // Now check if we have a 3D map
        if (bvMaterial3D != nullptr) {
          // type is binned
          smj[m_cfg.typekey] = "interpolated3D";
          smj[m_cfg.mapkey] = true;
          bUtility = &(bvMaterial3D->binUtility());
          // convert the data
          if (m_cfg.writeData) {
            json mmat = json::array();
            MaterialGrid3D grid = bvMaterial3D->getMapper().getGrid();
            for (size_t bin = 0; bin < grid.size(); bin++) {
              mmat.push_back(encodeMaterial(grid.at(bin)));
            }
            smj[m_cfg.datakey] = std::move(mmat);
          }
        }
      }
    }
  }
  // add the bin utility
  if (bUtility != nullptr && !bUtility->binningData().empty()) {
    std::vector<std::string> binkeys = {m_cfg.bin0key, m_cfg.bin1key,
                                        m_cfg.bin2key};
    // loop over dimensions and write
    auto& binningData = bUtility->binningData();
    // loop over the dimensions
    for (size_t ibin = 0; ibin < binningData.size(); ++ibin) {
      json binj;
      auto cbData = binningData[ibin];
      binj.push_back(Acts::binningValueNames[cbData.binvalue]);
      if (cbData.option == Acts::closed) {
        binj.push_back("closed");
      } else {
        binj.push_back("open");
      }
      binj.push_back(cbData.bins());
      // If protoMaterial has a non uniform binning (non default) then it is
      // used by default in the mapping
      if (smj[m_cfg.typekey] == "proto" && cbData.bins() > 1)
        smj[m_cfg.mapkey] = true;
      // If it's not a proto map, write min / max
      if (smj[m_cfg.typekey] != "proto") {
        std::pair<double, double> minMax = {cbData.min, cbData.max};
        binj.push_back(minMax);
      }
      smj[binkeys[ibin]] = binj;
    }
    std::vector<double> transfo;
    Acts::Transform3D transfo_matrix = bUtility->transform();
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        transfo.push_back(transfo_matrix(j, i));
      }
    }
    smj[m_cfg.transfokeys] = transfo;
  }
  return smj;
}

void Acts::JsonGeometryConverter::addSurfaceToJson(json& sjson,
                                                   const Surface* surface) {
  // Get the ID of the surface (redundant but help readability)
  std::ostringstream SurfaceID;
  SurfaceID << surface->geometryId();
  sjson[m_cfg.surfacegeometryidkey] = SurfaceID.str();

  // Cast the surface bound to both disk and cylinder
  const Acts::SurfaceBounds& surfaceBounds = surface->bounds();
  auto sTransform = surface->transform(GeometryContext());

  const Acts::RadialBounds* radialBounds =
      dynamic_cast<const Acts::RadialBounds*>(&surfaceBounds);
  const Acts::CylinderBounds* cylinderBounds =
      dynamic_cast<const Acts::CylinderBounds*>(&surfaceBounds);
  const Acts::AnnulusBounds* annulusBounds =
      dynamic_cast<const Acts::AnnulusBounds*>(&surfaceBounds);

  if (radialBounds != nullptr) {
    sjson[m_cfg.surfacetypekey] = "Disk";
    sjson[m_cfg.surfacepositionkey] = sTransform.translation().z();
    sjson[m_cfg.surfacerangekey] = {radialBounds->rMin(), radialBounds->rMax()};
  }
  if (cylinderBounds != nullptr) {
    sjson[m_cfg.surfacetypekey] = "Cylinder";
    sjson[m_cfg.surfacepositionkey] = cylinderBounds->get(CylinderBounds::eR);
    sjson[m_cfg.surfacerangekey] = {
        -1 * cylinderBounds->get(CylinderBounds::eHalfLengthZ),
        cylinderBounds->get(CylinderBounds::eHalfLengthZ)};
  }
  if (annulusBounds != nullptr) {
    sjson[m_cfg.surfacetypekey] = "Annulus";
    sjson[m_cfg.surfacepositionkey] = sTransform.translation().z();
    sjson[m_cfg.surfacerangekey] = {
        {annulusBounds->rMin(), annulusBounds->rMax()},
        {annulusBounds->phiMin(), annulusBounds->phiMax()}};
  }
}

/// Create the Material Matrix
Acts::MaterialSlabMatrix Acts::JsonGeometryConverter::jsonToMaterialMatrix(
    const json& data) {
  Acts::MaterialSlabMatrix mpMatrix;
  // the input data must be array[array[object]]
  for (auto& outer : data) {
    Acts::MaterialSlabVector mpVector;
    for (auto& inner : outer) {
      mpVector.emplace_back(decodeMaterialSlab(inner));
    }
    mpMatrix.push_back(std::move(mpVector));
  }
  return mpMatrix;
}

/// Create the BinUtility for this
Acts::BinUtility Acts::JsonGeometryConverter::jsonToBinUtility(
    const json& bin) {
  if (bin.size() >= 3) {
    // finding the iterator position to determine the binning value
    auto bit = std::find(Acts::binningValueNames.begin(),
                         Acts::binningValueNames.end(), bin[0]);
    size_t indx = std::distance(Acts::binningValueNames.begin(), bit);
    Acts::BinningValue bval = Acts::BinningValue(indx);
    Acts::BinningOption bopt = bin[1] == "open" ? Acts::open : Acts::closed;
    unsigned int bins = bin[2];
    float min = 0;
    float max = 0;
    if (bin.size() >= 4 && bin[3].size() == 2) {
      min = bin[3][0];
      max = bin[3][1];
    }
    return Acts::BinUtility(bins, min, max, bopt, bval);
  }
  return Acts::BinUtility();
}

/// Create the local to global transform
Acts::Transform3D Acts::JsonGeometryConverter::jsonToTransform(
    const json& transfo) {
  Transform3D transform;
  int i = 0;
  int j = 0;
  for (auto& element : transfo) {
    transform(j, i) = element;
    j++;
    if (j == 4) {
      i++;
      j = 0;
    }
  }
  return transform;
}

Acts::BinUtility Acts::JsonGeometryConverter::DefaultBin(
    const Acts::Surface& surface) {
  Acts::BinUtility bUtility;

  const Acts::SurfaceBounds& surfaceBounds = surface.bounds();
  const Acts::RadialBounds* radialBounds =
      dynamic_cast<const Acts::RadialBounds*>(&surfaceBounds);
  const Acts::CylinderBounds* cylinderBounds =
      dynamic_cast<const Acts::CylinderBounds*>(&surfaceBounds);
  const Acts::AnnulusBounds* annulusBounds =
      dynamic_cast<const Acts::AnnulusBounds*>(&surfaceBounds);
  const Acts::RectangleBounds* rectangleBounds =
      dynamic_cast<const Acts::RectangleBounds*>(&surfaceBounds);

  if (radialBounds != nullptr) {
    bUtility += BinUtility(
        1,
        radialBounds->get(RadialBounds::eAveragePhi) -
            radialBounds->get(RadialBounds::eHalfPhiSector),
        radialBounds->get(RadialBounds::eAveragePhi) +
            radialBounds->get(RadialBounds::eHalfPhiSector),
        (radialBounds->get(RadialBounds::eHalfPhiSector) - M_PI) < s_epsilon
            ? Acts::closed
            : Acts::open,
        Acts::binPhi);
    bUtility += BinUtility(1, radialBounds->rMin(), radialBounds->rMax(),
                           Acts::open, Acts::binR);
    return bUtility;
  }
  if (cylinderBounds != nullptr) {
    bUtility += BinUtility(
        1,
        cylinderBounds->get(CylinderBounds::eAveragePhi) -
            cylinderBounds->get(CylinderBounds::eHalfPhiSector),
        cylinderBounds->get(CylinderBounds::eAveragePhi) +
            cylinderBounds->get(CylinderBounds::eHalfPhiSector),
        (cylinderBounds->get(CylinderBounds::eHalfPhiSector) - M_PI) < s_epsilon
            ? Acts::closed
            : Acts::open,
        Acts::binPhi);
    bUtility +=
        BinUtility(1, -1 * cylinderBounds->get(CylinderBounds::eHalfLengthZ),
                   cylinderBounds->get(CylinderBounds::eHalfLengthZ),
                   Acts::open, Acts::binZ);
    return bUtility;
  }
  if (annulusBounds != nullptr) {
    bUtility += BinUtility(1, annulusBounds->get(AnnulusBounds::eMinPhiRel),
                           annulusBounds->get(AnnulusBounds::eMaxPhiRel),
                           Acts::open, Acts::binPhi);
    bUtility += BinUtility(1, annulusBounds->rMin(), annulusBounds->rMax(),
                           Acts::open, Acts::binR);
    return bUtility;
  }
  if (rectangleBounds != nullptr) {
    bUtility += BinUtility(1, rectangleBounds->get(RectangleBounds::eMinX),
                           rectangleBounds->get(RectangleBounds::eMaxX),
                           Acts::open, Acts::binX);
    bUtility += BinUtility(1, rectangleBounds->get(RectangleBounds::eMinY),
                           rectangleBounds->get(RectangleBounds::eMaxY),
                           Acts::open, Acts::binY);
    return bUtility;
  }
  ACTS_INFO(
      "No corresponding bound found for the surface : " << surface.name());
  return bUtility;
}

Acts::BinUtility Acts::JsonGeometryConverter::DefaultBin(
    const Acts::TrackingVolume& volume) {
  Acts::BinUtility bUtility;

  auto cyBounds =
      dynamic_cast<const CylinderVolumeBounds*>(&(volume.volumeBounds()));
  auto cutcylBounds =
      dynamic_cast<const CutoutCylinderVolumeBounds*>(&(volume.volumeBounds()));
  auto cuBounds =
      dynamic_cast<const CuboidVolumeBounds*>(&(volume.volumeBounds()));

  if (cyBounds != nullptr) {
    bUtility += BinUtility(1, cyBounds->get(CylinderVolumeBounds::eMinR),
                           cyBounds->get(CylinderVolumeBounds::eMaxR),
                           Acts::open, Acts::binR);
    bUtility += BinUtility(
        1, -cyBounds->get(CylinderVolumeBounds::eHalfPhiSector),
        cyBounds->get(CylinderVolumeBounds::eHalfPhiSector),
        (cyBounds->get(CylinderVolumeBounds::eHalfPhiSector) - M_PI) < s_epsilon
            ? Acts::closed
            : Acts::open,
        Acts::binPhi);
    bUtility +=
        BinUtility(1, -cyBounds->get(CylinderVolumeBounds::eHalfLengthZ),
                   cyBounds->get(CylinderVolumeBounds::eHalfLengthZ),
                   Acts::open, Acts::binZ);
    return bUtility;
  }
  if (cutcylBounds != nullptr) {
    bUtility +=
        BinUtility(1, cutcylBounds->get(CutoutCylinderVolumeBounds::eMinR),
                   cutcylBounds->get(CutoutCylinderVolumeBounds::eMaxR),
                   Acts::open, Acts::binR);
    bUtility += BinUtility(1, -M_PI, M_PI, Acts::closed, Acts::binPhi);
    bUtility += BinUtility(
        1, -cutcylBounds->get(CutoutCylinderVolumeBounds::eHalfLengthZ),
        cutcylBounds->get(CutoutCylinderVolumeBounds::eHalfLengthZ), Acts::open,
        Acts::binZ);
    return bUtility;
  } else if (cuBounds != nullptr) {
    bUtility += BinUtility(1, -cuBounds->get(CuboidVolumeBounds::eHalfLengthX),
                           cuBounds->get(CuboidVolumeBounds::eHalfLengthX),
                           Acts::open, Acts::binX);
    bUtility += BinUtility(1, -cuBounds->get(CuboidVolumeBounds::eHalfLengthY),
                           cuBounds->get(CuboidVolumeBounds::eHalfLengthY),
                           Acts::open, Acts::binY);
    bUtility += BinUtility(1, -cuBounds->get(CuboidVolumeBounds::eHalfLengthZ),
                           cuBounds->get(CuboidVolumeBounds::eHalfLengthZ),
                           Acts::open, Acts::binZ);
    return bUtility;
  }
  ACTS_INFO(
      "No corresponding bound found for the volume : " << volume.volumeName());
  return bUtility;
}
