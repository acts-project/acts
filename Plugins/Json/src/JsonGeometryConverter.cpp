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
using jsonKey = Acts::JsonGeometryConverter::jsonKey;

using volumeMaterialPointer = const Acts::IVolumeMaterial*;
using surfaceMaterialPointer = const Acts::ISurfaceMaterial*;

using volumePointer = const Acts::TrackingVolume*;
using surfacePointer = const Acts::Surface*;

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

/// Create the Material Matrix
Acts::MaterialSlabMatrix decodeMaterialMatrix(
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

// Add the bin utility to the json representation
void encodeBinUtility(json& smj, const Acts::BinUtility* bUtility) {

  if (bUtility != nullptr && !bUtility->binningData().empty()) {
    std::vector<std::string> binkeys = {jsonKey().bin0key, jsonKey().bin1key,
                                        jsonKey().bin2key};
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
      if (smj[jsonKey().typekey] == "proto" && cbData.bins() > 1)
        smj[jsonKey().mapkey] = true;
      // If it's not a proto map, write min / max
      if (smj[jsonKey().typekey] != "proto") {
        std::pair<double, double> minMax = {cbData.min, cbData.max};
        binj.push_back(minMax);
      }
      smj[binkeys[ibin]] = binj;
    }
    // write the transform matrix
    std::vector<double> transfo;
    Acts::Transform3 transfo_matrix = bUtility->transform();
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        transfo.push_back(transfo_matrix(j, i));
      }
    }
    smj[jsonKey().transfokeys] = transfo;
  }
}

/// Create the BinUtility for this
Acts::BinUtility decodeBinUtility(
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
Acts::Transform3 decodeTransform(
    const json& transfo) {
  Acts::Transform3 transform;
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

Acts::BinUtility defaultSurfaceBin(
    const Acts::Surface* surface) {
  Acts::BinUtility bUtility;

  // Check which type of bounds is associated to the surface 
  const Acts::SurfaceBounds& surfaceBounds = surface->bounds();
  const Acts::RadialBounds* radialBounds =
      dynamic_cast<const Acts::RadialBounds*>(&surfaceBounds);
  const Acts::CylinderBounds* cylinderBounds =
      dynamic_cast<const Acts::CylinderBounds*>(&surfaceBounds);
  const Acts::AnnulusBounds* annulusBounds =
      dynamic_cast<const Acts::AnnulusBounds*>(&surfaceBounds);
  const Acts::RectangleBounds* rectangleBounds =
      dynamic_cast<const Acts::RectangleBounds*>(&surfaceBounds);

  if (radialBounds != nullptr) {
    bUtility += Acts::BinUtility(
        1,
        radialBounds->get(Acts::RadialBounds::eAveragePhi) -
            radialBounds->get(Acts::RadialBounds::eHalfPhiSector),
        radialBounds->get(Acts::RadialBounds::eAveragePhi) +
            radialBounds->get(Acts::RadialBounds::eHalfPhiSector),
        (radialBounds->get(Acts::RadialBounds::eHalfPhiSector) - M_PI) < Acts::s_epsilon
            ? Acts::closed
            : Acts::open,
        Acts::binPhi);
    bUtility += Acts::BinUtility(1, radialBounds->rMin(), radialBounds->rMax(),
                           Acts::open, Acts::binR);
    return bUtility;
  }
  if (cylinderBounds != nullptr) {
    bUtility += Acts::BinUtility(
        1,
        cylinderBounds->get(Acts::CylinderBounds::eAveragePhi) -
            cylinderBounds->get(Acts::CylinderBounds::eHalfPhiSector),
        cylinderBounds->get(Acts::CylinderBounds::eAveragePhi) +
            cylinderBounds->get(Acts::CylinderBounds::eHalfPhiSector),
        (cylinderBounds->get(Acts::CylinderBounds::eHalfPhiSector) - M_PI) < Acts::s_epsilon
            ? Acts::closed
            : Acts::open,
        Acts::binPhi);
    bUtility +=
        Acts::BinUtility(1, -1 * cylinderBounds->get(Acts::CylinderBounds::eHalfLengthZ),
                   cylinderBounds->get(Acts::CylinderBounds::eHalfLengthZ),
                   Acts::open, Acts::binZ);
    return bUtility;
  }
  if (annulusBounds != nullptr) {
    bUtility += Acts::BinUtility(1, annulusBounds->get(Acts::AnnulusBounds::eMinPhiRel),
                           annulusBounds->get(Acts::AnnulusBounds::eMaxPhiRel),
                           Acts::open, Acts::binPhi);
    bUtility += Acts::BinUtility(1, annulusBounds->rMin(), annulusBounds->rMax(),
                           Acts::open, Acts::binR);
    return bUtility;
  }
  if (rectangleBounds != nullptr) {
    bUtility += Acts::BinUtility(1, rectangleBounds->get(Acts::RectangleBounds::eMinX),
                           rectangleBounds->get(Acts::RectangleBounds::eMaxX),
                           Acts::open, Acts::binX);
    bUtility += Acts::BinUtility(1, rectangleBounds->get(Acts::RectangleBounds::eMinY),
                           rectangleBounds->get(Acts::RectangleBounds::eMaxY),
                           Acts::open, Acts::binY);
    return bUtility;
  }
  return bUtility;
}

Acts::BinUtility defaultVolumeBin(
    const Acts::TrackingVolume* volume) {
  Acts::BinUtility bUtility;

  // Check which type of bound is associated to the volume 
  auto cyBounds =
      dynamic_cast<const Acts::CylinderVolumeBounds*>(&(volume->volumeBounds()));
  auto cutcylBounds =
      dynamic_cast<const Acts::CutoutCylinderVolumeBounds*>(&(volume->volumeBounds()));
  auto cuBounds =
      dynamic_cast<const Acts::CuboidVolumeBounds*>(&(volume->volumeBounds()));

  if (cyBounds != nullptr) {
    bUtility += Acts::BinUtility(1, cyBounds->get(Acts::CylinderVolumeBounds::eMinR),
                           cyBounds->get(Acts::CylinderVolumeBounds::eMaxR),
                           Acts::open, Acts::binR);
    bUtility += Acts::BinUtility(
        1, -cyBounds->get(Acts::CylinderVolumeBounds::eHalfPhiSector),
        cyBounds->get(Acts::CylinderVolumeBounds::eHalfPhiSector),
        (cyBounds->get(Acts::CylinderVolumeBounds::eHalfPhiSector) - M_PI) < Acts::s_epsilon
            ? Acts::closed
            : Acts::open,
        Acts::binPhi);
    bUtility +=
        Acts::BinUtility(1, -cyBounds->get(Acts::CylinderVolumeBounds::eHalfLengthZ),
                   cyBounds->get(Acts::CylinderVolumeBounds::eHalfLengthZ),
                   Acts::open, Acts::binZ);
    return bUtility;
  }
  if (cutcylBounds != nullptr) {
    bUtility +=
        Acts::BinUtility(1, cutcylBounds->get(Acts::CutoutCylinderVolumeBounds::eMinR),
                   cutcylBounds->get(Acts::CutoutCylinderVolumeBounds::eMaxR),
                   Acts::open, Acts::binR);
    bUtility += Acts::BinUtility(1, -M_PI, M_PI, Acts::closed, Acts::binPhi);
    bUtility += Acts::BinUtility(
        1, -cutcylBounds->get(Acts::CutoutCylinderVolumeBounds::eHalfLengthZ),
        cutcylBounds->get(Acts::CutoutCylinderVolumeBounds::eHalfLengthZ), Acts::open,
        Acts::binZ);
    return bUtility;
  } else if (cuBounds != nullptr) {
    bUtility += Acts::BinUtility(1, -cuBounds->get(Acts::CuboidVolumeBounds::eHalfLengthX),
                           cuBounds->get(Acts::CuboidVolumeBounds::eHalfLengthX),
                           Acts::open, Acts::binX);
    bUtility += Acts::BinUtility(1, -cuBounds->get(Acts::CuboidVolumeBounds::eHalfLengthY),
                           cuBounds->get(Acts::CuboidVolumeBounds::eHalfLengthY),
                           Acts::open, Acts::binY);
    bUtility += Acts::BinUtility(1, -cuBounds->get(Acts::CuboidVolumeBounds::eHalfLengthZ),
                           cuBounds->get(Acts::CuboidVolumeBounds::eHalfLengthZ),
                           Acts::open, Acts::binZ);
    return bUtility;
  }
  return bUtility;
}

}  // namespace

namespace Acts{

void to_json(json& j, const surfaceMaterialPointer& material) {
  // A bin utility needs to be written
  const Acts::BinUtility* bUtility = nullptr;
  // Check if we have a proto material
  auto psMaterial = dynamic_cast<const Acts::ProtoSurfaceMaterial*>(material);
  if (psMaterial != nullptr) {
    // Type is proto material
    j[jsonKey().typekey] = "proto";
    // by default the protoMaterial is not used for mapping
    j[jsonKey().mapkey] = false;
    // write the bin utility
    bUtility = &(psMaterial->binUtility());
    encodeBinUtility(j, bUtility);
    return;
  }
  // Now check if we have a homogeneous material
  auto hsMaterial =
      dynamic_cast<const Acts::HomogeneousSurfaceMaterial*>(material);
  if (hsMaterial != nullptr) {
    // type is homogeneous
    j[jsonKey().typekey] = "homogeneous";
    j[jsonKey().mapkey] = true;
    j[jsonKey().datakey] = json::array({
        json::array({
            encodeMaterialSlab(hsMaterial->materialSlab(0, 0)),
        }),
    });
    // write the bin utility
    encodeBinUtility(j, bUtility);
    return; 
  } 
  // Only option remaining: BinnedSurface material
  auto bsMaterial =
      dynamic_cast<const Acts::BinnedSurfaceMaterial*>(material);
  if (bsMaterial != nullptr) {
    // type is binned
    j[jsonKey().typekey] = "binned";
    j[jsonKey().mapkey] = true;
    bUtility = &(bsMaterial->binUtility());
    // convert the data
    // get the material matrix
    json mmat = json::array();
    for (const auto& mpVector : bsMaterial->fullMaterial()) {
      json mvec = json::array();
      for (const auto& mp : mpVector) {
        mvec.push_back(encodeMaterialSlab(mp));
      }
      mmat.push_back(std::move(mvec));
    }
    j[jsonKey().datakey] = std::move(mmat);
    // write the bin utility
    encodeBinUtility(j, bUtility);
    return;     
  }
  // No material the json object is left empty.
  return;
}

void from_json(const json& j, surfaceMaterialPointer& material) {
  // By default no material is return.
  material = nullptr;
  // The bin utility for deescribing the data is retrive.
  Acts::BinUtility bUtility;
  for (auto& [key, value] : j.items()) {
    if (key == jsonKey().transfokeys and not value.empty()) {
      bUtility = Acts::BinUtility(decodeTransform(value));
      break;
    }
  }
  // Convert the material
  Acts::MaterialSlabMatrix mpMatrix;
  // Structured binding
  for (auto& [key, value] : j.items()) {
    // Check json keys
    if (key == jsonKey().bin0key and not value.empty()) {
      bUtility += decodeBinUtility(value);
    } else if (key == jsonKey().bin1key and not value.empty()) {
      bUtility += decodeBinUtility(value);
    }
    if (key == jsonKey().datakey and not value.empty()) {
      mpMatrix = decodeMaterialMatrix(value);
    }
  }

  // Return the appropriate typr of material
  if (mpMatrix.empty()) {
    material = new Acts::ProtoSurfaceMaterial(bUtility);
  } else if (bUtility.bins() == 1) {
    material = new Acts::HomogeneousSurfaceMaterial(mpMatrix[0][0]);
  } else {
    material = new Acts::BinnedSurfaceMaterial (bUtility, mpMatrix);
  }
}

void to_json(json& j, const volumeMaterialPointer& material) {
  // A bin utility needs to be written
  const Acts::BinUtility* bUtility = nullptr;
  // Check if we have a proto material
  auto pvMaterial = dynamic_cast<const Acts::ProtoVolumeMaterial*>(material);
  if (pvMaterial != nullptr) {
    // Type is proto material
    j[jsonKey().typekey] = "proto";
    // By default the protoMaterial is not used for mapping
    j[jsonKey().mapkey] = false;
    bUtility = &(pvMaterial->binUtility());
    // Write the bin utility
    encodeBinUtility(j, bUtility);
    return;  
  }
  // Now check if we have a homogeneous material
  auto hvMaterial =
      dynamic_cast<const Acts::HomogeneousVolumeMaterial*>(material);
  if (hvMaterial != nullptr) {
    // type is homogeneous
    j[jsonKey().typekey] = "homogeneous";
    j[jsonKey().mapkey] = true;
    // array of encoded materials w/ one entry
    j[jsonKey().datakey] = json::array({
        encodeMaterial(hvMaterial->material({0, 0, 0})),
    });
    // Write the bin utility
    encodeBinUtility(j, bUtility);
    return;    
  }
  // Only option remaining: material map
  auto bvMaterial2D = dynamic_cast<
      const Acts::InterpolatedMaterialMap<Acts::MaterialMapper<Acts::MaterialGrid2D>>*>(
      material);
      // Now check if we have a 2D map
  if (bvMaterial2D != nullptr) {
    // type is binned
    j[jsonKey().typekey] = "interpolated2D";
    j[jsonKey().mapkey] = true;
    bUtility = &(bvMaterial2D->binUtility());
    // convert the data
    json mmat = json::array();
    Acts::MaterialGrid2D grid = bvMaterial2D->getMapper().getGrid();
    for (size_t bin = 0; bin < grid.size(); bin++) {
      mmat.push_back(encodeMaterial(grid.at(bin)));
    }
    j[jsonKey().datakey] = std::move(mmat);
    // Write the bin utility
    encodeBinUtility(j, bUtility);
    return;  
  } 
  // Only option remaining: material map
  auto bvMaterial3D = dynamic_cast<const Acts::InterpolatedMaterialMap<
      Acts::MaterialMapper<Acts::MaterialGrid3D>>*>(material);
  // Now check if we have a 3D map
  if (bvMaterial3D != nullptr) {
    // type is binned
    j[jsonKey().typekey] = "interpolated3D";
    j[jsonKey().mapkey] = true;
    bUtility = &(bvMaterial3D->binUtility());
    // convert the data
    json mmat = json::array();
    Acts::MaterialGrid3D grid = bvMaterial3D->getMapper().getGrid();
    for (size_t bin = 0; bin < grid.size(); bin++) {
      mmat.push_back(encodeMaterial(grid.at(bin)));
    }
    j[jsonKey().datakey] = std::move(mmat);
    // Write the bin utility
    encodeBinUtility(j, bUtility);
    return;    
  }
}

void from_json(const json& j, volumeMaterialPointer& material) {
  material = nullptr;
  // The bin utility for deescribing the data
  Acts::BinUtility bUtility;
  for (auto& [key, value] : j.items()) {
    if (key == jsonKey().transfokeys and not value.empty()) {
      bUtility = Acts::BinUtility(decodeTransform(value));
      break;
    }
  }
  // Convert the material
  std::vector<Acts::Material> mmat;
  // Structured binding
  for (auto& [key, value] : j.items()) {
    // Check json keys
    if (key == jsonKey().bin0key and not value.empty()) {
      bUtility += decodeBinUtility(value);
    } else if (key == jsonKey().bin1key and not value.empty()) {
      bUtility += decodeBinUtility(value);
    } else if (key == jsonKey().bin2key and not value.empty()) {
      bUtility += decodeBinUtility(value);
    }
    if (key == jsonKey().datakey and not value.empty()) {
      for (const auto& bin : value) {
        mmat.push_back(decodeMaterial(bin));
      }
    }
  }

  // We have protoMaterial
  if (mmat.empty()) {
    material = new Acts::ProtoVolumeMaterial(bUtility);
    return;
  } 
  if (mmat.size() == 1) {
    material = new Acts::HomogeneousVolumeMaterial(mmat[0]);
    return;
  } 
  if (bUtility.dimensions() == 2) {
    std::function<Acts::Vector2(Acts::Vector3)> transfoGlobalToLocal;
    Acts::Grid2D grid = createGrid2D(bUtility, transfoGlobalToLocal);

    Acts::Grid2D::point_t min = grid.minPosition();
    Acts::Grid2D::point_t max = grid.maxPosition();
    Acts::Grid2D::index_t nBins = grid.numLocalBins();

    Acts::EAxis axis1(min[0], max[0], nBins[0]);
    Acts::EAxis axis2(min[1], max[1], nBins[1]);

    // Build the grid and fill it with data
    Acts::MaterialGrid2D mGrid(std::make_tuple(axis1, axis2));

    for (size_t bin = 0; bin < mmat.size(); bin++) {
      mGrid.at(bin) = mmat[bin].parameters();
    }
    Acts::MaterialMapper<Acts::MaterialGrid2D> matMap(transfoGlobalToLocal, mGrid);
    material =
        new Acts::InterpolatedMaterialMap<Acts::MaterialMapper<Acts::MaterialGrid2D>>(
            std::move(matMap), bUtility);
    return;
  }
  if (bUtility.dimensions() == 3) {
    std::function<Acts::Vector3(Acts::Vector3)> transfoGlobalToLocal;
    Acts::Grid3D grid = createGrid3D(bUtility, transfoGlobalToLocal);

    Acts::Grid3D::point_t min = grid.minPosition();
    Acts::Grid3D::point_t max = grid.maxPosition();
    Acts::Grid3D::index_t nBins = grid.numLocalBins();

    Acts::EAxis axis1(min[0], max[0], nBins[0]);
    Acts::EAxis axis2(min[1], max[1], nBins[1]);
    Acts::EAxis axis3(min[2], max[2], nBins[2]);

    // Build the grid and fill it with data
    Acts::MaterialGrid3D mGrid(std::make_tuple(axis1, axis2, axis3));

    for (size_t bin = 0; bin < mmat.size(); bin++) {
      mGrid.at(bin) = mmat[bin].parameters();
    }
    Acts::MaterialMapper<Acts::MaterialGrid3D> matMap(transfoGlobalToLocal, mGrid);
    material =
        new Acts::InterpolatedMaterialMap<Acts::MaterialMapper<Acts::MaterialGrid3D>>(
            std::move(matMap), bUtility);
    return;
  }
}

void to_json(json& j, const surfacePointer& surface) {
  // Get the ID of the surface (redundant but help readability)
  std::ostringstream SurfaceID;
  SurfaceID << surface->geometryId();
  j[jsonKey().geometryidkey] = SurfaceID.str();
  if(surface->surfaceMaterial() == nullptr){
    // Cast the surface bound to both disk and cylinder
    const Acts::SurfaceBounds& surfaceBounds = surface->bounds();
    auto sTransform = surface->transform(Acts::GeometryContext());
    j[jsonKey().typekey] = "proto";
    // by default the protoMaterial is not used for mapping
    j[jsonKey().mapkey] = false;
    const Acts::BinUtility bUtility = defaultSurfaceBin(surface);
    encodeBinUtility(j, &bUtility);

    const Acts::RadialBounds* radialBounds =
        dynamic_cast<const Acts::RadialBounds*>(&surfaceBounds);
    if (radialBounds != nullptr) {
      j[jsonKey().surfacetypekey] = "Disk";
      j[jsonKey().surfacepositionkey] = sTransform.translation().z();
      j[jsonKey().surfacerangekey] = {radialBounds->rMin(), radialBounds->rMax()};
      return;
    }
    const Acts::CylinderBounds* cylinderBounds =
      dynamic_cast<const Acts::CylinderBounds*>(&surfaceBounds);
    if (cylinderBounds != nullptr) {
      j[jsonKey().surfacetypekey] = "Cylinder";
      j[jsonKey().surfacepositionkey] = cylinderBounds->get(Acts::CylinderBounds::eR);
      j[jsonKey().surfacerangekey] = {
          -1 * cylinderBounds->get(Acts::CylinderBounds::eHalfLengthZ),
          cylinderBounds->get(Acts::CylinderBounds::eHalfLengthZ)};
      return;
    }
    const Acts::AnnulusBounds* annulusBounds =
      dynamic_cast<const Acts::AnnulusBounds*>(&surfaceBounds);
    if (annulusBounds != nullptr) {
      j[jsonKey().surfacetypekey] = "Annulus";
      j[jsonKey().surfacepositionkey] = sTransform.translation().z();
      j[jsonKey().surfacerangekey] = {
          {annulusBounds->rMin(), annulusBounds->rMax()},
          {annulusBounds->phiMin(), annulusBounds->phiMax()}};
      return;
    }
  }
  else{
    to_json(j, surface->surfaceMaterial());
  }
}

// Not intended for use, here for complitness sake
// to_json on surface is used as input for the material mapping
void from_json(const json& /*j*/, surfacePointer& /*surface*/) {
  return;
}

void to_json(json& j, const volumePointer& volume) {
  std::ostringstream volumeID;
  volumeID << volume->geometryId();
  j[jsonKey().geometryidkey] = volumeID.str();
  j[jsonKey().namekey]= volume->volumeName();
  if(volume->volumeMaterial() == nullptr){
    j[jsonKey().typekey] = "proto";
    // by default the protoMaterial is not used for mapping
    j[jsonKey().mapkey] = false;
    const Acts::BinUtility bUtility = defaultVolumeBin(volume);
    encodeBinUtility(j, &bUtility);
  }
  else{
    to_json(j, volume->volumeMaterial());
  }
  return;
}

// Not intended for use, here for complitness sake
// to_json on volume is used as input for the material mapping
void from_json(const json& /*j*/, volumePointer& /*volume*/) {
  return;
}  

}

Acts::JsonGeometryConverter::JsonGeometryConverter(
    const Acts::JsonGeometryConverter::Config& cfg)
    : m_cfg(std::move(cfg)),
      m_volumeMaterialConverter(m_volumeName),
      m_volumeConverter(m_volumeName),
      m_surfaceMaterialConverter(m_surfaceName),
      m_surfaceConverter(m_surfaceName){

  // Validate the configuration
  if (!m_cfg.logger) {
    throw std::invalid_argument("Missing logger");
  }
}

/// Convert method
///
json Acts::JsonGeometryConverter::materialMapsToJson(
    const DetectorMaterialMaps& maps) {

  VolumeMaterialMap volumeMap = maps.second;
  std::vector<std::pair<GeometryIdentifier, const IVolumeMaterial*>> mapVolumeInit;
  for (auto it=volumeMap.begin(); it!=volumeMap.end(); it++){
    mapVolumeInit.push_back({it->first,it->second.get()});
  }

  GeometryHierarchyMap<const IVolumeMaterial*> HierarchyVolumeMap(
        mapVolumeInit);
  json materialVolume = m_volumeMaterialConverter.toJson(HierarchyVolumeMap);

  SurfaceMaterialMap surfaceMap = maps.first;
  std::vector<std::pair<GeometryIdentifier, const ISurfaceMaterial*>> mapSurfaceInit;
  for (auto it=surfaceMap.begin(); it!=surfaceMap.end(); it++){
    mapSurfaceInit.push_back({it->first,it->second.get()});
  }

  GeometryHierarchyMap<const ISurfaceMaterial*>
      HierarchySurfaceMap(mapSurfaceInit);

  json materialSurface = m_surfaceMaterialConverter.toJson(HierarchySurfaceMap);

  json materialMap;
  materialMap["Volumes"] = materialVolume;
  materialMap["Surfaces"] = materialSurface;
  return materialMap;
}

Acts::JsonGeometryConverter::DetectorMaterialMaps Acts::JsonGeometryConverter::jsonToMaterialMaps(const json& materialmap) {

  json materialVolume = materialmap["Volumes"];
  GeometryHierarchyMap<const IVolumeMaterial*>
      HierarchyVolumeMap = m_volumeMaterialConverter.fromJson(materialVolume);

  VolumeMaterialMap volumeMap;
  for (size_t i = 0; i < HierarchyVolumeMap.size(); i++) {
    std::shared_ptr<const IVolumeMaterial> volumePointer(HierarchyVolumeMap.valueAt(i));
    volumeMap.insert(
        {HierarchyVolumeMap.idAt(i), std::move(volumePointer)});
  }

    json materialSurface = materialmap["Surfaces"];
  GeometryHierarchyMap<const ISurfaceMaterial*>
      HierarchySurfaceMap = m_surfaceMaterialConverter.fromJson(materialSurface);
  
  SurfaceMaterialMap surfaceMap;
  for (size_t i = 0; i < HierarchySurfaceMap.size(); i++) {
    std::shared_ptr<const ISurfaceMaterial> volumePointer(HierarchySurfaceMap.valueAt(i));
    surfaceMap.insert({HierarchySurfaceMap.idAt(i), std::move(volumePointer)});
  }

  Acts::JsonGeometryConverter::DetectorMaterialMaps maps = {surfaceMap, volumeMap};

  // Return the filled maps
  return maps;
}

json Acts::JsonGeometryConverter::trackingGeometryToJson(
    const Acts::TrackingGeometry& tGeometry) {

  std::vector <
      std::pair<GeometryIdentifier, const TrackingVolume*>> volumeHierarchy;
  std::vector < std::pair<GeometryIdentifier, const Surface*>> surfaceHierarchy;
  convertToHierarchy(volumeHierarchy, surfaceHierarchy,
                     tGeometry.highestTrackingVolume());
  GeometryHierarchyMap<const TrackingVolume*> HierarchyVolumeMap(volumeHierarchy);
  json jsonVolumes = m_volumeConverter.toJson(HierarchyVolumeMap);
  GeometryHierarchyMap<const Surface*> HierarchySurfaceMap(surfaceHierarchy);
  json jsonSurfaces = m_surfaceConverter.toJson(HierarchySurfaceMap);
  json hierarchyMap;
  hierarchyMap["1.Volumes"] = jsonVolumes;
  hierarchyMap["2.Surfaces"] = jsonSurfaces;
  return hierarchyMap;
}

void Acts::JsonGeometryConverter::convertToHierarchy(
    std::vector <std::pair<GeometryIdentifier, const TrackingVolume*>> &volumeHierarchy,
    std::vector <std::pair<GeometryIdentifier, const Surface*>> &surfaceHierarchy,
    const Acts::TrackingVolume* tVolume){

  if ( (tVolume->volumeMaterial() != nullptr || m_cfg.processNonMaterial == true) && m_cfg.processVolumes == true) {
    volumeHierarchy.push_back({tVolume->geometryId(), tVolume});
  }
  // there are confined volumes
  if (tVolume->confinedVolumes() != nullptr) {
    // get through the volumes
    auto& volumes = tVolume->confinedVolumes()->arrayObjects();
    // loop over the volumes
    for (auto& vol : volumes) {
      // recursive call
      convertToHierarchy(volumeHierarchy, surfaceHierarchy, vol.get());
    }
  }
  // there are dense volumes
  if (m_cfg.processDenseVolumes && !tVolume->denseVolumes().empty()) {
    // loop over the volumes
    for (auto& vol : tVolume->denseVolumes()) {
      // recursive call
      convertToHierarchy(volumeHierarchy, surfaceHierarchy, vol.get());
    }
  }
  if (tVolume->confinedLayers() != nullptr) {
    // get the layers
    auto& layers = tVolume->confinedLayers()->arrayObjects();
    // loop over the layers
    for (auto& lay : layers) {
      if(lay->geometryId() == m_cfg.processRepresenting == true){
        auto& layRep = lay->surfaceRepresentation();
        if ( (layRep.surfaceMaterial() != nullptr || m_cfg.processNonMaterial == true) && layRep.geometryId() != GeometryIdentifier()){
          surfaceHierarchy.push_back({layRep.geometryId(), &layRep});
        }
      }
      if (lay->approachDescriptor() != nullptr && m_cfg.processApproaches == true) {
        for (auto& asf : lay->approachDescriptor()->containedSurfaces()) {
          if (asf->surfaceMaterial() != nullptr || m_cfg.processNonMaterial == true) {  
            surfaceHierarchy.push_back({asf->geometryId(), asf});
          }
        }
      }
      if (lay->surfaceArray() != nullptr && m_cfg.processSensitives == true){
       for (auto& ssf : lay->surfaceArray()->surfaces()) {
          if (ssf->surfaceMaterial() != nullptr || m_cfg.processNonMaterial == true) {   
            surfaceHierarchy.push_back({ssf->geometryId(), ssf});
          }
        }
      }
    }
  }
  // Let's finally check the boundaries
  for (auto& bsurf : tVolume->boundarySurfaces()) {
    // the surface representation
    auto& bssfRep = bsurf->surfaceRepresentation();
    if (bssfRep.geometryId().volume() == tVolume->geometryId().volume() && m_cfg.processBoundaries == true) {
      if (bssfRep.surfaceMaterial() != nullptr || m_cfg.processNonMaterial == true) {    
        surfaceHierarchy.push_back({bssfRep.geometryId(), &bssfRep});
      }
    }
  }

}
