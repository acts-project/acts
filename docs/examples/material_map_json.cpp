// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "ActsPlugins/Json/TrackingGeometryMaterialJsonConverter.hpp"

/// Read material, apply it to an existing geometry, and write it as JSON.
void exampleMaterialMapJson(Acts::TrackingGeometry& geometry) {
  //! [Read and write material map]
  Acts::TrackingGeometryMaterialJsonConverter converter;
  auto material = converter.fromFile("material.cbor.zst");
  material.apply(geometry);
  converter.toFile(material, "material.json");
  //! [Read and write material map]
}
