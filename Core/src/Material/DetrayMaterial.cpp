// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/DetrayExceptions.hpp"
#include "Acts/Material/BinnedSurfaceMaterial.hpp"
#include "Acts/Material/GridSurfaceMaterial.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/ProtoSurfaceMaterial.hpp"

#include <detray/io/frontend/payloads.hpp>

namespace Acts {
std::unique_ptr<DetraySurfaceMaterial> BinnedSurfaceMaterial::toDetrayPayload(
    const detray::io::volume_payload& /*volume*/) const {
  throw DetrayUnsupportedMaterialException("BinnedSurfaceMaterial");
}

template <>
std::unique_ptr<DetraySurfaceMaterial>
ProtoSurfaceMaterialT<Acts::BinUtility>::toDetrayPayload(
    const detray::io::volume_payload& /*volume*/) const {
  // Does not apply to detray
  return nullptr;
}

template <>
std::unique_ptr<DetraySurfaceMaterial>
ProtoSurfaceMaterialT<std::vector<DirectedProtoAxis>>::toDetrayPayload(
    const detray::io::volume_payload& /*volume*/) const {
  // Does not apply to detray
  return nullptr;
}

std::unique_ptr<DetraySurfaceMaterial>
detail::IGridSurfaceMaterialBase::toDetrayPayload(
    const detray::io::volume_payload& /*volume*/) const {
  throw DetrayUnsupportedMaterialException("detail::IGridSurfaceMaterialBase");
}

namespace {
detray::io::material_slab_payload convertMaterialSlab(
    const MaterialSlab& slab) {
  detray::io::material_slab_payload payload;
  // Fill the material parameters and the thickness
  const auto& material = slab.material();
  payload.thickness = slab.thickness();
  payload.mat = detray::io::material_payload{
      {material.X0(), material.L0(), material.Ar(), material.Z(),
       material.massDensity(), material.molarDensity(), 0.}};
  payload.type = detray::io::material_id::slab;
  return payload;
}
}  // namespace

std::unique_ptr<DetraySurfaceMaterial>
HomogeneousSurfaceMaterial::toDetrayPayload(
    const detray::io::volume_payload& /*volume*/) const {
  return std::make_unique<DetraySurfaceMaterial>(
      convertMaterialSlab(materialSlab()));
}

}  // namespace Acts
