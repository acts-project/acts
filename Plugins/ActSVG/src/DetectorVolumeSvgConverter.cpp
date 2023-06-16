// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/ActSVG/DetectorVolumeSvgConverter.hpp"

#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Detector/Portal.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Plugins/ActSVG/SurfaceGridSvgConverter.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <set>

Acts::Svg::ProtoVolume Acts::Svg::DetectorVolumeConverter::convert(
    const GeometryContext& gctx, const Experimental::DetectorVolume& dVolume,
    const DetectorVolumeConverter::Options& volumeOptions) {
  // The local logger
  ACTS_LOCAL_LOGGER(
      getDefaultLogger("DetectorVolumeSvgConverter", volumeOptions.logLevel));

  ProtoVolume pVolume;
  pVolume._name = dVolume.name();

  auto volumeTransform = dVolume.transform(gctx);

  // The detector volume is of cylindrical shape
  const auto& boundValues = dVolume.volumeBounds().values();
  if (dVolume.volumeBounds().type() == Acts::VolumeBounds::eCylinder) {
    // we keep 6 for the moment
    for (unsigned int ib = 0; ib < 2u; ++ib) {
      pVolume._bound_values.push_back(
          static_cast<actsvg::scalar>(boundValues[ib]));
    }
    pVolume._bound_values.push_back(
        static_cast<actsvg::scalar>(volumeTransform.translation().z()));
    for (unsigned int ib = 2u; ib < 5u; ++ib) {
      pVolume._bound_values.push_back(
          static_cast<actsvg::scalar>(boundValues[ib]));
    }
  }

  // Adding the portals to the lot
  for (auto [ip, p] : enumerate(dVolume.portals())) {
    auto pPortal =
        PortalConverter::convert(gctx, *p, volumeOptions.portalOptions);
    auto pPortalCandidate = volumeOptions.portalIndices.find(p);
    if (pPortalCandidate != volumeOptions.portalIndices.end()) {
      pPortal._name = "portal_" + std::to_string(pPortalCandidate->second);
    } else {
      pPortal._name = dVolume.name() + "_portal_" + std::to_string(ip);
    }
    pVolume._portals.push_back(pPortal);
  }

  SurfaceGridConverter::Options sgOptions;
  auto [surfaces, grid, associations] =
      SurfaceGridConverter::convert(gctx, dVolume, sgOptions);
  pVolume._surfaces = surfaces;
  pVolume._surface_grid = grid;
  pVolume._grid_associations = associations;

  // Adding the surfaces
  std::vector<ProtoSurface> pSurfaces;
  for (auto [is, s] : enumerate(dVolume.surfaces())) {
    auto pSurface =
        SurfaceConverter::convert(gctx, *s, volumeOptions.surfaceOptions);
    pSurface._name = dVolume.name() + "_surface_" + std::to_string(is);
    pSurfaces.push_back(pSurface);
  }
  pVolume._v_surfaces = pSurfaces;

  return pVolume;
}

std::array<actsvg::svg::object, 2u> Acts::Svg::View::layer(const ProtoVolume& volume) {
  // The sheets
  actsvg::svg::object module_sheet;
  actsvg::svg::object grid_sheet;

  if (volume._surface_grid._type == ProtoGrid::e_z_phi) {
    module_sheet = actsvg::display::barrel_sheet(
        volume._name + "_modules", volume, {800, 800},
        actsvg::display::e_module_info);
    grid_sheet =
        actsvg::display::barrel_sheet(volume._name + "_grid", volume,
                                      {800, 800}, actsvg::display::e_grid_info);

  } else if (volume._surface_grid._type == ProtoGrid::e_r_phi) {
    module_sheet = actsvg::display::endcap_sheet(
        volume._name + "_modules", volume, {800, 800},
        actsvg::display::e_module_info);
    grid_sheet =
        actsvg::display::endcap_sheet(volume._name + "_grid", volume,
                                      {800, 800}, actsvg::display::e_grid_info);
  }

  return {module_sheet, grid_sheet};
}
