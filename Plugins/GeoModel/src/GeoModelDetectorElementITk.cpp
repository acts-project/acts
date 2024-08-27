// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/GeoModel/GeoModelDetectorElementITk.hpp"

#include "Acts/Surfaces/AnnulusBounds.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"

namespace Acts {

std::tuple<std::shared_ptr<Acts::GeoModelDetectorElementITk>,
           std::shared_ptr<Acts::Surface>>
Acts::GeoModelDetectorElementITk::convertFromGeomodel(
    std::shared_ptr<GeoModelDetectorElement> detEl,
    std::shared_ptr<Surface> srf, const GeometryContext &gctx, int hardware,
    int barrelEndcap, int layerWheel, int etaModule, int phiModule, int side) {
  auto helper = [&]<typename surface_t, typename bounds_t>() {
    auto bounds = std::make_shared<bounds_t>(
        dynamic_cast<const bounds_t &>(srf->bounds()));

    auto itkEl = std::make_shared<GeoModelDetectorElementITk>(
        detEl->physicalVolume(), nullptr, detEl->transform(gctx),
        detEl->thickness());
    auto surface = Surface::makeShared<surface_t>(bounds, *itkEl.get());

    itkEl->m_hardware = hardware;
    itkEl->m_barrelEndcap = barrelEndcap;
    itkEl->m_layerWheel = layerWheel;
    itkEl->m_etaModule = etaModule;
    itkEl->m_etaModule = phiModule;
    itkEl->m_side = side;

    itkEl->attachSurface(surface);
    itkEl->setDatabaseEntryName(detEl->databaseEntryName());
    return std::pair{itkEl, surface};
  };

  if (srf->type() == Acts::Surface::Plane &&
      srf->bounds().type() == Acts::SurfaceBounds::eRectangle) {
    return helper.operator()<Acts::PlaneSurface, Acts::RectangleBounds>();
  }
  if (srf->type() == Acts::Surface::Disc &&
      srf->bounds().type() == Acts::SurfaceBounds::eAnnulus) {
    return helper.operator()<Acts::DiscSurface, Acts::AnnulusBounds>();
  }

  throw std::runtime_error(
      "Only Plane+Rectangle and Disc+Annulus are converted for the ITk");
}

}  // namespace Acts
