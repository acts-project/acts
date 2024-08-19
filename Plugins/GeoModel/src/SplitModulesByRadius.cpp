// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/GeoModel/SplitModulesByRadius.hpp"

#include "Acts/Surfaces/AnnulusBounds.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <iostream>

namespace Acts {

std::optional<std::vector<std::pair<std::shared_ptr<Surface>,
                                    std::shared_ptr<GeoModelDetectorElement>>>>
ModuleByRadiusSplitter::split(
    std::shared_ptr<Surface> surface,
    std::shared_ptr<GeoModelDetectorElement> detElement,
    const Acts::GeometryContext &gctx) const {
  const auto &logger = *m_logger;

  auto annulusBounds = dynamic_cast<const AnnulusBounds *>(&surface->bounds());

  if (annulusBounds == nullptr) {
    ACTS_DEBUG("Not annulus bounds");
    return std::nullopt;
  }

  std::optional<std::vector<std::pair<
      std::shared_ptr<Surface>, std::shared_ptr<GeoModelDetectorElement>>>>
      result;

  for (auto [design, radii] : m_designs) {
    if ((std::abs(radii.front() - annulusBounds->rMin()) > m_tolerance) ||
        (std::abs(radii.back() - annulusBounds->rMax()) > m_tolerance)) {
      ACTS_VERBOSE("Skip design '" << design << "' for element '"
                                   << detElement->databaseEntryName());
      continue;
    }

    radii.front() = annulusBounds->rMin();
    radii.back() = annulusBounds->rMax();

    ACTS_DEBUG("Accept design '" << design << "' for element '"
                                 << detElement->databaseEntryName());

    result.emplace();
    const auto origValues = annulusBounds->values();

    for (auto i = 0ul; i < radii.size() - 1; ++i) {
      std::array<double, AnnulusBounds::eSize> values;
      std::copy(origValues.begin(), origValues.end(), values.begin());
      values[AnnulusBounds::eMinR] = radii[i];
      values[AnnulusBounds::eMaxR] = radii[i + 1];

      auto bounds = std::make_shared<AnnulusBounds>(values);

      auto newDetElement =
          GeoModelDetectorElement::createDetectorElement<DiscSurface>(
              detElement->physicalVolume(), bounds, detElement->transform(gctx),
              detElement->thickness());
      auto newSurface = newDetElement->surface().getSharedPtr();

      newDetElement->setDatabaseEntryName(detElement->databaseEntryName());

      result->emplace_back(newSurface, newDetElement);
    }

    return result;
  }

  ACTS_WARNING("Could not split '" << detElement->databaseEntryName()
                                   << "' (rmin: " << annulusBounds->rMin()
                                   << ", rmax: " << annulusBounds->rMax()
                                   << ")");
  return result;
}

}  // namespace Acts
