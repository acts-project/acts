// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/GeoModel/GeoModelModuleSplitter.hpp"

#include "Acts/Surfaces/AnnulusBounds.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <iostream>

namespace Acts {

std::vector<std::shared_ptr<GeoModelDetectorElement>>
GeoModelModuleSplitter::split(
    std::shared_ptr<GeoModelDetectorElement> detElement,
    const Acts::GeometryContext &gctx) const {
  const auto &logger = *m_logger;

  const auto surface = detElement->surface().getSharedPtr();
  const auto annulusBounds =
      dynamic_cast<const AnnulusBounds *>(&surface->bounds());

  if (annulusBounds == nullptr) {
    throw std::runtime_error("Surfaces does not has annulus bounds");
  }

  std::vector<std::shared_ptr<GeoModelDetectorElement>> result;

  for (const auto& [patternName, radii] : m_splitPatterns) {
    if ((std::abs(radii.front() - annulusBounds->rMin()) > m_tolerance) ||
        (std::abs(radii.back() - annulusBounds->rMax()) > m_tolerance)) {
      ACTS_VERBOSE("Skip pattern '" << patternName << "' for element '"
                                    << detElement->databaseEntryName());
      continue;
    }

    radii.front() = annulusBounds->rMin();
    radii.back() = annulusBounds->rMax();

    ACTS_DEBUG("Accept pattern '" << patternName << "' for element '"
                                  << detElement->databaseEntryName() << "'");

    result.reserve(radii.size() - 1);
    for (std::size_t i = 0ul; i < radii.size() - 1; ++i) {
      ACTS_VERBOSE("Make new annulus bounds: " << [&]() {
        std::stringstream ss;
        for (auto v : annulusBounds->values()) {
          ss << v << " ";
        }
        return ss.str();
      }());

      auto bounds = std::make_shared<AnnulusBounds>(
          radii[i], radii[i + 1], annulusBounds->get(AnnulusBounds::eMinPhiRel),
          annulusBounds->get(AnnulusBounds::eMaxPhiRel),
          Vector2{annulusBounds->get(AnnulusBounds::eOriginX),
                  annulusBounds->get(AnnulusBounds::eOriginY)},
          annulusBounds->get(AnnulusBounds::eAveragePhi));

      auto newDetElement =
          GeoModelDetectorElement::createDetectorElement<DiscSurface>(
              detElement->physicalVolume(), bounds, detElement->transform(gctx),
              detElement->thickness());
      newDetElement->setDatabaseEntryName(detElement->databaseEntryName());

      result.push_back(newDetElement);
    }

    return result;
  }

  throw std::runtime_error([&]() {
    std::stringstream ss;
    ss << "Could not split '" << detElement->databaseEntryName()
       << "' (rmin: " << annulusBounds->rMin()
       << ", rmax: " << annulusBounds->rMax() << ")";
    return ss.str();
  }());
}

}  // namespace Acts
