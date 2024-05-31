// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/GeoModel/detail/GeoModelExtentHelper.hpp"

#include "Acts/Plugins/GeoModel/detail/GeoModelDbHelper.hpp"
#include "Acts/Utilities/BinningData.hpp"
#include "Acts/Utilities/Enumerate.hpp"

std::vector<Acts::BinningValue>
Acts::detail::GeoModelExentHelper::readConstaints(
    const std::vector<std::string>& boundsEntry, std::size_t boundsPos,
    const std::string& ctype) {
  // Check the bounds entry
  if (boundsEntry.size() < 2u) {
    throw std::invalid_argument(
        "GeoModelBlueprintCreater: Bounds entry has to have at least 2 "
        "entries (type, values)");
  }
  std::set<Acts::BinningValue> constraints;
  // Switch on the bounds type
  if (boundsEntry[boundsPos] == "cyl") {
    // Capture the values
    std::vector<std::string> valuesEntry =
        GeoModelDbHelper::splitString(boundsEntry[boundsPos + 1u], ",");
    if (valuesEntry.size() < 4u) {
      throw std::invalid_argument(
          "GeoModelBlueprintCreater: Cylinder bounds entry has to have at "
          "least 4 entries (rmin, rmax, zmin, zmax)");
    }
    // Raw database values to extent entries
    constexpr std::array<BinningValue, 6u> bvCyl = {binR, binR,   binZ,
                                                    binZ, binPhi, binPhi};

    for (auto [iv, value] : enumerate(valuesEntry)) {
      if (value == ctype || value[0u] == ctype[0u]) {
        constraints.insert(bvCyl.at(iv));
      }
    }
  }
  return {constraints.begin(), constraints.end()};
}

std::tuple<Acts::VolumeBounds::BoundsType, Acts::Extent>
Acts::detail::GeoModelExentHelper::extentFromTable(
    const std::vector<std::string>& boundsEntry, std::size_t boundsPos,
    const Acts::Extent& externalExtent, const Acts::Extent& internalExtent) {
  // Check the bounds entry
  if (boundsEntry.size() < 2u) {
    throw std::invalid_argument(
        "GeoModelBlueprintCreater: Bounds entry has to have at least 2 "
        "entries (type, values)");
  }

  // Start with the mother extent
  // -> and shrink it to size
  VolumeBounds::BoundsType boundsType = VolumeBounds::BoundsType::eOther;
  Extent extent;
  // Switch on the bounds type
  if (boundsEntry[boundsPos] == "cyl") {
    // Set the bounds type
    boundsType = VolumeBounds::BoundsType::eCylinder;
    // Capture the values
    std::vector<std::string> valuesEntry =
        GeoModelDbHelper::splitString(boundsEntry[boundsPos + 1u], ",");
    if (valuesEntry.size() < 4u) {
      throw std::invalid_argument(
          "GeoModelBlueprintCreater: Cylinder bounds entry has to have at "
          "least 4 entries (rmin, rmax, zmin, zmax)");
    }
    // Raw database values to extent entries
    constexpr std::array<BinningValue, 6u> bvCyl = {binR, binR,   binZ,
                                                    binZ, binPhi, binPhi};
    for (auto [iv, value] : enumerate(valuesEntry)) {
      // Get the binning value
      BinningValue bValue = bvCyl.at(iv);
      ActsScalar val = std::numeric_limits<ActsScalar>::max();
      bool isMin = (iv % 2 == 0);
      // Case "e" : exxternal extent
      if (value == "e") {
        // External parameters do not constrain it
        if (!externalExtent.constrains(bValue)) {
          throw std::invalid_argument(
              "GeoModelExtentHelper: External extent does not constrain. ");
        }
        val = isMin ? externalExtent.min(bValue) : externalExtent.max(bValue);
      } else if (value == "i" || value[0u] == 'i') {
        // Add the envelope
        ActsScalar envelope = 0.;
        if (value.size() > 2u) {
          std::vector<std::string> valEntry =
              GeoModelDbHelper::splitString(value, "+");
          envelope = std::stod(valEntry[1]);
        }

        // Internals do not constrain it
        if (!internalExtent.constrains(bValue)) {
          throw std::invalid_argument(
              "GeoModelExtentHelper: Internals do not constrain.");
        }
        val = isMin ? internalExtent.min(bValue) - envelope
                    : internalExtent.max(bValue) + envelope;
        // Case "i" : Inherited from mother
      } else {
        val = std::stod(value);
      }
      // Case value is a number -> shrink mother to it or set it
      if (isMin) {
        extent.setMin(bValue, val);
      } else {
        extent.setMax(bValue, val);
      }
    }
  }
  return std::make_tuple(boundsType, extent);
}
