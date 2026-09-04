// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/GbtsBinning.hpp"
#include "Acts/Seeding/GbtsLayerDescription.hpp"

#include <cstdint>
#include <vector>

namespace Acts::Experimental::detail {

/// Layer helper with eta-bin access, built and owned by GbtsGeometry.
class GbtsLayer final {
 public:
  /// @param layerDescription Layer description for the layer
  /// @param etaBinWidth Eta bin width
  /// @param bin0 Starting bin index
  GbtsLayer(const GbtsLayerDescription& layerDescription, float etaBinWidth,
            std::uint32_t bin0);

  /// Get eta bin for given z and r coordinates
  /// @param zh Z coordinate
  /// @param rh Radius coordinate
  /// @return Eta bin index
  std::uint32_t getEtaBin(float zh, float rh) const;

  /// Get the eta binning the geometry gave this layer
  /// @return The layer's binning
  const GbtsLayerBinning& binning() const { return m_binning; }

  /// Get the layer description
  /// @return Reference to the layer description
  const GbtsLayerDescription& layerDescription() const {
    return m_layerDescription;
  }

  /// Verify bin compatibility
  /// @param otherLayer Other layer to compare with
  /// @param b1 First bin index
  /// @param b2 Second bin index
  /// @param minZ0 Minimum z0 coordinate
  /// @param maxZ0 Maximum z0 coordinate
  /// @return True if bins are compatible
  bool checkCompatibility(const GbtsLayer& otherLayer, std::uint32_t b1,
                          std::uint32_t b2, float minZ0, float maxZ0) const;

 private:
  GbtsLayerDescription m_layerDescription;
  GbtsLayerBinning m_binning;

  std::vector<float> m_minRadius;
  std::vector<float> m_maxRadius;
  std::vector<float> m_minBinCoord;
  std::vector<float> m_maxBinCoord;
};

}  // namespace Acts::Experimental::detail
