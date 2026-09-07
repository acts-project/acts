// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

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

  /// Get number of bins
  /// @return Number of bins
  std::uint32_t numOfBins() const {
    return static_cast<std::uint32_t>(m_bins.size());
  }

  /// Get bins
  /// @return Vector of bin indices
  const std::vector<std::uint32_t>& bins() const { return m_bins; }

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

  std::vector<std::uint32_t> m_bins;
  std::vector<float> m_minRadius;
  std::vector<float> m_maxRadius;
  std::vector<float> m_minBinCoord;
  std::vector<float> m_maxBinCoord;

  float m_minEta{};
  float m_maxEta{};
  /// Width of one eta bin.
  float m_etaBin{};
  float m_r1{};
  float m_z1{};
  float m_r2{};
  float m_z2{};
  std::uint32_t m_nBins{};
};

}  // namespace Acts::Experimental::detail
