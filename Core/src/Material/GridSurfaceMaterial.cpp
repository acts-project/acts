// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Material/GridSurfaceMaterial.hpp"

#include <ostream>
#include <stdexcept>
#include <type_traits>

namespace Acts {

namespace {

/// @brief Resolve the number of entries held by a storage alternative
std::size_t storageSize(const GridSurfaceMaterial::Storage& storage) {
  return std::visit(
      [](const auto& s) -> std::size_t {
        using T = std::decay_t<decltype(s)>;
        if constexpr (std::is_same_v<T, GridSurfaceMaterial::Direct>) {
          return s.size();
        } else {
          return s.indices.size();
        }
      },
      storage);
}

/// @brief Look up the material slab for a given global bin index
struct SlabAt {
  std::size_t bin;

  const MaterialSlab& operator()(const GridSurfaceMaterial::Direct& s) const {
    return s.at(bin);
  }
  const MaterialSlab& operator()(const GridSurfaceMaterial::Indexed& s) const {
    return s.material.at(s.indices.at(bin));
  }
  const MaterialSlab& operator()(
      const GridSurfaceMaterial::GloballyIndexed& s) const {
    return s.material->at(s.indices.at(bin));
  }
};

}  // namespace

GridSurfaceMaterial::GridSurfaceMaterial(MultiAxisSpec2D binning,
                                         Storage storage, double splitFactor,
                                         MappingType mappingType)
    : ISurfaceMaterial(splitFactor, mappingType),
      m_binning(std::move(binning)),
      m_multiAxis(m_binning.buildMultiAxis()),
      m_storage(std::move(storage)) {
  std::size_t nBins = m_multiAxis->getNTotalBins(true);
  if (storageSize(m_storage) != nBins) {
    throw std::invalid_argument(
        "GridSurfaceMaterial: storage size does not match the number of "
        "bins (including under-/overflow) implied by the binning.");
  }
}

const MaterialSlab& GridSurfaceMaterial::materialSlab(const Vector2& lp) const {
  std::size_t bin = m_multiAxis->getGlobalBinFromPoint({lp[0], lp[1]});
  return std::visit(SlabAt{bin}, m_storage);
}

const MaterialSlab& GridSurfaceMaterial::materialSlab(
    const Vector3& /*gp*/) const {
  throw std::logic_error(
      "GridSurfaceMaterial: global (Vector3) material lookup is not "
      "supported, use materialSlab(const Vector2&) instead.");
}

ISurfaceMaterial& GridSurfaceMaterial::scale(double factor) {
  std::visit(
      [factor](auto& s) {
        using T = std::decay_t<decltype(s)>;
        if constexpr (std::is_same_v<T, Direct>) {
          for (auto& msl : s) {
            msl.scaleThickness(static_cast<float>(factor));
          }
        } else if constexpr (std::is_same_v<T, Indexed>) {
          for (auto& msl : s.material) {
            msl.scaleThickness(static_cast<float>(factor));
          }
        } else {
          if (s.sharedEntries) {
            throw std::invalid_argument(
                "GridSurfaceMaterial: scaling a globally indexed material "
                "with shared entries is not supported.");
          }
          for (std::size_t index : s.indices) {
            (*s.material)[index].scaleThickness(static_cast<float>(factor));
          }
        }
      },
      m_storage);
  return *this;
}

std::ostream& GridSurfaceMaterial::toStream(std::ostream& sl) const {
  sl << "Acts::GridSurfaceMaterial : " << std::endl;
  sl << m_binning << std::endl;
  return sl;
}

}  // namespace Acts
