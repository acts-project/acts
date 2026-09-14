// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Material/GridSurfaceMaterial.hpp"

#include "Acts/Utilities/AxisSpec.hpp"
#include "Acts/Utilities/Helpers.hpp"

#include <array>
#include <optional>
#include <ostream>
#include <stdexcept>

namespace Acts {

namespace {

/// @brief Resolve the number of entries held by a storage alternative
std::size_t storageSize(const GridSurfaceMaterial::Storage& storage) {
  return std::visit(
      overloaded{
          [](const GridSurfaceMaterial::Direct& s) { return s.size(); },
          [](const GridSurfaceMaterial::Indexed& s) {
            return s.indices.size();
          },
          [](const GridSurfaceMaterial::GloballyIndexed& s) {
            return s.indices.size();
          },
      },
      storage);
}

/// @brief Capture a pair of axes as a fully specified 2D binning spec
MultiAxisSpec2D multiAxisSpecFromAxes(const IAxis& axis0, const IAxis& axis1) {
  std::array<AxisSpec, 2> specs{AxisSpec::FromAxis(axis0),
                                AxisSpec::FromAxis(axis1)};
  return MultiAxisSpec2D(specs);
}

/// @brief Flatten a column-major (regular-bin only) 2D payload into a
/// per-global-bin vector, including under-/overflow bins
template <typename value_t>
std::vector<value_t> flattenPayload2D(
    const IMultiAxis2D& multiAxis,
    const std::vector<std::vector<value_t>>& payload) {
  std::vector<value_t> flat(multiAxis.getNTotalBins(true));
  IMultiAxis2D::LocalBins nBins = multiAxis.getNBins();
  for (std::size_t i0 = 0; i0 < nBins[0]; ++i0) {
    for (std::size_t i1 = 0; i1 < nBins[1]; ++i1) {
      IMultiAxis2D::LocalBins lbin{i0 + 1, i1 + 1};
      flat[multiAxis.getGlobalBinFromLocalBins(lbin)] = payload[i0][i1];
    }
  }
  return flat;
}

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

std::unique_ptr<GridSurfaceMaterial> GridSurfaceMaterial::createDirect(
    const IAxis& axis0, const IAxis& axis1,
    const std::vector<std::vector<MaterialSlab>>& payload) {
  MultiAxisSpec2D binning = multiAxisSpecFromAxes(axis0, axis1);
  auto multiAxis = binning.buildMultiAxis();
  Direct storage = flattenPayload2D(*multiAxis, payload);
  return std::make_unique<GridSurfaceMaterial>(std::move(binning),
                                               std::move(storage));
}

std::unique_ptr<GridSurfaceMaterial> GridSurfaceMaterial::createDirect(
    const MultiAxisSpec2D& binning, const Surface& surface,
    const std::vector<std::vector<MaterialSlab>>& payload) {
  auto axes = resolveMultiAxis(binning, surface);
  return createDirect(axes->getAxis(0), axes->getAxis(1), payload);
}

std::unique_ptr<GridSurfaceMaterial> GridSurfaceMaterial::createIndexed(
    const IAxis& axis0, const IAxis& axis1, std::vector<MaterialSlab> material,
    const std::vector<std::vector<std::size_t>>& payload) {
  MultiAxisSpec2D binning = multiAxisSpecFromAxes(axis0, axis1);
  auto multiAxis = binning.buildMultiAxis();
  Indexed storage{flattenPayload2D(*multiAxis, payload), std::move(material)};
  return std::make_unique<GridSurfaceMaterial>(std::move(binning),
                                               std::move(storage));
}

std::unique_ptr<GridSurfaceMaterial> GridSurfaceMaterial::createIndexed(
    const MultiAxisSpec2D& binning, const Surface& surface,
    std::vector<MaterialSlab> material,
    const std::vector<std::vector<std::size_t>>& payload) {
  auto axes = resolveMultiAxis(binning, surface);
  return createIndexed(axes->getAxis(0), axes->getAxis(1), std::move(material),
                       payload);
}

std::unique_ptr<GridSurfaceMaterial> GridSurfaceMaterial::createGloballyIndexed(
    const IAxis& axis0, const IAxis& axis1,
    std::shared_ptr<std::vector<MaterialSlab>> material,
    const std::vector<std::vector<std::size_t>>& payload) {
  MultiAxisSpec2D binning = multiAxisSpecFromAxes(axis0, axis1);
  auto multiAxis = binning.buildMultiAxis();
  GloballyIndexed storage{flattenPayload2D(*multiAxis, payload),
                          std::move(material)};
  return std::make_unique<GridSurfaceMaterial>(std::move(binning),
                                               std::move(storage));
}

std::unique_ptr<GridSurfaceMaterial> GridSurfaceMaterial::createGloballyIndexed(
    const MultiAxisSpec2D& binning, const Surface& surface,
    std::shared_ptr<std::vector<MaterialSlab>> material,
    const std::vector<std::vector<std::size_t>>& payload) {
  auto axes = resolveMultiAxis(binning, surface);
  return createGloballyIndexed(axes->getAxis(0), axes->getAxis(1),
                               std::move(material), payload);
}

std::vector<AxisDirection> GridSurfaceMaterial::localAxisDirections() const {
  std::optional<AxisDirection> dir0 = m_binning.axisSpec(0).direction();
  std::optional<AxisDirection> dir1 = m_binning.axisSpec(1).direction();
  if (!dir0.has_value() || !dir1.has_value()) {
    return {};
  }
  return {*dir0, *dir1};
}

const MaterialSlab& GridSurfaceMaterial::materialSlab(const Vector2& lp) const {
  std::size_t bin = m_multiAxis->getGlobalBinFromPoint({lp[0], lp[1]});
  return std::visit(
      overloaded{
          [bin](const Direct& s) -> const MaterialSlab& { return s.at(bin); },
          [bin](const Indexed& s) -> const MaterialSlab& {
            return s.material.at(s.indices.at(bin));
          },
          [bin](const GloballyIndexed& s) -> const MaterialSlab& {
            return s.material->at(s.indices.at(bin));
          },
      },
      m_storage);
}

const MaterialSlab& GridSurfaceMaterial::materialSlab(
    const Vector3& /*gp*/) const {
  throw std::logic_error(
      "GridSurfaceMaterial: global (Vector3) material lookup is not "
      "supported, use materialSlab(const Vector2&) instead.");
}

ISurfaceMaterial& GridSurfaceMaterial::scale(double factor) {
  std::visit(
      overloaded{
          [factor](Direct& s) {
            for (auto& msl : s) {
              msl.scaleThickness(static_cast<float>(factor));
            }
          },
          [factor](Indexed& s) {
            for (auto& msl : s.material) {
              msl.scaleThickness(static_cast<float>(factor));
            }
          },
          [factor](GloballyIndexed& s) {
            for (std::size_t index : s.indices) {
              (*s.material)[index].scaleThickness(static_cast<float>(factor));
            }
          },
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
