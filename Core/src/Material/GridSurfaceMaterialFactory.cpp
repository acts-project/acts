// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Material/GridSurfaceMaterialFactory.hpp"

#include "Acts/Utilities/AxisSpec.hpp"

#include <array>

namespace Acts {

namespace {

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

namespace GridSurfaceMaterialFactory {

std::unique_ptr<GridSurfaceMaterial> createDirect(
    const IAxis& axis0, const IAxis& axis1,
    const std::vector<std::vector<MaterialSlab>>& payload) {
  MultiAxisSpec2D binning = multiAxisSpecFromAxes(axis0, axis1);
  auto multiAxis = binning.buildMultiAxis();
  GridSurfaceMaterial::Direct storage = flattenPayload2D(*multiAxis, payload);
  return std::make_unique<GridSurfaceMaterial>(std::move(binning),
                                               std::move(storage));
}

std::unique_ptr<GridSurfaceMaterial> createDirect(
    const MultiAxisSpec2D& binning, const Surface& surface,
    const std::vector<std::vector<MaterialSlab>>& payload) {
  auto axes = resolveMultiAxis(binning, surface);
  return createDirect(axes->getAxis(0), axes->getAxis(1), payload);
}

std::unique_ptr<GridSurfaceMaterial> createIndexed(
    const IAxis& axis0, const IAxis& axis1, std::vector<MaterialSlab> material,
    const std::vector<std::vector<std::size_t>>& payload) {
  MultiAxisSpec2D binning = multiAxisSpecFromAxes(axis0, axis1);
  auto multiAxis = binning.buildMultiAxis();
  GridSurfaceMaterial::Indexed storage{flattenPayload2D(*multiAxis, payload),
                                       std::move(material)};
  return std::make_unique<GridSurfaceMaterial>(std::move(binning),
                                               std::move(storage));
}

std::unique_ptr<GridSurfaceMaterial> createIndexed(
    const MultiAxisSpec2D& binning, const Surface& surface,
    std::vector<MaterialSlab> material,
    const std::vector<std::vector<std::size_t>>& payload) {
  auto axes = resolveMultiAxis(binning, surface);
  return createIndexed(axes->getAxis(0), axes->getAxis(1), std::move(material),
                       payload);
}

std::unique_ptr<GridSurfaceMaterial> createGloballyIndexed(
    const IAxis& axis0, const IAxis& axis1,
    std::shared_ptr<std::vector<MaterialSlab>> material,
    const std::vector<std::vector<std::size_t>>& payload) {
  MultiAxisSpec2D binning = multiAxisSpecFromAxes(axis0, axis1);
  auto multiAxis = binning.buildMultiAxis();
  GridSurfaceMaterial::GloballyIndexed storage{
      flattenPayload2D(*multiAxis, payload), std::move(material)};
  return std::make_unique<GridSurfaceMaterial>(std::move(binning),
                                               std::move(storage));
}

std::unique_ptr<GridSurfaceMaterial> createGloballyIndexed(
    const MultiAxisSpec2D& binning, const Surface& surface,
    std::shared_ptr<std::vector<MaterialSlab>> material,
    const std::vector<std::vector<std::size_t>>& payload) {
  auto axes = resolveMultiAxis(binning, surface);
  return createGloballyIndexed(axes->getAxis(0), axes->getAxis(1),
                               std::move(material), payload);
}

}  // namespace GridSurfaceMaterialFactory

}  // namespace Acts
