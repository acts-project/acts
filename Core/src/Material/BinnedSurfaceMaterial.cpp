// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Material/BinnedSurfaceMaterial.hpp"

#include "Acts/Material/MaterialSlab.hpp"

#include <array>
#include <ostream>
#include <stdexcept>
#include <utility>
#include <vector>

Acts::BinnedSurfaceMaterial::BinnedSurfaceMaterial(
    const BinUtility& binUtility, MaterialSlabVector fullProperties,
    double splitFactor, MappingType mappingType)
    : ISurfaceMaterial(splitFactor, mappingType) {
  auto converted = convertBinUtility(binUtility);
  m_directedProtoAxes = std::move(converted.first);
  m_globalToLocalTransform = std::move(converted.second);
  // fill the material with deep copy
  m_fullMaterial.push_back(std::move(fullProperties));
}

Acts::BinnedSurfaceMaterial::BinnedSurfaceMaterial(
    const BinUtility& binUtility, MaterialSlabMatrix fullProperties,
    double splitFactor, MappingType mappingType)
    : ISurfaceMaterial(splitFactor, mappingType),
      m_fullMaterial(std::move(fullProperties)) {
  auto converted = convertBinUtility(binUtility);
  m_directedProtoAxes = std::move(converted.first);
  m_globalToLocalTransform = std::move(converted.second);
}

Acts::BinnedSurfaceMaterial::BinnedSurfaceMaterial(
    DirectedProtoAxis directedProtoAxis, Transform3 globalToLocalTransform,
    MaterialSlabVector fullProperties, double splitFactor,
    MappingType mappingType)
    : ISurfaceMaterial(splitFactor, mappingType),
      m_directedProtoAxes{std::move(directedProtoAxis)},
      m_globalToLocalTransform(std::move(globalToLocalTransform)) {
  m_fullMaterial.push_back(std::move(fullProperties));
}

Acts::BinnedSurfaceMaterial::BinnedSurfaceMaterial(
    std::array<DirectedProtoAxis, 2> directedProtoAxes,
    Transform3 globalToLocalTransform, MaterialSlabMatrix fullProperties,
    double splitFactor, MappingType mappingType)
    : ISurfaceMaterial(splitFactor, mappingType),
      m_directedProtoAxes{std::move(directedProtoAxes[0u]),
                          std::move(directedProtoAxes[1u])},
      m_globalToLocalTransform(std::move(globalToLocalTransform)),
      m_fullMaterial(std::move(fullProperties)) {}

Acts::BinnedSurfaceMaterial& Acts::BinnedSurfaceMaterial::scale(double factor) {
  for (auto& materialVector : m_fullMaterial) {
    for (auto& materialBin : materialVector) {
      materialBin.scaleThickness(factor);
    }
  }
  return (*this);
}

const Acts::MaterialSlab& Acts::BinnedSurfaceMaterial::materialSlab(
    const Vector2& lp) const {
  if (m_directedProtoAxes.empty()) {
    throw std::logic_error(
        "BinnedSurfaceMaterial has no DirectedProtoAxis configured.");
  }
  const Vector3 localPosition(lp[0], lp[1], 0.);
  std::size_t ibin0 = 0u;
  std::size_t ibin1 = 0u;
  if (m_directedProtoAxes.size() == 1u) {
    std::array<double, 1u> casted{};
    const std::array<AxisDirection, 1u> castDirs{
        m_directedProtoAxes[0u].getAxisDirection()};
    GridAccessHelpers::fillCasts(localPosition, castDirs, casted,
                                 std::make_index_sequence<1u>{});
    ibin0 = correctedBinIndex(m_directedProtoAxes[0u].getAxis(), casted[0u]);
  } else {
    std::array<double, 2u> casted{};
    const std::array<AxisDirection, 2u> castDirs{
        m_directedProtoAxes[0u].getAxisDirection(),
        m_directedProtoAxes[1u].getAxisDirection()};
    GridAccessHelpers::fillCasts(localPosition, castDirs, casted,
                                 std::make_index_sequence<2u>{});
    ibin0 = correctedBinIndex(m_directedProtoAxes[0u].getAxis(), casted[0u]);
    ibin1 = correctedBinIndex(m_directedProtoAxes[1u].getAxis(), casted[1u]);
  }
  return m_fullMaterial[ibin1][ibin0];
}

const Acts::MaterialSlab& Acts::BinnedSurfaceMaterial::materialSlab(
    const Acts::Vector3& gp) const {
  if (m_directedProtoAxes.empty()) {
    throw std::logic_error(
        "BinnedSurfaceMaterial has no DirectedProtoAxis configured.");
  }
  const Vector3 localPosition = m_globalToLocalTransform * gp;
  std::size_t ibin0 = 0u;
  std::size_t ibin1 = 0u;
  if (m_directedProtoAxes.size() == 1u) {
    std::array<double, 1u> casted{};
    const std::array<AxisDirection, 1u> castDirs{
        m_directedProtoAxes[0u].getAxisDirection()};
    GridAccessHelpers::fillCasts(localPosition, castDirs, casted,
                                 std::make_index_sequence<1u>{});
    ibin0 = correctedBinIndex(m_directedProtoAxes[0u].getAxis(), casted[0u]);
  } else {
    std::array<double, 2u> casted{};
    const std::array<AxisDirection, 2u> castDirs{
        m_directedProtoAxes[0u].getAxisDirection(),
        m_directedProtoAxes[1u].getAxisDirection()};
    GridAccessHelpers::fillCasts(localPosition, castDirs, casted,
                                 std::make_index_sequence<2u>{});
    ibin0 = correctedBinIndex(m_directedProtoAxes[0u].getAxis(), casted[0u]);
    ibin1 = correctedBinIndex(m_directedProtoAxes[1u].getAxis(), casted[1u]);
  }
  return m_fullMaterial[ibin1][ibin0];
}

std::ostream& Acts::BinnedSurfaceMaterial::toStream(std::ostream& sl) const {
  const std::size_t bins0 =
      !m_fullMaterial.empty() ? m_fullMaterial[0u].size() : 0u;
  const std::size_t bins1 = m_fullMaterial.size();
  sl << "Acts::BinnedSurfaceMaterial : " << std::endl;
  sl << "   - Number of Material bins [0,1] : " << bins0 << " / " << bins1
     << std::endl;
  sl << "   - Number of DirectedProtoAxis   : " << m_directedProtoAxes.size()
     << std::endl;
  sl << "   - Parse full update material    : " << std::endl;  //
  // output  the full material
  unsigned int imat1 = 0;
  for (auto& materialVector : m_fullMaterial) {
    unsigned int imat0 = 0;
    // the vector iterator
    for (auto& materialBin : materialVector) {
      sl << " Bin [" << imat1 << "][" << imat0 << "] - " << (materialBin);
      ++imat0;
    }
    ++imat1;
  }
  sl << "  - DirectedProtoAxes: " << m_directedProtoAxes << std::endl;
  sl << "  - GlobalToLocalTransform:\n"
     << m_globalToLocalTransform.matrix() << std::endl;
  return sl;
}
