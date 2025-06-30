// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/EventData/SpacePointContainer2.hpp"

namespace Acts::Experimental {

SpacePointContainer2::SpacePointContainer2(const SpacePointContainer2 &other)
    : m_x(other.m_x),
      m_y(other.m_y),
      m_z(other.m_z),
      m_sourceLinkOffsets(other.m_sourceLinkOffsets),
      m_sourceLinkCounts(other.m_sourceLinkCounts),
      m_sourceLinks(other.m_sourceLinks),
      m_rColumn(other.m_rColumn),
      m_phiColumn(other.m_phiColumn),
      m_timeColumn(other.m_timeColumn),
      m_varianceZColumn(other.m_varianceZColumn),
      m_varianceRColumn(other.m_varianceRColumn),
      m_topStripVectorColumn(other.m_topStripVectorColumn),
      m_bottomStripVectorColumn(other.m_bottomStripVectorColumn),
      m_stripCenterDistanceColumn(other.m_stripCenterDistanceColumn),
      m_topStripCenterColumn(other.m_topStripCenterColumn),
      m_copyFromIndexColumn(other.m_copyFromIndexColumn) {
  for (const auto &[name, column] : other.m_namedExtraColumns) {
    m_namedExtraColumns.emplace(name, column->copy());
  }
}

SpacePointContainer2 &SpacePointContainer2::operator=(
    const SpacePointContainer2 &other) {
  if (this == &other) {
    return *this;
  }

  m_x = other.m_x;
  m_y = other.m_y;
  m_z = other.m_z;
  m_sourceLinkOffsets = other.m_sourceLinkOffsets;
  m_sourceLinkCounts = other.m_sourceLinkCounts;
  m_sourceLinks = other.m_sourceLinks;
  m_rColumn = other.m_rColumn;
  m_phiColumn = other.m_phiColumn;
  m_timeColumn = other.m_timeColumn;
  m_varianceZColumn = other.m_varianceZColumn;
  m_varianceRColumn = other.m_varianceRColumn;
  m_topStripVectorColumn = other.m_topStripVectorColumn;
  m_bottomStripVectorColumn = other.m_bottomStripVectorColumn;
  m_stripCenterDistanceColumn = other.m_stripCenterDistanceColumn;
  m_topStripCenterColumn = other.m_topStripCenterColumn;
  m_copyFromIndexColumn = other.m_copyFromIndexColumn;

  m_extraColumns.clear();
  for (const auto &[name, column] : other.m_namedExtraColumns) {
    m_namedExtraColumns.emplace(name, column->copy());
  }
  return *this;
}

void SpacePointContainer2::reserve(std::size_t size, float averageSourceLinks) {
  m_x.reserve(size);
  m_y.reserve(size);
  m_z.reserve(size);
  m_sourceLinkOffsets.reserve(size);
  m_sourceLinkCounts.reserve(size);
  m_sourceLinks.reserve(static_cast<std::size_t>(size * averageSourceLinks));

  for (auto &column : m_extraColumns) {
    column->reserve(size);
  }
}

void SpacePointContainer2::clear() {
  m_x.clear();
  m_y.clear();
  m_z.clear();
  m_sourceLinkOffsets.clear();
  m_sourceLinkCounts.clear();
  m_sourceLinks.clear();

  for (auto &column : m_extraColumns) {
    column->clear();
  }
}

void SpacePointContainer2::createExtraColumns(
    SpacePointKnownExtraColumn columns) {
  using enum SpacePointKnownExtraColumn;

  if ((columns & R) != None && !m_rColumn.has_value()) {
    m_rColumn = SpacePointExtraColumnHolder<float>();
    m_rColumn->resize(size());
    m_extraColumns.push_back(&*m_rColumn);
  }
  if ((columns & Phi) != None && !m_phiColumn.has_value()) {
    m_phiColumn = SpacePointExtraColumnHolder<float>();
    m_phiColumn->resize(size());
    m_extraColumns.push_back(&*m_phiColumn);
  }
  if ((columns & Time) != None && !m_timeColumn.has_value()) {
    m_timeColumn = SpacePointExtraColumnHolder<float>(NoTime);
    m_timeColumn->resize(size());
    m_extraColumns.push_back(&*m_timeColumn);
  }
  if ((columns & VarianceZ) != None && !m_varianceZColumn.has_value()) {
    m_varianceZColumn = SpacePointExtraColumnHolder<float>();
    m_varianceZColumn->resize(size());
    m_extraColumns.push_back(&*m_varianceZColumn);
  }
  if ((columns & VarianceR) != None && !m_varianceRColumn.has_value()) {
    m_varianceRColumn = SpacePointExtraColumnHolder<float>();
    m_varianceRColumn->resize(size());
    m_extraColumns.push_back(&*m_varianceRColumn);
  }
  if ((columns & TopStripVector) != None &&
      !m_topStripVectorColumn.has_value()) {
    m_topStripVectorColumn = SpacePointExtraColumnHolder<Eigen::Vector3f>();
    m_topStripVectorColumn->resize(size());
    m_extraColumns.push_back(&*m_topStripVectorColumn);
  }
  if ((columns & BottomStripVector) != None &&
      !m_bottomStripVectorColumn.has_value()) {
    m_bottomStripVectorColumn = SpacePointExtraColumnHolder<Eigen::Vector3f>();
    m_bottomStripVectorColumn->resize(size());
    m_extraColumns.push_back(&*m_bottomStripVectorColumn);
  }
  if ((columns & StripCenterDistance) != None &&
      !m_stripCenterDistanceColumn.has_value()) {
    m_stripCenterDistanceColumn =
        SpacePointExtraColumnHolder<Eigen::Vector3f>();
    m_stripCenterDistanceColumn->resize(size());
    m_extraColumns.push_back(&*m_stripCenterDistanceColumn);
  }
  if ((columns & TopStripCenter) != None &&
      !m_topStripCenterColumn.has_value()) {
    m_topStripCenterColumn = SpacePointExtraColumnHolder<Eigen::Vector3f>();
    m_topStripCenterColumn->resize(size());
    m_extraColumns.push_back(&*m_topStripCenterColumn);
  }
  if ((columns & CopyFromIndex) != None && !m_copyFromIndexColumn.has_value()) {
    m_copyFromIndexColumn = SpacePointExtraColumnHolder<std::size_t>();
    m_copyFromIndexColumn->resize(size());
    m_extraColumns.push_back(&*m_copyFromIndexColumn);
  }
}

bool SpacePointContainer2::hasExtraColumns(
    SpacePointKnownExtraColumn columns) const {
  using enum SpacePointKnownExtraColumn;

  if ((columns & R) != None && !m_rColumn.has_value()) {
    return false;
  }
  if ((columns & Phi) != None && !m_phiColumn.has_value()) {
    return false;
  }
  if ((columns & Time) != None && !m_timeColumn.has_value()) {
    return false;
  }
  if ((columns & VarianceZ) != None && !m_varianceZColumn.has_value()) {
    return false;
  }
  if ((columns & VarianceR) != None && !m_varianceRColumn.has_value()) {
    return false;
  }
  if ((columns & TopStripVector) != None &&
      !m_topStripVectorColumn.has_value()) {
    return false;
  }
  if ((columns & BottomStripVector) != None &&
      !m_bottomStripVectorColumn.has_value()) {
    return false;
  }
  if ((columns & StripCenterDistance) != None &&
      !m_stripCenterDistanceColumn.has_value()) {
    return false;
  }
  if ((columns & TopStripCenter) != None &&
      !m_topStripCenterColumn.has_value()) {
    return false;
  }
  if ((columns & CopyFromIndex) != None && !m_copyFromIndexColumn.has_value()) {
    return false;
  }
  return true;
}

}  // namespace Acts::Experimental
