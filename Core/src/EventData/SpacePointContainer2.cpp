// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/EventData/SpacePointContainer2.hpp"

namespace Acts::Experimental {

SpacePointContainer2::SpacePointContainer2(SpacePointColumns columns) noexcept {
  createColumns(columns);
}

SpacePointContainer2::SpacePointContainer2(
    const SpacePointContainer2 &other) noexcept
    : m_sourceLinkOffsets(other.m_sourceLinkOffsets),
      m_sourceLinkCounts(other.m_sourceLinkCounts),
      m_sourceLinks(other.m_sourceLinks) {
  copyColumns(other);
}

SpacePointContainer2::SpacePointContainer2(
    SpacePointContainer2 &&other) noexcept
    : m_sourceLinkOffsets(std::move(other.m_sourceLinkOffsets)),
      m_sourceLinkCounts(std::move(other.m_sourceLinkCounts)),
      m_sourceLinks(std::move(other.m_sourceLinks)) {
  moveColumns(other);
}

SpacePointContainer2 &SpacePointContainer2::operator=(
    const SpacePointContainer2 &other) noexcept {
  if (this == &other) {
    return *this;
  }

  m_sourceLinkOffsets = other.m_sourceLinkOffsets;
  m_sourceLinkCounts = other.m_sourceLinkCounts;
  m_sourceLinks = other.m_sourceLinks;

  copyColumns(other);

  return *this;
}

SpacePointContainer2 &SpacePointContainer2::operator=(
    SpacePointContainer2 &&other) noexcept {
  if (this == &other) {
    return *this;
  }

  m_sourceLinkOffsets = std::move(other.m_sourceLinkOffsets);
  m_sourceLinkCounts = std::move(other.m_sourceLinkCounts);
  m_sourceLinks = std::move(other.m_sourceLinks);

  moveColumns(other);

  return *this;
}

void SpacePointContainer2::copyColumns(const SpacePointContainer2 &other) {
  m_allColumns.reserve(other.m_allColumns.size());

  for (const auto &[name, column] : other.m_namedColumns) {
    m_namedColumns.try_emplace(name, column->copy());
  }

  m_allocatedColumns = other.m_allocatedColumns;
  knownColumns() = other.knownColumns();

  initializeColumns();
}

void SpacePointContainer2::moveColumns(SpacePointContainer2 &other) noexcept {
  m_allColumns.reserve(other.m_allColumns.size());

  for (auto &[name, column] : other.m_namedColumns) {
    m_namedColumns.try_emplace(name, std::move(column));
  }

  other.m_namedColumns.clear();
  other.m_allColumns.clear();

  m_allocatedColumns = other.m_allocatedColumns;
  knownColumns() = std::move(other).knownColumns();

  initializeColumns();
}

void SpacePointContainer2::initializeColumns() noexcept {
  m_allColumns.clear();

  for (const auto &[name, column] : m_namedColumns) {
    m_allColumns.push_back(column.get());
  }

  const auto appendExtraColumn = [this]<typename T>(std::optional<T> &column) {
    if (column.has_value()) {
      m_allColumns.push_back(&*column);
    }
  };
  std::apply([&](auto &...args) { ((appendExtraColumn(args)), ...); },
             knownColumns());
}

void SpacePointContainer2::reserve(std::size_t size,
                                   float averageSourceLinks) noexcept {
  m_sourceLinkOffsets.reserve(size);
  m_sourceLinkCounts.reserve(size);
  m_sourceLinks.reserve(static_cast<std::size_t>(size * averageSourceLinks));

  for (auto &column : m_allColumns) {
    column->reserve(size);
  }
}

void SpacePointContainer2::clear() noexcept {
  m_sourceLinkOffsets.clear();
  m_sourceLinkCounts.clear();
  m_sourceLinks.clear();

  for (auto &column : m_allColumns) {
    column->clear();
  }
}

void SpacePointContainer2::createColumns(SpacePointColumns columns) noexcept {
  using enum SpacePointColumns;

  if ((columns & X) != None && !m_xColumn.has_value()) {
    m_xColumn = SpacePointColumnHolder<float>();
    m_xColumn->resize(size());
    m_allColumns.push_back(&*m_xColumn);
  }
  if ((columns & Y) != None && !m_yColumn.has_value()) {
    m_yColumn = SpacePointColumnHolder<float>();
    m_yColumn->resize(size());
    m_allColumns.push_back(&*m_yColumn);
  }
  if ((columns & Z) != None && !m_zColumn.has_value()) {
    m_zColumn = SpacePointColumnHolder<float>();
    m_zColumn->resize(size());
    m_allColumns.push_back(&*m_zColumn);
  }
  if ((columns & R) != None && !m_rColumn.has_value()) {
    m_rColumn = SpacePointColumnHolder<float>();
    m_rColumn->resize(size());
    m_allColumns.push_back(&*m_rColumn);
  }
  if ((columns & Phi) != None && !m_phiColumn.has_value()) {
    m_phiColumn = SpacePointColumnHolder<float>();
    m_phiColumn->resize(size());
    m_allColumns.push_back(&*m_phiColumn);
  }
  if ((columns & Time) != None && !m_timeColumn.has_value()) {
    m_timeColumn = SpacePointColumnHolder<float>(NoTime);
    m_timeColumn->resize(size());
    m_allColumns.push_back(&*m_timeColumn);
  }
  if ((columns & VarianceZ) != None && !m_varianceZColumn.has_value()) {
    m_varianceZColumn = SpacePointColumnHolder<float>();
    m_varianceZColumn->resize(size());
    m_allColumns.push_back(&*m_varianceZColumn);
  }
  if ((columns & VarianceR) != None && !m_varianceRColumn.has_value()) {
    m_varianceRColumn = SpacePointColumnHolder<float>();
    m_varianceRColumn->resize(size());
    m_allColumns.push_back(&*m_varianceRColumn);
  }
  if ((columns & TopStripVector) != None &&
      !m_topStripVectorColumn.has_value()) {
    m_topStripVectorColumn = SpacePointColumnHolder<Eigen::Vector3f>();
    m_topStripVectorColumn->resize(size());
    m_allColumns.push_back(&*m_topStripVectorColumn);
  }
  if ((columns & BottomStripVector) != None &&
      !m_bottomStripVectorColumn.has_value()) {
    m_bottomStripVectorColumn = SpacePointColumnHolder<Eigen::Vector3f>();
    m_bottomStripVectorColumn->resize(size());
    m_allColumns.push_back(&*m_bottomStripVectorColumn);
  }
  if ((columns & StripCenterDistance) != None &&
      !m_stripCenterDistanceColumn.has_value()) {
    m_stripCenterDistanceColumn = SpacePointColumnHolder<Eigen::Vector3f>();
    m_stripCenterDistanceColumn->resize(size());
    m_allColumns.push_back(&*m_stripCenterDistanceColumn);
  }
  if ((columns & TopStripCenter) != None &&
      !m_topStripCenterColumn.has_value()) {
    m_topStripCenterColumn = SpacePointColumnHolder<Eigen::Vector3f>();
    m_topStripCenterColumn->resize(size());
    m_allColumns.push_back(&*m_topStripCenterColumn);
  }
  if ((columns & CopyFromIndex) != None && !m_copyFromIndexColumn.has_value()) {
    m_copyFromIndexColumn = SpacePointColumnHolder<std::size_t>();
    m_copyFromIndexColumn->resize(size());
    m_allColumns.push_back(&*m_copyFromIndexColumn);
  }

  m_allocatedColumns = m_allocatedColumns | columns;
}

}  // namespace Acts::Experimental
