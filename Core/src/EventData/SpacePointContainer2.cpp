// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/EventData/SpacePointContainer2.hpp"

#include "Acts/Utilities/Helpers.hpp"

namespace Acts::Experimental {

SpacePointContainer2::SpacePointContainer2(
    const SpacePointContainer2 &other) noexcept
    : m_x(other.m_x),
      m_y(other.m_y),
      m_z(other.m_z),
      m_sourceLinkOffsets(other.m_sourceLinkOffsets),
      m_sourceLinkCounts(other.m_sourceLinkCounts),
      m_sourceLinks(other.m_sourceLinks) {
  copyExtraColumns(other);
}

SpacePointContainer2::SpacePointContainer2(
    SpacePointContainer2 &&other) noexcept
    : m_x(std::move(other.m_x)),
      m_y(std::move(other.m_y)),
      m_z(std::move(other.m_z)),
      m_sourceLinkOffsets(std::move(other.m_sourceLinkOffsets)),
      m_sourceLinkCounts(std::move(other.m_sourceLinkCounts)),
      m_sourceLinks(std::move(other.m_sourceLinks)) {
  moveExtraColumns(other);
}

SpacePointContainer2 &SpacePointContainer2::operator=(
    const SpacePointContainer2 &other) noexcept {
  if (this == &other) {
    return *this;
  }

  m_x = other.m_x;
  m_y = other.m_y;
  m_z = other.m_z;
  m_sourceLinkOffsets = other.m_sourceLinkOffsets;
  m_sourceLinkCounts = other.m_sourceLinkCounts;
  m_sourceLinks = other.m_sourceLinks;

  copyExtraColumns(other);

  return *this;
}

SpacePointContainer2 &SpacePointContainer2::operator=(
    SpacePointContainer2 &&other) noexcept {
  if (this == &other) {
    return *this;
  }

  m_x = std::move(other.m_x);
  m_y = std::move(other.m_y);
  m_z = std::move(other.m_z);
  m_sourceLinkOffsets = std::move(other.m_sourceLinkOffsets);
  m_sourceLinkCounts = std::move(other.m_sourceLinkCounts);
  m_sourceLinks = std::move(other.m_sourceLinks);

  moveExtraColumns(other);

  return *this;
}

void SpacePointContainer2::copyExtraColumns(const SpacePointContainer2 &other) {
  m_extraColumns.reserve(other.m_extraColumns.size());

  for (const auto &[name, column] : other.m_namedExtraColumns) {
    m_namedExtraColumns.try_emplace(name, column->copy());
  }

  m_knownExtraColumns = other.m_knownExtraColumns;
  knownExtraColumns() = other.knownExtraColumns();

  initializeExtraColumns();
}

void SpacePointContainer2::moveExtraColumns(
    SpacePointContainer2 &other) noexcept {
  m_extraColumns.reserve(other.m_extraColumns.size());

  for (auto &[name, column] : other.m_namedExtraColumns) {
    m_namedExtraColumns.try_emplace(name, std::move(column));
  }

  other.m_namedExtraColumns.clear();
  other.m_extraColumns.clear();

  m_knownExtraColumns = other.m_knownExtraColumns;
  knownExtraColumns() = std::move(other).knownExtraColumns();

  initializeExtraColumns();
}

void SpacePointContainer2::initializeExtraColumns() noexcept {
  m_extraColumns.clear();

  for (const auto &[name, column] : m_namedExtraColumns) {
    m_extraColumns.push_back(column.get());
  }

  const auto appendExtraColumn = [this]<typename T>(std::optional<T> &column) {
    if (column.has_value()) {
      m_extraColumns.push_back(&*column);
    }
  };
  std::apply([&](auto &...args) { ((appendExtraColumn(args)), ...); },
             knownExtraColumns());
}

void SpacePointContainer2::reserve(std::size_t size,
                                   float averageSourceLinks) noexcept {
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

void SpacePointContainer2::clear() noexcept {
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
    SpacePointKnownExtraColumn columns) noexcept {
  using enum SpacePointKnownExtraColumn;

  if (ACTS_CHECK_BIT(columns, R) && !m_rColumn.has_value()) {
    m_rColumn = SpacePointExtraColumnHolder<float>();
    m_rColumn->resize(size());
  }
  if (ACTS_CHECK_BIT(columns, Phi) && !m_phiColumn.has_value()) {
    m_phiColumn = SpacePointExtraColumnHolder<float>();
    m_phiColumn->resize(size());
  }
  if (ACTS_CHECK_BIT(columns, Time) && !m_timeColumn.has_value()) {
    m_timeColumn = SpacePointExtraColumnHolder<float>(NoTime);
    m_timeColumn->resize(size());
  }
  if (ACTS_CHECK_BIT(columns, VarianceZ) && !m_varianceZColumn.has_value()) {
    m_varianceZColumn = SpacePointExtraColumnHolder<float>();
    m_varianceZColumn->resize(size());
  }
  if (ACTS_CHECK_BIT(columns, VarianceR) && !m_varianceRColumn.has_value()) {
    m_varianceRColumn = SpacePointExtraColumnHolder<float>();
    m_varianceRColumn->resize(size());
  }
  if (ACTS_CHECK_BIT(columns, TopStripVector) &&
      !m_topStripVectorColumn.has_value()) {
    m_topStripVectorColumn = SpacePointExtraColumnHolder<Eigen::Vector3f>();
    m_topStripVectorColumn->resize(size());
  }
  if (ACTS_CHECK_BIT(columns, BottomStripVector) &&
      !m_bottomStripVectorColumn.has_value()) {
    m_bottomStripVectorColumn = SpacePointExtraColumnHolder<Eigen::Vector3f>();
    m_bottomStripVectorColumn->resize(size());
  }
  if (ACTS_CHECK_BIT(columns, StripCenterDistance) &&
      !m_stripCenterDistanceColumn.has_value()) {
    m_stripCenterDistanceColumn =
        SpacePointExtraColumnHolder<Eigen::Vector3f>();
    m_stripCenterDistanceColumn->resize(size());
  }
  if (ACTS_CHECK_BIT(columns, TopStripCenter) &&
      !m_topStripCenterColumn.has_value()) {
    m_topStripCenterColumn = SpacePointExtraColumnHolder<Eigen::Vector3f>();
    m_topStripCenterColumn->resize(size());
  }
  if (ACTS_CHECK_BIT(columns, CopyFromIndex) &&
      !m_copyFromIndexColumn.has_value()) {
    m_copyFromIndexColumn = SpacePointExtraColumnHolder<std::size_t>();
    m_copyFromIndexColumn->resize(size());
  }

  m_knownExtraColumns = m_knownExtraColumns | columns;

  initializeExtraColumns();
}

void SpacePointContainer2::dropExtraColumns(
    SpacePointKnownExtraColumn columns) noexcept {
  using enum SpacePointKnownExtraColumn;

  if (ACTS_CHECK_BIT(columns, R)) {
    m_rColumn.reset();
  }
  if (ACTS_CHECK_BIT(columns, Phi)) {
    m_phiColumn.reset();
  }
  if (ACTS_CHECK_BIT(columns, Time)) {
    m_timeColumn.reset();
  }
  if (ACTS_CHECK_BIT(columns, VarianceZ)) {
    m_varianceZColumn.reset();
  }
  if (ACTS_CHECK_BIT(columns, VarianceR)) {
    m_varianceRColumn.reset();
  }
  if (ACTS_CHECK_BIT(columns, TopStripVector)) {
    m_topStripVectorColumn.reset();
  }
  if (ACTS_CHECK_BIT(columns, BottomStripVector)) {
    m_bottomStripVectorColumn.reset();
  }
  if (ACTS_CHECK_BIT(columns, StripCenterDistance)) {
    m_stripCenterDistanceColumn.reset();
  }
  if (ACTS_CHECK_BIT(columns, TopStripCenter)) {
    m_topStripCenterColumn.reset();
  }
  if (ACTS_CHECK_BIT(columns, CopyFromIndex)) {
    m_copyFromIndexColumn.reset();
  }

  m_knownExtraColumns = m_knownExtraColumns & ~columns;

  initializeExtraColumns();
}

void SpacePointContainer2::dropExtraColumn(const std::string &name) {
  auto it = m_namedExtraColumns.find(name);

  if (it == m_namedExtraColumns.end()) {
    throw std::runtime_error("Extra column '" + name + "' does not exist");
  }

  m_namedExtraColumns.erase(it);

  initializeExtraColumns();
}

bool SpacePointContainer2::hasExtraColumn(
    const std::string &name) const noexcept {
  return m_namedExtraColumns.contains(name);
}

}  // namespace Acts::Experimental
