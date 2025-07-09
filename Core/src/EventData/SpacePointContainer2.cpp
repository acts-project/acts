// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/EventData/SpacePointContainer2.hpp"

#include "Acts/Utilities/Helpers.hpp"

#include <unordered_set>

namespace Acts::Experimental {

SpacePointContainer2::SpacePointContainer2(SpacePointColumns columns) noexcept {
  createColumns(columns);
}

SpacePointContainer2::SpacePointContainer2(
    const SpacePointContainer2 &other) noexcept
    : m_size(other.m_size), m_sourceLinks(other.m_sourceLinks) {
  copyColumns(other);
}

SpacePointContainer2::SpacePointContainer2(
    SpacePointContainer2 &&other) noexcept
    : m_size(other.m_size), m_sourceLinks(std::move(other.m_sourceLinks)) {
  moveColumns(other);
}

SpacePointContainer2 &SpacePointContainer2::operator=(
    const SpacePointContainer2 &other) noexcept {
  if (this == &other) {
    return *this;
  }

  m_sourceLinks = other.m_sourceLinks;

  copyColumns(other);

  return *this;
}

SpacePointContainer2 &SpacePointContainer2::operator=(
    SpacePointContainer2 &&other) noexcept {
  if (this == &other) {
    return *this;
  }

  m_sourceLinks = std::move(other.m_sourceLinks);

  moveColumns(other);

  return *this;
}

void SpacePointContainer2::copyColumns(const SpacePointContainer2 &other) {
  m_namedColumns.reserve(other.m_namedColumns.size());

  for (const auto &[name, column] : other.m_namedColumns) {
    std::unique_ptr<SpacePointColumnHolderBase> holder =
        column.second != nullptr ? column.second->copy() : nullptr;
    m_namedColumns.try_emplace(name,
                               std::pair{holder.get(), std::move(holder)});
  }

  m_knownColumns = other.m_knownColumns;
  knownColumns() = other.knownColumns();
}

void SpacePointContainer2::moveColumns(SpacePointContainer2 &other) noexcept {
  m_namedColumns.reserve(other.m_namedColumns.size());

  for (auto &[name, column] : other.m_namedColumns) {
    m_namedColumns.try_emplace(name, std::move(column));
  }

  other.m_namedColumns.clear();

  m_knownColumns = other.m_knownColumns;
  knownColumns() = std::move(other).knownColumns();
}

void SpacePointContainer2::reserve(std::uint32_t size,
                                   float averageSourceLinks) noexcept {
  for (const auto &[name, column] : m_namedColumns) {
    column.first->reserve(size);
  }

  if (hasColumns(SpacePointColumns::SourceLinks)) {
    m_sourceLinks.reserve(
        static_cast<std::uint32_t>(size * averageSourceLinks));
  }
}

void SpacePointContainer2::clear() noexcept {
  for (const auto &[name, column] : m_namedColumns) {
    column.first->clear();
  }

  m_sourceLinks.clear();
}

void SpacePointContainer2::assignSourceLinks(
    IndexType index, std::span<const SourceLink> sourceLinks) {
  if (index >= size()) {
    throw std::out_of_range("Index out of range in SpacePointContainer2");
  }
  if (!m_sourceLinkOffsetColumn.has_value() ||
      !m_sourceLinkCountColumn.has_value()) {
    throw std::logic_error("No source links column available");
  }
  if (entry(m_sourceLinkCountColumn->proxy(), index) != 0) {
    throw std::logic_error("Source links already assigned to the space point");
  }

  m_sourceLinkOffsetColumn->proxy().data().push_back(
      static_cast<SpacePointIndex2>(m_sourceLinks.size()));
  m_sourceLinkCountColumn->proxy().data().push_back(
      static_cast<std::uint8_t>(sourceLinks.size()));
  m_sourceLinks.insert(m_sourceLinks.end(), sourceLinks.begin(),
                       sourceLinks.end());
}

void SpacePointContainer2::createColumns(SpacePointColumns columns) noexcept {
  using enum SpacePointColumns;

  const auto createColumn =
      [&]<typename T>(SpacePointColumns mask, std::string_view name,
                      T defaultValue,
                      std::optional<SpacePointColumnHolder<T>> &column) {
        if (ACTS_CHECK_BIT(columns, mask) && !column.has_value()) {
          column = SpacePointColumnHolder<T>(std::move(defaultValue));
          column->resize(size());
          m_namedColumns.try_emplace(std::string(name),
                                     std::pair{&column.value(), nullptr});
          m_knownColumns = m_knownColumns | mask;
        }
      };

  [&]<std::size_t... Is>(std::index_sequence<Is...>) {
    ((createColumn(
         std::get<Is>(knownColumnMaks()), std::get<Is>(knownColumnNames()),
         std::get<Is>(knownColumnDefaults()), std::get<Is>(knownColumns()))),
     ...);
  }(std::index_sequence_for<decltype(knownColumns())>{});
}

void SpacePointContainer2::dropColumns(SpacePointColumns columns) noexcept {
  using enum SpacePointColumns;

  const auto dropColumn =
      [&]<typename T>(SpacePointColumns mask, std::string_view name,
                      std::optional<SpacePointColumnHolder<T>> &column) {
        if (ACTS_CHECK_BIT(columns, mask) && column.has_value()) {
          m_namedColumns.erase(std::string(name));
          column.reset();
          m_knownColumns = m_knownColumns & ~mask;
        }
      };

  [&]<std::size_t... Is>(std::index_sequence<Is...>) {
    ((dropColumn(std::get<Is>(knownColumnMaks()),
                 std::get<Is>(knownColumnNames()),
                 std::get<Is>(knownColumns()))),
     ...);
  }(std::index_sequence_for<decltype(knownColumns())>{});
}

void SpacePointContainer2::dropColumn(const std::string &name) {
  if (reservedColumn(name)) {
    throw std::runtime_error("Cannot drop reserved column: " + name);
  }

  auto it = m_namedColumns.find(name);
  if (it == m_namedColumns.end()) {
    throw std::runtime_error("Column does not exist: " + name);
  }

  m_namedColumns.erase(it);
}

bool SpacePointContainer2::reservedColumn(const std::string &name) noexcept {
  static const std::unordered_set<std::string> reservedColumns = std::apply(
      [](auto... reservedNames) {
        return std::unordered_set<std::string>({reservedNames...});
      },
      knownColumnNames());

  return reservedColumns.find(name) != reservedColumns.end();
}

}  // namespace Acts::Experimental
