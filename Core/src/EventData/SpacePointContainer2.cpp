// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/EventData/SpacePointContainer2.hpp"

#include "Acts/Utilities/Helpers.hpp"

#include <string_view>
#include <unordered_set>

namespace {

template <typename Tuple>
using tuple_indices =
    std::make_index_sequence<std::tuple_size_v<std::remove_reference_t<Tuple>>>;

}

namespace Acts {

static_assert(std::random_access_iterator<SpacePointContainer2::iterator>);
static_assert(
    std::random_access_iterator<SpacePointContainer2::const_iterator>);
static_assert(
    std::random_access_iterator<SpacePointContainer2::MutableSubset::Iterator>);
static_assert(
    std::random_access_iterator<SpacePointContainer2::ConstSubset::Iterator>);

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

  other.m_size = 0;
}

SpacePointContainer2 &SpacePointContainer2::operator=(
    const SpacePointContainer2 &other) noexcept {
  if (this == &other) {
    return *this;
  }

  copyColumns(other);
  m_sourceLinks = other.m_sourceLinks;
  m_size = other.m_size;

  return *this;
}

SpacePointContainer2 &SpacePointContainer2::operator=(
    SpacePointContainer2 &&other) noexcept {
  if (this == &other) {
    return *this;
  }

  moveColumns(other);
  m_sourceLinks = std::move(other.m_sourceLinks);
  m_size = other.m_size;

  other.m_size = 0;

  return *this;
}

void SpacePointContainer2::copyColumns(const SpacePointContainer2 &other) {
  m_allColumns.reserve(other.m_allColumns.size());
  m_dynamicColumns.reserve(other.m_dynamicColumns.size());

  const auto copyKnownColumns =
      [&]<typename T>(std::string_view name,
                      std::optional<ColumnHolder<T>> &column,
                      const std::optional<ColumnHolder<T>> &otherColumn) {
        column = otherColumn;
        if (column.has_value()) {
          m_allColumns.try_emplace(std::string(name), &column.value());
        }
      };

  [&]<std::size_t... Is>(std::index_sequence<Is...>) {
    ((copyKnownColumns(std::get<Is>(knownColumnNames()),
                       std::get<Is>(knownColumns()),
                       std::get<Is>(other.knownColumns()))),
     ...);
  }(tuple_indices<decltype(knownColumns())>{});
  m_knownColumns = other.m_knownColumns;

  for (auto &[name, column] : other.m_dynamicColumns) {
    std::unique_ptr<ColumnHolderBase> columnCopy = column->copy();
    m_allColumns.try_emplace(name, columnCopy.get());
    m_dynamicColumns.try_emplace(name, std::move(columnCopy));
  }
}

void SpacePointContainer2::moveColumns(SpacePointContainer2 &other) noexcept {
  m_allColumns.reserve(other.m_allColumns.size());
  m_dynamicColumns.reserve(other.m_dynamicColumns.size());

  const auto moveKnownColumns =
      [&]<typename T>(std::string_view name,
                      std::optional<ColumnHolder<T>> &column,
                      std::optional<ColumnHolder<T>> &otherColumn) {
        column = std::move(otherColumn);
        if (column.has_value()) {
          m_allColumns.try_emplace(std::string(name), &column.value());
        }
      };

  [&]<std::size_t... Is>(std::index_sequence<Is...>) {
    ((moveKnownColumns(std::get<Is>(knownColumnNames()),
                       std::get<Is>(knownColumns()),
                       std::get<Is>(other.knownColumns()))),
     ...);
  }(tuple_indices<decltype(knownColumns())>{});
  m_knownColumns = other.m_knownColumns;

  for (auto &[name, column] : other.m_dynamicColumns) {
    m_allColumns.try_emplace(name, column.get());
    m_dynamicColumns.try_emplace(name, std::move(column));
  }

  other.m_allColumns.clear();
  other.m_knownColumns = SpacePointColumns::None;
  other.m_dynamicColumns.clear();
}

void SpacePointContainer2::reserve(std::uint32_t size,
                                   float averageSourceLinks) noexcept {
  if (hasColumns(SpacePointColumns::SourceLinks)) {
    m_sourceLinks.reserve(
        static_cast<std::uint32_t>(size * averageSourceLinks));
  }

  for (const auto &[name, column] : m_allColumns) {
    column->reserve(size);
  }
}

void SpacePointContainer2::clear() noexcept {
  m_size = 0;
  m_sourceLinks.clear();

  for (const auto &[name, column] : m_allColumns) {
    column->clear();
  }
}

MutableSpacePointProxy2 SpacePointContainer2::createSpacePoint() noexcept {
  ++m_size;

  for (const auto &[name, column] : m_allColumns) {
    column->emplace_back();
  }

  return MutableProxy(*this, size() - 1);
}

void SpacePointContainer2::copyFrom(Index index,
                                    const SpacePointContainer2 &sourceContainer,
                                    Index sourceIndex,
                                    SpacePointColumns columnsToCopy) {
  if (index >= size() || sourceIndex >= sourceContainer.size()) {
    throw std::out_of_range(
        "Index out of range in SpacePointContainer2::copyFrom");
  }
  if ((columnsToCopy & sourceContainer.m_knownColumns) != columnsToCopy) {
    throw std::logic_error(
        "Source container does not have all columns to copy");
  }
  if ((columnsToCopy & m_knownColumns) != columnsToCopy) {
    throw std::logic_error(
        "Destination container does not have all columns to copy");
  }

  const auto copyColumn =
      [&]<typename T>(SpacePointColumns mask,
                      std::optional<ColumnHolder<T>> &destinationColumn,
                      const std::optional<ColumnHolder<T>> &sourceColumn) {
        if (ACTS_CHECK_BIT(columnsToCopy, mask)) {
          assert(destinationColumn.has_value() &&
                 "Column is not available in destination container");
          assert(sourceColumn.has_value() &&
                 "Column is not available in source container");
          destinationColumn->proxy(*this)[index] =
              sourceColumn->proxy(sourceContainer)[sourceIndex];
        }
      };

  [&]<std::size_t... Is>(std::index_sequence<Is...>) {
    ((copyColumn(std::get<Is>(knownColumnMasks()), std::get<Is>(knownColumns()),
                 std::get<Is>(sourceContainer.knownColumns()))),
     ...);
  }(tuple_indices<decltype(knownColumns())>{});

  if (ACTS_CHECK_BIT(columnsToCopy, SpacePointColumns::SourceLinks)) {
    assignSourceLinks(index, sourceContainer.sourceLinks(sourceIndex));
  }
}

void SpacePointContainer2::assignSourceLinks(
    Index index, std::span<const SourceLink> sourceLinks) {
  if (index >= size()) {
    throw std::out_of_range("Index out of range in SpacePointContainer2");
  }
  if (!m_sourceLinkOffsetColumn.has_value() ||
      !m_sourceLinkCountColumn.has_value()) {
    throw std::logic_error("No source links column available");
  }
  if (m_sourceLinkCountColumn->proxy(*this)[index] != 0) {
    throw std::logic_error("Source links already assigned to the space point");
  }

  m_sourceLinkOffsetColumn->proxy(*this)[index] =
      static_cast<SpacePointIndex2>(m_sourceLinks.size());
  m_sourceLinkCountColumn->proxy(*this)[index] =
      static_cast<std::uint8_t>(sourceLinks.size());
  m_sourceLinks.insert(m_sourceLinks.end(), sourceLinks.begin(),
                       sourceLinks.end());
}

void SpacePointContainer2::createColumns(SpacePointColumns columns) noexcept {
  using enum SpacePointColumns;

  const auto createColumn =
      [&]<typename T>(SpacePointColumns mask, std::string_view name,
                      T defaultValue, std::optional<ColumnHolder<T>> &column) {
        if (ACTS_CHECK_BIT(columns, mask) && !column.has_value()) {
          column = ColumnHolder<T>(std::move(defaultValue));
          column->resize(size());
          m_allColumns.try_emplace(std::string(name), &column.value());
          m_knownColumns = m_knownColumns | mask;
        }
      };

  [&]<std::size_t... Is>(std::index_sequence<Is...>) {
    ((createColumn(
         std::get<Is>(knownColumnMasks()), std::get<Is>(knownColumnNames()),
         std::get<Is>(knownColumnDefaults()), std::get<Is>(knownColumns()))),
     ...);
  }(tuple_indices<decltype(knownColumns())>{});
}

void SpacePointContainer2::dropColumns(SpacePointColumns columns) noexcept {
  using enum SpacePointColumns;

  const auto dropColumn = [&]<typename T>(
                              SpacePointColumns mask, std::string_view name,
                              std::optional<ColumnHolder<T>> &column) {
    if (ACTS_CHECK_BIT(columns, mask) && column.has_value()) {
      m_allColumns.erase(std::string(name));
      column.reset();
      m_knownColumns = m_knownColumns & ~mask;
    }
  };

  [&]<std::size_t... Is>(std::index_sequence<Is...>) {
    ((dropColumn(std::get<Is>(knownColumnMasks()),
                 std::get<Is>(knownColumnNames()),
                 std::get<Is>(knownColumns()))),
     ...);
  }(tuple_indices<decltype(knownColumns())>{});

  if (ACTS_CHECK_BIT(columns, SpacePointColumns::SourceLinks)) {
    m_sourceLinks.clear();
  }
}

void SpacePointContainer2::dropColumn(const std::string &name) {
  if (reservedColumn(name)) {
    throw std::runtime_error("Cannot drop reserved column: " + name);
  }

  auto it = m_allColumns.find(name);
  if (it == m_allColumns.end()) {
    throw std::runtime_error("Column does not exist: " + name);
  }

  m_allColumns.erase(it);
  m_dynamicColumns.erase(name);
}

bool SpacePointContainer2::reservedColumn(const std::string &name) noexcept {
  static const auto reservedColumns = std::apply(
      [](auto... reservedNames) {
        return std::unordered_set<std::string, std::hash<std::string_view>,
                                  std::equal_to<>>({reservedNames...});
      },
      knownColumnNames());

  return reservedColumns.contains(name);
}

}  // namespace Acts
