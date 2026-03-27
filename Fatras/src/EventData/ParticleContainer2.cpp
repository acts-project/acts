// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsFatras/EventData/ParticleContainer2.hpp"

#include "Acts/Utilities/Helpers.hpp"

#include <unordered_set>

namespace {

template <typename Tuple>
using tuple_indices =
    std::make_index_sequence<std::tuple_size_v<std::remove_reference_t<Tuple>>>;

}

namespace ActsFatras {

static_assert(std::random_access_iterator<ParticleContainer2::iterator>);
static_assert(std::random_access_iterator<ParticleContainer2::const_iterator>);
static_assert(
    std::random_access_iterator<ParticleContainer2::MutableSubset::iterator>);
static_assert(
    std::random_access_iterator<ParticleContainer2::ConstSubset::iterator>);

ParticleContainer2::ParticleContainer2(ParticleColumns2 columns) noexcept {
  createColumns(columns);
}

ParticleContainer2::ParticleContainer2(const ParticleContainer2 &other) noexcept
    : m_size(other.m_size), m_parentIndices(other.m_parentIndices) {
  copyColumns(other);
}

ParticleContainer2::ParticleContainer2(ParticleContainer2 &&other) noexcept
    : m_size(other.m_size), m_parentIndices(std::move(other.m_parentIndices)) {
  moveColumns(other);

  other.m_size = 0;
}

ParticleContainer2 &ParticleContainer2::operator=(
    const ParticleContainer2 &other) noexcept {
  if (this == &other) {
    return *this;
  }

  copyColumns(other);
  m_parentIndices = other.m_parentIndices;
  m_size = other.m_size;

  return *this;
}

ParticleContainer2 &ParticleContainer2::operator=(
    ParticleContainer2 &&other) noexcept {
  if (this == &other) {
    return *this;
  }

  moveColumns(other);
  m_parentIndices = std::move(other.m_parentIndices);
  m_size = other.m_size;

  other.m_size = 0;

  return *this;
}

void ParticleContainer2::copyColumns(const ParticleContainer2 &other) {
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

void ParticleContainer2::moveColumns(ParticleContainer2 &other) noexcept {
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
  other.m_knownColumns = ParticleColumns2::None;
  other.m_dynamicColumns.clear();
}

void ParticleContainer2::reserve(std::uint32_t size,
                                 float averageParentIndices) noexcept {
  if (hasColumns(ParticleColumns2::Parents)) {
    m_parentIndices.reserve(
        static_cast<std::uint32_t>(size * averageParentIndices));
  }

  for (const auto &[name, column] : m_allColumns) {
    column->reserve(size);
  }
}

void ParticleContainer2::clear() noexcept {
  m_size = 0;
  m_parentIndices.clear();

  for (const auto &[name, column] : m_allColumns) {
    column->clear();
  }
}

MutableParticleProxy2 ParticleContainer2::createParticle() noexcept {
  ++m_size;

  for (const auto &[name, column] : m_allColumns) {
    column->emplace_back();
  }

  return MutableProxy(*this, size() - 1);
}

void ParticleContainer2::createColumns(ParticleColumns2 columns) noexcept {
  using enum ParticleColumns2;

  const auto createColumn =
      [&]<typename T>(ParticleColumns2 mask, std::string_view name,
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

void ParticleContainer2::dropColumns(ParticleColumns2 columns) noexcept {
  using enum ParticleColumns2;

  const auto dropColumn = [&]<typename T>(
                              ParticleColumns2 mask, std::string_view name,
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

  if (ACTS_CHECK_BIT(columns, ParticleColumns2::Parents)) {
    m_parentIndices.clear();
  }
}

void ParticleContainer2::dropColumn(const std::string &name) {
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

bool ParticleContainer2::reservedColumn(const std::string &name) noexcept {
  static const auto reservedColumns = std::apply(
      [](auto... reservedNames) {
        return std::unordered_set<std::string, std::hash<std::string_view>,
                                  std::equal_to<>>({reservedNames...});
      },
      knownColumnNames());

  return reservedColumns.contains(name);
}

}  // namespace ActsFatras
