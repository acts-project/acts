// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Utilities to clusterize a collection of cells on a regular grid
// by sorting them in such a way that adjacent cells with the same
// label belong to the same cluster, where a regular grid is a grid
// on which cells have exactly 8 neighbors (except for edge cells),
// and the indices of adjacent cells differ by one in the corresponding
// direction.

#pragma once
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <concepts>
#include <limits>
#include <ranges>
#include <span>
#include <type_traits>

namespace Acts::InPlaceClusterization {

// For cells that should be clusterized the following functions
// need to be implemented:
namespace traits {
// @brief Get the number of dimensions of the grid the cell lives on
template <typename cell_t>
constexpr auto getCellDimension() {
  return std::tuple_size<decltype(cell_t::m_coordinates)>();
}

/// @brief Get the coordinates of a cell
/// @param a the cell
/// @param axis_i the index of the axis (0 is the first axis).
/// @return the coordinate of the given cell of the given axis
/// The result is expected to be a signed integer. Unsigned integers should also
/// work. Float types may work but due to limited precision connections may not
/// be identified correctly in all cases. To correctly handle float types the
/// functions @ref isConnected and @ref canAbortSearch may need to be specialized.
template <typename cell_t>
auto getCellCoordinate(const cell_t &a, unsigned int axis_i) {
  return a.m_coordinates[axis_i];
}

/// @brief Set a label for a given cell
/// the label type must fit numbers as high as the number of cells which
/// are clustered.
template <typename cell_t, std::unsigned_integral index_t>
void setLabel(cell_t &a, index_t label) {
  a.m_label = label;
}

/// @brief Get the label associated to the given cell
template <typename cell_t>
auto getLabel(const cell_t &a) {
  return a.m_label;
}
}  // namespace traits

/// @brief A possible Cell type, that can be clustered
/// coordinates_t are expected to be signed integers
/// index_t should be an unsigned integer
/// Unsigned should work, float coordinate types may cause problems. See
/// getCellCoordinate.
template <typename coordinates_t, std::size_t N, std::unsigned_integral index_t>
struct Cell {
  Cell(const std::array<coordinates_t, N> &coordinates, index_t src_index)
      : m_coordinates(coordinates), m_label(index_t{}), m_srcIndex(src_index) {}

  std::array<coordinates_t, N>
      m_coordinates;   ///< the coordinates of the cell on the regular grid
  index_t m_label;     ///< a label which will be assigned by the clusterization
  index_t m_srcIndex;  ///< the index to find the source cell
};

/// @brief concept of a call object that can be clustered.
template <typename cell_t>
concept CellWithLabel = requires(const cell_t &a, const cell_t &b,
                                 decltype(traits::getLabel(a)) idx,
                                 unsigned int axis_i) {
  { traits::getLabel(a) } -> std::unsigned_integral;
  traits::setLabel(a, idx);
  {
    traits::getCellCoordinate(a, axis_i) < traits::getCellCoordinate(b, axis_i)
  } -> std::convertible_to<bool>;
  { traits::getCellDimension<cell_t>() } -> std::convertible_to<std::size_t>;
};

/// @brief base concept for the container that can be used for a cell collection
template <typename container_t>
concept SequenceContainer = requires(container_t &cont, unsigned int index) {
  { cont[index] } -> std::convertible_to<typename container_t::value_type>;
  { cont.empty() } -> std::convertible_to<bool>;
  { cont.size() } -> std::convertible_to<std::size_t>;
} && std::permutable<typename container_t::iterator>;

/// @brief concept of a cell container that can be clustered
template <typename cell_collection_t>
concept CellCollection = SequenceContainer<cell_collection_t> &&
                         CellWithLabel<typename cell_collection_t::value_type>;

/// @brief test whether cells are connected considering common corners and edges.
/// @param coordinates_diff absolute values of the differences of the cell coordinates
/// Test whether two cells with the given coordinate differences are connected
/// either horizontally, vertically, diagonally or are equal.
template <typename coordinates_t>
bool isConnectedCommonEdgeOrCorner(const coordinates_t &coordinates_diff) {
  bool connected = true;
  for (const auto &a_coordinate_diff : coordinates_diff) {
    connected &= a_coordinate_diff <= 1;
  }
  return connected;
}

/// @brief test whether cells are connected considering common edges only.
/// @param coordinates_diff absolute values of the differences of the cell coordinates
/// Test whether two cells with the given coordinate differences are connected
/// either horizontally, vertically or are equal.
template <typename coordinates_t>
bool isConnectedCommonEdge(const coordinates_t &coordinates_diff) {
  int connections = 0;
  for (const auto &a_coordinate_diff : coordinates_diff) {
    connections += a_coordinate_diff;
  }
  return connections <= 1;
}

/// @brief test whether cells are connected considering common edges or common corners
enum class EConnectionType { CommonEdgeOrCorner, CommonEdge };
template <EConnectionType CONNECTION_TYPE, typename coordinates_t>
bool isConnected(const coordinates_t &coordinates_diff) {
  if constexpr (CONNECTION_TYPE == EConnectionType::CommonEdgeOrCorner) {
    return isConnectedCommonEdgeOrCorner(coordinates_diff);
  } else {
    return isConnectedCommonEdge(coordinates_diff);
  }
}

/// @brief test whether the search for connections can be aborted.
/// the search is performed in descending order in one coordinate, once the
/// difference in this coordinate is too large the search can be aborted because
/// there cannot be any further connections.
template <typename coordinates_t>
bool canAbortSearch(const coordinates_t &coordinates_diff,
                    unsigned int sort_axis_i) {
  assert(sort_axis_i < coordinates_diff.size());
  return coordinates_diff[sort_axis_i] > 1;
}

/// @brief compute the absolute difference of two coordinates
template <typename coordinate_t>
auto absDifference(coordinate_t a, coordinate_t b) {
  if constexpr (std::is_signed_v<coordinate_t>) {
    return std::abs(a - b);
  } else {
    return (a > b ? a - b : b - a);
  }
}

/// @brief Label cells which are in ascending order of the coordinate of the given axis
/// @param cells the sorted cell collection that should be labeled.
template <unsigned int SORT_AXIS, CellCollection cell_collection_t,
          std::unsigned_integral index_t = unsigned int,
          EConnectionType CONNECTION_TYPE>
void labelSortedCells(cell_collection_t &cells) {
  using cell_t = std::remove_cvref_t<decltype(cells[0])>;
  static constexpr std::size_t NDim = traits::getCellDimension<cell_t>();
  static_assert(NDim > 0);
  using label_t = std::remove_cvref_t<decltype(traits::getLabel(cells[0]))>;
  using coordinate_t =
      std::remove_cvref_t<decltype(traits::getCellCoordinate(cells[0], 0u))>;
  static_assert(std::numeric_limits<index_t>::max() >=
                std::numeric_limits<label_t>::max());
  // cells are sorted in the first coordinate
  // thus can stop searching for adjacent cells
  // if distance in the sorted coordinate is too large
  for (index_t idx_a = 0; idx_a < cells.size(); ++idx_a) {
    traits::setLabel(cells[idx_a], idx_a);
    for (index_t idx_b = idx_a; idx_b-- > 0;) {
      // NOLINTBEGIN(cppcoreguidelines-pro-type-member-init)
      // all elements of diff will be set to a value no need to initialize
      // with a dummy value.
      std::array<coordinate_t, NDim> diff;
      for (unsigned int axis_i = 0; axis_i < NDim; ++axis_i) {
        diff[axis_i] =
            absDifference(traits::getCellCoordinate(cells[idx_a], axis_i),
                          traits::getCellCoordinate(cells[idx_b], axis_i));
      }
      // NOLINTEND(cppcoreguidelines-pro-type-member-init)

      if (isConnected<CONNECTION_TYPE>(diff)) {
        if (cells[idx_a].m_label < idx_a) {
          // if the label is not its index, the cell was merged to a cluster
          // thus the two clusters need to be merged rather than just merging
          // in this one cell
          // NOLINTBEGIN(cppcoreguidelines-init-variables)
          // min_label and max_label will be set before use, no reason to set
          // to unused dummy value before.
          unsigned int min_label;
          unsigned int max_label;
          if (cells[idx_a].m_label < cells[idx_b].m_label) {
            min_label = traits::getLabel(cells[idx_a]);
            max_label = traits::getLabel(cells[idx_b]);
            traits::setLabel(cells[idx_b], min_label);
          } else {
            max_label = traits::getLabel(cells[idx_a]);
            min_label = traits::getLabel(cells[idx_b]);
            traits::setLabel(cells[idx_a], min_label);
          }
          // NOLINTEND(cppcoreguidelines-init-variables)
          // nothing will be done if the min and the max label are identical
          if (min_label != max_label) {
            // can only encounter cells with label max_label down to index
            // max_label
            for (index_t idx = idx_a; idx-- > max_label;) {
              if (traits::getLabel(cells[idx]) == max_label) {
                traits::setLabel(cells[idx], min_label);
              }
            }
          }
        } else {
          traits::setLabel(cells[idx_a], traits::getLabel(cells[idx_b]));
        }
      } else if (canAbortSearch(diff, SORT_AXIS)) {
        // difference too large in sorted coordinate, there won't be any more
        // candidates for merging. Can abort search for cell idx_a
        idx_b = 0;
        break;
      }
    }
  }
}

/// @brief Sort the cells in ascending order of the coordinate of the specified axis
template <unsigned int AXIS, CellCollection cell_collection_t,
          std::size_t NDim = 2, std::unsigned_integral index_t = unsigned int>
  requires(NDim > 0 && AXIS < NDim)
void sortCellsByCoordinate(cell_collection_t &cells) {
  using cell_t = std::remove_cvref_t<decltype(cells[0])>;
  std::sort(cells.begin(), cells.end(), [](const cell_t &a, const cell_t &b) {
    return traits::getCellCoordinate(a, AXIS) <
           traits::getCellCoordinate(b, AXIS);
  });
}

/// @brief Sort the cells in ascending order of the asscoiated label.
template <CellCollection cell_collection_t>
void groupCellsByLabel(cell_collection_t &cells) {
  using cell_t = std::remove_cvref_t<decltype(cells[0])>;
  // stable sort to retain coordinate ordering
  std::stable_sort(cells.begin(), cells.end(),
                   [](const cell_t &a, const cell_t &b) {
                     return traits::getLabel(a) < traits::getLabel(b);
                   });
}

/// @brief Sort the cell collection in such a way that cells of a cluster are adjacent.
/// @param cells the cell collection which will be sorted.
/// Will sort the cells first by the coordinate of the specified axis to
/// accelerate connection tests. Then associate the same label to connected
/// cells where the label corresponds to the smallest cell index (in the sorted
/// cell collection) of a cluster, where the cells are in ascending order of one
/// of the coordinates Finally will sort the cells by the label in ascending
/// order. As a consequence cells which belong to one cluster are adjacent and
/// have the same label. index_t must be large enough to fit the number of cells
/// in the collection.
template <unsigned int SORT_AXIS, CellCollection cell_collection_t,
          std::unsigned_integral index_t = unsigned int,
          EConnectionType CONNECTION_TYPE = EConnectionType::CommonEdgeOrCorner>
void clusterize(cell_collection_t &cells) {
  assert(cells.size() <= std::numeric_limits<index_t>::max());
  sortCellsByCoordinate<SORT_AXIS>(cells);
  labelSortedCells<SORT_AXIS, cell_collection_t, index_t, CONNECTION_TYPE>(
      cells);
  groupCellsByLabel(cells);
}

/// @brief determine the number of clusters.
/// It is assumed that @ref clusterize was called previously for the given
/// cell collection i.e. the cells are labeled and sorted by label.
/// It will count the number of different labels do determine the number of
/// clusters in the collection.
template <CellCollection cell_collection_t,
          std::unsigned_integral index_t = unsigned int>
std::size_t countLabels(const cell_collection_t &cells) {
  assert(cells.size() <= std::numeric_limits<index_t>::max());
  std::size_t nlabels = !cells.empty();
  for (index_t idx_a = 1; idx_a < cells.size(); ++idx_a) {
    nlabels +=
        traits::getLabel(cells[idx_a]) != traits::getLabel(cells[idx_a - 1]);
  }
  return nlabels;
}

/// @brief call the given function for each cluster of a label sorted cell collection.
/// It is assumed that @ref clusterize was called previously for the given
/// cell collection i.e. the cells are labeled and sorted by label.
/// Will call the given function for each cluster, where the function
/// gets a reference to the full cell collection and the begin and end
/// cell index which defines the cell range of a cluster.
template <CellCollection cell_collection_t, typename func_t>
void for_each_cluster(cell_collection_t &cells, func_t func) {
  if (!cells.empty()) {
    using index_t = decltype(traits::getLabel(cells[0]));
    index_t idx_begin = 0;
    index_t idx = 1;
    for (; idx < cells.size(); ++idx) {
      if (traits::getLabel(cells[idx]) != traits::getLabel(cells[idx - 1])) {
        func(cells, idx_begin, idx);
        idx_begin = idx;
      }
    }
    if (idx_begin < idx) {
      func(cells, idx_begin, idx);
    }
  }
}

/// @brief Add element ranges for each cluster in the cell collection to the given range container
/// @param cells labeled and label sorted cell collection
/// @param ranges a range container
/// It is assumed that @ref clusterize was called previously for the given
/// cell collection i.e. the cells are labeled and sorted by label.
/// Will add one cell range per cluster to the given range container.
/// The value type of the range container has to be constructible by
/// two indices the index of the first cell in a cluster and the index
/// of the cell after the last cell of the cluster.
template <CellCollection cell_collection_t, typename range_collection_t>
void addCellRanges(const cell_collection_t &cells, range_collection_t &ranges) {
  for_each_cluster(
      cells, [&ranges]([[maybe_unused]] const cell_collection_t &all_cells,
                       unsigned int idx_begin, unsigned int idx_end) {
        ranges.emplace_back(idx_begin, idx_end);
      });
}

/// @brief Sort the cells of each cluster by one coordinate and add cell ranges to the given range container
/// @param cells labeled and label sorted cell collection
/// @param ranges a range container
/// It is assumed that @ref clusterize was called previously for the given
/// cell collection i.e. the cells are labeled and sorted by label.
/// Will first sort the cells of each cluster by the coordinate of the specified
/// axis and add one cell range per cluster to the given range container. The
/// value type of the range container has to be constructible by two indices the
/// index of the first cell in a cluster and the index of the cell after the
/// last cell of the cluster.
template <unsigned int AXIS, CellCollection cell_collection_t,
          typename range_collection_t>
void addCellRangesAndSort(cell_collection_t &cells,
                          range_collection_t &ranges) {
  for_each_cluster(cells, [&ranges](cell_collection_t &all_cells,
                                    unsigned int idx_begin,
                                    unsigned int idx_end) {
    auto cluster_range =
        std::span(all_cells.begin() + idx_begin, all_cells.begin() + idx_end);
    using cell_t = std::remove_cvref_t<decltype(cells.front())>;
    std::ranges::sort(cluster_range, [](const cell_t &a, const cell_t &b) {
      return traits::getCellCoordinate(a, AXIS) <
             traits::getCellCoordinate(b, AXIS);
    });
    ranges.emplace_back(idx_begin, idx_end);
  });
}

}  // namespace Acts::InPlaceClusterization
