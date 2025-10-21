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
#include <limits>
#include <ranges>
#include <span>
#include <type_traits>

namespace InPlaceClusterization {

// For cells that should be clusterized the following functions
// need to be implemented:
namespace traits {
// @brief Get the number of dimensions of the grid the cell lives on
template <typename T_Cell>
constexpr auto getCellDimension() {
  return std::tuple_size<decltype(T_Cell::m_coordinates)>();
}

/// @brief Get the coordinates of a cell
/// @return An array of the cell dimension containing the coordinates of the cell for all grid axes.
template <typename T_Cell, typename T_Axis>
inline auto getCellCoordinate(const T_Cell &a, T_Axis axis_i) {
  return a.m_coordinates[axis_i];
}

/// @brief Set a label for a given cell
/// the label type must fit numbers as high as the number of cells which
/// are clustered.
template <typename T_Cell, typename T_Index>
inline void setLabel(T_Cell &a, T_Index label) {
  a.m_label = label;
}

/// @brief Get the label associated to the given cell
template <typename T_Cell>
inline auto getLabel(const T_Cell &a) {
  return a.m_label;
}
}  // namespace traits

/// @brief A possible Cell type, that can be clustered
/// T_Coordinates should be a signed or unsigned integer
/// T_Index should be an unsigned integer
template <typename T_Coordinates, std::size_t N, typename T_Index>
struct Cell {
  Cell(const std::array<T_Coordinates, N> &coordinates, T_Index src_index)
      : m_coordinates(coordinates), m_label(T_Index{}), m_srcIndex(src_index) {}

  std::array<T_Coordinates, N>
      m_coordinates;   ///< the coordinates of the cell on the regular grid
  T_Index m_label;     ///< a label which will be assigned by the clusterization
  T_Index m_srcIndex;  ///< the index to find the source cell
};

/// @brief concept of a call object that can be clustered.
// @TODO also require decltype(getLabel(cell)) is unsigned integer
template <typename T_Cell>
concept CellWithLabel = requires(const T_Cell &a) {
  traits::getLabel(a);
} && requires(const T_Cell &a, decltype(traits::getLabel(a)) idx) {
  traits::setLabel(a, idx);
} && requires(const T_Cell &a, unsigned int axis_i) {
  traits::getCellCoordinate(a, axis_i);
} && requires { traits::getCellDimension<T_Cell>(); };

/// @brief base concept for the container that can be used for a cell collection
// @TODO require that container can be sorted.
template <typename T_Container>
concept SequenceContainer = requires(T_Container &cont, unsigned int index) {
  { cont[index] } -> std::convertible_to<typename T_Container::value_type>;
} && requires(T_Container &cont) {
  { cont.empty() } -> std::convertible_to<bool>;
} && requires(T_Container &cont) {
  { cont.size() } -> std::convertible_to<std::size_t>;
};

/// @brief concept of a cell container that can be clustered
template <typename T_CellCollection>
concept CellCollection = SequenceContainer<T_CellCollection> &&
                         CellWithLabel<typename T_CellCollection::value_type>;

/// @brief test whether cells are connected assuming 8 fold connectivity
/// @param coordinates_diff absolute values of the differences of the cell coordinates
/// Test whether two cells with the given coordinate differences are connected
// either horizontally, vertically, diagonally or are equal.
template <typename T_Coordinates>
inline bool isConnected8(const T_Coordinates &coordinates_diff) {
  bool connected = true;
  for (const auto &a_coordinate_diff : coordinates_diff) {
    connected &= a_coordinate_diff <= 1;
  }
  return connected;
}

/// @brief test whether cells are connected assuming 4 fold connectivity
/// @param coordinates_diff absolute values of the differences of the cell coordinates
/// Test whether two cells with the given coordinate differences are connected
// either horizontally, vertically or are equal.
template <typename T_Coordinates>
inline bool isConnected4(const T_Coordinates &coordinates_diff) {
  int connections = 0;
  for (const auto &a_coordinate_diff : coordinates_diff) {
    connections += a_coordinate_diff;
  }
  return connections <= 1;
}

/// @brief test whether cells are connected assuming 4 or 8 fold connectivity
enum EFold { k8FoldConnection, k4FoldConnection };
template <EFold FOLD, typename T_Coordinates>
inline bool isConnected(const T_Coordinates &coordinates_diff) {
  if constexpr (FOLD == k8FoldConnection) {
    return isConnected8(coordinates_diff);
  } else {
    return isConnected4(coordinates_diff);
  }
}

/// @brief Label cells which are in ascending order of the coordinate of the given axis
/// @param cells the sorted cell collection that should be labeled.
template <unsigned int SORT_AXIS, CellCollection T_CellCollection,
          typename T_Index = unsigned int, EFold FOLD>
inline void labelSortedCells(T_CellCollection &cells) {
  using T_Cell = std::remove_cvref_t<decltype(cells[0])>;
  static constexpr std::size_t NDim = traits::getCellDimension<T_Cell>();
  static_assert(NDim > 0);
  using T_Label = std::remove_cvref_t<decltype(traits::getLabel(cells[0]))>;
  using T_Coordinate =
      std::remove_cvref_t<decltype(traits::getCellCoordinate(cells[0], 0))>;
  static_assert(std::numeric_limits<T_Index>::max() >=
                std::numeric_limits<T_Label>::max());
  // cells are sorted in the first coordinate
  // thus can stop searching for adjacent cells
  // if distance in the sorted coordinate is too large
  for (T_Index idx_a = 0; idx_a < cells.size(); ++idx_a) {
    traits::setLabel(cells[idx_a], idx_a);
    for (T_Index idx_b = idx_a; idx_b-- > 0;) {
      std::array<T_Coordinate, NDim> diff;
      for (unsigned int axis_i = 0; axis_i < NDim; ++axis_i) {
        diff[axis_i] =
            std::abs(traits::getCellCoordinate(cells[idx_a], axis_i) -
                     traits::getCellCoordinate(cells[idx_b], axis_i));
      }

      if (isConnected<FOLD>(diff)) {
        if (cells[idx_a].m_label < idx_a) {
          // if the label is not its index, the cell was merged to a cluster
          // thus the two clusters need to be merged rather than just merging
          // in this one cell
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
          // nothing will be done if the min and the max label are identical
          if (min_label != max_label) {
            // can only encounter cells with label max_label down to index
            // max_label
            for (T_Index idx = idx_a; idx-- > max_label;) {
              if (traits::getLabel(cells[idx]) == max_label) {
                traits::setLabel(cells[idx], min_label);
              }
            }
          }
        } else {
          traits::setLabel(cells[idx_a], traits::getLabel(cells[idx_b]));
        }
      } else if (diff[SORT_AXIS] > 1) {
        // difference too large in sorted coordinate, there won't be any more
        // candidates for merging. Can abort search for cell idx_a
        idx_b = 0;
        break;
      }
    }
  }
}

/// @brief Sort the cells in ascending order of the coordinate of the specified axis
template <unsigned int AXIS, CellCollection T_CellCollection,
          std::size_t NDim = 2, typename T_Index = unsigned int>
  requires(NDim > 0 && AXIS < NDim)
void sortCellsByCoordinate(T_CellCollection &cells) {
  using T_Cell = std::remove_cvref_t<decltype(cells[0])>;
  std::sort(cells.begin(), cells.end(), [](const T_Cell &a, const T_Cell &b) {
    return traits::getCellCoordinate(a, AXIS) <
           traits::getCellCoordinate(b, AXIS);
  });
}

/// @brief Sort the cells in ascending order of the asscoiated label.
template <CellCollection T_CellCollection>
void groupCellsByLabel(T_CellCollection &cells) {
  using T_Cell = std::remove_cvref_t<decltype(cells[0])>;
  // stable sort to retain coordinate ordering
  std::stable_sort(cells.begin(), cells.end(),
                   [](const T_Cell &a, const T_Cell &b) {
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
/// have the same label. T_Index must be large enough to fit the number of cells
/// in the collection.
template <unsigned int SORT_AXIS, CellCollection T_CellCollection,
          typename T_Index = unsigned int, EFold FOLD = k8FoldConnection>
inline void clusterize(T_CellCollection &cells) {
  assert(cells.size() <= std::numeric_limits<T_Index>::max());
  sortCellsByCoordinate<SORT_AXIS>(cells);
  labelSortedCells<SORT_AXIS, T_CellCollection, T_Index, FOLD>(cells);
  groupCellsByLabel(cells);
}

/// @brief determine the number of clusters.
/// It is assumed that @ref clusterize was called previously for the given
/// cell collection i.e. the cells are labeled and sorted by label.
/// It will count the number of different labels do determine the number of
/// clusters in the collection.
template <CellCollection T_CellCollection, typename T_Index = unsigned int>
inline std::size_t countLabels(const T_CellCollection &cells) {
  assert(cells.size() <= std::numeric_limits<T_Index>::max());
  std::size_t nlabels = !cells.empty();
  for (T_Index idx_a = 1; idx_a < cells.size(); ++idx_a) {
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
template <CellCollection T_CellCollection, typename T_Func>
inline void for_each_cluster(T_CellCollection &cells, T_Func func) {
  if (!cells.empty()) {
    using T_Index = decltype(traits::getLabel(cells[0]));
    T_Index idx_begin = 0;
    T_Index idx = 1;
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
template <CellCollection T_CellCollection, typename T_RangeCollection>
inline void addCellRanges(const T_CellCollection &cells,
                          T_RangeCollection &ranges) {
  for_each_cluster(cells,
                   [&ranges]([[maybe_unused]] const T_CellCollection &all_cells,
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
template <unsigned int AXIS, CellCollection T_CellCollection,
          typename T_RangeCollection>
inline void addCellRangesAndSort(T_CellCollection &cells,
                                 T_RangeCollection &ranges) {
  for_each_cluster(cells, [&ranges](T_CellCollection &all_cells,
                                    unsigned int idx_begin,
                                    unsigned int idx_end) {
    auto cluster_range =
        std::span(all_cells.begin() + idx_begin, all_cells.begin() + idx_end);
    using T_Cell = std::remove_cvref_t<decltype(cells.front())>;
    std::ranges::sort(cluster_range, [](const T_Cell &a, const T_Cell &b) {
      return traits::getCellCoordinate(a, AXIS) <
             traits::getCellCoordinate(b, AXIS);
    });
    ranges.emplace_back(idx_begin, idx_end);
  });
}

}  // namespace InPlaceClusterization
