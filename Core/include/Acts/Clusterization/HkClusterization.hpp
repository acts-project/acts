#pragma once

#include <memory>
#include <vector>

namespace Acts {

typedef int Label;
#define NO_LABEL 0

// Simple wrapper around a generic cell type
// The Cell type should have the following accessor functions defined:
//   int get_cell_row(const Cell&)
//   int get_cell_column(const Cell&)
template <typename Cell>
struct LabeledCell {
  Cell const* ptr;
  mutable Label lbl;
  explicit LabeledCell(const Cell& cell) : ptr{std::addressof(cell)}, lbl{0} {}
};

// Definition of the Cell api for the LabeledCell<Cell> wrapper
template <typename Cell>
int get_cell_row(const LabeledCell<Cell>& cell) {
  return get_cell_row(*cell.ptr);
}

template <typename Cell>
int get_cell_column(const LabeledCell<Cell>& cell) {
  return get_cell_column(*cell.ptr);
}

template <typename Cell>
Label& get_cell_label(LabeledCell<Cell>& lcell) {
  return lcell.lbl;
}

/// @brief labelClusters
///
/// In-place connected component labelling using the Hoshen-Kopelman algorithm.
/// The `Cell` type must have the following functions defined:
///   int  get_cell_row(const Cell&),
///   int  get_cell_column(const Cell&)
///   int& get_cell_label(Cell&)
/// If a particular Cell type does not have a label slot, use the
/// provided LabeledCell<Cell> wrapper type.
///
/// @param [in] cells the cell collection to be labeled
/// @param [in] commonCorner flag indicating if cells sharing a common corner should be merged into one cluster
/// @return nothing
template <typename Cell, typename CellCollection>
void labelClusters(CellCollection& cells, bool commonCorner);

/// @brief mergeClusters
///
/// Merge a set of cells previously labeled (for instance with `labelClusters`)
/// into actual clusters. The Cluster type must have the following function
/// defined:
///   void cluster_add_cell(Cluster&, const Cell&)
///
/// @param [in] cells the labeled cell collection
/// @return nothing
template <typename Cell, typename Cluster, typename CellCollection,
          typename ClusterCollection = std::vector<Cluster>>
ClusterCollection mergeClusters(CellCollection& cells);

/// @brief createClusters
/// Conveniance function which runs both labelClusters and createClusters.
template <typename Cell, typename Cluster, typename CellCollection,
          typename ClusterCollection = std::vector<Cluster>>
ClusterCollection createClusters(CellCollection& cells, bool commonCorner) {
  labelClusters<Cell>(cells, commonCorner);
  return mergeClusters<Cell, Cluster>(cells);
}

}  // namespace Acts

#include "Acts/Clusterization/HkClusterization.ipp"
