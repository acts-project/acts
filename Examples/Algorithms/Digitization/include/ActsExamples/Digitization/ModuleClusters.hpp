// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Clusterization/Clusterization.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "ActsExamples/Digitization/MeasurementCreation.hpp"
#include "ActsExamples/EventData/Cluster.hpp"
#include "ActsExamples/EventData/SimHit.hpp"

#include <cstddef>
#include <set>
#include <utility>
#include <variant>
#include <vector>

namespace ActsExamples {

struct DigitizedParameters;

/// Intermediate, mutable representation of a measurement inside a module.
///
/// Carries the bound parameters as a struct of arrays (parallel to the layout
/// of @c DigitizedParameters) together with the cluster or single cell it
/// originates from. It is the working unit during the clustering and merging
/// performed by @c ModuleClusters.
struct ModuleValue {
  /// Bound indices addressed by this value.
  std::vector<Acts::BoundIndices> paramIndices;
  /// Parameter values, one per entry in @c paramIndices.
  std::vector<double> paramValues;
  /// Parameter variances, one per entry in @c paramIndices.
  std::vector<double> paramVariances;
  /// The payload: either a full cluster or a single, not-yet-merged cell.
  std::variant<Cluster, Cluster::Cell> value;
  /// Simulated hits contributing to this value.
  std::set<SimHitIndex> sources;
  /// Connected-component label assigned during geometric clustering.
  Acts::Ccl::Label label = Acts::Ccl::NO_LABEL;
};

/// Collects the digitized values of a single module and turns them into
/// measurements, optionally merging cells into clusters.
///
/// Values are accumulated via @c add and finalized via @c digitizedParameters.
/// When merging is enabled, cells are first grouped by a connected-component
/// analysis on the module segmentation and then split again where the
/// non-geometric parameters (e.g. time) are incompatible within @c nsigma.
class ModuleClusters {
 public:
  /// Construct an empty collection for one module.
  ///
  /// @param segmentation the module segmentation used to position cells
  /// @param geoIndices the bound indices resolved from the cell geometry
  /// @param merge whether to merge cells into clusters
  /// @param nsigma compatibility window (in sigma) for merging non-geometric
  ///        parameters
  /// @param commonCorner whether cells touching only at a corner are connected
  ModuleClusters(std::vector<Acts::DirectedProtoAxis> segmentation,
                 std::vector<Acts::BoundIndices> geoIndices, bool merge,
                 double nsigma, bool commonCorner)
      : m_segmentation(std::move(segmentation)),
        m_geoIndices(std::move(geoIndices)),
        m_merge(merge),
        m_nsigma(nsigma),
        m_commonCorner(commonCorner) {}

  /// Add a digitized measurement to the module.
  ///
  /// In merging mode the cluster is broken up into its individual cells;
  /// otherwise the parameters are stored as-is (pass-through mode).
  ///
  /// @param params the digitized parameters to add
  /// @param simhit the index of the simulated hit they originate from
  void add(DigitizedParameters params, SimHitIndex simhit);

  /// Finalize the module and return the digitized measurements.
  ///
  /// Triggers the merging step if enabled and converts the internal values
  /// back into @c DigitizedParameters.
  ///
  /// @return the digitized parameters paired with their contributing sim hits
  std::vector<std::pair<DigitizedParameters, std::set<SimHitIndex>>>
  digitizedParameters();

 private:
  /// Module segmentation used to position cells.
  std::vector<Acts::DirectedProtoAxis> m_segmentation;
  /// Bound indices resolved from the cell geometry.
  std::vector<Acts::BoundIndices> m_geoIndices;
  /// Accumulated values, holding cells before and clusters after merging.
  std::vector<ModuleValue> m_moduleValues;
  /// Whether cells are merged into clusters.
  bool m_merge;
  /// Compatibility window (in sigma) for merging non-geometric parameters.
  double m_nsigma;
  /// Whether cells touching only at a corner are considered connected.
  bool m_commonCorner;

  /// Collect the unique cells, summing the activation of duplicates.
  ///
  /// @return the deduplicated cell values
  std::vector<ModuleValue> createCellCollection() const;

  /// Merge the accumulated cells into clusters in @c m_moduleValues.
  ///
  /// Performs a connected-component analysis on the geometry and then splits
  /// groups whose non-geometric parameters are incompatible.
  void merge();

  /// Squash a group of values into a single, weighted cluster value.
  ///
  /// Cell positions are averaged by activation weight and the geometric
  /// parameters are filled before the non-geometric ones to keep slot
  /// positions aligned with the source layout.
  ///
  /// @param values the values forming one cluster
  ///
  /// @throws std::runtime_error if @p values is empty
  ///
  /// @return the merged cluster value
  ModuleValue squash(std::vector<ModuleValue>& values) const;

  /// Find the entries of @p indices that are not geometric.
  ///
  /// @param indices the bound indices to inspect
  ///
  /// @return the positions into @p indices of the non-geometric entries
  std::vector<std::size_t> nonGeoEntries(
      std::vector<Acts::BoundIndices>& indices) const;

  /// Group values whose non-geometric parameters are compatible within
  /// @c nsigma.
  ///
  /// @param values the values to group
  ///
  /// @return the values partitioned into compatible groups
  std::vector<std::vector<ModuleValue>> mergeParameters(
      std::vector<ModuleValue> values) const;
};

}  // namespace ActsExamples
