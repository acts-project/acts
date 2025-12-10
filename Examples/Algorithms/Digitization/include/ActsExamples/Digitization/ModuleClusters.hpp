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
#include "Acts/Utilities/BinUtility.hpp"
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

struct ModuleValue {
  std::vector<Acts::BoundIndices> paramIndices = {};
  std::vector<double> paramValues = {};
  std::vector<double> paramVariances = {};
  std::variant<Cluster, Cluster::Cell> value;
  std::set<SimHitContainer::size_type> sources = {};
  Acts::Ccl::Label label = {Acts::Ccl::NO_LABEL};
};

class ModuleClusters {
 public:
  using simhit_t = SimHitContainer::size_type;

  ModuleClusters(Acts::BinUtility segmentation,
                 std::vector<Acts::BoundIndices> geoIndices, bool merge,
                 double nsigma, bool commonCorner, bool alt = false)
      : m_segmentation(std::move(segmentation)),
        m_geoIndices(std::move(geoIndices)),
        m_merge(merge),
        m_nsigma(nsigma),
        m_commonCorner(commonCorner),
        m_useInPlaceClusterization(alt) {}

  void add(DigitizedParameters params, simhit_t simhit);
  std::vector<std::pair<DigitizedParameters, std::set<simhit_t>>>
  digitizedParameters();

 private:
  Acts::BinUtility m_segmentation;
  std::vector<Acts::BoundIndices> m_geoIndices;
  std::vector<ModuleValue> m_moduleValues;
  bool m_merge;
  double m_nsigma;
  bool m_commonCorner;
  bool m_useInPlaceClusterization;

  std::vector<ModuleValue> createCellCollection();
  void merge();
  void mergeWithInPlaceClusterization();
  template <typename index_t, unsigned int AXIS>
  void mergeWithInPlaceClusterizationImpl();

  template <class cell_list_t>
  ModuleValue squash(cell_list_t& values) const;
  std::vector<std::size_t> nonGeoEntries(
      const std::vector<Acts::BoundIndices>& indices) const;
  template <class cell_list_t>
  std::vector<std::vector<ModuleValue>> mergeParameters(
      cell_list_t& values) const;
};

}  // namespace ActsExamples
