// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <vector>

#include "Acts/Utilities/BinUtility.hpp"
#include "ActsExamples/Digitization/MeasurementCreation.hpp"
#include "ActsExamples/EventData/SimHit.hpp"


namespace ActsExamples {
class ModuleClusters {
public:
  using simhit_t = SimHitContainer::size_type;

  ModuleClusters(Acts::BinUtility segmentation, std::vector<Acts::BoundIndices> geoIndices, bool merge = false) :
    m_segmentation(std::move(segmentation)),
    m_geoIndices(std::move(geoIndices)),
    m_merge(merge)
    {}

  void add(DigitizedParameters params, simhit_t simhit);
  void merge() { std::cerr << "NOT IMPLEMENTED!" << std::endl;}

  std::vector<std::pair<DigitizedParameters, std::vector<simhit_t>>> digitizedParameters();

private:

  struct ModuleValue {
    std::vector<Acts::BoundIndices> paramIndices = {};
    std::vector<Acts::ActsScalar> paramValues = {};
    std::vector<Acts::ActsScalar> paramVariances = {};
    std::variant<Cluster, Cluster::Cell> value;
    std::vector<simhit_t> sources = {};
  };


  Acts::BinUtility m_segmentation;
  std::vector<Acts::BoundIndices> m_geoIndices;
  bool m_merge;
  std::vector<ModuleValue> m_moduleValues;

};
}
