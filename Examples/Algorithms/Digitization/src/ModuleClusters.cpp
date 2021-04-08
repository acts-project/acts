// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Digitization/ModuleClusters.hpp"

namespace ActsExamples {

void ModuleClusters::add(DigitizedParameters params, simhit_t simhit)
{
  ModuleValue mval;
  mval.paramIndices = std::move(params.indices);
  mval.paramValues =  std::move(params.values);
  mval.paramVariances = std::move(params.variances);
  mval.sources = {simhit};

  if (m_merge) { // Break-up the cluster
    for (auto cell : params.cluster.channels) {
      ModuleValue mval_cell = mval;
      mval_cell.value = cell;
      m_moduleValues.push_back(std::move(mval_cell));
    }
  } else { // pass-through mode
    mval.value = std::move(params.cluster);
    m_moduleValues.push_back(std::move(mval));
  }
}

std::vector<std::pair<DigitizedParameters, std::vector<ModuleClusters::simhit_t>>>
ModuleClusters::digitizedParameters()
{
  if (m_merge) { // (re-)build the clusters
    merge();
  }
  std::vector<std::pair<DigitizedParameters, std::vector<simhit_t>>> retv;
  for (ModuleValue& mval : m_moduleValues) {
    if (std::holds_alternative<Cluster::Cell>(mval.value)) {
      // Should never happen! Either the cluster should have
      // passed-through or the cells merged into clusters in the
      // merge() step. Treat this as a bug.
      throw std::runtime_error("Invalid cluster!");
    }
    DigitizedParameters dpars;
    dpars.indices = mval.paramIndices;
    dpars.values = mval.paramValues;
    dpars.variances = mval.paramVariances;
    dpars.cluster = std::get<Cluster>(mval.value);
    retv.emplace_back(std::move(dpars), mval.sources);
  }
  return retv;
}


} // namespace ActsExamples
