// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsFatras/Kernel/detail/RandomNumberDistributions.hpp"

ActsFatras::LandauDist::param_type::param_type(double mean_, double scale_)
    : mean(mean_), scale(scale_) {}

bool ActsFatras::LandauDist::param_type::operator==(
    const param_type &other) const {
  return (mean == other.mean) && (scale == other.scale);
}

ActsFatras::LandauDist::LandauDist(double mean, double scale)
    : m_cfg(mean, scale) {}

ActsFatras::LandauDist::LandauDist(const param_type &cfg) : m_cfg(cfg) {}

ActsFatras::LandauDist::result_type ActsFatras::LandauDist::min() const {
  return -std::numeric_limits<double>::infinity();
}

ActsFatras::LandauDist::result_type ActsFatras::LandauDist::max() const {
  return std::numeric_limits<double>::infinity();
}

bool ActsFatras::LandauDist::operator==(const LandauDist &other) const {
  return (m_cfg == other.m_cfg);
}
