// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Digitization/Segmentation.hpp"

#include <cmath>
#include <vector>
namespace Acts {
namespace Test {

template <typename measurement_t>
class TestCluster {
 public:
  TestCluster(measurement_t meas) : m_measurement(meas) {}
  TestCluster(measurement_t meas,
              std::shared_ptr<const Acts::Segmentation> moduleSegmentation)
      : m_measurement(std::move(meas)),
        m_segmentation(std::move(moduleSegmentation)) {}
  measurement_t measurement() const { return m_measurement; }
  const Acts::Segmentation& segmentation() const {
    return (*m_segmentation.get());
  }

 private:
  measurement_t m_measurement;
  std::shared_ptr<const Acts::Segmentation> m_segmentation;
};

}  // namespace Test
}  // namespace Acts