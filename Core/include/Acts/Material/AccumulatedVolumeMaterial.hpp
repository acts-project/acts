// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <limits>

#include "Acts/Material/Material.hpp"

namespace Acts {

/// @class AccumulatedVolumeMaterial
///
/// This class is used by the VolumeMaterialMapper in order to
/// accumulate/collect material information during the mapping process. The
/// object represents the collection of a single grid point.
///
/// It calculates the average of the material parameters when called, and
/// returns a material with the averaged properties.
class AccumulatedVolumeMaterial {
 public:
  AccumulatedVolumeMaterial() = default;
  ~AccumulatedVolumeMaterial() = default;

  /// @brief This function collects the classification values of a material in
  /// the container
  ///
  /// @param [in] mat Material that will be collected
  void accumulate(const Material& mat);

  /// @brief Total average of each material classification value stored in the
  /// object independently
  ///
  /// @return Material consisting of the averaged values
  Material average();

 private:
  float m_totalX0{std::numeric_limits<float>::infinity()};
  float m_totalL0{std::numeric_limits<float>::infinity()};
  float m_totalAr{0.};
  float m_totalZ{0.};
  float m_totalRho{0.};
  unsigned int m_materialEntries{0};
  unsigned int m_vacuumEntries{0};
};

}  // namespace Acts
