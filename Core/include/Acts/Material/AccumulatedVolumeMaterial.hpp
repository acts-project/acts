// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <limits>

#include "Acts/Material/MaterialProperties.hpp"

namespace Acts {

/// Accumulate and average volume-based material properties.
///
/// This class is intended to be used during the mapping process.
class AccumulatedVolumeMaterial {
 public:
  /// Add one entry with the given material properties.
  void accumulate(const MaterialProperties& mat);

  /// Compute the average material collected so far.
  ///
  /// @returns Vacuum properties if no matter has been accumulated yet.
  Material average();

 private:
  double m_totalX0{std::numeric_limits<double>::infinity()};
  double m_totalL0{std::numeric_limits<double>::infinity()};
  double m_totalAr{0.};
  double m_totalZ{0.};
  double m_totalRho{0.};
  double m_thickness{1.};
  unsigned int m_materialEntries{0};
};

}  // namespace Acts
