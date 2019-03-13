// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// AccumulatedVolumeMaterial.h, Acts project
///////////////////////////////////////////////////////////////////

#pragma once

#include "Acts/Material/Material.hpp"

namespace Acts {
/// @class AccumulatedVolumeMaterial
///
/// This class is used by the VolumeMaterialMapper in order to
/// accumulate/collect material information during the mapping process. The
/// object represents the collection of a single grid point.
///
/// It calculates the average of the material classification values when called,
/// and returns a material with the averaged properties.
class AccumulatedVolumeMaterial
{
public:
  /// @brief Default constructor
  AccumulatedVolumeMaterial() = default;

  /// @brief Default destructor
  ~AccumulatedVolumeMaterial() = default;

  /// @brief This function collects the classification values of a material in
  /// the container
  ///
  /// @param [in] mat Material that will be collected
  void
  accumulate(const Material& mat);

  /// @brief Total average of each material classification value stored in the
  /// object independently
  ///
  /// @return Material consisting of the averaged values
  Material
  average();

private:
  float m_totalX0{std::numeric_limits<float>::infinity()};  //!< accumulate the
                                                            //! contribution to
  //! X0
  float m_totalL0{std::numeric_limits<float>::infinity()};  //!< accumulate the
                                                            //! contribution to
  //! L0
  float m_totalA{0.};    //!< accumulate the contribution to A
  float m_totalZ{0.};    //!< accumulate the contribution to Z
  float m_totalRho{0.};  //!< accumulate the contribution to rho

  unsigned int m_materialEntries{0};  //!< the number of material points
  unsigned int m_vacuumEntries{0};    //!< the number of vacuum points
};
}  // namespace Acts