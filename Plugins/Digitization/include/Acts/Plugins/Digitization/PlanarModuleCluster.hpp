// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

/// Set the identifier PLUGIN
#ifdef ACTS_CORE_IDENTIFIER_PLUGIN
#include ACTS_CORE_IDENTIFIER_PLUGIN
#else
using Identifier = unsigned long long;
#endif

#include "Acts/EventData/Measurement.hpp"
#include "Acts/Plugins/Digitization/DigitizationCell.hpp"
#include "Acts/Plugins/Digitization/DigitizationModule.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"

namespace Acts {

template <ParID_t... params>
using Measurement_t = Measurement<Identifier, params...>;

class PlanarModuleCluster : public Measurement_t<ParDef::eLOC_0, ParDef::eLOC_1>
{
public:
  /// Constructor from DigitizationCells
  ///
  /// @param [in] mSurface is the module surface
  /// @param [in] cIdentifier is the channel identifier of the local position
  /// @param [in] cov is the covariance matrix
  /// @param [in] loc0 is the local position in the first coordinate
  /// @param [in] loc1 is the local position in the second coordinate
  /// @param [in] dCells is the vector of digitization cells
  PlanarModuleCluster(std::shared_ptr<const Surface> mSurface,
                      const Identifier&              cIdentifier,
                      ActsSymMatrixD<2>              cov,
                      double                         loc0,
                      double                         loc1,
                      std::vector<DigitizationCell>  dCells,
                      const DigitizationModule*      dModule = nullptr)
    : Measurement_t<ParDef::eLOC_0, ParDef::eLOC_1>(std::move(mSurface),
                                                    cIdentifier,
                                                    std::move(cov),
                                                    loc0,
                                                    loc1)
    , m_digitizationCells(std::move(dCells))
    , m_digitizationModule(dModule)
  {
  }

  /// access to the digitization cells
  ///
  /// @return the vector to the digitization cells
  const std::vector<DigitizationCell>&
  digitizationCells() const;

  /// access to the digitization module
  ///
  /// @return the pointer to the digitization module
  const DigitizationModule*
  digitizationModule() const;

private:
  std::vector<DigitizationCell> m_digitizationCells;  /// the digitization cells
  const DigitizationModule* m_digitizationModule;  /// the digitization module
};

inline const std::vector<DigitizationCell>&
PlanarModuleCluster::digitizationCells() const
{
  return m_digitizationCells;
}

inline const DigitizationModule*
PlanarModuleCluster::digitizationModule() const
{
  return m_digitizationModule;
}
}
