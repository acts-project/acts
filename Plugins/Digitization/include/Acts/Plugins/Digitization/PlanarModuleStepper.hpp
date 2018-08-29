// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include <memory>
#include <vector>
#include "Acts/Plugins/Digitization/DigitizationCell.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Logger.hpp"

namespace Acts {

class DigitizationModule;

/// @class PlanarModuleStepper
///
/// Module for fast, geometric digitization
/// this is a planar module stepper that calculates the step length
/// in given segmentations and retrunes digitisation steps

class PlanarModuleStepper
{
public:
  /// Constructor
  ///
  /// @param pmsConfig is the configuration
  /// @param mlogger is the logging istance
  PlanarModuleStepper(std::unique_ptr<const Logger> mlogger
                      = getDefaultLogger("PlanarModuleStepper", Logging::INFO));

  /// Destructor
  ~PlanarModuleStepper() = default;

  /// Calculate the steps caused by this track - full simulation interface
  ///
  /// @param dmodule is the digitization module
  /// @param startPoint is the starting position of the stepping
  /// @param endPoint is the end postion of the stepping
  ///
  /// @return is the vector of digitization steps
  std::vector<DigitizationStep>
  cellSteps(const DigitizationModule& dmodule,
            const Vector3D&           startPoint,
            const Vector3D&           endPoint) const;

  /// Calculate the steps caused by this track - fast simulation interface
  ///
  /// @param dmodule is the digitization module
  /// @param moduleIntersection is the 2d intersection at the module surface
  /// @param trackDirection is the track direction at the instersection
  ///
  /// @return is the vector of digitization steps
  std::vector<DigitizationStep>
  cellSteps(const DigitizationModule& dmodule,
            const Vector2D&           moduleIntersection,
            const Vector3D&           trackDirection) const;

  /// Set logging instance
  ///
  /// @param logger is the logging instance to be set
  void
  setLogger(std::unique_ptr<const Logger> logger)
  {
    m_logger = std::move(logger);
  }

private:
  /// Private access method to the logging instance
  const Logger&
  logger() const
  {
    return *m_logger;
  }

  /// logging instance
  std::unique_ptr<const Logger> m_logger;
};

}  // namespace
