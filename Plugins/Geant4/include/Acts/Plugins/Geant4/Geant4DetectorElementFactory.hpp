// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Common.hpp"
#include "Acts/Geometry/Extent.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Delegate.hpp"

#include <tuple>
#include <vector>

class G4VPhysicalVolume;

namespace Acts {

class Geant4DetectorElement;

/// A factory to convert Geant4 physical volumes
/// into Geant4 detector elements
///
class Geant4DetectorElementFactory {
 public:
  /// Nested configuration struct that holds
  /// global lifetime configuration
  struct Config {
    ActsScalar scaleConversion = 1.;
  };

  /// Nested cache that collects the current
  /// 
  struct Cache {};

  /// Nested option struct that allows
  /// per call changable configuration
  struct Options {};

  /// The Geant4 detector element factory
  ///
  /// @param cfg the configuration struct
  Geant4DetectorElementFactory(Config cfg);

 private:
  Config m_cfg;
};

}  // namespace Acts
