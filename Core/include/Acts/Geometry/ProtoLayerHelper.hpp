// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/ProtoLayer.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Units.hpp"

namespace Acts {

/// @class ProtoLayerHelper
///
/// This class is designed to parse a vector of Surfaces and sort them
/// into corresponding proto layers.
///
/// @todo write more documentation on how this is done
class ProtoLayerHelper {
 public:
  /// Nested configuration struct
  struct Config {};

  /// Constructor with explicit config
  ///
  /// @param cfg Explicit config struct
  /// @param logger logging instance
  ProtoLayerHelper(const Config& cfg,
                   std::unique_ptr<const Logger> logger =
                       getDefaultLogger("ProtoLayerHelper", Logging::INFO))
      : m_cfg(cfg), m_logger(std::move(logger)) {}

  /// Destructor
  virtual ~ProtoLayerHelper() = default;

  /// Sort the surfaces into ProtoLayers
  ///
  /// @param gctx The geometry context (usually building context at this stage)
  /// @param surfaces The surfaces to be sorted into arrays
  /// @param bValue The binning value for the sorting
  /// @param joinTolerance The tolerance for which bins are joined
  ///
  /// @return A vector of ProtoLayers
  std::vector<ProtoLayer> protoLayers(
      const GeometryContext& gctx, const std::vector<const Surface*>& surfaces,
      BinningValue bValue, double joinTolerance) const;

 private:
  /// Configuration struct
  Config m_cfg;

  /// Logging instance
  std::unique_ptr<const Logger> m_logger;
};

}  // namespace Acts