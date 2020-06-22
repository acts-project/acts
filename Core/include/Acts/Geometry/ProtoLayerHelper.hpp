// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <map>

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
  /// Prescription how to sort the surfaces onto proto layers + tolerance
  using SortingConfig = std::pair<BinningValue, double>;
  /// Prescription how to bin surfaces on one proto layer + tolerance
  using BinningConfig = std::pair<BinningValue, double>;

  using BinnedRange = std::pair<double, double>;

  struct Config {};

  /// Constructor with explicit config
  ///
  /// @param cfg Explicit config struct
  /// @param logger logging instance
  ProtoLayerHelper(const Config& cfg,
                   std::unique_ptr<const Logger> logger =
                       getDefaultLogger("ProtoLayerHelper", Logging::INFO))
      : m_cfg(cfg), m_logger(std::move(logger)) {}

  virtual ~ProtoLayerHelper() = default;

  /// Sort the surfaces into ProtoLayers
  ///
  /// @param gctx The geometry context (usually building context at this stage)
  /// @param surfaces The surfaces to be sorted into arrays
  /// @param sorting The sorting setup, one single sorting
  /// @param binning The binning setup, if emtpy, no binning attempts done
  ///
  /// @return A vector of ProtoLayers
  std::vector<ProtoLayer> protoLayers(
      const GeometryContext& gctx, const std::vector<const Surface*>& surfaces,
      const SortingConfig& sorting,
      const std::vector<BinningConfig>& binning = {}) const;

  /// Sort the surfaces into ProtoLayers, sequential sorting
  ///
  /// @param gctx The geometry context (usually building context at this stage)
  /// @param surfaces The surfaces to be sorted into arrays
  /// @param sortings The sequential sorting setup, one single sorting
  /// @param binning The binning setup, if emtpy, no binning attempts done
  ///
  /// @return A vector of ProtoLayers
  std::vector<ProtoLayer> protoLayers(
      const GeometryContext& gctx, const std::vector<const Surface*>& surfaces,
      const std::vector<SortingConfig>& sortings,
      const std::vector<BinningConfig>& binning = {}) const;

  /// Method to estimate the binnings of the ProtoLayers
  ///
  /// @param gctx The geometry context (usually building context at this stage)
  /// @param surfaces The sorted surfaces to be put into a ProtoLayer
  /// @param binning The binning setup
  std::map<BinningValue, ProtoLayer::BinningRange> gridBinning(
      const GeometryContext& gctx, const std::vector<const Surface*> surfaces,
      const std::vector<BinningConfig>& binning) const;

 private:
  /// Configuration struct
  Config m_cfg;

  /// Logging instance
  std::unique_ptr<const Logger> m_logger;

  /// Private access to logger
  const Logger& logger() const { return *m_logger; }
};

}  // namespace Acts
