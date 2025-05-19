// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Plugins/DD4hep/DD4hepAlignmentStore.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Framework/IContextDecorator.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <array>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

namespace ActsExamples {
struct AlgorithmContext;

/// @brief A simple alignment decorator for the DD4hep geometry
/// allowing to load a single static alignment onto the geometry
///
/// The alignments are stored in a HierarchyMap in a hierarchical way
class DD4hepAlignmentDecorator : public IContextDecorator {
 public:
  /// @brief nested configuration struct
  struct Config {
    /// The alignment store map higher bound IOV (i.e. event number)
    std::vector<std::tuple<std::array<std::size_t, 2u>,
                           std::shared_ptr<Acts::IDD4hepAlignmentStore>>>
        alignmentStores;
    /// The nominal alignment store (before first bound, after last bound)
    std::shared_ptr<Acts::IDD4hepAlignmentStore> nominalStore = nullptr;
  };

  /// Constructor
  ///
  /// @param cfg Configuration struct
  /// @param logger The logging framework
  explicit DD4hepAlignmentDecorator(
      const Config& cfg,
      std::unique_ptr<const Acts::Logger> logger = Acts::getDefaultLogger(
          "DD4hepAlignmentDecorator", Acts::Logging::INFO));

  /// Virtual destructor
  ~DD4hepAlignmentDecorator() override = default;

  /// @brief decorates (adds, modifies) the AlgorithmContext
  /// with a geometric rotation per event
  ///
  /// @note If decorators depend on each other, they have to be
  /// added in order.
  ///
  /// @param context the bare (or at least non-const) Event context
  ProcessCode decorate(AlgorithmContext& context) override;

  /// @brief decorator name() for screen output
  const std::string& name() const override { return m_name; }

 private:
  Config m_cfg;                                  ///< the configuration class
  std::unique_ptr<const Acts::Logger> m_logger;  ///!< the logging instance
  const std::string m_name = "DD4hepAlignmentDecorator";

  /// Private access to the logging instance
  const Acts::Logger& logger() const { return *m_logger; }
};

}  // namespace ActsExamples
