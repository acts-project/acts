// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/TransformStore.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Framework/IContextDecorator.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <array>
#include <functional>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

namespace ActsExamples {
struct AlgorithmContext;

/// @brief A simple alignment decorator for the  geometry
/// showcasing an IOV based alignment possibility
///
class AlignmentDecorator : public IContextDecorator {
 public:
  /// @brief Nested configuration struct
  ///
  /// It can either be loaded with a transform store()s e.g. through I/O
  /// or with alignment generators that create the transform stores on the fly
  /// to demosntrate time dependent alignment,
  struct Config {
    /// The alignment store map higher bound IOV (i.e. event number)
    std::vector<std::tuple<std::array<std::size_t, 2u>,
                           std::shared_ptr<Acts::ITransformStore>>>
        alignmentStores;
    /// The nominal alignment store (before first bound, after last bound)
    std::shared_ptr<Acts::ITransformStore> nominalStore = nullptr;

    /// Create a nominal alignment on the fly, it receives a caller function
    /// that could provide random numbers
    std::function<std::shared_ptr<Acts::ITransformStore>(
        std::function<double()>&)>
        alignmentGemerator = nullptr;
  };

  /// Constructor
  ///
  /// @param cfg Configuration struct
  /// @param logger The logging framework
  explicit AlignmentDecorator(const Config& cfg,
                              std::unique_ptr<const Acts::Logger> logger =
                                  Acts::getDefaultLogger("AlignmentDecorator",
                                                         Acts::Logging::INFO));

  /// Virtual destructor
  ~AlignmentDecorator() override = default;

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
  const std::string m_name = "AlignmentDecorator";

  /// Private access to the logging instance
  const Acts::Logger& logger() const { return *m_logger; }
};

}  // namespace ActsExamples
