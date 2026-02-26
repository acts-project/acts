// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"

namespace podio {
class Frame;
}

namespace ActsExamples {

/// This is an abstract base class for PODIO input converters, that
/// centralizes the retrieval of the input @c podio::Frame object to convert
/// from.
class PodioInputConverter : public IAlgorithm {
 public:
  /// Constructor for the PODIO input converter.
  ///
  /// @param name The name of the algorithm
  /// @param inputFrame The input frame to convert
  /// @param logger The logger instance
  PodioInputConverter(const std::string& name, const std::string& inputFrame,
                      std::unique_ptr<const Acts::Logger> logger = nullptr);

  /// Destructor for the PODIO input converter.
  ~PodioInputConverter() override;

  /// Execute the algorithm. Subclasses do not implement this method.
  ///
  /// @param ctx The algorithm context
  /// @return The process code
  ProcessCode execute(const AlgorithmContext& ctx) const final;

  /// Convert the input @c podio::Frame object to the internal EDM format.
  ///
  /// @param ctx The algorithm context
  /// @param frame The input @c podio::Frame object
  /// @return The process code
  virtual ProcessCode convert(const AlgorithmContext& ctx,
                              const podio::Frame& frame) const = 0;

 private:
  class Impl;

  std::unique_ptr<Impl> m_impl;
};

}  // namespace ActsExamples
