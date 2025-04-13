// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IWriter.hpp"

#include <memory>
#include <string>

namespace ActsExamples {

namespace detail {
class PodioWriterImpl;
}  // namespace detail

class PodioWriter final : public IWriter {
 public:
  struct Config {
    /// The path to the output file.
    std::string outputPath;

    /// Retrieve a @c podio::Frame from the event store using this name.
    /// @note If not set, a new frame will be created.
    std::optional<std::string> inputFrame = std::nullopt;

    /// The podio `category` name to write the frame to
    std::string category;

    /// The collection names to write to the output file.
    std::vector<std::string> collections;
  };

  PodioWriter(const Config& config, Acts::Logging::Level level);
  ~PodioWriter() override;

  std::string name() const override;

  /// Readonly access to the config
  const Config& config() const;

  ProcessCode finalize() override;

  ProcessCode write(const AlgorithmContext& ctx) override;

 private:
  const Acts::Logger& logger() const { return *m_logger; }

  std::unique_ptr<const Acts::Logger> m_logger;

  std::unique_ptr<detail::PodioWriterImpl> m_impl;
};

}  // namespace ActsExamples
