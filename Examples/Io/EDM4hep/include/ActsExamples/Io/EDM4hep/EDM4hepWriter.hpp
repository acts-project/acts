// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Plugins/Podio/PodioUtil.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IWriter.hpp"

#include <memory>
#include <string>

#include <DD4hep/DetElement.h>

namespace ActsExamples {

namespace detail {
class EDM4hepWriterImpl;
}  // namespace detail

class EDM4hepWriter final : public IWriter {
 public:
  struct Config {
    std::string outputPath;

    /// Retrieve a @c podio::Frame from the event store using this name.
    std::string inputFrame;

    /// The podio `category` name to write the frame to
    std::string category;
  };

  EDM4hepWriter(const Config& config, Acts::Logging::Level level);
  ~EDM4hepWriter() override;

  std::string name() const final;

  /// Readonly access to the config
  const Config& config() const;

  ProcessCode finalize() final;

  ProcessCode write(const AlgorithmContext& ctx) final;

 private:
  const Acts::Logger& logger() const { return *m_logger; }

  std::unique_ptr<const Acts::Logger> m_logger;

  std::unique_ptr<detail::EDM4hepWriterImpl> m_impl;
};

}  // namespace ActsExamples
