// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/IReader.hpp"

#include <filesystem>
#include <memory>
namespace ActsExamples {

namespace detail {
class EDM4hepReaderImpl;
}
class EDM4hepReader : public IReader {
 public:
  struct Config {
    std::filesystem::path inputPath;
    std::string outputFrame;
    /// The podio `category` name to read the frame from
    std::string category;
  };

  EDM4hepReader(const Config& config, Acts::Logging::Level level);
  ~EDM4hepReader() override;

  std::string name() const final;
  std::pair<std::size_t, std::size_t> availableEvents() const final;
  ProcessCode read(const ActsExamples::AlgorithmContext& context) final;

  const Config& config() const;

 private:
  std::unique_ptr<detail::EDM4hepReaderImpl> m_impl;
  std::unique_ptr<const Acts::Logger> m_logger;

  const Acts::Logger& logger() const { return *m_logger; }
};
}  // namespace ActsExamples
