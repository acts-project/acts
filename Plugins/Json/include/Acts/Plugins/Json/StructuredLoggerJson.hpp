// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"

#include <filesystem>
#include <mutex>

#include <nlohmann/json.hpp>

namespace Acts {

class StructuredLoggerJson : public StructuredLoggerBase {
 public:
  StructuredLoggerJson(std::filesystem::path outputFile);
  ~StructuredLoggerJson();

  void log(const Logger &logger, std::string_view msg) override;
  void log(std::string_view tag, std::string_view name,
           const Eigen::MatrixXd &m) override;
  void log(std::string_view tag, std::string_view name,
           const std::span<double> &s) override;

 private:
  std::mutex m_mutex;
  nlohmann::json m_json;
  std::filesystem::path m_path;
};

}  // namespace Acts
