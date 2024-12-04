// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <Acts/Plugins/Json/StructuredLoggerJson.hpp>

#include <fstream>

namespace Acts {

StructuredLoggerJson::StructuredLoggerJson(std::filesystem::path outputFile)
    : m_path(outputFile) {
  m_json = nlohmann::json::array();
}

StructuredLoggerJson::~StructuredLoggerJson() {
  std::cout << "Dump json to " << m_path << std::endl;
  std::ofstream os(m_path);
  os << m_json.dump(2);
}

void StructuredLoggerJson::log(const Logger &logger, std::string_view msg) {
  nlohmann::json j;
  j["tag"] = logger.name();
  j["level"] = static_cast<int>(logger.level());
  j["type"] = "message";
  j["key"] = "";
  j["val"] = msg;

  std::lock_guard<std::mutex> guard(m_mutex);
  m_json.push_back(j);
  std::cout << "log msg to json" << std::endl;
}

void StructuredLoggerJson::log(std::string_view tag, std::string_view name,
                               const Eigen::MatrixXd &m) {
  nlohmann::json j;
  j["tag"] = tag;
  j["level"] = "STRUCTURED";
  j["key"] = name;

  if (m.rows() == 1) {
    j["type"] = "vector";
    j["val"] = m.row(0);
  } else {
    j["type"] = "matrix";
    j["val"] = nlohmann::json::array();
    for (const auto &row : m.rowwise()) {
      j["val"].push_back(row);
    }
  }

  std::lock_guard<std::mutex> guard(m_mutex);
  m_json.push_back(j);
}

void StructuredLoggerJson::log(std::string_view tag, std::string_view name,
                               const std::span<double> &s) {
  nlohmann::json j;
  j["tag"] = tag;
  j["level"] = "STRUCTURED";
  j["type"] = "vector";
  j["key"] = name;
  j["val"] = s;

  std::lock_guard<std::mutex> guard(m_mutex);
  m_json.push_back(j);
}

}  // namespace Acts
