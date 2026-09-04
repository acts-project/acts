// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <boost/test/framework.hpp>
#include <boost/test/tree/test_unit.hpp>

#include <filesystem>

namespace ActsTests {

class TemporaryDirectory {
  static bool skipTemporary() {
    static const bool flag = [] {
      const char* envvar = std::getenv("ACTS_NO_TEST_TEMP");
      if (envvar != nullptr) {
        std::string val{envvar};
        if (!val.empty()) {
          return true;
        }
      }
      return false;
    }();

    return flag;
  }

 public:
  TemporaryDirectory() {
    if (skipTemporary()) {
      m_dir = std::filesystem::current_path();
      return;
    }

    auto base = std::filesystem::temp_directory_path();

    const auto& testName =
        boost::unit_test::framework::current_test_case().p_name.value;

    m_dir = base / ("acts_unit_test_" + testName);

    std::filesystem::create_directory(m_dir);
  }

  ~TemporaryDirectory() {
    if (!skipTemporary()) {
      std::filesystem::remove_all(m_dir);
    }
  }

  const std::filesystem::path& path() const { return m_dir; }

 private:
  std::filesystem::path m_dir;
};

}  // namespace ActsTests
