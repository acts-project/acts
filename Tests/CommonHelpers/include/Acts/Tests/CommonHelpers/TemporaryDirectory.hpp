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

namespace Acts::Test {
class TemporaryDirectory {
 public:
  TemporaryDirectory() {
    auto base = std::filesystem::temp_directory_path();

    const auto& testName =
        boost::unit_test::framework::current_test_case().p_name.value;

    m_dir = base / ("acts_unit_test_" + testName);

    std::filesystem::create_directory(m_dir);
  }

  ~TemporaryDirectory() { std::filesystem::remove_all(m_dir); }

  const std::filesystem::path& path() const { return m_dir; }

 private:
  std::filesystem::path m_dir;
};
}  // namespace Acts::Test
