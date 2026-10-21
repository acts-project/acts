// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Local include(s).
#include "traccc/options/performance.hpp"

#include "traccc/examples/utils/printable.hpp"

// System include(s).
#include <format>

namespace traccc::opts {

performance::performance() : interface("Performance Measurement Options") {
  m_desc.add_options()("check-performance",
                       boost::program_options::bool_switch(&run),
                       "Run performance checks");
}

std::unique_ptr<configuration_printable> performance::as_printable() const {
  auto cat = std::make_unique<configuration_category>(m_description);

  cat->add_child(std::make_unique<configuration_kv_pair>(
      "Run performance checks", std::format("{}", run)));

  return cat;
}
}  // namespace traccc::opts
