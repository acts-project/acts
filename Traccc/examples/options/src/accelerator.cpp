// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Local include(s).
#include "traccc/options/accelerator.hpp"

#include "traccc/examples/utils/printable.hpp"

// System include(s).
#include <format>

namespace traccc::opts {

accelerator::accelerator() : interface("Accelerator Options") {
  m_desc.add_options()(
      "use-gpu-texture-memory",
      boost::program_options::bool_switch(&use_gpu_texture_memory),
      "Use GPU texture memory on the accelerator");
}

std::unique_ptr<configuration_printable> accelerator::as_printable() const {
  auto cat = std::make_unique<configuration_category>(m_description);

  cat->add_child(std::make_unique<configuration_kv_pair>(
      "Use GPU texture memory", std::format("{}", use_gpu_texture_memory)));

  return cat;
}

}  // namespace traccc::opts
