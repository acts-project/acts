// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Local include(s).
#include "traccc/options/details/interface.hpp"

namespace traccc::opts {

interface::interface(std::string_view group_desc)
    : m_desc(std::string{group_desc}), m_description(group_desc) {}

void interface::read(const boost::program_options::variables_map&) {}

const boost::program_options::options_description& interface::options() const {
  return m_desc;
}
}  // namespace traccc::opts
