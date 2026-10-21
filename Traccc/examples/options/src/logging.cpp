// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Local include(s).
#include "traccc/options/logging.hpp"

#include <Acts/Utilities/Logger.hpp>

#include "traccc/examples/utils/printable.hpp"

// System include(s).
#include <format>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace {
class verbosity_counter : public boost::program_options::typed_value<int> {
 public:
  verbosity_counter(int* store)
      : boost::program_options::typed_value<int>(store) {
    default_value(0);
    zero_tokens();
  }

  virtual ~verbosity_counter() = default;

  virtual void xparse(boost::any& store,
                      const std::vector<std::string>&) const {
    store = boost::any(++m_value);
  }

 private:
  mutable int m_value{0};
};
}  // namespace

namespace traccc::opts {

logging::logging() : interface("Logging Options") {
  m_desc.add_options()("verbose,v", new verbosity_counter(&m_verbosity_decr),
                       "Increase verbosity (can be specified multiple times)")(
      "quiet,q", new verbosity_counter(&m_verbosity_incr),
      "Decrease verbosity (can be specified multiple times)");
}

logging::operator traccc::Logging::Level() const {
  int verb = m_verbosity_incr - m_verbosity_decr;

  if (verb <= -2) {
    return traccc::Logging::Level::VERBOSE;
  } else if (verb == -1) {
    return traccc::Logging::Level::DEBUG;
  } else if (verb == 0) {
    return traccc::Logging::Level::INFO;
  } else if (verb == 1) {
    return traccc::Logging::Level::WARNING;
  } else if (verb == 2) {
    return traccc::Logging::Level::ERROR;
  } else {
    return traccc::Logging::Level::FATAL;
  }
}

std::unique_ptr<configuration_printable> logging::as_printable() const {
  auto cat = std::make_unique<configuration_category>(m_description);

  traccc::Logging::Level lvl(*this);

  std::string verbosity_string{Acts::Logging::levelName(lvl)};

  cat->add_child(std::make_unique<configuration_kv_pair>("Verbosity level",
                                                         verbosity_string));

  return cat;
}

}  // namespace traccc::opts
