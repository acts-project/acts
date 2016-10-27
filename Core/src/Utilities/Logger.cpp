// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTS/Utilities/Logger.hpp"

namespace Acts {

std::unique_ptr<Logger>
getDefaultLogger(const std::string&    name,
                 const Logging::Level& lvl,
                 std::ostream*         log_stream)
{
  using namespace Logging;
  auto output = std::make_unique<LevelOutputDecorator>(
      std::make_unique<NamedOutputDecorator>(
          std::make_unique<TimedOutputDecorator>(
              std::make_unique<DefaultOutputPolicy>(log_stream)),
          name));
  auto print = std::make_unique<DefaultPrintPolicy>(lvl);
  return std::make_unique<Logger>(std::move(output), std::move(print));
}
}  // end of namespace Acts
