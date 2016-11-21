// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_LOGFILE_ERASURE_FIXTURE_HPP
#define ACTS_LOGFILE_ERASURE_FIXTURE_HPP 1

#include <fstream>
#include <string>

namespace Acts {

namespace Test {

  struct logfile_erasure_fixture
  {
    logfile_erasure_fixture(const std::string& name)
    {
      std::ofstream f(name);
      f.close();
    }
  };

}  // namespace Test

}  // namespace Acts
#endif  // ACTS_LOGFILE_ERASURE_FIXTURE_HPP
