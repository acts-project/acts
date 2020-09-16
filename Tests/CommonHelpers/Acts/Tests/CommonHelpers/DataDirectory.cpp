// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Tests/CommonHelpers/DataDirectory.hpp"

#include <boost/filesystem.hpp>

std::string Acts::Test::getDataPath(const std::string& relativePath) {
  using boost::filesystem::path;

  path dataDir(ACTS_TEST_DATA_DIR);
  path absPath = dataDir / path(relativePath);
  // ensure absolute, consistent path
  // weakly_canonical allows non-existing trailing components
  return weakly_canonical(absPath).string();
}
