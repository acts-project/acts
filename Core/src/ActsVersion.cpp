// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/ActsVersion.hpp"

#include <ostream>
#include <string_view>

namespace Acts {

VersionInfo::VersionInfo(unsigned int major, unsigned int minor,
                         unsigned int patch, const char* const commitHash)
    : major(major), minor(minor), patch(patch), commitHash(commitHash) {}

VersionInfo VersionInfo::fromLibrary() {
  // this is filled by the Core shared library
  // while the constants below depend on the include
  return VersionInfo{VersionMajor, VersionMinor, VersionPatch, CommitHash};
}

bool VersionInfo::operator==(const VersionInfo& other) const {
  return major == other.major && minor == other.minor && patch == other.patch &&
         std::string_view{commitHash} == std::string_view{other.commitHash};
}

std::ostream& operator<<(std::ostream& os, const VersionInfo& vi) {
  os << vi.major << "." << vi.minor << "." << vi.patch << " (commit "
     << vi.commitHash << ")";
  return os;
}
}  // namespace Acts
