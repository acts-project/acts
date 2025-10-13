// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/ActsVersion.hpp"

#include <ostream>
#include <string_view>

namespace Acts {

VersionInfo::VersionInfo(unsigned int majorIn, unsigned int minorIn,
                         unsigned int patchIn,
                         std::optional<std::string_view> commitHashIn)
    : versionMajor(majorIn),
      versionMinor(minorIn),
      versionPatch(patchIn),
      commitHash(commitHashIn) {}

VersionInfo VersionInfo::fromLibrary() {
  // this is filled by the Core shared library
  // while the constants below depend on the include
  return VersionInfo{VersionMajor, VersionMinor, VersionPatch, std::nullopt};
}

bool VersionInfo::operator==(const VersionInfo& other) const {
  return versionMajor == other.versionMajor &&
         versionMinor == other.versionMinor &&
         versionPatch == other.versionPatch && commitHash == other.commitHash;
}

VersionInfo VersionInfo::withoutCommit() const {
  auto copy = *this;
  copy.commitHash = std::nullopt;
  return copy;
}

std::ostream& operator<<(std::ostream& os, const VersionInfo& vi) {
  os << vi.versionMajor << "." << vi.versionMinor << "." << vi.versionPatch;
  constexpr static std::string_view s_unknown = "UNKNOWN";
  if (vi.commitHash) {
    os << " (commit " << vi.commitHash.value_or(s_unknown) << ")";
  }
  return os;
}
}  // namespace Acts
