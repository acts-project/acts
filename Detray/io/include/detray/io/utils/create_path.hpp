// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// System include(s).
#include <filesystem>
#include <string>
#include <system_error>

namespace detray::io {

/// Check if a file exists
inline bool file_exists(const std::string& outdir) {
  auto path = std::filesystem::path(outdir);

  return std::filesystem::exists(path);
}

/// @returns an alternative file name by counting up, if file exists
inline std::string alt_file_name(const std::string& name) {
  auto path = std::filesystem::path(name);
  auto parent_path = path.parent_path();
  std::string stem = path.stem();
  std::string extension = path.extension();

  /// @returns alternate file stem upon collision
  auto get_alternate_file_stem = [](std::string& file_stem,
                                    const std::size_t n) {
    const std::string delim{"_"};

    // File stem already comes with a number, simply update it
    if (n > 2u) {
      std::size_t pos{file_stem.rfind(delim)};
      if (pos == std::string::npos || (pos + 1 == file_stem.size())) {
        throw std::runtime_error("Malformed file name");
      }

      // Leave the delimiter where it is and replace the file number
      ++pos;
      file_stem.replace(pos, file_stem.size() - pos, std::to_string(n));

      return file_stem;
    }

    return file_stem + delim + std::to_string(n);
  };

  std::size_t n_trials{2u};
  while (std::filesystem::exists(path)) {
    stem = get_alternate_file_stem(stem, n_trials);
    path = parent_path / std::filesystem::path{stem + extension};
    ++n_trials;

    // The maximum here is arbitrary
    if (n_trials >= 10000u) {
      throw std::runtime_error("Too many versions of file exist: " + stem +
                               extension);
    }
  }

  path = parent_path / std::filesystem::path{stem + extension};
  return path.string();
}

/// Check if a given file path exists and generate it if not
inline auto create_path(const std::string& outdir) {
  auto path = std::filesystem::path(outdir);

  if (!path.empty() && !std::filesystem::exists(path)) {
    if (std::error_code err; !std::filesystem::create_directories(path, err)) {
      throw std::runtime_error(err.message());
    }
  }

  return path;
}

}  // namespace detray::io
