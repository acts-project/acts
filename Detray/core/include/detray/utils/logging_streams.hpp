// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2025-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#if DETRAY_LOG_LVL >= 0

// Project include(s)
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/utils/type_list.hpp"

// System include(s)
#include <cstring>
#include <iomanip>
#include <iostream>
#include <regex>
#include <source_location>
#include <string>

namespace detray::log::detail {

/// @returns the name of the current source file without the full path
/// @see
/// https://stackoverflow.com/questions/31050113/how-to-extract-the-source-filename-without-path-and-suffix-at-compile-time
constexpr const char *source_file_name(const char *path) {
  const char *file = path;
  while (*path) {
    if (*path++ == '/') {
      file = path;
    }
  }
  return file;
}

/// Print a type name for type @tparam T in the logs
template <typename T>
inline std::string_view process_typename() {
  static const std::string type_name = [] {
    std::string s{""};
    try {
      s = detray::types::demangle_type_name<T>();
    } catch (...) {
      return std::string{"unknown"};
    }

    if (s.empty()) {
      return std::string{"unknown"};
    }

    std::regex re{"detray::"};
    s = std::regex_replace(s, re, "");

    // Special case type-list for readability
    // Only attempt if the full type is a list
    std::regex re_list{R"(^types::list<(.*)>$)"};

    if (std::smatch match; std::regex_match(s, match, re_list)) {
      s = "[" + std::string{match[1]} + "]";
    }

    return s;
  }();

  return type_name;
}

}  // namespace detray::log::detail

// Enable device-side logging
#if defined(__CUDACC__) || defined(__HIP__) || \
    defined(CL_SYCL_LANGUAGE_VERSION) || defined(SYCL_LANGUAGE_VERSION)
#define __DEVICE_LOGGING__
#endif

#ifdef __DEVICE_LOGGING__
#define DETRAY_TYPENAME(type) "unknown type"
#else
#define DETRAY_TYPENAME(type) detray::log::detail::process_typename<type>()
#endif

#ifdef __DEVICE_LOGGING__
#define DETRAY_LOG_VECTOR(x) ""
#else
#define DETRAY_LOG_VECTOR(x)          \
  [&]() {                             \
    std::stringstream _vec_os;        \
    std::size_t _vec_i = 0;           \
    for (const auto &_vec_elem : x) { \
      if (_vec_i > 0) {               \
        _vec_os << ", ";              \
      }                               \
      _vec_os << _vec_elem;           \
      _vec_i++;                       \
    }                                 \
    return _vec_os.str();             \
  }()
#endif

// String that represents the backend that emitted the log message
#if defined(__CUDACC__)
#define __BACKEND__ "CUDA"
#elif defined(__HIP__)
#define __BACKEND__ "HIP"
#elif defined(CL_SYCL_LANGUAGE_VERSION) || defined(SYCL_LANGUAGE_VERSION)
#define __BACKEND__ "SYCL"
#else
#define __BACKEND__ "HOST"
#endif

#ifdef __DEVICE_LOGGING__
#define __FILENAME__ detray::log::detail::source_file_name(__FILE__)
#else
#define __FILENAME__ \
  (std::strrchr(__FILE__, '/') ? std::strrchr(__FILE__, '/') + 1 : __FILE__)
#endif

// Print 'x' in the host logs only
#ifndef __DEVICE_LOGGING__

#define DETRAY_LOG_STREAM(stream, lib, lvl, x)                               \
  std::stream << lib << " " << std::left << std::setw(7) << lvl << std::left \
              << std::setw(9) << " (HOST):" << std::left << std::setw(29)    \
              << __FILENAME__ << "l." << std::left << std::setw(5)           \
              << std::source_location::current().line() << x << std::endl

// Log output to stdout
#define DETRAY_LOG_STREAM_STDOUT(lib, lvl, x) \
  DETRAY_LOG_STREAM(cout, lib, lvl, x)

// Log output to stderr
#define DETRAY_LOG_STREAM_STDERR(lib, lvl, x) \
  DETRAY_LOG_STREAM(clog, lib, lvl, x)

#else  // ifndef __DEVICE_LOGGING__
#define DETRAY_LOG_STREAM_STDOUT(lib, lvl, x)
#define DETRAY_LOG_STREAM_STDERR(lib, lvl, x)
#endif

// Print 'x' in printf
// @note 'x' is the format string and the variadic arguments the corresponding
// values
// @note logging is currently disabled for SYCL builds
#if !defined(CL_SYCL_LANGUAGE_VERSION) && !defined(SYCL_LANGUAGE_VERSION)
// TODO: Move to device specific header for e.g. rate limiting by thread idx
#ifdef __DEVICE_LOGGING__

#define DETRAY_LOG_PRINTF(lib, lvl, x, ...)                         \
  printf("%s %-7s (%s): %-29sl.%-5d" x "\n", lib, lvl, __BACKEND__, \
         __FILENAME__, __LINE__ __VA_OPT__(, ) __VA_ARGS__)

#define DETRAY_LOG_PRINTF_STDOUT(lib, lvl, x, ...) \
  DETRAY_LOG_PRINTF(lib, lvl, x, __VA_ARGS__)

#define DETRAY_LOG_PRINTF_STDERR(lib, lvl, x, ...) \
  DETRAY_LOG_PRINTF(lib, lvl, x, __VA_ARGS__)

#else  // host-side printf logging

#define DETRAY_LOG_PRINTF(stream, lib, lvl, x, ...)                          \
  fprintf(stream, "%s %-7s (%s): %-29sl.%-5d" x "\n", lib, lvl, __BACKEND__, \
          __FILENAME__, __LINE__ __VA_OPT__(, ) __VA_ARGS__)

#define DETRAY_LOG_PRINTF_STDOUT(lib, lvl, x, ...) \
  DETRAY_LOG_PRINTF(stdout, lib, lvl, x, __VA_ARGS__)

#define DETRAY_LOG_PRINTF_STDERR(lib, lvl, x, ...) \
  DETRAY_LOG_PRINTF(stderr, lib, lvl, x, __VA_ARGS__)

#endif  // ifdef __DEVICE_LOGGING__

#else  // ifndef SYCL
#define DETRAY_LOG_PRINTF_STDOUT(lib, lvl, x, ...)
#define DETRAY_LOG_PRINTF_STDERR(lib, lvl, x, ...)
#endif

#else  // DETRAY_LOG_LVL < 0
#define DETRAY_LOG_STREAM_STDOUT(lib, lvl, x)
#define DETRAY_LOG_STREAM_STDERR(lib, lvl, x)
#define DETRAY_LOG_PRINTF_STDOUT(lib, lvl, x, ...)
#define DETRAY_LOG_PRINTF_STDERR(lib, lvl, x, ...)
#define DETRAY_TYPENAME(type)
#endif

//
// Define log levels
//

// Errors and warnings
#define DETRAY_FATAL_STREAM(lib, x) DETRAY_LOG_STREAM_STDERR(lib, "FATAL", x)
#define DETRAY_ERROR_STREAM(lib, x) DETRAY_LOG_STREAM_STDERR(lib, "ERROR", x)
#define DETRAY_WARN_STREAM(lib, x) DETRAY_LOG_STREAM_STDERR(lib, "WARNING", x)

#define DETRAY_FATAL_PRINTF(lib, x, ...) \
  DETRAY_LOG_PRINTF_STDERR(lib, "FATAL", x, __VA_ARGS__)
#define DETRAY_ERROR_PRINTF(lib, x, ...) \
  DETRAY_LOG_PRINTF_STDERR(lib, "ERROR", x, __VA_ARGS__)
#define DETRAY_WARN_PRINTF(lib, x, ...) \
  DETRAY_LOG_PRINTF_STDERR(lib, "WARNING", x, __VA_ARGS__)

// Info
#if DETRAY_LOG_LVL > 0
#define DETRAY_INFO_STREAM(lib, x) DETRAY_LOG_STREAM_STDOUT(lib, "INFO", x)
#define DETRAY_INFO_PRINTF(lib, x, ...) \
  DETRAY_LOG_PRINTF_STDOUT(lib, "INFO", x, __VA_ARGS__)
#else
#define DETRAY_INFO_STREAM(lib, x)
#define DETRAY_INFO_PRINTF(lib, x, ...)
#endif

// Verbose
#if DETRAY_LOG_LVL > 1
#define DETRAY_VERBOSE_STREAM(lib, x) \
  DETRAY_LOG_STREAM_STDERR(lib, "VERBOSE", x)
#define DETRAY_VERBOSE_PRINTF(lib, x, ...) \
  DETRAY_LOG_PRINTF_STDERR(lib, "VERBOSE", x, __VA_ARGS__)
#else
#define DETRAY_VERBOSE_STREAM(lib, x)
#define DETRAY_VERBOSE_PRINTF(lib, x, ...)
#endif

// Debug
#if DETRAY_LOG_LVL > 2
#define DETRAY_DEBUG_STREAM(lib, x) DETRAY_LOG_STREAM_STDERR(lib, "DEBUG", x)
#define DETRAY_DEBUG_PRINTF(lib, x, ...) \
  DETRAY_LOG_PRINTF_STDERR(lib, "DEBUG", x, __VA_ARGS__)
#else
#define DETRAY_DEBUG_STREAM(lib, x)
#define DETRAY_DEBUG_PRINTF(lib, x, ...)
#endif
