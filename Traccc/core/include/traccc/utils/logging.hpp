/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Acts include(s).
#include <Acts/Utilities/Logger.hpp>

// detray include(s)
#include <detray/utils/logging_streams.hpp>

namespace traccc {

/// Use the @c Acts::Logging namespace.
namespace Logging = ::Acts::Logging;

/// Use the @c Acts::Logger type.
using Logger = ::Acts::Logger;

/// Construct a logger with default settings for this project
///
/// @param name the name of the log writer
/// @param lvl the log level
/// @param log_stream the stream to write the log to
///
/// @return a unique pointer to the logger
///
std::unique_ptr<const Logger> getDefaultLogger(
    const std::string& name, const Logging::Level& lvl = Logging::INFO,
    std::ostream* log_stream = &std::cout);

/// Construct a dummy logger that does nothing
///
/// @return a reference to the dummy logger
///
const Logger& getDummyLogger();

}  // namespace traccc

#define TRACCC_LOCAL_LOGGER(x) ACTS_LOCAL_LOGGER(x)
#define TRACCC_LOG(level, x) ACTS_LOG(level, log)
#define TRACCC_VERBOSE(x) ACTS_VERBOSE(x)
#define TRACCC_DEBUG(x) ACTS_DEBUG(x)
#define TRACCC_INFO(x) ACTS_INFO(x)
#define TRACCC_WARNING(x) ACTS_WARNING(x)
#define TRACCC_ERROR(x) ACTS_ERROR(x)
#define TRACCC_FATAL(x) ACTS_FATAL(x)

// Define traccc logging macros that are safe to use in device compiled code:

// Printed only when compiled for host
// TODO: Harmonize with ACTS logger
#define TRACCC_FATAL_HOST(x) DETRAY_FATAL_STREAM("TRACCC", x)
#define TRACCC_ERROR_HOST(x) DETRAY_ERROR_STREAM("TRACCC", x)
#define TRACCC_WARNING_HOST(x) DETRAY_WARN_STREAM("TRACCC", x)
#define TRACCC_INFO_HOST(x) DETRAY_INFO_STREAM("TRACCC", x)
#define TRACCC_VERBOSE_HOST(x) DETRAY_VERBOSE_STREAM("TRACCC", x)
#define TRACCC_DEBUG_HOST(x) DETRAY_DEBUG_STREAM("TRACCC", x)

// Printed in both host and device execution
// TODO: Implement rate limiting
#define TRACCC_FATAL_HOST_DEVICE(x, ...) \
    DETRAY_FATAL_PRINTF("TRACCC", x, __VA_ARGS__)
#define TRACCC_ERROR_HOST_DEVICE(x, ...) \
    DETRAY_ERROR_PRINTF("TRACCC", x, __VA_ARGS__)
#define TRACCC_WARNING_HOST_DEVICE(x, ...) \
    DETRAY_WARN_PRINTF("TRACCC", x, __VA_ARGS__)
#define TRACCC_INFO_HOST_DEVICE(x, ...) \
    DETRAY_INFO_PRINTF("TRACCC", x, __VA_ARGS__)
#define TRACCC_VERBOSE_HOST_DEVICE(x, ...) \
    DETRAY_VERBOSE_PRINTF("TRACCC", x, __VA_ARGS__)
#define TRACCC_DEBUG_HOST_DEVICE(x, ...) \
    DETRAY_DEBUG_PRINTF("TRACCC", x, __VA_ARGS__)

// Printed only when compiled for device
#ifdef __DEVICE_LOGGING__

#define TRACCC_FATAL_DEVICE(x, ...) \
    DETRAY_FATAL_PRINTF("TRACCC", x, __VA_ARGS__)
#define TRACCC_ERROR_DEVICE(x, ...) \
    DETRAY_ERROR_PRINTF("TRACCC", x, __VA_ARGS__)
#define TRACCC_WARNING_DEVICE(x, ...) \
    DETRAY_WARN_PRINTF("TRACCC", x, __VA_ARGS__)
#define TRACCC_INFO_DEVICE(x, ...) DETRAY_INFO_PRINTF("TRACCC", x, __VA_ARGS__)
#define TRACCC_VERBOSE_DEVICE(x, ...) \
    DETRAY_VERBOSE_PRINTF("TRACCC", x, __VA_ARGS__)
#define TRACCC_DEBUG_DEVICE(x, ...) \
    DETRAY_DEBUG_PRINTF("TRACCC", x, __VA_ARGS__)

#else

#define TRACCC_FATAL_DEVICE(x, ...)
#define TRACCC_ERROR_DEVICE(x, ...)
#define TRACCC_WARNING_DEVICE(x, ...)
#define TRACCC_INFO_DEVICE(x, ...)
#define TRACCC_VERBOSE_DEVICE(x, ...)
#define TRACCC_DEBUG_DEVICE(x, ...)

#endif
