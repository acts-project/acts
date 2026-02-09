// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <array>
#include <atomic>
#include <csignal>
#include <cstddef>
#include <limits>
#include <memory>
#include <mutex>
#include <stack>
#include <vector>

#include <boost/container/static_vector.hpp>
#include <boost/stacktrace/stacktrace_fwd.hpp>

namespace ActsPlugins {
/// @addtogroup fpemonitoring_plugin
/// @{

/// Floating-point exception types
enum class FpeType : std::uint32_t {
  INTDIV = FPE_INTDIV,
  INTOVF = FPE_INTOVF,
  FLTDIV = FPE_FLTDIV,
  FLTOVF = FPE_FLTOVF,
  FLTUND = FPE_FLTUND,
  FLTRES = FPE_FLTRES,
  FLTINV = FPE_FLTINV,
  FLTSUB = FPE_FLTSUB,
};

/// Output stream operator for FpeType
/// @param os The output stream
/// @param type The FPE type to output
/// @return The output stream
std::ostream &operator<<(std::ostream &os, FpeType type);

/// @brief Monitor for floating-point exceptions with stack trace capture
class FpeMonitor {
 public:
  /// @brief Buffer for storing stack traces
  struct Buffer {
    /// Constructor
    /// @param bufferSize Size of buffer to allocate
    explicit Buffer(std::size_t bufferSize)
        : m_data{std::make_unique<std::byte[]>(bufferSize)},
          m_size{bufferSize} {}

    Buffer(const Buffer &) = delete;
    /// Move constructor
    /// @param other Buffer to move from
    Buffer(Buffer &&other) noexcept
        : m_data(std::move(other.m_data)),
          m_size(other.m_size),
          m_offset(other.m_offset) {
      other.m_size = 0;
      other.m_offset = 0;
    }

    /// Get pointer and size of remaining buffer space
    /// @return Pair of pointer to next available byte and remaining size
    std::pair<void *, std::size_t> next() const {
      return {m_data.get() + m_offset, m_size - m_offset};
    }

    /// Advance buffer offset
    /// @param offset Number of bytes to advance
    void pushOffset(std::size_t offset) {
      assert(m_offset + offset < m_size);
      m_offset = offset;
    }

    /// Reset buffer offset to beginning
    void reset() { m_offset = 0; }

    /// Get total buffer size
    /// @return Total buffer size in bytes
    std::size_t size() const { return m_size; }
    /// Get current buffer offset
    /// @return Current offset in bytes
    std::size_t offset() const { return m_offset; }

    /// Get raw pointer to buffer data
    /// @return Pointer to buffer data
    std::byte *data() { return m_data.get(); }

   private:
    std::unique_ptr<std::byte[]> m_data;
    std::size_t m_size{};
    std::size_t m_offset{};
  };

  /// @brief Result of FPE monitoring containing counts and stack traces
  struct Result {
    /// @brief Information about a floating-point exception occurrence
    struct FpeInfo {
      /// Number of times this exception occurred
      std::size_t count;
      /// Type of floating-point exception
      FpeType type;
      /// Stack trace where the exception occurred
      std::shared_ptr<const boost::stacktrace::stacktrace> st;

      /// Constructor
      /// @param countIn Number of occurrences
      /// @param typeIn Exception type
      /// @param stIn Stack trace
      FpeInfo(std::size_t countIn, FpeType typeIn,
              std::shared_ptr<const boost::stacktrace::stacktrace> stIn);
      ~FpeInfo();
    };

    /// Merge with another result and return new result
    /// @param with Result to merge with
    /// @return Merged result
    Result merged(const Result &with) const;
    /// Merge another result into this one
    /// @param with Result to merge
    void merge(const Result &with);

    /// Check if an exception type was encountered
    /// @param type Exception type to check
    /// @return True if encountered
    bool encountered(FpeType type) const;
    /// Get count of exceptions of a specific type
    /// @param type Exception type
    /// @return Number of occurrences
    unsigned int count(FpeType type) const;

    /// Get all stack traces
    /// @return Vector of FPE information
    const std::vector<FpeInfo> &stackTraces() const;
    /// Get number of stack traces
    /// @return Number of stack traces
    unsigned int numStackTraces() const;

    /// Remove duplicate stack traces
    void deduplicate();

    /// Check if result contains a specific exception type and stack trace
    /// @param type Exception type
    /// @param st Stack trace to check
    /// @return True if contained
    bool contains(FpeType type, const boost::stacktrace::stacktrace &st) const;

    /// Print summary of exceptions
    /// @param os Output stream
    /// @param depth Maximum stack trace depth
    void summary(
        std::ostream &os,
        std::size_t depth = std::numeric_limits<std::size_t>::max()) const;

    Result() = default;

    /// Check if there are any stack traces
    /// @return True if stack traces are present
    bool hasStackTraces() const { return !m_stackTraces.empty(); }

    /// Add an exception occurrence from raw stack data
    /// @param type Exception type
    /// @param stackPtr Pointer to stack data
    /// @param bufferSize Size of stack buffer
    void add(FpeType type, void *stackPtr, std::size_t bufferSize);

   private:
    std::vector<FpeInfo> m_stackTraces;
    std::array<unsigned int, 32> m_counts{};

    friend FpeMonitor;
  };

  FpeMonitor();
  /// Constructor with exception mask
  /// @param excepts Mask of exceptions to monitor
  explicit FpeMonitor(int excepts);
  /// Move constructor
  /// @param other Monitor to move from
  FpeMonitor(FpeMonitor &&other) = default;
  ~FpeMonitor();

  /// Get monitoring result
  /// @return Reference to result object
  Result &result();

  /// Process recorded exceptions
  void consumeRecorded();

  /// Re-enable exception monitoring
  void rearm();

  /// Convert stack trace to string
  /// @param st Stack trace to convert
  /// @param depth Maximum depth of stack trace to include
  /// @return String representation of stack trace
  static std::string stackTraceToString(const boost::stacktrace::stacktrace &st,
                                        std::size_t depth);
  /// Get source location from stack frame
  /// @param frame Stack frame to extract location from
  /// @return String representation of source location
  static std::string getSourceLocation(const boost::stacktrace::frame &frame);

  /// Check if stack trace symbolization is available
  /// @return True if symbolization is available
  static bool canSymbolize();

 private:
  void enable();
  void disable();

  static void ensureSignalHandlerInstalled();
  static void signalHandler(int signal, siginfo_t *si, void *ctx);

  struct GlobalState {
    std::atomic_bool isSignalHandlerInstalled{false};
    std::mutex mutex{};
  };

  static std::stack<FpeMonitor *> &stack();
  static GlobalState &globalState();

  int m_excepts = 0;

  Result m_result;

  Buffer m_buffer{65536};

  boost::container::static_vector<std::tuple<FpeType, void *, std::size_t>, 128>
      m_recorded;
};

/// @}
}  // namespace ActsPlugins
