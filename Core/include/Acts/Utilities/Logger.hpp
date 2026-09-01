// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// STL include(s)
#include <ctime>
#include <iomanip>
#include <iostream>
#include <memory>
#include <mutex>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <utility>

/// @addtogroup logging
/// @{

/// @defgroup logging_macros Logging Macros
/// @ingroup logging
/// @brief Helper macros for logging with @ref Acts::Logger
///
/// When a logger accessible via the `logger()` method, see @ref logging_patterns,
/// use these macros to perform the actual logging:
///
/// @snippet{trimleft} examples/logging.cpp Logging Macros
///
/// The macros support stream-style formatting with `<<` operators.
/// @{

/// @brief Macro to use a local Acts::Logger object
///
/// @param log_object logger instance of type
//         `std::unique_ptr<const Acts::Logger>`
///
/// @pre In the current scope, the symbol @c logger is not yet defined.
/// @post The ownership of the given @c log_object is transferred and
///       @c log_object should not be used directly any more.
///
/// This macro allows to use a locally defined logging object with the ACTS_*
/// logging macros. The envisaged usage is the following:
///
/// @snippet{trimleft} examples/logging.cpp Local logger macro
#define ACTS_LOCAL_LOGGER(log_object)                                          \
  struct __local_acts_logger {                                                 \
    explicit __local_acts_logger(std::unique_ptr<const ::Acts::Logger> logger) \
        : m_logger(std::move(logger)) {}                                       \
                                                                               \
    const ::Acts::Logger& operator()() const { return *m_logger; }             \
                                                                               \
    std::unique_ptr<const ::Acts::Logger> m_logger;                            \
  };                                                                           \
  __local_acts_logger logger(log_object);

/// Log a message at the specified level with an explicit logger instance
/// @param lgr The logger instance (must be a Acts::Logger reference)
/// @param level The logging level
/// @param x The message to log
#define ACTS_LOG_WITH_LOGGER(lgr, level, x) \
  do {                                      \
    if ((lgr).doPrint(level)) {             \
      std::ostringstream os;                \
      os << x;                              \
      (lgr).log(level, os.str());           \
    }                                       \
  } while (0)

/// Log a message at the specified level
/// @param level The logging level
/// @param x The message to log
#define ACTS_LOG(level, x) ACTS_LOG_WITH_LOGGER(logger(), level, x)

/// @brief macro for verbose debug output
///
/// @param x debug message
///
/// @pre @c logger() must be a valid expression in the scope where this
///      macro is used and it must return a Acts::Logger object.
///
/// The debug message is printed if the current Acts::Logging::Level <=
/// Acts::Logging::VERBOSE.
#define ACTS_VERBOSE(x) ACTS_LOG(Acts::Logging::VERBOSE, x)

/// @brief macro for debug debug output
///
/// @param x debug message
///
/// @pre @c logger() must be a valid expression in the scope where this
///      macro is used and it must return a Acts::Logger object.
///
/// The debug message is printed if the current Acts::Logging::Level <=
/// Acts::Logging::DEBUG.
#define ACTS_DEBUG(x) ACTS_LOG(Acts::Logging::DEBUG, x)

/// @brief macro for info debug output
///
/// @param x debug message
///
/// @pre @c logger() must be a valid expression in the scope where this
///      macro is used and it must return a Acts::Logger object.
///
/// The debug message is printed if the current Acts::Logging::Level <=
/// Acts::Logging::INFO.
#define ACTS_INFO(x) ACTS_LOG(Acts::Logging::INFO, x)

/// @brief macro for warning debug output
///
/// @param x debug message
///
/// @pre @c logger() must be a valid expression in the scope where this
///      macro is used and it must return a Acts::Logger object.
///
/// The debug message is printed if the current Acts::Logging::Level <=
/// Acts::Logging::WARNING.
#define ACTS_WARNING(x) ACTS_LOG(Acts::Logging::WARNING, x)

/// @brief macro for error debug output
///
/// @param x debug message
///
/// @pre @c logger() must be a valid expression in the scope where this
///      macro is used and it must return a Acts::Logger object.
///
/// The debug message is printed if the current Acts::Logging::Level <=
/// Acts::Logging::ERROR.
#define ACTS_ERROR(x) ACTS_LOG(Acts::Logging::ERROR, x)

/// @brief macro for fatal debug output
///
/// @param x debug message
///
/// @pre @c logger() must be a valid expression in the scope where this
///      macro is used and it must return a Acts::Logger object.
///
/// The debug message is printed if the current Acts::Logging::Level <=
/// Acts::Logging::FATAL.
#define ACTS_FATAL(x) ACTS_LOG(Acts::Logging::FATAL, x)

/// @}
/// @}

namespace Acts {

namespace Logging {

/// @addtogroup logging
/// @{

/// @brief constants steering the debug output
///
/// All messages with a debug level equal or higher than the currently set
/// debug output level will be printed.
enum Level {
  VERBOSE = 0,  ///< Detailed diagnostic trace information
  DEBUG,        ///< Debug information during development
  INFO,         ///< General information messages
  WARNING,      ///< Non-critical error conditions
  ERROR,        ///< Error conditions which require follow-up
  FATAL,        ///< Unrecoverable error conditions
  MAX           ///< Filler level
};

/// @brief Get the string name for a logging level
/// @param level The logging level
/// @return String representation of the logging level
inline std::string_view levelName(Level level) {
  switch (level) {
    case Level::VERBOSE:
      return "VERBOSE";
    case Level::DEBUG:
      return "DEBUG";
    case Level::INFO:
      return "INFO";
    case Level::WARNING:
      return "WARNING";
    case Level::ERROR:
      return "ERROR";
    case Level::FATAL:
      return "FATAL";
    case Level::MAX:
      return "MAX";
    default:
      throw std::invalid_argument{"Unknown level"};
  }
}

/// @defgroup logging_thresholds Logging Thresholds
/// @ingroup logging
/// @brief Functions and classes to manage logging failure thresholds
///
/// Generally, log levels in ACTS are only of informative value: even
/// @ref Acts::Logging::Level::ERROR and @ref Acts::Logging::Level::FATAL will only print
/// messages, **and not terminate execution**.
///
/// This is desirable in an experiment context, where jobs should not
/// immediately terminate when ACTS encounters something that is logged as an
/// error. In a test context, however, this behavior is not optimal: the tests
/// should ensure in known configurations errors do not occur, or only in
/// specific circumstances. To solve this, ACTS implements an optional log
/// *threshold* mechanism.
///
/// The failure threshold is a property of the @ref Acts::Logger: a message at
/// or above it is printed and then raises @ref Acts::Logging::ThresholdFailure.
/// @ref Acts::Logging::Level::MAX, the default, never fails.
///
/// A logger carries the threshold it was built with; nothing is consulted when
/// a message is logged.
///
/// @ref Acts::getDefaultLogger arms the loggers it builds at
/// @ref Acts::Logging::detail::getDefaultFailureThreshold. A @ref Acts::Logger
/// built directly from a print and a filter policy is not armed.
///
/// Core does not read the environment; the `ACTS_LOG_FAILURE_THRESHOLD`
/// environment variable is applied by the Python bindings on import.
///
/// @note The check happens before the level filter, so a coarse filter level
///       cannot hide a message that should fail the job, and it applies to
///       every @ref Acts::OutputPrintPolicy.
///
/// @{

/// Custom exception class so threshold failures can be caught
class ThresholdFailure : public std::runtime_error {
  using std::runtime_error::runtime_error;
};

namespace detail {

/// @brief Get the threshold @ref Acts::getDefaultLogger arms new loggers at
/// @return the default, or @ref Level::MAX for "do not arm"
Level getDefaultFailureThreshold();

/// @brief Set the threshold @ref Acts::getDefaultLogger arms new loggers at
///
/// @warning Global state, and only reaches loggers built after the call.
/// @param level The new default, or @ref Level::MAX for "do not arm"
void setDefaultFailureThreshold(Level level);

}  // namespace detail

/// @}

/// @brief abstract base class for printing debug output
///
/// Implementations of this interface need to define how and where to @a print
/// debug messages (e.g. to a file, to a stream into a database etc).
class OutputPrintPolicy {
 public:
  /// virtual default destructor
  virtual ~OutputPrintPolicy() = default;

  /// @brief handle output of debug message
  ///
  /// @param [in] lvl   debug output level of message
  /// @param [in] input text of debug message
  virtual void flush(const Level& lvl, const std::string& input) = 0;

  /// Return the name of the print policy
  /// @return the name
  virtual const std::string& name() const = 0;

  /// Make a copy of this print policy with a new name
  /// @param name the new name
  /// @return the copy
  virtual std::unique_ptr<OutputPrintPolicy> clone(
      const std::string& name) const = 0;
};

/// @brief abstract base class for filtering debug output
///
/// Implementations of this interface need to define whether a debug message
/// with a certain debug level is processed or filtered out.
class OutputFilterPolicy {
 public:
  /// virtual default destructor
  virtual ~OutputFilterPolicy() = default;

  /// @brief decide whether a debug message should be processed
  ///
  /// @param [in] lvl debug level of debug message
  ///
  /// @return @c true of debug message should be processed, @c false if debug
  ///         message should be skipped
  virtual bool doPrint(const Level& lvl) const = 0;

  /// Get the level of this filter policy
  /// @return the levele
  virtual Level level() const = 0;

  /// Make a copy of this filter policy with a new level
  /// @param level the new level
  /// @return the new copy
  virtual std::unique_ptr<OutputFilterPolicy> clone(Level level) const = 0;
};

/// @brief default filter policy for debug messages
///
/// All debug messages with a debug level equal or larger to the specified
/// threshold level are processed.
class DefaultFilterPolicy final : public OutputFilterPolicy {
 public:
  /// @brief constructor
  ///
  /// @param [in] lvl threshold debug level
  explicit DefaultFilterPolicy(Level lvl) : m_level(lvl) {}

  /// virtual default destructor
  ~DefaultFilterPolicy() override = default;

  /// @brief decide whether a debug message should be processed
  ///
  /// @param [in] lvl debug level of debug message
  ///
  /// @return @c true if @p lvl >= #m_level, otherwise @c false
  bool doPrint(const Level& lvl) const override { return m_level <= lvl; }

  /// Get the level of this filter policy
  /// @return the levele
  Level level() const override { return m_level; }

  /// Make a copy of this filter policy with a new level
  /// @param level the new level
  /// @return the new copy
  std::unique_ptr<OutputFilterPolicy> clone(Level level) const override {
    return std::make_unique<DefaultFilterPolicy>(level);
  }

 private:
  /// threshold debug level for messages to be processed
  Level m_level;
};

/// @brief base class for decorating the debug output
///
/// Derived classes may augment the debug message with additional information.
/// Chaining different decorators is possible to customize the output to your
/// needs.
class OutputDecorator : public OutputPrintPolicy {
 public:
  /// @brief constructor wrapping actual output print policy
  ///
  /// @param [in] wrappee output print policy object which is wrapped by this
  ///        decorator object
  explicit OutputDecorator(std::unique_ptr<OutputPrintPolicy> wrappee)
      : m_wrappee(std::move(wrappee)) {}

  /// @brief flush the debug message to the destination stream
  ///
  /// @param [in] lvl   debug level of debug message
  /// @param [in] input text of debug message
  ///
  /// This function delegates the flushing of the debug message to its wrapped
  /// object.
  void flush(const Level& lvl, const std::string& input) override {
    m_wrappee->flush(lvl, input);
  }

  /// Return the name of the output decorator (forwards to wrappee)
  /// @return the name
  const std::string& name() const override { return m_wrappee->name(); }

 protected:
  /// wrapped object for printing the debug message
  std::unique_ptr<OutputPrintPolicy> m_wrappee;
};

/// @brief decorate debug message with a name
///
/// The debug message is complemented with a name.
class NamedOutputDecorator final : public OutputDecorator {
 public:
  /// @brief constructor
  ///
  /// @param [in] wrappee  output print policy object to be wrapped
  /// @param [in] name     name to be added to debug message
  /// @param [in] maxWidth maximum width of field used for name
  NamedOutputDecorator(std::unique_ptr<OutputPrintPolicy> wrappee,
                       const std::string& name, unsigned int maxWidth = 15)
      : OutputDecorator(std::move(wrappee)),
        m_name(name),
        m_maxWidth(maxWidth) {}

  /// @brief flush the debug message to the destination stream
  ///
  /// @param [in] lvl   debug level of debug message
  /// @param [in] input text of debug message
  ///
  /// This function prepends the given name to the debug message and then
  /// delegates the flushing of the whole message to its wrapped object.
  void flush(const Level& lvl, const std::string& input) override {
    std::ostringstream os;
    os << std::left << std::setw(static_cast<int>(m_maxWidth))
       << m_name.substr(0, m_maxWidth - 3) << input;
    OutputDecorator::flush(lvl, os.str());
  }

  /// Make a copy of this print policy with a new name
  /// @param name the new name
  /// @return the copy
  std::unique_ptr<OutputPrintPolicy> clone(
      const std::string& name) const override {
    return std::make_unique<NamedOutputDecorator>(m_wrappee->clone(name), name,
                                                  m_maxWidth);
  }

  /// Get this named output decorators name
  /// @return the name
  const std::string& name() const override { return m_name; }

 private:
  /// name to be prepended
  std::string m_name;

  /// maximum width of field for printing the name
  unsigned int m_maxWidth;
};

/// @brief decorate debug message with a time stamp
///
/// The debug message is complemented with a time stamp.
class TimedOutputDecorator final : public OutputDecorator {
 public:
  /// @brief constructor
  ///
  /// @param [in] wrappee output print policy object to be wrapped
  /// @param [in] format  format of time stamp (see std::strftime)
  explicit TimedOutputDecorator(std::unique_ptr<OutputPrintPolicy> wrappee,
                                const std::string& format = "%X")
      : OutputDecorator(std::move(wrappee)), m_format(format) {}

  /// @brief flush the debug message to the destination stream
  ///
  /// @param [in] lvl   debug level of debug message
  /// @param [in] input text of debug message
  ///
  /// This function prepends a time stamp to the debug message and then
  /// delegates the flushing of the whole message to its wrapped object.
  void flush(const Level& lvl, const std::string& input) override {
    std::ostringstream os;
    os << std::left << std::setw(12) << now() << input;
    OutputDecorator::flush(lvl, os.str());
  }

  /// Make a copy of this print policy with a new name
  /// @param name the new name
  /// @return the copy
  std::unique_ptr<OutputPrintPolicy> clone(
      const std::string& name) const override {
    return std::make_unique<TimedOutputDecorator>(m_wrappee->clone(name),
                                                  m_format);
  }

 private:
  /// @brief get current time stamp
  ///
  /// @return current time stamp as string
  std::string now() const {
    char buffer[20];
    time_t t{};
    std::time(&t);
    struct tm tbuf{};
    std::strftime(buffer, sizeof(buffer), m_format.c_str(),
                  localtime_r(&t, &tbuf));
    return buffer;
  }

  /// format of the time stamp (see std::strftime for details)
  std::string m_format;
};

/// @brief decorate debug message with a thread ID
///
/// The debug message is complemented with a thread ID.
class ThreadOutputDecorator final : public OutputDecorator {
 public:
  /// @brief constructor
  ///
  /// @param [in] wrappee output print policy object to be wrapped
  explicit ThreadOutputDecorator(std::unique_ptr<OutputPrintPolicy> wrappee)
      : OutputDecorator(std::move(wrappee)) {}

  /// @brief flush the debug message to the destination stream
  ///
  /// @param [in] lvl   debug level of debug message
  /// @param [in] input text of debug message
  ///
  /// This function prepends the thread ID to the debug message and then
  /// delegates the flushing of the whole message to its wrapped object.
  void flush(const Level& lvl, const std::string& input) override {
    std::ostringstream os;
    os << std::left << std::setw(20) << std::this_thread::get_id() << input;
    OutputDecorator::flush(lvl, os.str());
  }

  /// Make a copy of this print policy with a new name
  /// @param name the new name
  /// @return the copy
  std::unique_ptr<OutputPrintPolicy> clone(
      const std::string& name) const override {
    return std::make_unique<ThreadOutputDecorator>(m_wrappee->clone(name));
  }
};

/// @brief decorate debug message with its debug level
///
/// The debug message is complemented with its debug level.
class LevelOutputDecorator final : public OutputDecorator {
 public:
  /// @brief constructor
  ///
  /// @param [in] wrappee output print policy object to be wrapped
  explicit LevelOutputDecorator(std::unique_ptr<OutputPrintPolicy> wrappee)
      : OutputDecorator(std::move(wrappee)) {}

  /// @brief flush the debug message to the destination stream
  ///
  /// @param [in] lvl   debug level of debug message
  /// @param [in] input text of debug message
  ///
  /// This function prepends the debug level to the debug message and then
  /// delegates the flushing of the whole message to its wrapped object.
  void flush(const Level& lvl, const std::string& input) override {
    std::ostringstream os;
    os << std::left << std::setw(10) << toString(lvl) << input;
    OutputDecorator::flush(lvl, os.str());
  }

  /// Make a copy of this print policy with a new name
  /// @param name the new name
  /// @return the copy
  std::unique_ptr<OutputPrintPolicy> clone(
      const std::string& name) const override {
    return std::make_unique<LevelOutputDecorator>(m_wrappee->clone(name));
  }

 private:
  /// @brief convert debug level to string
  ///
  /// @param [in] lvl debug level
  ///
  /// @return string representation of debug level
  std::string toString(const Level& lvl) const {
    static const char* const buffer[] = {"VERBOSE", "DEBUG", "INFO",
                                         "WARNING", "ERROR", "FATAL"};
    return buffer[lvl];
  }
};

/// @brief default print policy for debug messages
///
/// This class allows to print debug messages without further modifications to
/// a specified output stream.
class DefaultPrintPolicy final : public OutputPrintPolicy {
 public:
  /// @brief constructor
  ///
  /// @param [in] out pointer to output stream object
  ///
  /// @pre @p out is non-zero
  explicit DefaultPrintPolicy(std::ostream* out = &std::cout) : m_out(out) {}

  /// @brief flush the debug message to the destination stream
  ///
  /// @param [in] lvl   debug level of debug message
  /// @param [in] input text of debug message
  void flush(const Level& /*lvl*/, const std::string& input) final {
    // Mutex to serialize access to std::cout
    static std::mutex s_stdoutMutex;
    std::unique_lock lock{s_stdoutMutex,
                          std::defer_lock};  // prep empty, we might not need it

    if (m_out == &std::cout) {
      lock.lock();  // lock only if we are printing to std::cout
    }

    (*m_out) << input << std::endl;
  }

  /// Fulfill @c OutputPrintPolicy interface. This policy doesn't actually have a
  /// name, so the assumption is that somewhere in the decorator hierarchy,
  /// there is something that returns a name without delegating to a wrappee,
  /// before reaching this overload.
  /// @note This method will throw an exception
  /// @return the name, but it never returns
  const std::string& name() const override {
    throw std::runtime_error{
        "Default print policy doesn't have a name. Is there no named output in "
        "the decorator chain?"};
  };

  /// Make a copy of this print policy with a new name
  /// @return the copy
  std::unique_ptr<OutputPrintPolicy> clone(
      const std::string& /*name*/) const override {
    return std::make_unique<DefaultPrintPolicy>(m_out);
  };

 private:
  /// pointer to destination output stream
  std::ostream* m_out;
};

/// @}

}  // namespace Logging

/// @brief class for printing debug output
/// @ingroup logging
///
/// This class provides the user interface for printing debug messages with
/// different levels of severity.
///
class Logger {
 public:
  /// @brief construct from output print and filter policy
  ///
  /// @param [in] pPrint  policy for printing debug messages
  /// @param [in] pFilter policy for filtering debug messages
  /// @param [in] failureThreshold level at or above which a message raises
  ///                               @ref Logging::ThresholdFailure.
  ///                               @ref Logging::Level::MAX, the default,
  ///                               never fails.
  Logger(std::unique_ptr<Logging::OutputPrintPolicy> pPrint,
         std::unique_ptr<Logging::OutputFilterPolicy> pFilter,
         Logging::Level failureThreshold = Logging::Level::MAX)
      : m_printPolicy(std::move(pPrint)),
        m_filterPolicy(std::move(pFilter)),
        m_failureThreshold(failureThreshold) {}

  /// @brief decide whether a message with a given debug level has to be printed
  ///
  /// A message that exceeds the failure threshold is always printed, so that
  /// the message which fails the job is visible regardless of the filter level.
  ///
  /// @param [in] lvl debug level of debug message
  ///
  /// @return @c true if debug message should be printed, otherwise @c false
  bool doPrint(const Logging::Level& lvl) const {
    return m_filterPolicy->doPrint(lvl) || exceedsFailureThreshold(lvl);
  }

  /// @brief log a debug message
  ///
  /// @param [in] lvl debug level of debug message
  /// @param [in] input text of debug message
  /// @throws Logging::ThresholdFailure if @p lvl is at or above the failure
  ///         threshold of this logger
  void log(const Logging::Level& lvl, const std::string& input) const {
    const bool fail = exceedsFailureThreshold(lvl);
    if (fail || m_filterPolicy->doPrint(lvl)) {
      m_printPolicy->flush(lvl, input);
    }
    if (fail) {
      throw Logging::ThresholdFailure(
          "Log message at level " + std::string{Logging::levelName(lvl)} +
          " exceeds the ACTS_LOG_FAILURE_THRESHOLD of logger '" +
          diagnosticName() + "', bailing out. See " +
          "https://cern.ch/acts-log-thresh");
    }
  }

  /// @brief log a message without the failure-threshold check
  ///
  /// For call sites that must not throw, in particular destructors, where a
  /// @ref Logging::ThresholdFailure would terminate the process. The message is
  /// emitted exactly as @ref log would emit it, but it can never fail the job.
  ///
  /// @param [in] lvl debug level of debug message
  /// @param [in] input text of debug message
  void logWithoutFailure(const Logging::Level& lvl,
                         const std::string& input) const noexcept {
    try {
      if (m_filterPolicy->doPrint(lvl)) {
        m_printPolicy->flush(lvl, input);
      }
    } catch (const std::exception&) {
      // a print policy that throws must not take the process down either
    }
  }

  /// @brief The failure threshold of this logger
  ///
  /// @return the level at or above which a message raises
  ///         @ref Logging::ThresholdFailure; @ref Logging::Level::MAX if this
  ///         logger never fails
  Logging::Level failureThreshold() const { return m_failureThreshold; }

  /// Return the print policy for this logger
  /// @return the print policy
  const Logging::OutputPrintPolicy& printPolicy() const {
    return *m_printPolicy;
  }

  /// Return the filter policy for this logger
  /// @return the filter policy
  const Logging::OutputFilterPolicy& filterPolicy() const {
    return *m_filterPolicy;
  }

  /// Return the level of the filter policy of this logger
  /// @return the level
  Logging::Level level() const { return m_filterPolicy->level(); }

  /// Return the name of the print policy of this logger
  /// @return the name
  const std::string& name() const { return m_printPolicy->name(); }

  /// Make a copy of this logger, optionally changing the name or the level
  /// @param _name the optional new name
  /// @param _level the optional new level
  /// @return Unique pointer to a cloned logger
  std::unique_ptr<Logger> clone(
      const std::optional<std::string>& _name = std::nullopt,
      const std::optional<Logging::Level>& _level = std::nullopt) const {
    return std::make_unique<Logger>(
        m_printPolicy->clone(_name.has_value() ? *_name : name()),
        m_filterPolicy->clone(_level.has_value() ? *_level : level()),
        m_failureThreshold);
  }

  /// Make a copy of this logger with a different failure threshold.
  ///
  /// @note Like @ref clone this copies the policy objects, so a stateful print
  ///       policy starts over in the copy.
  ///
  /// @param _failureThreshold the level at or above which messages should raise
  ///        @ref Logging::ThresholdFailure, or @ref Logging::Level::MAX to
  ///        never fail
  /// @return Unique pointer to a cloned logger
  std::unique_ptr<Logger> withFailureThreshold(
      Logging::Level _failureThreshold) const {
    return std::make_unique<Logger>(m_printPolicy->clone(diagnosticName()),
                                    m_filterPolicy->clone(level()),
                                    _failureThreshold);
  }

  /// Make a copy of this logger that never raises @ref
  /// Logging::ThresholdFailure. Use this for code paths that log an error on
  /// purpose, such as a unit test exercising an error path.
  ///
  /// @note Like @ref clone this copies the policy objects, so a stateful print
  ///       policy starts over in the copy.
  ///
  /// @return Unique pointer to a cloned logger
  std::unique_ptr<Logger> withoutFailureThreshold() const {
    return withFailureThreshold(Logging::Level::MAX);
  }

  /// Make a copy of the logger, with a new level. Convenience function for
  /// if you only want to change the level but not the name.
  /// @param _level the new level
  /// @return the new logger
  std::unique_ptr<Logger> clone(Logging::Level _level) const {
    return clone(std::nullopt, _level);
  }

  /// Make a copy of the logger, with a suffix added to the end of it's
  /// name. You can also optionally supply a new level
  /// @param suffix the suffix to add to the end of the name
  /// @param _level the optional new level
  /// @return Unique pointer to a cloned logger with modified name
  std::unique_ptr<Logger> cloneWithSuffix(
      const std::string& suffix,
      std::optional<Logging::Level> _level = std::nullopt) const {
    return clone(name() + suffix, _level.value_or(level()));
  }

  /// Helper function so a logger reference can be used as is with the logging
  /// macros
  /// @return Reference to this logger
  const Logger& operator()() const { return *this; }

 private:
  /// @brief whether a message at @p lvl trips this logger's failure threshold
  bool exceedsFailureThreshold(Logging::Level lvl) const {
    return lvl >= m_failureThreshold;
  }

  /// @brief best-effort name for diagnostics
  ///
  /// @ref name is allowed to throw when the print policy chain carries no name,
  /// so neither the failure message nor a clone may depend on it succeeding.
  /// @return the logger name, or a placeholder
  std::string diagnosticName() const {
    try {
      return name();
    } catch (const std::exception&) {
      return "<unnamed>";
    }
  }

  /// policy object for printing debug messages
  std::unique_ptr<Logging::OutputPrintPolicy> m_printPolicy;

  /// policy object for filtering debug messages
  std::unique_ptr<Logging::OutputFilterPolicy> m_filterPolicy;

  /// level at or above which a message raises Logging::ThresholdFailure;
  /// Logging::Level::MAX means this logger never fails
  Logging::Level m_failureThreshold;
};

/// @brief get default debug output logger
///
/// @param [in] name       name of the logger instance
/// @param [in] lvl        debug threshold level
/// @param [in] log_stream output stream used for printing debug messages
///
/// This function returns a pointer to a Logger instance with the following
/// decorations enabled:
/// - time stamps
/// - name of logging instance
/// - debug level
///
/// @param [in] failureThreshold level at or above which a message raises
///                               @ref Acts::Logging::ThresholdFailure, or
///                               `std::nullopt` to take the process-wide
///                               default from
///                               @ref Acts::Logging::detail::getDefaultFailureThreshold
///
/// @return pointer to logging instance
std::unique_ptr<const Logger> getDefaultLogger(
    const std::string& name, const Logging::Level& lvl,
    std::ostream* log_stream = &std::cout,
    std::optional<Logging::Level> failureThreshold = std::nullopt);

/// Get a dummy logger that discards all output
/// @return Reference to dummy logger instance
const Logger& getDummyLogger();

}  // namespace Acts
