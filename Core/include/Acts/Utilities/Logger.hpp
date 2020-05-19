// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
// STL include(s)
#include <ctime>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <thread>

// clang-format off
/// @brief macro to use a local Acts::Logger object
/// @ingroup Logging
///
/// @param log_object logger instance of type
//         <tt>std::unique_ptr<const Acts::Logger></tt>
///
/// @pre In the current scope, the symbol @c logger is not yet defined.
/// @post The ownership of the given @c log_object is transferred and
///       @c log_object should not be used directly any more.
///
/// This macro allows to use a locally defined logging object with the ACTS_*
/// logging macros. The envisaged usage is the following:
///
/// @code{.cpp}
/// void myFunction() {
///    std::unique_ptr<const Acts::Logger> myLogger
///        = /* .. your initialization .. */;
///    ACTS_LOCAL_LOGGER(std::move(myLogger));
///
///    ACTS_VERBOSE("hello world!");
/// }
/// @endcode
#define ACTS_LOCAL_LOGGER(log_object)                                          \
  struct __local_acts_logger                                                   \
  {                                                                            \
    __local_acts_logger(std::unique_ptr<const ::Acts::Logger> logger):         \
      m_logger(std::move(logger))                                              \
    {}                                                                         \
                                                                               \
    const ::Acts::Logger& operator()() const                                   \
    {                                                                          \
      return *m_logger;                                                        \
    }                                                                          \
                                                                               \
    std::unique_ptr<const ::Acts::Logger> m_logger;                            \
  };                                                                           \
  __local_acts_logger logger(log_object);

/// @brief macro for verbose debug output
/// @ingroup Logging
///
/// @param x debug message
///
/// @pre @c logger() must be a valid expression in the scope where this
///      macro is used and it must return a Acts::Logger object.
///
/// The debug message is printed if the current Acts::Logging::Level <=
/// Acts::Logging::VERBOSE.
#define ACTS_VERBOSE(x)                                                        \
  if (logger().doPrint(Acts::Logging::VERBOSE))                                \
    logger().log(Acts::Logging::VERBOSE) << x;

/// @brief macro for debug debug output
/// @ingroup Logging
///
/// @param x debug message
///
/// @pre @c logger() must be a valid expression in the scope where this
///      macro is used and it must return a Acts::Logger object.
///
/// The debug message is printed if the current Acts::Logging::Level <=
/// Acts::Logging::DEBUG.
#define ACTS_DEBUG(x)                                                          \
  if (logger().doPrint(Acts::Logging::DEBUG))                                  \
    logger().log(Acts::Logging::DEBUG) << x;

/// @brief macro for info debug output
/// @ingroup Logging
///
/// @param x debug message
///
/// @pre @c logger() must be a valid expression in the scope where this
///      macro is used and it must return a Acts::Logger object.
///
/// The debug message is printed if the current Acts::Logging::Level <=
/// Acts::Logging::INFO.
#define ACTS_INFO(x)                                                           \
  if (logger().doPrint(Acts::Logging::INFO))                                   \
    logger().log(Acts::Logging::INFO) << x;

/// @brief macro for warning debug output
/// @ingroup Logging
///
/// @param x debug message
///
/// @pre @c logger() must be a valid expression in the scope where this
///      macro is used and it must return a Acts::Logger object.
///
/// The debug message is printed if the current Acts::Logging::Level <=
/// Acts::Logging::WARNING.
#define ACTS_WARNING(x)                                                        \
  if (logger().doPrint(Acts::Logging::WARNING))                                \
    logger().log(Acts::Logging::WARNING) << x;

/// @brief macro for error debug output
/// @ingroup Logging
///
/// @param x debug message
///
/// @pre @c logger() must be a valid expression in the scope where this
///      macro is used and it must return a Acts::Logger object.
///
/// The debug message is printed if the current Acts::Logging::Level <=
/// Acts::Logging::ERROR.
#define ACTS_ERROR(x)                                                          \
  if (logger().doPrint(Acts::Logging::ERROR))                                  \
    logger().log(Acts::Logging::ERROR) << x;

/// @brief macro for fatal debug output
/// @ingroup Logging
///
/// @param x debug message
///
/// @pre @c logger() must be a valid expression in the scope where this
///      macro is used and it must return a Acts::Logger object.
///
/// The debug message is printed if the current Acts::Logging::Level <=
/// Acts::Logging::FATAL.
#define ACTS_FATAL(x)                                                          \
  if (logger().doPrint(Acts::Logging::FATAL))                                  \
    logger().log(Acts::Logging::FATAL) << x;
// clang-format on

namespace Acts {

/// @brief debug output related helper classes and functions
/// @ingroup Logging
namespace Logging {
/// @brief constants steering the debug output
///
/// All messages with a debug level equal or higher than the currently set
/// debug output level will be printed.
enum Level {
  VERBOSE = 0,  ///< VERBOSE level
  DEBUG,        ///< DEBUG level
  INFO,         ///< INFO level
  WARNING,      ///< WARNING level
  ERROR,        ///< ERROR level
  FATAL         ///< FATAL level
};

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
  virtual void flush(const Level& lvl, const std::ostringstream& input) = 0;
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
};

/// @brief thread-safe output stream
///
/// This classes caches the output internally and only flushes it to the
/// destination stream once it is destroyed. Using local instances of this
/// class therefore provides a thread-safe way for printing debug messages.
class OutStream final {
  /// function type for output flushing
  using OutputFunc = std::function<void(const std::ostringstream&)>;

 public:
  /// @brief construct stream object
  ///
  /// @param [in] output function object called for flushing the internal
  ///        cache
  explicit OutStream(OutputFunc output)
      : m_stream(), m_outputFunctor(std::move(output)) {}

  /// @brief copy constructor
  ///
  /// @param [in] copy stream object to copy
  OutStream(const OutStream& copy)
      : m_stream(), m_outputFunctor(copy.m_outputFunctor) {
    m_stream << copy.m_stream.str();
  }

  /// @brief destructor
  ///
  /// When calling the destructor, the internal cache is flushed using the
  /// function provided during construction.
  ~OutStream() { m_outputFunctor(m_stream); }

  /// @brief stream input operator forwarded to internal cache
  ///
  /// @tparam T input type
  ///
  /// @param [in] input content added to the stream
  template <typename T>
  OutStream& operator<<(T&& input) {
    m_stream << std::forward<T>(input);
    return *this;
  }

  /// @brief forward stream modifiers to internal cache
  ///
  /// @tparam T stream type
  ///
  /// @param [in] f stream modifier
  template <typename T>
  OutStream& operator<<(T& (*f)(T&)) {
    f(m_stream);
    return *this;
  }

 private:
  /// internal cache of stream
  std::ostringstream m_stream;

  /// output function called for flushing cache upon destruction
  OutputFunc m_outputFunctor;
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
  explicit DefaultFilterPolicy(const Level& lvl) : m_level(lvl) {}

  /// virtual default destructor
  ~DefaultFilterPolicy() override = default;

  /// @brief decide whether a debug message should be processed
  ///
  /// @param [in] lvl debug level of debug message
  ///
  /// @return @c true if @p lvl >= #m_level, otherwise @c false
  bool doPrint(const Level& lvl) const override { return m_level <= lvl; }

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
  void flush(const Level& lvl, const std::ostringstream& input) override {
    m_wrappee->flush(lvl, input);
  }

 private:
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
  void flush(const Level& lvl, const std::ostringstream& input) override {
    std::ostringstream os;
    os << std::left << std::setw(m_maxWidth) << m_name.substr(0, m_maxWidth - 3)
       << input.str();
    OutputDecorator::flush(lvl, os);
  }

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
  TimedOutputDecorator(std::unique_ptr<OutputPrintPolicy> wrappee,
                       const std::string& format = "%X")
      : OutputDecorator(std::move(wrappee)), m_format(format) {}

  /// @brief flush the debug message to the destination stream
  ///
  /// @param [in] lvl   debug level of debug message
  /// @param [in] input text of debug message
  ///
  /// This function prepends a time stamp to the debug message and then
  /// delegates the flushing of the whole message to its wrapped object.
  void flush(const Level& lvl, const std::ostringstream& input) override {
    std::ostringstream os;
    os << std::left << std::setw(12) << now() << input.str();
    OutputDecorator::flush(lvl, os);
  }

 private:
  /// @brief get current time stamp
  ///
  /// @return current time stamp as string
  std::string now() const {
    char buffer[20];
    time_t t;
    std::time(&t);
    std::strftime(buffer, sizeof(buffer), m_format.c_str(), localtime(&t));
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
  void flush(const Level& lvl, const std::ostringstream& input) override {
    std::ostringstream os;
    os << std::left << std::setw(20) << std::this_thread::get_id()
       << input.str();
    OutputDecorator::flush(lvl, os);
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
  void flush(const Level& lvl, const std::ostringstream& input) override {
    std::ostringstream os;
    os << std::left << std::setw(10) << toString(lvl) << input.str();
    OutputDecorator::flush(lvl, os);
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
  void flush(const Level& /*lvl*/, const std::ostringstream& input) final {
    (*m_out) << input.str() << std::endl;
  }

 private:
  /// pointer to destination output stream
  std::ostream* m_out;
};
}  // namespace Logging

/// @brief class for printing debug output
///
/// This class provides the user interface for printing debug messages with
/// different levels of severity.
///
/// @ingroup Logging
class Logger {
 public:
  /// @brief construct from output print and filter policy
  ///
  /// @param [in] pPrint  policy for printing debug messages
  /// @param [in] pFilter policy for filtering debug messages
  Logger(std::unique_ptr<Logging::OutputPrintPolicy> pPrint,
         std::unique_ptr<Logging::OutputFilterPolicy> pFilter)
      : m_printPolicy(std::move(pPrint)), m_filterPolicy(std::move(pFilter)) {}

  /// @brief decide whether a message with a given debug level has to be printed
  ///
  /// @param [in] lvl debug level of debug message
  ///
  /// @return @c true if debug message should be printed, otherwise @c false
  bool doPrint(const Logging::Level& lvl) const {
    return m_filterPolicy->doPrint(lvl);
  }

  /// @brief create output stream object with internal cache
  ///
  /// @param [in] lvl debug level of debug message
  ///
  /// This function creates and returns a stream object which behaves like a
  /// std::ostream and internally caches the debug message. The message will
  /// only be written to the destination stream once this stream object goes out
  /// of scope.
  ///
  /// @return output stream object with internal cache for debug message
  Logging::OutStream log(const Logging::Level& lvl) const {
    return Logging::OutStream(std::bind(&Logging::OutputPrintPolicy::flush,
                                        m_printPolicy.get(), lvl,
                                        std::placeholders::_1));
  }

 private:
  /// policy object for printing debug messages
  std::unique_ptr<Logging::OutputPrintPolicy> m_printPolicy;

  /// policy object for filtering debug messages
  std::unique_ptr<Logging::OutputFilterPolicy> m_filterPolicy;
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
/// @return pointer to logging instance
std::unique_ptr<const Logger> getDefaultLogger(
    const std::string& name, const Logging::Level& lvl,
    std::ostream* log_stream = &std::cout);

}  // namespace Acts
