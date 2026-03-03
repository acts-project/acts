@defgroup logging Logging
@brief Subsystem for logging and message handling

The ACTS logging system provides flexible, hierarchical logging with multiple
severity levels and decorators. Logger objects can easily be created using the
@ref Acts::getDefaultLogger function which should be sufficient to get you
started. In case you need more customized debug output, you can make use of the
output decorators defined in `Acts::Logging` or even write your own
implementation of @ref Acts::Logging::OutputDecorator.

## Logging Levels

The logging system supports the following severity levels (from lowest to highest):

- [`VERBOSE`](@ref Acts::Logging::VERBOSE): @copybrief Acts::Logging::VERBOSE
- [`DEBUG`](@ref Acts::Logging::DEBUG): @copybrief Acts::Logging::DEBUG
- [`INFO`](@ref Acts::Logging::INFO): @copybrief Acts::Logging::INFO
- [`WARNING`](@ref Acts::Logging::WARNING): @copybrief Acts::Logging::WARNING
- [`ERROR`](@ref Acts::Logging::ERROR): @copybrief Acts::Logging::ERROR
- [`FATAL`](@ref Acts::Logging::FATAL): @copybrief Acts::Logging::FATAL

## Using Logging Macros

@copydetails logging_macros

## Logging Patterns {#logging_patterns}

ACTS provides several patterns for integrating logging into your code, each
suited for different use cases.

### Member Logger Pattern

Use this pattern when a class needs persistent logging throughout its lifetime.
The logger is stored as a member variable and passed to the constructor.

**Best for:**

- Classes that perform ongoing work with logging needs
- Components that need consistent logger identity
- Objects with well-defined lifetimes
- Allowing caller to control logging configuration

@snippet{trimleft} examples/logging.cpp Member Logger Pattern

**Key characteristics:**

- Logger stored as `std::unique_ptr<const Acts::Logger>`
- Provides `logger()` accessor method returning const reference
- Constructor accepts `std::unique_ptr<const Acts::Logger>` by value and moves it
- Caller creates logger with `getDefaultLogger()` or other logger factory
- ACTS logging macros work directly in member functions
- Ownership transferred to the class instance

### Const Reference Argument Pattern

Use this pattern when a function or algorithm should accept a logger from the
caller, providing maximum flexibility.

**Best for:**

- Standalone functions and algorithms
- Places where it's inconvenient to store a member logger

@snippet{trimleft} examples/logging.cpp Const Ref Argument Pattern

**Key characteristics:**

- Logger passed as `const Acts::Logger&` parameter
- Often has default value using `getDummyLogger()`
- Allows caller to choose logging behavior
- Supports dependency injection for testing

### `getDummyLogger` Pattern

Use this pattern when logging is optional or should be disabled by default. In
contract to @ref getDefaultLogger, this function returns a logger that discards
all messages, resulting in negligible runtime overhead. It is useful in cases
where you don't want to default to logging to `std::cout`, but not all
call-sites give you easy access to a logger.

**Best for:**

- Optional logging that can be enabled when debugging
- Functions where most callers don't have a logger

@snippet{trimleft} examples/logging.cpp getDummyLogger Pattern

**Key characteristics:**

- Returns a global dummy logger that discards all output
- negligible runtime overhead when used (messages not processed)
- Common default for function parameters
- Can be overridden with real logger when needed

### `getDefaultLogger` Factory Function

Use this function to create standalone loggers with standard formatting
(timestamp, component name, log level). Note that the loggers returned by this
function will always look to standard output streams (e.g. `std::cout`).

## Advanced Topics

### Custom Output Streams

@ref Acts::getDefaultLogger accepts an optional output stream parameter:

@snippet{trimleft} examples/logging.cpp Custom Output Streams

### Logger Cloning

Loggers can be cloned to create independent instances with optional
modifications to the name and/or log level. This is particularly useful when
creating sub-component loggers or when you need multiple loggers with similar
configurations.

@snippet{trimleft} examples/logging.cpp Logger Cloning

#### Common Use Cases

**Sub-component loggers:** When a class has multiple internal components that need separate logging:

@snippet{trimleft} examples/logging.cpp Logger Cloning Sub-component

**Per-component log levels:** When building detectors or systems with different verbosity needs:

@snippet{trimleft} examples/logging.cpp Logger Cloning Per-component Levels

**Testing:** Creating test loggers with specific configurations:

@snippet{trimleft} examples/logging.cpp Logger Cloning Testing

## Logger Integration

In case you are using ACTS in another framework which comes with its own
logging facility (e.g. Gaudi) you can pipe the logging output from ACTS
tools and algorithms to your framework's logging system by supplying different
implementations of:

- @ref Acts::Logging::OutputFilterPolicy (for mapping logging levels)
- @ref Acts::Logging::OutputPrintPolicy (for passing the Acts output
  to your internal logging system)

There are two approaches to logger integration:

1. **Overriding [`getDefaultLogger`](@ref Acts::getDefaultLogger)** (now **discouraged**).
   This has the downside that log levels cannot be controlled from top-level
   experiment specific code. This means that it is non-trivial to steer the log
   level of an e.g. Gaudi algorithm via the `OutputLevel` property, and have
   the ACTS code respect this log level.

   @warning ACTS code has iteratively moved to not construct loggers via
   @ref Acts::getDefaultLogger as much as possible, in favor of using a
   const-reference to @ref Acts::Logger. The latter can be defaulted to a
   dummy logger using @ref Acts::getDummyLogger. It is more suitable to
   pass into functions that might be called from other ACTS functions (rather
   than constructing a local logger via `getDefaultLogger`, or creating logger
   instances on the fly).

   Since ACTS makes extensive use of @ref Acts::getDefaultLogger to provide
   sufficient information for debugging, you might want to provide a modified
   implementation of this function (using your output filter and printing
   policies) to also pipe this output to your framework. You can use the
   `LD_PRELOAD` approach by providing an appropriate implementation for a
   function of the following signature into a separate source file and compile
   it in a shared library:

   @snippet{trimleft} examples/logging.cpp Logger preload

   Then you can run your executable with:

   ```console
   LD_PRELOAD=<YOUR_SHARED_LIBRARY> path/to/your/executable
   ```

2. **Passing logger instances to high level components** (recommended), and rely on ACTS code to
   pass them into lower level classes / functions. This is the preferred approach
   as it allows proper control of log levels from top-level experiment code.

## Logging Thresholds

@copydetails logging_thresholds

Two main functions exist to interact with the failure threshold:

- @ref Acts::Logging::getFailureThreshold
- @ref Acts::Logging::setFailureThreshold
