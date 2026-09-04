// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/Logger.hpp"

#include <fstream>
#include <memory>

// Helper functions that take logger as const ref argument
void processSomething(int value,
                      const Acts::Logger &logger = Acts::getDummyLogger()) {
  ACTS_DEBUG("Processing value: " << value);
  ACTS_INFO("Process completed successfully");
}

void optionalLogging(const Acts::Logger &logger = Acts::getDummyLogger()) {
  ACTS_VERBOSE("This message will be discarded with dummy logger");
  // Dummy logger discards all output - useful for performance-critical code
}

int main([[maybe_unused]] int argc, [[maybe_unused]] char **argv) {
  //! [Member Logger Pattern]
  struct MyClass {
    std::unique_ptr<const Acts::Logger> m_logger;
    const Acts::Logger &logger() const { return *m_logger; }

    explicit MyClass(std::unique_ptr<const Acts::Logger> logger)
        : m_logger(std::move(logger)) {}

    void doWork() { ACTS_INFO("Doing work in MyClass"); }
  };

  // Usage
  MyClass obj(Acts::getDefaultLogger("MyClass", Acts::Logging::INFO));
  obj.doWork();
  //! [Member Logger Pattern]

  //! [Const Ref Argument Pattern]
  // Usage with custom logger
  auto customLogger = Acts::getDefaultLogger("Processor", Acts::Logging::DEBUG);
  processSomething(42, *customLogger);

  // Usage with default (dummy) logger
  processSomething(100);
  //! [Const Ref Argument Pattern]

  //! [getDummyLogger Pattern]
  // Call without logging overhead
  optionalLogging();

  // Or provide a real logger when needed
  auto debugLogger = Acts::getDefaultLogger("Debug", Acts::Logging::VERBOSE);
  optionalLogging(*debugLogger);
  //! [getDummyLogger Pattern]

  {
    auto verboseLogger =
        Acts::getDefaultLogger("ComponentB", Acts::Logging::VERBOSE);
    const auto &logger = *verboseLogger;

    int variable = 10;
    std::string errorMsg = "File not found";

    //! [Logging Macros]
    ACTS_VERBOSE("Detailed trace information");
    ACTS_DEBUG("Debug info: " << variable);
    ACTS_INFO("Operation completed");
    ACTS_WARNING("Potential issue detected");
    ACTS_ERROR("Operation failed: " << errorMsg);
    ACTS_FATAL("Critical failure");
    //! [Logging Macros]
  }

  //! [Logger Cloning]
  // Create a base logger
  auto baseLogger = Acts::getDefaultLogger("Fitter", Acts::Logging::INFO);

  // Clone with same name and level
  auto clonedLogger = baseLogger->clone();

  // Clone with a different name
  auto renamedLogger = baseLogger->clone("NewFitter");

  // Clone with a different log level (keeps same name)
  auto verboseClone = baseLogger->clone(Acts::Logging::VERBOSE);

  // Clone with both new name and new level
  auto customClone = baseLogger->clone("CustomFitter", Acts::Logging::DEBUG);

  // Clone with a suffix appended to the name
  auto actorLogger = baseLogger->cloneWithSuffix("Actor");
  // Result: logger named "FitterActor" with INFO level

  // Clone with a suffix and different level
  auto debugActorLogger =
      baseLogger->cloneWithSuffix("Actor", Acts::Logging::DEBUG);
  // Result: logger named "FitterActor" with DEBUG level

  // Common pattern: Create sub-component loggers
  auto updaterLogger = baseLogger->cloneWithSuffix("Updater");
  auto smootherLogger = baseLogger->cloneWithSuffix("Smoother");
  //! [Logger Cloning]

  //! [Logger Cloning Sub-component]
  // Example: Class with multiple internal components that need separate logging
  struct TrackFitter {
    std::unique_ptr<const Acts::Logger> m_logger;
    std::unique_ptr<const Acts::Logger> m_actorLogger;
    std::unique_ptr<const Acts::Logger> m_updaterLogger;

    explicit TrackFitter(std::unique_ptr<const Acts::Logger> logger)
        : m_logger(std::move(logger)),
          m_actorLogger(m_logger->cloneWithSuffix("Actor")),
          m_updaterLogger(m_logger->cloneWithSuffix("Updater")) {}
  };
  // Creates loggers: "Fitter", "FitterActor", "FitterUpdater"
  //! [Logger Cloning Sub-component]

  //! [Logger Cloning Per-component Levels]
  // Creating loggers with different verbosity levels for different components
  auto baseDetectorLogger =
      Acts::getDefaultLogger("Detector", Acts::Logging::INFO);
  auto surfaceLogger =
      baseDetectorLogger->clone("Surface", Acts::Logging::DEBUG);
  auto volumeLogger =
      baseDetectorLogger->clone("Volume", Acts::Logging::WARNING);
  //! [Logger Cloning Per-component Levels]

  //! [Logger Cloning Testing]
  // Creating test loggers with specific configurations
  auto productionLogger =
      Acts::getDefaultLogger("Production", Acts::Logging::WARNING);
  auto testLogger = productionLogger->clone("Test", Acts::Logging::VERBOSE);
  //! [Logger Cloning Testing]

  //! [Custom Output Streams]
  // Create a logger that writes to a file instead of stdout
  std::ofstream logFile("mylog.txt");
  auto fileLogger =
      Acts::getDefaultLogger("Component", Acts::Logging::INFO, &logFile);
  //! [Custom Output Streams]
}

//! [Local logger macro]
void myFunction() {
  auto myLogger = Acts::getDefaultLogger("Production", Acts::Logging::WARNING);
  ACTS_LOCAL_LOGGER(std::move(myLogger));

  ACTS_VERBOSE("hello world!");
}
//! [Local logger macro]

//! [Logger preload]
namespace Acts {
std::unique_ptr<const Logger> getDefaultLogger(const std::string &,
                                               const Logging::Level &,
                                               std::ostream *);
}
//! [Logger preload]
