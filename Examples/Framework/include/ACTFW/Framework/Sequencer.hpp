// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <Acts/Utilities/Logger.hpp>
#include <cstddef>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "ACTFW/Framework/IAlgorithm.hpp"
#include "ACTFW/Framework/IContextDecorator.hpp"
#include "ACTFW/Framework/IReader.hpp"
#include "ACTFW/Framework/IService.hpp"
#include "ACTFW/Framework/IWriter.hpp"

namespace FW {

/// A simple algorithm sequencer for event processing.
///
/// This is the backbone of the framework. It reads events from file,
/// runs the configured algorithms for each event, and writes selected data
/// back to a file.
class Sequencer {
 public:
  struct Config {
    /// number of events to skip at the beginning
    size_t skip = 0;
    /// number of events to process, SIZE_MAX to process all available events
    size_t events = SIZE_MAX;
    /// logging level
    Acts::Logging::Level logLevel = Acts::Logging::INFO;
    /// number of parallel threads to run, negative for automatic determination
    int numThreads = -1;
    /// output directory for timing information, empty for working directory
    std::string outputDir;
  };

  Sequencer(const Config& cfg);

  /// Add a service to the set of services.
  ///
  /// @throws std::invalid_argument if the service is NULL.
  void addService(std::shared_ptr<IService> service);
  /// Add a context decorator to the set of context decorators.
  ///
  /// @throws std::invalid_argument if the decorator is NULL.
  void addContextDecorator(std::shared_ptr<IContextDecorator> decorator);
  /// Add a reader to the set of readers.
  ///
  /// @throws std::invalid_argument if the reader is NULL.
  void addReader(std::shared_ptr<IReader> reader);
  /// Append an algorithm to the sequence of algorithms.
  ///
  /// @throws std::invalid_argument if the algorithm is NULL.
  void addAlgorithm(std::shared_ptr<IAlgorithm> algorithm);
  /// Add a writer to the set of writers.
  ///
  /// @throws std::invalid_argument if the writer is NULL.
  void addWriter(std::shared_ptr<IWriter> writer);

  /// Run the event loop.
  ///
  /// @return status code compatible with the `main()` return code
  /// @returns EXIT_SUCCESS when everying worked without problems
  /// @returns EXIT_FAILURE if something went wrong
  ///
  /// @note If the number of events to process is undefined, the sequencer
  /// will process events until the first reader signals the end-of-file. If
  /// given, it sets an upper bound.
  ///
  /// This function is intended to be run as the last thing in the tool
  /// main function and its return value can be used directly as the program
  /// return value, i.e.
  ///
  ///     int main(int argc, char* argv[])
  ///     {
  ///         Sequencer::Config cfg;
  ///         ... // configure the sequencer
  ///         Sequencer seq;
  ///         ... // set up the algorithms
  ///         return seq.run();
  ///     }
  ///
  /// This will run the start-of-run hook for all configured services, run all
  /// configured readers, algorithms, and writers for each event, then invoke
  /// the end-of-run hook for all configured writers.
  int run();

 private:
  /// List of all configured algorithm names.
  std::vector<std::string> listAlgorithmNames() const;
  /// Determine range of (requested) events; [SIZE_MAX, SIZE_MAX) for error.
  std::pair<size_t, size_t> determineEventsRange() const;

  Config m_cfg;
  std::vector<std::shared_ptr<IService>> m_services;
  std::vector<std::shared_ptr<IContextDecorator>> m_decorators;
  std::vector<std::shared_ptr<IReader>> m_readers;
  std::vector<std::shared_ptr<IAlgorithm>> m_algorithms;
  std::vector<std::shared_ptr<IWriter>> m_writers;
  std::unique_ptr<const Acts::Logger> m_logger;

  const Acts::Logger& logger() const { return *m_logger; }
};

}  // namespace FW
