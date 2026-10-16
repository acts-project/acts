// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Traccc/TracccCopyTool.hpp"

#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>

#include <traccc/edm/track_container.hpp>

namespace ActsExamples {

namespace detail {
template <typename T>
struct IsTracccTrackContainer : std::false_type {};
template <typename algebra_t>
struct IsTracccTrackContainer<traccc::edm::track_container<algebra_t>>
    : std::true_type {};
}  // namespace detail

/// Copy a traccc collection from the host to the device
///
/// @tparam collection_t the traccc collection types struct, see
///         @c TracccCopyTool
template <typename collection_t>
class TracccHostToDeviceCopyAlg final : public IAlgorithm {
 public:
  using HostCollection = typename collection_t::host;
  using DeviceCollection = typename collection_t::buffer;

  struct Config {
    /// The tool doing the copies
    std::shared_ptr<const TracccCopyTool> copyTool;
    /// Input host collection
    std::string input;
    /// Output device collection
    std::string output;
  };

  explicit TracccHostToDeviceCopyAlg(
      const Config& cfg, std::unique_ptr<const Acts::Logger> logger = nullptr)
      : IAlgorithm("TracccHostToDeviceCopyAlg", std::move(logger)),
        m_cfg(cfg) {
    if (m_cfg.copyTool == nullptr) {
      throw std::invalid_argument(
          "TracccHostToDeviceCopyAlg: copyTool is not set");
    }
    m_input.initialize(m_cfg.input);
    m_output.initialize(m_cfg.output);
  }

  ProcessCode execute(const AlgorithmContext& ctx) const override {
    const HostCollection& host = m_input(ctx);
    DeviceCollection device = m_cfg.copyTool->template toDevice<collection_t>(host);
    ACTS_DEBUG("Copied " << host.size() << " elements from '" << m_cfg.input
                         << "' to the " << m_cfg.copyTool->backend().name()
                         << " device");
    m_output(ctx, std::move(device));
    return ProcessCode::SUCCESS;
  }

  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;

  ReadDataHandle<HostCollection> m_input{this, "Input"};
  WriteDataHandle<DeviceCollection> m_output{this, "Output"};
};

/// Copy a traccc collection (or track container) from the device to the host
///
/// @tparam collection_t the traccc collection types struct, see
///         @c TracccCopyTool, or a traccc track container
template <typename collection_t>
class TracccDeviceToHostCopyAlg final : public IAlgorithm {
 public:
  static constexpr bool kIsTrackContainer =
      detail::IsTracccTrackContainer<collection_t>::value;

  using HostCollection = typename collection_t::host;
  using DeviceCollection = typename collection_t::buffer;

  struct Config {
    /// The tool doing the copies
    std::shared_ptr<const TracccCopyTool> copyTool;
    /// Input device collection
    std::string input;
    /// Output host collection
    std::string output;
    /// Track containers only: optional host measurements the tracks refer
    /// to, they are attached to the output track container
    std::string inputMeasurements;
  };

  explicit TracccDeviceToHostCopyAlg(
      const Config& cfg, std::unique_ptr<const Acts::Logger> logger = nullptr)
      : IAlgorithm("TracccDeviceToHostCopyAlg", std::move(logger)),
        m_cfg(cfg) {
    if (m_cfg.copyTool == nullptr) {
      throw std::invalid_argument(
          "TracccDeviceToHostCopyAlg: copyTool is not set");
    }
    if (!kIsTrackContainer && !m_cfg.inputMeasurements.empty()) {
      throw std::invalid_argument(
          "TracccDeviceToHostCopyAlg: inputMeasurements is only used for "
          "track containers");
    }
    m_input.initialize(m_cfg.input);
    m_output.initialize(m_cfg.output);
    m_inputMeasurements.maybeInitialize(m_cfg.inputMeasurements);
  }

  ProcessCode execute(const AlgorithmContext& ctx) const override {
    const DeviceCollection& device = m_input(ctx);
    const TracccCopyTool& tool = *m_cfg.copyTool;

    if constexpr (kIsTrackContainer) {
      traccc::edm::measurement_collection::const_view measurements;
      if (m_inputMeasurements.isInitialized()) {
        measurements = vecmem::get_data(m_inputMeasurements(ctx));
      }
      HostCollection host = tool.tracksToHost(device, measurements);
      ACTS_DEBUG("Copied " << host.tracks.size() << " tracks with "
                           << host.states.size() << " states from '"
                           << m_cfg.input << "' to the host");
      m_output(ctx, std::move(host));
    } else {
      HostCollection host = tool.template toHost<collection_t>(device);
      ACTS_DEBUG("Copied " << host.size() << " elements from '" << m_cfg.input
                           << "' to the host");
      m_output(ctx, std::move(host));
    }
    return ProcessCode::SUCCESS;
  }

  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;

  ReadDataHandle<DeviceCollection> m_input{this, "Input"};
  ReadDataHandle<traccc::edm::measurement_collection::host>
      m_inputMeasurements{this, "InputMeasurements"};
  WriteDataHandle<HostCollection> m_output{this, "Output"};
};

}  // namespace ActsExamples
