// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Traccc/TracccBackend.hpp"

#include <memory>
#include <stdexcept>
#include <type_traits>

#include <traccc/edm/measurement_collection.hpp>
#include <traccc/edm/track_container.hpp>
#include <vecmem/containers/data/buffer_type.hpp>
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/memory/memory_resource.hpp>
#include <vecmem/utils/copy.hpp>

namespace ActsExamples {

/// Host <-> device copies of traccc event data through a traccc backend
///
/// The tool is generic in the traccc collection that is copied. A collection
/// is described by its "types" struct providing the @c host and @c buffer
/// types, e.g. @c traccc::edm::measurement_collection,
/// @c traccc::edm::spacepoint_collection or
/// @c traccc::bound_track_parameters_collection_types. Track containers,
/// which hold separate track and track state collections, are supported as
/// well, see @c tracksToHost.
///
/// All copies are done on the backend stream/queue and are synchronous, i.e.
/// the result can be used right away.
class TracccCopyTool {
 public:
  explicit TracccCopyTool(std::shared_ptr<const TracccBackend> backend)
      : m_backend(std::move(backend)) {
    if (m_backend == nullptr) {
      throw std::invalid_argument("TracccCopyTool: backend is not set");
    }
  }

  const TracccBackend& backend() const { return *m_backend; }

  /// Memory resource used for collections copied to the host. This is the
  /// host-accessible (e.g. pinned) memory of the backend if it provides one.
  vecmem::memory_resource& hostMemoryResource() const {
    vecmem::memory_resource* host = m_backend->memoryResource().host;
    return host != nullptr ? *host : m_hostMr;
  }

  /// Copy a host collection into a new device buffer
  ///
  /// @tparam collection_t the traccc collection types struct
  /// @param host the host collection to copy
  /// @param type type of the device buffer. Resizable buffers carry their
  ///        size on the device, which some traccc algorithms rely on
  ///        (e.g. the measurement sorting).
  template <typename collection_t>
  typename collection_t::buffer toDevice(
      const typename collection_t::host& host,
      vecmem::data::buffer_type type =
          vecmem::data::buffer_type::resizable) const {
    vecmem::copy& copy = m_backend->copy();
    typename collection_t::buffer buffer(
        static_cast<unsigned int>(host.size()),
        m_backend->memoryResource().main, type);
    copy.setup(buffer)->wait();
    copy(vecmem::get_data(host), buffer, vecmem::copy::type::host_to_device)
        ->wait();
    return buffer;
  }

  /// Copy a device buffer into a new host collection
  ///
  /// @tparam collection_t the traccc collection types struct
  template <typename collection_t>
  typename collection_t::host toHost(
      const typename collection_t::buffer& buffer) const {
    using host_t = typename collection_t::host;
    // SoA host collections take the memory resource by reference, plain
    // vecmem vectors (std::pmr::vector) by pointer
    auto makeHost = [this]() {
      if constexpr (std::is_constructible_v<host_t, vecmem::memory_resource&>) {
        return host_t(hostMemoryResource());
      } else {
        return host_t(&hostMemoryResource());
      }
    };
    host_t host = makeHost();
    m_backend->copy()(buffer, host, vecmem::copy::type::device_to_host)
        ->wait();
    return host;
  }

  /// Copy a device track container into a new host track container
  ///
  /// A track container holds separate track and track state collections, so
  /// it can not go through the generic @c toHost.
  ///
  /// @param buffer the device tracks and track states
  /// @param measurements the (host) measurements the tracks refer to
  template <typename algebra_t = traccc::default_algebra>
  typename traccc::edm::track_container<algebra_t>::host tracksToHost(
      const typename traccc::edm::track_container<algebra_t>::buffer& buffer,
      const traccc::edm::measurement_collection::const_view& measurements =
          {}) const {
    typename traccc::edm::track_container<algebra_t>::host host{
        hostMemoryResource(), measurements};
    vecmem::copy& copy = m_backend->copy();
    copy(buffer.tracks, host.tracks, vecmem::copy::type::device_to_host)
        ->wait();
    copy(buffer.states, host.states, vecmem::copy::type::device_to_host)
        ->wait();
    return host;
  }

 private:
  std::shared_ptr<const TracccBackend> m_backend;
  mutable vecmem::host_memory_resource m_hostMr;
};

}  // namespace ActsExamples
