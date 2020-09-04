// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Sycl/Utilities/QueueWrapper.hpp"

#include "Acts/Plugins/Sycl/Utilities/DeviceSelector.hpp"

#include <string>

#include <CL/sycl.hpp>

namespace Acts::Sycl {

QueueWrapper::QueueWrapper(const std::string& deviceNameSubtring) {
  initialize(deviceNameSubtring);
}

QueueWrapper::QueueWrapper(QueueWrapper&& parent) noexcept
    : m_queue(parent.m_queue), m_ownsQueue(parent.m_ownsQueue) {
  parent.m_queue = nullptr;
  parent.m_ownsQueue = false;
}

QueueWrapper::QueueWrapper(const QueueWrapper& other) {
  m_queue = other.m_queue;
  m_ownsQueue = false;
}

QueueWrapper::~QueueWrapper() {
  if (m_ownsQueue && (m_queue != nullptr)) {
    delete m_queue;
  }
};

QueueWrapper& QueueWrapper::operator=(QueueWrapper&& rhs) noexcept {
  // Check whether we have to do anything
  if (this == &rhs) {
    return *this;
  }

  // Destroy this queue,
  if (m_ownsQueue && (m_queue != nullptr)) {
    delete m_queue;
  }

  // Perform the move
  m_queue = rhs.m_queue;
  m_ownsQueue = rhs.m_ownsQueue;
  rhs.m_queue = nullptr;
  rhs.m_ownsQueue = false;

  // Return this object.
  return *this;
};

QueueWrapper& QueueWrapper::operator=(const QueueWrapper& other) {
  // Check whether we have to do anything
  if (this == &other) {
    return *this;
  }

  m_queue = other.m_queue;
  m_ownsQueue = false;
  return *this;
};

cl::sycl::queue* QueueWrapper::getQueue() const {
  return m_queue;
}

void QueueWrapper::initialize(const std::string& deviceNameSubstring) {
  // SYCL kernel exceptions are asynchronous
  auto exception_handler = [](cl::sycl::exception_list exceptions) {
    for (std::exception_ptr const& e : exceptions) {
      try {
        std::rethrow_exception(e);
      } catch (cl::sycl::exception const& e) {
        std::cerr << "Caught asynchronous SYCL exception:\n" << e.what();
        exit(0);
      }
    }
  };

  // Create queue with custom device selector
  m_queue = new cl::sycl::queue(DeviceSelector(deviceNameSubstring),
                                exception_handler);

  m_ownsQueue = true;

  // See which device we are running on.
  std::cerr << "Running on: "
            << m_queue->get_device().get_info<cl::sycl::info::device::name>()
            << std::endl;
};

}  // namespace Acts::Sycl
