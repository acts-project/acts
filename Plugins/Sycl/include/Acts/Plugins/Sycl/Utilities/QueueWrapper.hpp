// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include <string>

/// Forward declaration of incomplete type cl::sycl::queue
inline namespace cl {
namespace sycl {
class queue;
}
};  // namespace cl

namespace Acts::Sycl {

class QueueWrapper {
 public:
  /// Create queue with default selector or given name
  QueueWrapper(const std::string& = "");
  /// Move constructor
  QueueWrapper(QueueWrapper&& parent) noexcept;
  /// Copy constructor
  QueueWrapper(const QueueWrapper& other);
  /// Destructor
  ~QueueWrapper();

  /// Move assignment operator
  QueueWrapper& operator=(QueueWrapper&& rhs) noexcept;
  /// Copy assignment operator
  QueueWrapper& operator=(const QueueWrapper& other);

  /// Get stored pointer
  cl::sycl::queue* getQueue() const;

 private:
  /// Raw pointer to SYCL queue object
  cl::sycl::queue* m_queue;
  /// Owns queue
  bool m_ownsQueue;
  /// @brief This function creates the SYCL queue object.
  void initialize(const std::string& = "");

};  // class QueueWrapper
}  // namespace Acts::Sycl
