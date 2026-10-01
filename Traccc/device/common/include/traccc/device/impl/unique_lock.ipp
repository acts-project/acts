// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

namespace traccc::device {
template <typename Mutex>
TRACCC_HOST_DEVICE unique_lock<Mutex>::unique_lock(mutex_type& m,
                                                   std::defer_lock_t) {
  m_mutex_ptr = &m;
  m_owns_lock = false;
}

template <typename Mutex>
TRACCC_HOST_DEVICE unique_lock<Mutex>::unique_lock(mutex_type& m,
                                                   std::try_to_lock_t) {
  m_mutex_ptr = &m;
  m_owns_lock = m_mutex_ptr->try_lock();
}

template <typename Mutex>
TRACCC_HOST_DEVICE unique_lock<Mutex>::unique_lock(mutex_type& m,
                                                   std::adopt_lock_t) {
  m_mutex_ptr = &m;
  m_owns_lock = true;
}

template <typename Mutex>
TRACCC_HOST_DEVICE unique_lock<Mutex>::~unique_lock() {
  if (m_owns_lock) {
    m_mutex_ptr->unlock();
  }
}

template <typename Mutex>
TRACCC_HOST_DEVICE void unique_lock<Mutex>::lock() {
  assert(!m_owns_lock);
  m_mutex_ptr->lock();
  m_owns_lock = true;
}

template <typename Mutex>
TRACCC_HOST_DEVICE bool unique_lock<Mutex>::try_lock() {
  assert(!m_owns_lock);
  m_owns_lock = m_mutex_ptr->try_lock();
  return m_owns_lock;
}

template <typename Mutex>
TRACCC_HOST_DEVICE void unique_lock<Mutex>::unlock() {
  assert(m_owns_lock);
  m_mutex_ptr->unlock();
  m_owns_lock = false;
}

template <typename Mutex>
TRACCC_HOST_DEVICE bool unique_lock<Mutex>::owns_lock() const {
  return m_owns_lock;
}
}  // namespace traccc::device
