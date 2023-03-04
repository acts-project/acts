// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

namespace Acts::detail_tc {

/// Internal holder type for referencing a backend without ownership
template <typename T>
struct RefHolder {
  T* ptr;

  RefHolder(T* _ptr) : ptr{_ptr} {}
  RefHolder(T& ref) : ptr{&ref} {}

  const T& operator*() const { return *ptr; }
  T& operator*() { return *ptr; }

  const T* operator->() const { return ptr; }
  T* operator->() { return ptr; }
};

/// Internal holder type holding a backend container by value
template <typename T>
struct ValueHolder {
  T val;

  // Let's be clear with the user that we take the ownership
  // Only require rvalues and avoid hidden copies
  ValueHolder(T& _val) = delete;
  ValueHolder(T&& _val) : val{std::move(_val)} {}

  // Does it makes sense to allow copy operations?

  const T& operator*() const { return val; }
  T& operator*() { return val; }

  const T* operator->() const { return &val; }
  T* operator->() { return &val; }
};

}  // namespace Acts::detail_tc
