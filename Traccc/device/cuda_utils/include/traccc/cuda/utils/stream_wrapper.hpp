// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

namespace traccc::cuda {

/// Wrapper class around @c cudaStream_t
///
/// It is necessary for passing around CUDA stream objects in code that should
/// not be directly exposed to the CUDA header(s).
///
/// Note that unlike @c vecmem::cuda::stream_wrapper, this type can not own a
/// stream of its own. It can only view a stream that is owned by
/// "somebody else".
///
class stream_wrapper {
 public:
  /// Wrap an existing @c cudaStream_t object, without taking ownership
  explicit stream_wrapper(void* stream);

  /// Copy constructor
  stream_wrapper(const stream_wrapper& parent) = default;
  /// Move constructor
  stream_wrapper(stream_wrapper&& parent) = default;

  /// Copy assignment
  stream_wrapper& operator=(const stream_wrapper& rhs) = default;
  /// Move assignment
  stream_wrapper& operator=(stream_wrapper&& rhs) = default;

  /// Device that the stream is associated to
  int device() const;

  /// Access a typeless pointer to the managed @c cudaStream_t object
  void* cudaStream() const;

  /// Wait for all queued tasks from the stream to complete
  void synchronize() const;

 private:
  /// Bare pointer to the wrapped @c cudaStream_t object
  void* m_stream;

};  // class stream_wrapper

}  // namespace traccc::cuda
