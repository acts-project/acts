// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// CUDA plugin include(s).
#include "Acts/Plugins/Cuda/Utilities/StreamWrapper.hpp"
#include "Acts/Plugins/Cuda/Utilities/ErrorCheck.cuh"
#include "StreamHandlers.cuh"

// CUDA include(s).
#include <cuda_runtime.h>

namespace Acts {
namespace Cuda {

  StreamWrapper::StreamWrapper( void* stream )
  : m_stream( stream ) {

  }

  StreamWrapper::StreamWrapper( StreamWrapper&& parent )
  : m_stream( parent.m_stream ) {

    parent.m_stream = nullptr;
  }

  StreamWrapper::~StreamWrapper() {

    // Destroy the stream, if we still hold it.
    if( m_stream ) {
      //ACTS_CUDA_ERROR_CHECK( cudaStreamDestroy( getStreamFrom( *this ) ) );
    }
  }

  StreamWrapper& StreamWrapper::operator=( StreamWrapper&& rhs ) {

    // Check whether anything needs to be done.
    if( this == &rhs ) {
      return *this;
    }

    // Destroy the current stream, if we hold one.
    if( m_stream ) {
      //ACTS_CUDA_ERROR_CHECK( cudaStreamDestroy( getStreamFrom( *this ) ) );
    }

    // Perform the move.
    m_stream = rhs.m_stream;
    rhs.m_stream = nullptr;

    // Return this object.
    return *this;
  }

} // namespace Cuda
} // namespace Acts
