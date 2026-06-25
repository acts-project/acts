// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/CompositeSpacePoint.hpp"
#include "Acts/Utilities/detail/ContainerIterator.hpp"

#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <utility>

#include <cuda_runtime.h>

namespace ActsExamples {

struct cuda_composite_space_point_arrays {
  double* local_x = nullptr;
  double* local_y = nullptr;
  double* local_z = nullptr;

  double* sensor_x = nullptr;
  double* sensor_y = nullptr;
  double* sensor_z = nullptr;

  double* next_sensor_x = nullptr;
  double* next_sensor_y = nullptr;
  double* next_sensor_z = nullptr;

  double* plane_normal_x = nullptr;
  double* plane_normal_y = nullptr;
  double* plane_normal_z = nullptr;

  double* drift_radius = nullptr;
  double* time = nullptr;

  double* cov_0 = nullptr;
  double* cov_1 = nullptr;
  double* cov_2 = nullptr;

  std::uint8_t* is_straw = nullptr;
  std::uint8_t* has_time = nullptr;
  std::uint8_t* measures_loc0 = nullptr;
  std::uint8_t* measures_loc1 = nullptr;
};

inline int cuda_test() {
  int device_count = 0;

  if (cudaGetDeviceCount(&device_count) != cudaSuccess) {
    return 1;
  }

  if (device_count == 0) {
    return 2;
  }

  double* value = nullptr;

  if (cudaMallocManaged(&value, sizeof(double)) != cudaSuccess) {
    return 3;
  }

  *value = 41.0;

  if (cudaDeviceSynchronize() != cudaSuccess) {
    cudaFree(value);
    return 4;
  }

  if (cudaFree(value) != cudaSuccess) {
    return 6;
  }

  return 0;
}

}


