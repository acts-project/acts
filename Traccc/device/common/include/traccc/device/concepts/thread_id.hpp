// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <concepts>

namespace traccc::device::concepts {
/**
 * @brief Concept to ensure that a type behaves like a thread identification
 * type which allows us to access thread and block IDs. This concept assumes
 * one-dimensional grids.
 *
 * @tparam T The thread identifier-like type.
 */
template <typename T>
concept thread_id1 = requires(T& i) {
  /*
   * This function should return the local thread identifier in a *flat* way,
   * e.g. compressing two or three dimensional blocks into one dimension.
   */
  { i.getLocalThreadId() } -> std::integral;

  /*
   * This function should return the local thread identifier in the X-axis.
   */
  { i.getLocalThreadIdX() } -> std::integral;

  /*
   * This function should return the global thread identifier in a *flat*
   * way, e.g. compressing two or three dimensional blocks into one
   * dimension.
   */
  { i.getGlobalThreadId() } -> std::integral;

  /*
   * This function should return the global thread identifier in the X-axis.
   */
  { i.getGlobalThreadIdX() } -> std::integral;

  /*
   * This function should return the block identifier in the X-axis.
   */
  { i.getBlockIdX() } -> std::integral;

  /*
   * This function should return the block size in the X-axis.
   */
  { i.getBlockIdX() } -> std::integral;

  /*
   * This function should return the grid identifier in the X-axis.
   */
  { i.getBlockIdX() } -> std::integral;
};
}  // namespace traccc::device::concepts
