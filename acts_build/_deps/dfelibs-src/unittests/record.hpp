// SPDX-License-Identifier: MIT

/// \file
/// \brief Example record shared among namedtuple tests

#pragma once

#include <cstdint>

#include "dfe/dfe_namedtuple.hpp"

struct Record {
  int16_t x = 0;
  int32_t y = 0;
  int64_t z = 0;
  uint64_t a = 0;
  bool THIS_IS_UNUSED = false;
  float b = 0;
  double c = 0;
  bool d = false;

  DFE_NAMEDTUPLE(Record, x, y, z, a, b, c, d)
};

Record
make_record(size_t i) {
  Record r;
  r.x = i;
  r.y = -2 * i;
  r.z = 4 * i;
  r.a = 8 * i;
  r.THIS_IS_UNUSED = ((i % 2) == 0);
  r.b = 0.23126121f * i;
  r.c = -42.53425 * i;
  r.d = ((i % 2) != 0);
  return r;
}
