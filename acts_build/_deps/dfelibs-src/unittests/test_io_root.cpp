// SPDX-License-Identifier: MIT

/// \file
/// \brief Unit tests for root i/o

#include <boost/test/unit_test.hpp>

#include <cstdint>

#include "dfe/dfe_io_root.hpp"
#include "dfe/dfe_namedtuple.hpp"

static constexpr size_t kNRecords = 4096;

struct Record {
  // regular types
  float f;
  double d;
  unsigned long long int ull;
  unsigned long int ul;
  unsigned int ui;
  unsigned short int us;
  long long int ll;
  long int l;
  int i;
  short int s;
  unsigned char uc;
  signed char sc;
  char c;
  bool b;
  // stdint types
  uint64_t u64;
  uint32_t u32;
  uint16_t u16;
  uint8_t u8;
  int64_t i64;
  int32_t i32;
  int16_t i16;
  int8_t i8;
  // ROOT types
  UChar_t ru8;
  UShort_t ru16;
  UInt_t ru32;
  ULong64_t ru64;
  Char_t ri8;
  Short_t ri16;
  Int_t ri32;
  Long64_t ri64;
  Bool_t rb;

  // order is intentionally not the same as the definition order
  DFE_NAMEDTUPLE(
    Record, u64, u32, u16, u8, i64, i32, i16, i8, ull, ul, ui, us, f, d, ll, l,
    i, s, uc, sc, c, b, ru8, ru16, ru32, ru64, ri8, ri16, ri32, ri64, rb);
};

Record
make_record(size_t i) {
  Record r;
  r.f = 0.234f * i;
  r.d = -0.234 * i;
  r.ull = 64 * i;
  r.ul = 32 * i;
  r.ui = 16 * i;
  r.us = 8 * i;
  r.ll = -64 * i;
  r.l = -32 * i;
  r.i = -16 * i;
  r.s = -8 * i;
  r.uc = 'x';
  r.sc = 'y';
  r.c = 'z';
  r.b = ((i % 2) == 0);
  r.u64 = 8 * i;
  r.u32 = 4 * i;
  r.u16 = 2 * i;
  r.u8 = i;
  r.i64 = INT64_MIN;
  r.i32 = -4 * i;
  r.i16 = -2 * i;
  r.i8 = -i;
  r.ru8 = 'a';
  r.ru16 = UINT16_MAX;
  r.ru32 = UINT32_MAX;
  r.ru64 = UINT64_MAX;
  r.ri8 = 'b';
  r.ri16 = INT16_MIN;
  r.ri32 = INT32_MIN;
  r.ri64 = INT64_MIN;
  r.rb = ((i % 4) == 0u);
  return r;
}

BOOST_TEST_DONT_PRINT_LOG_VALUE(Record::Tuple)

BOOST_AUTO_TEST_CASE(root_namedtuple_write_read) {
  // write some data
  {
    dfe::NamedTupleRootWriter<Record> writer("test.root", "records");

    for (size_t i = 0; i < kNRecords; ++i) {
      BOOST_CHECK_NO_THROW(writer.append(make_record(i)));
    }
  }
  // read the data back
  {
    dfe::NamedTupleRootReader<Record> reader("test.root", "records");

    Record record;
    size_t n = 0;
    while (reader.read(record)) {
      auto expected = make_record(n);
      BOOST_TEST(
        record.tuple() == expected.tuple(), "inconsistent record "
                                              << n << "\nseen: " << record
                                              << "\nexpected: " << expected);
      n += 1;
    }
    BOOST_TEST(n == kNRecords);
  }
}

// TODO failure tests
