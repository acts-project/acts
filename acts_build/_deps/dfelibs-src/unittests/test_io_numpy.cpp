// SPDX-License-Identifier: MIT

/// \file
/// \brief Unit tests for npy i/o

#include <boost/test/unit_test.hpp>

#include "dfe/dfe_io_numpy.hpp"
#include "dfe/dfe_namedtuple.hpp"
#include "record.hpp"

static constexpr size_t kNRecords = 1024;

BOOST_TEST_DONT_PRINT_LOG_VALUE(Record::Tuple)

BOOST_AUTO_TEST_CASE(numpy_namedtuple_write) {
  // write some data
  {
    dfe::NamedTupleNumpyWriter<Record> writer("test.npy");

    for (size_t i = 0; i < kNRecords; ++i) {
      BOOST_CHECK_NO_THROW(writer.append(make_record(i)));
    }
  }
}
