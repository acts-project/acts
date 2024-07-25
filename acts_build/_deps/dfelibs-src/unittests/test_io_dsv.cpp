// SPDX-License-Identifier: MIT

/// \file
/// \brief Unit tests for (d)elimiter-(s)eparated (v)alues i/o

#include <boost/test/unit_test.hpp>

#include "dfe/dfe_io_dsv.hpp"
#include "dfe/dfe_namedtuple.hpp"
#include "record.hpp"

BOOST_TEST_DONT_PRINT_LOG_VALUE(Record::Tuple)

#define TEST_READER_RECORDS(reader) \
  do { \
    Record record; \
    for (size_t i = 0; reader.read(record); ++i) { \
      auto expected = make_record(i); \
      BOOST_TEST( \
        record.tuple() == expected.tuple(), \
        "inconsistent record " << i << " expected=(" << expected << ") read=(" \
                               << record << ")"); \
      BOOST_TEST(record.THIS_IS_UNUSED == Record().THIS_IS_UNUSED); \
    } \
  } while (false)

// full write/read chain tests

constexpr size_t kNRecords = 1024;

BOOST_AUTO_TEST_CASE(csv_namedtuple_write_read) {
  // write some data
  {
    dfe::NamedTupleCsvWriter<Record> writer("test.csv");

    for (size_t i = 0; i < kNRecords; ++i) {
      BOOST_CHECK_NO_THROW(writer.append(make_record(i)));
    }
  }
  // read the data back
  {
    dfe::NamedTupleCsvReader<Record> reader("test.csv");

    TEST_READER_RECORDS(reader);
    BOOST_TEST(reader.num_records() == kNRecords);
  }
  // read the data back w/o verifying the header
  {
    dfe::NamedTupleCsvReader<Record> reader("test.csv", {}, false);

    TEST_READER_RECORDS(reader);
    BOOST_TEST(reader.num_records() == kNRecords);
  }
}

BOOST_AUTO_TEST_CASE(tsv_namedtuple_write_read) {
  // write some data
  {
    dfe::NamedTupleTsvWriter<Record> writer("test.tsv");

    for (size_t i = 0; i < kNRecords; ++i) {
      BOOST_CHECK_NO_THROW(writer.append(make_record(i)));
    }
  }
  // read the data back
  {
    dfe::NamedTupleTsvReader<Record> reader("test.tsv");

    TEST_READER_RECORDS(reader);
    BOOST_TEST(reader.num_records() == kNRecords);
  }
  // read the data back w/o verifying the header
  {
    dfe::NamedTupleTsvReader<Record> reader("test.tsv", {}, false);

    TEST_READER_RECORDS(reader);
    BOOST_TEST(reader.num_records() == kNRecords);
  }
}

BOOST_AUTO_TEST_CASE(tsv_untyped_write) {
  // open writer with 4 columns
  dfe::TsvWriter writer({"col0", "col1", "a", "z"}, "untyped.tsv");

  std::vector<double> values = {0.1, 2.3, 4.2};

  BOOST_CHECK_NO_THROW(writer.append(0.0, 1.0, 12u, "abc"));
  BOOST_CHECK_NO_THROW(writer.append(1, 2, "xy", "by"));
  // with vector unpacking
  BOOST_CHECK_NO_THROW(writer.append(23u, values));
  // with vector unpacking but too many entries
  values.push_back(-2.0);
  values.push_back(-34.2);
  BOOST_CHECK_THROW(writer.append(23u, values), std::invalid_argument);
  // not enough columns
  BOOST_CHECK_THROW(writer.append(1.0, 2.0, 12u), std::invalid_argument);
  // too many columns
  BOOST_CHECK_THROW(
    writer.append(1, 2, false, true, 123.2), std::invalid_argument);
}

// construct a path to an example data file
std::string
make_data_path(const char* filename) {
  std::string path = DFE_UNITTESTS_DIRECTORY;
  path += "/data/namedtuple-";
  path += filename;
  return path;
}

// read Pandas-generated files

constexpr size_t kNOnfile = 32;

// w/ reordered columns

BOOST_AUTO_TEST_CASE(csv_namedtuple_read_reordered) {
  std::string path = make_data_path("reordered_columns.csv");
  dfe::NamedTupleCsvReader<Record> reader(path);

  TEST_READER_RECORDS(reader);
  BOOST_TEST(reader.num_records() == kNOnfile);
}

BOOST_AUTO_TEST_CASE(tsv_namedtuple_read_reordered) {
  std::string path = make_data_path("reordered_columns.tsv");
  dfe::NamedTupleTsvReader<Record> reader(path);

  TEST_READER_RECORDS(reader);
  BOOST_TEST(reader.num_records() == kNOnfile);
}

// w/ reordered and additional columns

#define TEST_READER_RECORDS_EXTRA(reader, nextra) \
  do { \
    Record record; \
    std::vector<int> extra; \
    for (size_t i = 0; reader.read(record, extra); ++i) { \
      auto expected = make_record(i); \
      BOOST_TEST( \
        record.tuple() == expected.tuple(), \
        "inconsistent record " << i << " expected=(" << expected << ") read=(" \
                               << record << ")"); \
      BOOST_TEST(record.THIS_IS_UNUSED == Record().THIS_IS_UNUSED); \
      BOOST_TEST(extra.size() == nextra); \
      for (auto x : extra) { \
        BOOST_TEST(x == i); \
      } \
    } \
  } while (false)

BOOST_AUTO_TEST_CASE(csv_namedtuple_read_extra_columns) {
  dfe::NamedTupleCsvReader<Record> reader(make_data_path("extra_columns.csv"));

  TEST_READER_RECORDS_EXTRA(reader, 3);
  BOOST_TEST(reader.num_records() == kNOnfile);
  BOOST_TEST(reader.num_extra_columns() == 3);
}

BOOST_AUTO_TEST_CASE(tsv_namedtuple_read_extra_columns) {
  dfe::NamedTupleTsvReader<Record> reader(make_data_path("extra_columns.tsv"));

  TEST_READER_RECORDS_EXTRA(reader, 3);
  BOOST_TEST(reader.num_records() == kNOnfile);
  BOOST_TEST(reader.num_extra_columns() == 3);
}

// w/ optional columns

static void
test_read_with_optionals(const std::vector<std::string>& optional_columns) {
  using Reader = dfe::NamedTupleCsvReader<Record>;

  Reader reader(make_data_path("missing_columns.csv"), optional_columns);
  Record data;
  Record reference = data;

  for (std::size_t i = 0; reader.read(data); ++i) {
    Record expected = make_record(i);
    // existing columns should match expected
    BOOST_TEST(data.y == expected.y);
    BOOST_TEST(data.z == expected.z);
    BOOST_TEST(data.a == expected.a);
    BOOST_TEST(data.c == expected.c);
    BOOST_TEST(data.d == expected.d);
    // missing columns should be left untouched
    BOOST_TEST(data.x == reference.x);
    BOOST_TEST(data.b == reference.b);
  }
  BOOST_TEST(reader.num_records() == kNOnfile);
  BOOST_TEST(reader.num_extra_columns() == 0);
}

BOOST_AUTO_TEST_CASE(csv_namedtuple_read_missing_optionals) {
  BOOST_TEST_CONTEXT("optional columns b,x") {
    // b and x are actually missing in the file
    test_read_with_optionals({"b", "x"});
  }
  BOOST_TEST_CONTEXT("optional columns a,b,x,z") {
    // b and x are actually missing but a and z exist
    test_read_with_optionals({"a", "b", "x", "z"});
  }
}

// failure tests for readers

BOOST_AUTO_TEST_CASE(tsv_namedtuple_read_bad_files) {
  using Reader = dfe::NamedTupleCsvReader<Record>;

  Record r;
  BOOST_CHECK_THROW(Reader("does/not/exist.tsv"), std::runtime_error);
  // optional columns w/o header verification can not work
  BOOST_CHECK_THROW(
    Reader(make_data_path("missing_columns.tsv"), {"b", "x"}, false),
    std::runtime_error);
  // optional columns that are not part of th record definition
  BOOST_CHECK_THROW(
    Reader(make_data_path("extra_columns.tsv"), {"asfd", "xlkj", "xxdfd"}),
    std::runtime_error);
  BOOST_CHECK_THROW(
    Reader(make_data_path("bad_header1.tsv")), std::runtime_error);
  BOOST_CHECK_THROW(
    Reader(make_data_path("bad_header2.tsv")), std::runtime_error);
  BOOST_CHECK_THROW(
    Reader(make_data_path("bad_header3.tsv")), std::runtime_error);
  BOOST_CHECK_THROW(
    Reader(make_data_path("missing_columns.tsv")), std::runtime_error);
  BOOST_CHECK_THROW(
    Reader(make_data_path("too_few_columns.tsv")).read(r), std::runtime_error);
  BOOST_CHECK_THROW(
    Reader(make_data_path("too_many_columns.tsv")).read(r), std::runtime_error);
}
