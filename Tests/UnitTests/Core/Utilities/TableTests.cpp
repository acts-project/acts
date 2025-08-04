// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/Table.hpp"

#include <stdexcept>
#include <string>

namespace Acts::Test {

BOOST_AUTO_TEST_SUITE(Utilities)

BOOST_AUTO_TEST_CASE(TableBasicUsage) {
  Table table;

  // Add columns with different alignments and formats
  table.addColumn("Name", "{}", Table::Alignment::Left);
  table.addColumn("Value", "{:.2f}", Table::Alignment::Right);
  table.addColumn("Count", "{}", Table::Alignment::Center);

  // Add rows
  table.addRow("Item1", 1.234, 42);
  table.addRow("LongerItemName", 5.678, 7);

  std::string result = table.toString();

  // Expected output with exact formatting and alignment
  std::string expected =
      "| Name           | Value | Count |\n"
      "|:---------------|------:|:-----:|\n"
      "| Item1          |  1.23 |  42   |\n"
      "| LongerItemName |  5.68 |   7   |\n";

  BOOST_CHECK_EQUAL(result, expected);
}

BOOST_AUTO_TEST_CASE(TableMarkdownFormat) {
  Table table;
  table.addColumn("Algorithm", "{}", Table::Alignment::Left);
  table.addColumn("Time", "{}", Table::Alignment::Right);
  table.addColumn("Percentage", "{:.1f}%", Table::Alignment::Center);

  table.addRow("Reader", "1.2ms", 45.6);
  table.addRow("Processor", "800us", 30.2);

  std::string result = table.toString();

  // Expected output with exact alignment indicators
  std::string expected =
      "| Algorithm |  Time | Percentage |\n"
      "|:----------|------:|:----------:|\n"
      "| Reader    | 1.2ms |   45.6%    |\n"
      "| Processor | 800us |   30.2%    |\n";

  BOOST_CHECK_EQUAL(result, expected);
}

BOOST_AUTO_TEST_CASE(TableEmptyTable) {
  Table table;

  std::string result = table.toString();
  BOOST_CHECK(result.empty());
}

BOOST_AUTO_TEST_CASE(TableArgumentCountMismatch) {
  Table table;
  table.addColumn("Name", "{}", Table::Alignment::Left);
  table.addColumn("Value", "{}", Table::Alignment::Right);

  // Too many arguments
  BOOST_CHECK_THROW(table.addRow("Test", 123, 456), std::runtime_error);

  // Too few arguments
  BOOST_CHECK_THROW(table.addRow("Test"), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(TableFormatString) {
  Table table;
  table.addColumn("Integer", "{:04d}", Table::Alignment::Right);
  table.addColumn("Float", "{:.3f}", Table::Alignment::Left);
  table.addColumn("Scientific", "{:.2e}", Table::Alignment::Center);

  table.addRow(42, 3.14159, 1234567.89);

  std::string result = table.toString();

  // Expected output with exact formatting
  std::string expected =
      "| Integer | Float | Scientific |\n"
      "|--------:|:------|:----------:|\n"
      "|    0042 | 3.142 |  1.23e+06  |\n";

  BOOST_CHECK_EQUAL(result, expected);
}

BOOST_AUTO_TEST_CASE(TableSingleColumn) {
  Table table;
  table.addColumn("Items", "{}", Table::Alignment::Center);

  table.addRow("First");
  table.addRow("Second");
  table.addRow("Third");

  std::string result = table.toString();

  // Expected output for single column with center alignment
  std::string expected =
      "| Items  |\n"
      "|:------:|\n"
      "| First  |\n"
      "| Second |\n"
      "| Third  |\n";

  BOOST_CHECK_EQUAL(result, expected);
}

BOOST_AUTO_TEST_CASE(TableMixedTypes) {
  Table table;
  table.addColumn("String", "{}", Table::Alignment::Left);
  table.addColumn("Integer", "{}", Table::Alignment::Right);
  table.addColumn("Double", "{:.2f}", Table::Alignment::Right);
  table.addColumn("Boolean", "{}", Table::Alignment::Center);

  table.addRow("Test", 42, 3.14159, true);
  table.addRow("Another", -17, -2.718, false);

  std::string result = table.toString();

  // Expected output with mixed types
  std::string expected =
      "| String  | Integer | Double | Boolean |\n"
      "|:--------|--------:|-------:|:-------:|\n"
      "| Test    |      42 |   3.14 |  true   |\n"
      "| Another |     -17 |  -2.72 |  false  |\n";

  BOOST_CHECK_EQUAL(result, expected);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Acts::Test
