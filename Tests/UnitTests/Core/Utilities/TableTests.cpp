// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/Table.hpp"

#include <sstream>
#include <stdexcept>
#include <string>

using namespace Acts;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(UtilitiesSuite)

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

BOOST_AUTO_TEST_CASE(TableStringAlignment) {
  Table table;

  // Test string alignment overload with full names
  table.addColumn("Left", "{}", "left");
  table.addColumn("Right", "{}", "right");
  table.addColumn("Center", "{}", "center");

  table.addRow("Item1", "Value1", "Data1");
  table.addRow("LongerItem", "V2", "D2");

  std::string result = table.toString();

  // Expected output with string-based alignments
  std::string expected =
      "| Left       |  Right | Center |\n"
      "|:-----------|-------:|:------:|\n"
      "| Item1      | Value1 | Data1  |\n"
      "| LongerItem |     V2 |   D2   |\n";

  BOOST_CHECK_EQUAL(result, expected);
}

BOOST_AUTO_TEST_CASE(TableStringAlignmentShort) {
  Table table;

  // Test string alignment overload with short names (case insensitive)
  table.addColumn("Col1", "{}", "L");
  table.addColumn("Col2", "{}", "r");
  table.addColumn("Col3", "{}", "C");

  table.addRow("A", "B", "C");
  table.addRow("XX", "YY", "ZZ");

  std::string result = table.toString();

  // Expected output with short alignment strings
  std::string expected =
      "| Col1 | Col2 | Col3 |\n"
      "|:-----|-----:|:----:|\n"
      "| A    |    B |  C   |\n"
      "| XX   |   YY |  ZZ  |\n";

  BOOST_CHECK_EQUAL(result, expected);
}

BOOST_AUTO_TEST_CASE(TableInvalidAlignment) {
  Table table;

  // Test invalid alignment string
  BOOST_CHECK_THROW(table.addColumn("Test", "{}", "invalid"),
                    std::invalid_argument);
  BOOST_CHECK_THROW(table.addColumn("Test", "{}", "middle"),
                    std::invalid_argument);
  BOOST_CHECK_THROW(table.addColumn("Test", "{}", ""), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(TableMarkdownMode) {
  Table table;
  table.addColumn("Name", "{}", "left");
  table.addColumn("Value", "{}", "right");

  table.addRow("Item1", "A");
  table.addRow("Item2", "BB");

  // Test with markdown markers (default)
  std::string withMarkers = table.toString();
  std::string expectedWithMarkers =
      "| Name  | Value |\n"
      "|:------|------:|\n"
      "| Item1 |     A |\n"
      "| Item2 |    BB |\n";

  BOOST_CHECK_EQUAL(withMarkers, expectedWithMarkers);

  // Test without markdown markers
  table.setMarkdownMode(false);
  std::string withoutMarkers = table.toString();
  std::string expectedWithoutMarkers =
      "| Name  | Value |\n"
      "|-------|-------|\n"
      "| Item1 |     A |\n"
      "| Item2 |    BB |\n";

  BOOST_CHECK_EQUAL(withoutMarkers, expectedWithoutMarkers);
}

BOOST_AUTO_TEST_CASE(TableOstreamOperator) {
  Table table;
  table.addColumn("Col1", "{}", "left");
  table.addColumn("Col2", "{}", "center");

  table.addRow("A", "B");
  table.addRow("C", "D");

  // Test ostream operator
  std::ostringstream oss;
  oss << table;

  std::string expected =
      "| Col1 | Col2 |\n"
      "|:-----|:----:|\n"
      "| A    |  B   |\n"
      "| C    |  D   |\n";

  BOOST_CHECK_EQUAL(oss.str(), expected);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
