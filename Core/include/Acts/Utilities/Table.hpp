// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <algorithm>
#include <cctype>
#include <format>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace Acts {

/// A utility class for creating formatted markdown tables with automatic
/// column sizing and alignment.
///
/// Usage:
/// ```cpp
/// Table table;
/// table.addColumn("Name", "{}", "left");        // String alignment
/// table.addColumn("Value", "{:.2f}", "right");  // or "r" for short
/// table.addColumn("Count", "{}", Table::Alignment::Center);  // Enum also
/// supported table.addRow("Item1", 1.23, 42); table.addRow("Item2", 4.56, 7);
/// std::cout << table.toString();
/// ```
class Table {
 public:
  /// Text alignment options for table columns.
  enum class Alignment { Left, Right, Center };

  /// Column definition used for layout and formatting.
  struct Column {
    /// Column header text
    std::string header;
    /// Format string for column values
    std::string format;
    /// Column alignment
    Alignment alignment;
    /// Column width in characters
    std::size_t width = 0;
  };

  /// Add a column with header, format string, and alignment
  /// @param header Column header text
  /// @param format Format string for the column values
  /// @param alignment Column alignment (default: Left)
  void addColumn(const std::string& header, const std::string& format,
                 Alignment alignment = Alignment::Left) {
    if (!m_rows.empty()) {
      throw std::runtime_error("Cannot add columns after rows have been added");
    }
    m_columns.push_back({header, format, alignment, header.length()});
  }

  /// Add a column with header, format string, and alignment as string
  /// @param header Column header text
  /// @param format Format string for the column values
  /// @param alignment String alignment: "left"/"l", "right"/"r", "center"/"c" (case insensitive)
  /// @throws std::invalid_argument if alignment string is not recognized
  void addColumn(const std::string& header, const std::string& format,
                 const std::string& alignment = "left") {
    addColumn(header, format, parseAlignment(alignment));
  }

  /// Set whether to include markdown alignment markers in output
  /// @param includeMarkers If true (default), includes markdown alignment markers
  void setMarkdownMode(bool includeMarkers) {
    m_includeMarkdownMarkers = includeMarkers;
  }

  /// Add a row with variable arguments matching the number of columns
  /// @param args Arguments to add as row values
  /// @throws std::runtime_error if argument count doesn't match column count
  template <typename... Args>
  void addRow(Args&&... args) {
    if (sizeof...(Args) != m_columns.size()) {
      throw std::runtime_error("Number of arguments (" +
                               std::to_string(sizeof...(Args)) +
                               ") does not match number of columns (" +
                               std::to_string(m_columns.size()) + ")");
    }

    std::vector<std::string> row;
    std::size_t colIndex = 0;

    auto addCell = [&](auto&& arg) {
      std::string formatted =
          std::vformat(m_columns[colIndex].format, std::make_format_args(arg));
      row.push_back(formatted);
      m_columns[colIndex].width =
          std::max(m_columns[colIndex].width, formatted.length());
      ++colIndex;
    };

    (addCell(args), ...);

    m_rows.push_back(std::move(row));
  }

  /// Generate the formatted table as a markdown string
  /// @return Markdown formatted string representation of the table
  std::string toString() const {
    if (m_columns.empty()) {
      return "";
    }

    std::string result;

    // Build header row
    result += "|";
    for (std::size_t i = 0; i < m_columns.size(); ++i) {
      result += std::format(
          " {} |", formatAligned(m_columns[i].header, m_columns[i].width,
                                 m_columns[i].alignment));
    }
    result += "\n";

    // Build separator row
    result += "|";
    for (std::size_t i = 0; i < m_columns.size(); ++i) {
      std::size_t contentWidth = m_columns[i].width;

      if (m_includeMarkdownMarkers) {
        switch (m_columns[i].alignment) {
          case Alignment::Left:
            result += std::format(":{:-<{}}", "", contentWidth + 1);
            break;
          case Alignment::Right:
            result += std::format("{:-<{}}:", "", contentWidth + 1);
            break;
          case Alignment::Center:
            result += std::format(":{:-<{}}:", "", contentWidth);
            break;
        }
      } else {
        result += std::format("{:-<{}}", "", contentWidth + 2);
      }
      result += "|";
    }
    result += "\n";

    // Build data rows
    for (const auto& row : m_rows) {
      result += "|";
      for (std::size_t i = 0; i < row.size(); ++i) {
        result += std::format(" {} |", formatAligned(row[i], m_columns[i].width,
                                                     m_columns[i].alignment));
      }
      result += "\n";
    }

    return result;
  }

  /// Stream output operator for Table
  friend std::ostream& operator<<(std::ostream& os, const Table& table) {
    return os << table.toString();
  }

 private:
  /// Parse alignment string to enum
  /// @throws std::invalid_argument if alignment string is not recognized
  static Alignment parseAlignment(const std::string& alignment) {
    std::string lower = alignment;
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);

    if (lower == "left" || lower == "l") {
      return Alignment::Left;
    } else if (lower == "right" || lower == "r") {
      return Alignment::Right;
    } else if (lower == "center" || lower == "c") {
      return Alignment::Center;
    } else {
      throw std::invalid_argument(
          "Invalid alignment string: '" + alignment +
          "'. Valid options: 'left'/'l', 'right'/'r', 'center'/'c'");
    }
  }

  std::string formatAligned(const std::string& text, std::size_t width,
                            Alignment alignment) const {
    switch (alignment) {
      case Alignment::Left:
        return std::format("{:<{}}", text, width);
      case Alignment::Right:
        return std::format("{:>{}}", text, width);
      case Alignment::Center:
        return std::format("{:^{}}", text, width);
    }
    return text;
  }

  std::vector<Column> m_columns;
  std::vector<std::vector<std::string>> m_rows;
  bool m_includeMarkdownMarkers = true;
};

}  // namespace Acts
