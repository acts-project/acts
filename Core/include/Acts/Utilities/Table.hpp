// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <format>
#include <sstream>
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
/// table.addColumn("Name", "{}", Table::Alignment::Left);
/// table.addColumn("Value", "{:.2f}", Table::Alignment::Right);
/// table.addRow("Item1", 1.23);
/// table.addRow("Item2", 4.56);
/// std::cout << table.toString();
/// ```
class Table {
 public:
  enum class Alignment { Left, Right, Center };

  struct Column {
    std::string header;
    std::string format;
    Alignment alignment;
    std::size_t width = 0;
  };

  /// Add a column with header, format string, and alignment
  void addColumn(const std::string& header, const std::string& format,
                 Alignment alignment = Alignment::Left) {
    m_columns.push_back({header, format, alignment, header.length()});
  }

  /// Add a row with variable arguments matching the number of columns
  /// @throws std::runtime_error if argument count doesn't match column count
  template <typename... Args>
  void addRow(Args&&... args) {
    if (sizeof...(Args) != m_columns.size()) {
      throw std::runtime_error(
          "Number of arguments (" + std::to_string(sizeof...(Args)) +
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
  std::string toString() const {
    if (m_columns.empty()) {
      return "";
    }

    std::ostringstream oss;

    // Print header
    oss << "|";
    for (std::size_t i = 0; i < m_columns.size(); ++i) {
      oss << " " << alignText(m_columns[i].header, m_columns[i].width, m_columns[i].alignment) << " |";
    }
    oss << "\n";

    // Print separator with alignment indicators
    oss << "|";
    for (std::size_t i = 0; i < m_columns.size(); ++i) {
      std::size_t contentWidth = m_columns[i].width;
      
      switch (m_columns[i].alignment) {
        case Alignment::Left:
          oss << ":" << std::string(contentWidth + 1, '-');
          break;
        case Alignment::Right:
          oss << std::string(contentWidth + 1, '-') << ":";
          break;
        case Alignment::Center:
          oss << ":" << std::string(contentWidth, '-') << ":";
          break;
      }
      oss << "|";
    }
    oss << "\n";

    // Print rows
    for (const auto& row : m_rows) {
      oss << "|";
      for (std::size_t i = 0; i < row.size(); ++i) {
        oss << " " << alignText(row[i], m_columns[i].width, m_columns[i].alignment) << " |";
      }
      oss << "\n";
    }

    return oss.str();
  }

 private:
  std::string alignText(const std::string& text, std::size_t width,
                        Alignment alignment) const {
    if (text.length() >= width) {
      return text;
    }

    std::size_t padding = width - text.length();

    switch (alignment) {
      case Alignment::Left:
        return text + std::string(padding, ' ');
      case Alignment::Right:
        return std::string(padding, ' ') + text;
      case Alignment::Center: {
        std::size_t leftPad = padding / 2;
        std::size_t rightPad = padding - leftPad;
        return std::string(leftPad, ' ') + text + std::string(rightPad, ' ');
      }
    }
    return text;
  }

  std::vector<Column> m_columns;
  std::vector<std::vector<std::string>> m_rows;
};

}  // namespace Acts