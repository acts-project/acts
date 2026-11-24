// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsAlignment/Kernel/Alignment.hpp"

#include "Acts/Utilities/detail/OstreamStateGuard.hpp"

#include <fstream>
#include <stdexcept>

namespace ActsAlignment {
void dumpResultAsMP2(const ActsAlignment::AlignmentResult& result,
                     std::ostream& out) {
  // protect stream from leaking modifiers
  Acts::detail::OstreamStateGuard g(out);
  for (auto [surface, index] : result.idxedAlignSurfaces) {
    for (std::size_t i = 0; i < Acts::eAlignmentSize; ++i) {
      std::size_t row = Acts::eAlignmentSize * index + i;
      out << std::setw(8)
          << row + 1
          // column 2: parameter delta from fit
          << "   " << std::setw(12)
          << result.deltaAlignmentParameters(row)
          // column 3: optional pre-sigma for parameter delta
          << "   " << std::setw(4)
          << 0.
          // column 4: change of parameter delta w.r.t start value
          << "   " << std::setw(12)
          << result.deltaAlignmentParameters(row)
          // column 5: uncertainty of parameter delta from the fit
          << "   " << std::setw(12)
          << std::sqrt(result.alignmentCovariance(row, row))
          // column 6: count of measurements seeing this parameter.
          // Not available in ACTS solver, write a placeholder of 99
          << "    99 " << std::endl;
    }
  }
}
void dumpResultAsMP2(const ActsAlignment::AlignmentResult& result,
                     const std::filesystem::path& fname) {
  std::ofstream of(fname);
  if (!of.is_open()) {
    throw std::runtime_error("Unable to open the output file");
  }
  dumpResultAsMP2(result, of);
}
}  // namespace ActsAlignment
