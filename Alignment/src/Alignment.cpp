// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsAlignment/Kernel/Alignment.hpp"

#include "Acts/Utilities/detail/OstreamStateGuard.hpp"
#include "ActsAlignment/Kernel/AlignmentError.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

namespace ActsAlignment {
std::vector<std::pair<const Acts::Surface*, std::size_t>>
orderedAlignmentSurfaces(const AlignmentResult& result) {
  std::vector<std::pair<const Acts::Surface*, std::size_t>> indexed(
      result.idxedAlignSurfaces.begin(), result.idxedAlignSurfaces.end());
  std::ranges::sort(indexed, std::less{},
                    [](const auto& entry) { return entry.second; });
  return indexed;
}

void dumpResultAsMP2(const ActsAlignment::AlignmentResult& result,
                     std::ostream& out) {
  // protect stream from leaking modifiers
  Acts::detail::OstreamStateGuard g(out);
  for (const auto& indexedSurface : orderedAlignmentSurfaces(result)) {
    const std::size_t index = indexedSurface.second;
    for (std::size_t i = 0; i < Acts::eAlignmentSize; ++i) {
      std::size_t row = Acts::eAlignmentSize * index + i;
      const double variance = result.alignmentCovariance(row, row);
      // A masked degree of freedom can leave a non-positive or non-finite
      // variance. Write zero instead of taking the square root.
      const double error = (variance <= 0. || !std::isfinite(variance))
                               ? 0.
                               : std::sqrt(variance);
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
          // column 5: uncertainty of the last linearization, not of the
          // total correction
          << "   " << std::setw(12)
          << error
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

Acts::Transform3 applyLocalAlignmentDelta(const Acts::Transform3& transform,
                                          const Acts::AlignmentVector& delta) {
  Acts::Transform3 updated = transform;
  updated.translate(delta.segment<3>(Acts::eAlignmentCenter0));
  const auto angles = delta.segment<3>(Acts::eAlignmentRotation0);
  updated *= Acts::AngleAxis3(angles(2), Acts::Vector3::UnitZ());
  updated *= Acts::AngleAxis3(angles(1), Acts::Vector3::UnitY());
  updated *= Acts::AngleAxis3(angles(0), Acts::Vector3::UnitX());
  return updated;
}

namespace {

/// Local step taking @p from onto @p to.
/// Inverse of @c applyLocalAlignmentDelta.
Acts::AlignmentVector localAlignmentDelta(const Acts::Transform3& from,
                                          const Acts::Transform3& to) {
  // ry = asin(-R20). rx and rz come from atan2. When |cos(ry)| is very
  // small, rz stays 0. Eigen's eulerAngles uses a different range.
  const Acts::Transform3 relative = from.inverse() * to;
  const auto rotation = relative.rotation();
  const double r20 = rotation(2, 0);
  const double ry = std::asin(std::clamp(-r20, -1.0, 1.0));
  const double cy = std::cos(ry);
  double rx = 0.;
  double rz = 0.;
  if (std::abs(cy) > 1e-8) {
    rx = std::atan2(rotation(2, 1), rotation(2, 2));
    rz = std::atan2(rotation(1, 0), rotation(0, 0));
  } else {
    rx = std::atan2(-rotation(1, 2), rotation(1, 1));
  }

  Acts::AlignmentVector delta = Acts::AlignmentVector::Zero();
  delta.segment<3>(Acts::eAlignmentCenter0) = relative.translation();
  delta[Acts::eAlignmentRotation0] = rx;
  delta[Acts::eAlignmentRotation1] = ry;
  delta[Acts::eAlignmentRotation2] = rz;
  return delta;
}

}  // namespace

Acts::Result<void> setTotalAlignmentDelta(
    AlignmentResult& result,
    const std::unordered_map<const Acts::Surface*, Acts::Transform3>&
        initialTransforms,
    const std::unordered_map<const Acts::Surface*, Acts::Transform3>&
        currentTransforms,
    bool includeLastDelta) {
  const Eigen::Index rows = result.deltaAlignmentParameters.size();
  if (rows <= 0) {
    return Acts::Result<void>::success();
  }
  const Eigen::Index width = static_cast<Eigen::Index>(Acts::eAlignmentSize);
  Acts::DynamicVector total = Acts::DynamicVector::Zero(rows);
  for (const auto& [surface, index] : result.idxedAlignSurfaces) {
    const auto initial = initialTransforms.find(surface);
    const auto current = currentTransforms.find(surface);
    const Eigen::Index offset = static_cast<Eigen::Index>(index) * width;
    if (surface == nullptr || initial == initialTransforms.end() ||
        current == currentTransforms.end() || offset < 0 ||
        offset + width > rows) {
      return AlignmentError::TotalDeltaFailure;
    }
    Acts::Transform3 finalTransform = current->second;
    if (includeLastDelta) {
      const Acts::AlignmentVector step =
          result.deltaAlignmentParameters.segment(offset, width);
      finalTransform = applyLocalAlignmentDelta(finalTransform, step);
    }
    total.segment(offset, width) =
        localAlignmentDelta(initial->second, finalTransform);
  }
  result.deltaAlignmentParameters = std::move(total);
  return Acts::Result<void>::success();
}

}  // namespace ActsAlignment
