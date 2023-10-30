// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Detector/detail/BlueprintHelper.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"

#include <array>

namespace {

std::array<Acts::Vector3, 2u> cylEndpointsZ(
    const Acts::Experimental::Blueprint::Node& node) {
  Acts::Vector3 axisZ = node.transform.rotation().col(2);
  auto halfZ = node.boundaryValues[2];
  Acts::Vector3 center = node.transform.translation();
  Acts::Vector3 p0 = center - halfZ * axisZ;
  Acts::Vector3 p1 = center + halfZ * axisZ;
  return {p0, p1};
}

}  // namespace

void Acts::Experimental::detail::BlueprintHelper::sort(Blueprint::Node& node,
                                                       bool recursive) {
  if (node.children.size() < 2u) {
    return;
  }
  // Sort along x, y, z
  if (node.binning.size() == 1) {
    auto bVal = node.binning.front();
    // x,y,z binning along the axis
    if (bVal == binX or bVal == binY or bVal == binZ) {
      Vector3 nodeCenter = node.transform.translation();
      Vector3 nodeSortAxis = node.transform.rotation().col(bVal);
      std::sort(
          node.children.begin(), node.children.end(),
          [&](const auto& a, const auto& b) {
            return (a->transform.translation() - nodeCenter).dot(nodeSortAxis) <
                   (b->transform.translation() - nodeCenter).dot(nodeSortAxis);
          });
    } else if (bVal == binR and node.boundsType == VolumeBounds::eCylinder) {
      std::sort(node.children.begin(), node.children.end(),
                [](const auto& a, const auto& b) {
                  return 0.5 * (a->boundaryValues[0] + a->boundaryValues[1]) <
                         0.5 * (b->boundaryValues[0] + b->boundaryValues[1]);
                });
    }
  }

  // Sort the children
  if (recursive) {
    for (auto& child : node.children) {
      sort(*child, true);
    }
  }
}

void Acts::Experimental::detail::BlueprintHelper::fillGaps(
    Blueprint::Node& node, bool adjustToParent) {
  // Return if this is a leaf node
  if (node.isLeaf()) {
    return;
  }

  if (node.boundsType == VolumeBounds::eCylinder) {
    // Nodes must be sorted
    sort(node, false);

    // Container values
    auto cInnerR = node.boundaryValues[0];
    auto cOuterR = node.boundaryValues[1];
    auto cHalfZ = node.boundaryValues[2];

    std::vector<std::unique_ptr<Blueprint::Node>> gaps;
    // Only 1D binning implemented for the moment
    if (node.binning.size() == 1) {
      auto bVal = node.binning.front();
      if (bVal == binZ) {
        // adjust inner/outer radius
        if (adjustToParent) {
          std::for_each(node.children.begin(), node.children.end(),
                        [&](auto& child) {
                          child->boundaryValues[0] = cInnerR;
                          child->boundaryValues[1] = cOuterR;
                        });
        }
        auto [negC, posC] = cylEndpointsZ(node);
        // Assume sorted along the local z axis
        unsigned int igap = 0;
        for (auto& child : node.children) {
          auto [neg, pos] = cylEndpointsZ(*child);
          ActsScalar gapSpan = (neg - negC).norm();
          if (gapSpan > s_onSurfaceTolerance) {
            // Fill a gap node
            auto gapName = node.name + "_gap_" + std::to_string(igap);
            auto gapTransform = Transform3::Identity();
            gapTransform.rotate(node.transform.rotation());
            gapTransform.translate(0.5 * (neg + negC));
            auto gap = std::make_unique<Blueprint::Node>(
                gapName, gapTransform, VolumeBounds::eCylinder,
                std::vector<ActsScalar>{cInnerR, cOuterR, 0.5 * gapSpan});
            gaps.push_back(std::move(gap));
            ++igap;
          }
          // Set to new current negative value
          negC = pos;
        }
        // Check if a last one needs to be filled
        ActsScalar gapSpan = (negC - posC).norm();
        if (gapSpan > s_onSurfaceTolerance) {
          // Fill a gap node
          auto gapName = node.name + "_gap_" + std::to_string(igap);
          auto gapTransform = Transform3::Identity();
          gapTransform.rotate(node.transform.rotation());
          gapTransform.translate(0.5 * (negC + posC));
          auto gap = std::make_unique<Blueprint::Node>(
              gapName, gapTransform, VolumeBounds::eCylinder,
              std::vector<ActsScalar>{cInnerR, cOuterR, 0.5 * gapSpan});
          gaps.push_back(std::move(gap));
        }

      } else if (bVal == binR) {
        // We have binning in R present
        if (adjustToParent) {
          std::for_each(node.children.begin(), node.children.end(),
                        [&](auto& child) {
                          child->transform = node.transform;
                          child->boundaryValues[2] = cHalfZ;
                        });
        }
        // Fill the gaps in R
        unsigned int igap = 0;
        ActsScalar lastR = cInnerR;
        for (auto& child : node.children) {
          ActsScalar iR = child->boundaryValues[0];
          if (std::abs(iR - lastR) > s_onSurfaceTolerance) {
            auto gap = std::make_unique<Blueprint::Node>(
                node.name + "_gap_" + std::to_string(igap), node.transform,
                VolumeBounds::eCylinder,
                std::vector<ActsScalar>{lastR, iR, cHalfZ});
            gaps.push_back(std::move(gap));
            ++igap;
          }
          // Set to new outer radius
          lastR = child->boundaryValues[1];
        }
        // Check if a last one needs to be filled
        if (std::abs(lastR - cOuterR) > s_onSurfaceTolerance) {
          auto gap = std::make_unique<Blueprint::Node>(
              node.name + "_gap_" + std::to_string(igap), node.transform,
              VolumeBounds::eCylinder,
              std::vector<ActsScalar>{lastR, cOuterR, cHalfZ});
          gaps.push_back(std::move(gap));
        }
      } else {
        throw std::runtime_error(
            "BlueprintHelper: gap filling not implemented for "
            "cylinder and this binning type.");
      }
    }
    // Insert
    for (auto& gap : gaps) {
      node.add(std::move(gap));
    }

    // Sort again after inserting
    sort(node, false);
    // Fill the gaps recursively
    for (auto& child : node.children) {
      fillGaps(*child, adjustToParent);
    }

  } else {
    throw std::runtime_error(
        "BlueprintHelper: gap filling not implemented for "
        "this boundary type");
  }
}
