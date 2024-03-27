// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/CylinderVolumeStack.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <sstream>

namespace Acts {

struct CylinderVolumeStack::VolumeTuple {
  std::shared_ptr<Volume> volume;
  CylinderVolumeBounds* bounds;
  std::shared_ptr<CylinderVolumeBounds> updatedBounds;
  Transform3 localTransform;
};

CylinderVolumeStack::CylinderVolumeStack(
    std::vector<std::shared_ptr<Volume>>& volumes, BinningValue direction,
    AttachmentStrategy strategy, const Logger& logger)
    : Volume(createOuterVolume(volumes, direction, strategy, logger)),
      m_direction(direction),
      m_volumes(volumes) {}

Volume CylinderVolumeStack::createOuterVolume(
    std::vector<std::shared_ptr<Volume>>& volumes, BinningValue direction,
    AttachmentStrategy strategy, const Logger& logger) {
  // @TODO: Check no volume has avg phi etc, this is unsupported for now
  ACTS_DEBUG("Creating CylinderVolumeStack from " << volumes.size()
                                                  << " volumes");
  if (volumes.empty()) {
    throw std::invalid_argument(
        "CylinderVolumeStack requires at least one volume");
  }

  // For alignment check, we have to pick one of the volumes as the base
  const auto& groupTransform = volumes.front()->transform();
  ACTS_DEBUG("Group transform is:\n" << groupTransform.matrix());

  std::vector<VolumeTuple> volumeTuples;
  volumeTuples.reserve(volumes.size());

  for (const auto& volume : volumes) {
    auto* cylinderBounds =
        dynamic_cast<CylinderVolumeBounds*>(&volume->volumeBounds());
    if (cylinderBounds == nullptr) {
      throw std::invalid_argument{
          "CylinderVolumeStack requires all volumes to "
          "have CylinderVolumeBounds"};
    }
    volumeTuples.emplace_back(
        volume, cylinderBounds,
        std::make_shared<CylinderVolumeBounds>(*cylinderBounds),
        groupTransform.inverse() * volume->transform());
  }

  ACTS_DEBUG("*** Initial volume configuration:");
  printVolumeSequence(volumeTuples, logger);

  if (volumes.size() == 1) {
    ACTS_DEBUG("Only one volume, returning")
    return Volume(*volumes.front());
  }

  if (direction == Acts::binZ) {
    ACTS_DEBUG("Checking volume alignment");
    std::size_t n = 0;
    for (auto& vt : volumeTuples) {
      ACTS_DEBUG("Checking volume #"
                 << n << " at z: " << vt.localTransform.translation()[eZ]);
      ACTS_DEBUG("- Local transform is:\n" << vt.localTransform.matrix());

      // In the group coordinate system:

      // a) the volumes cannot have rotations meaningfully different from 1
      if (!vt.localTransform.rotation().isApprox(
              Transform3::Identity().rotation())) {
        ACTS_ERROR("Volumes are not aligned: rotation is different");
        throw std::invalid_argument(
            "Volumes are not aligned: rotation is different");
      }

      ACTS_DEBUG(" -> Rotation is ok!");

      // @TODO: What's a good tolerance here?
      constexpr auto tolerance = s_onSurfaceTolerance;

      // b) the volumes cannot have translation in x or y
      Vector2 translation = vt.localTransform.translation().head<2>();
      if (std::abs(translation[0]) > tolerance ||  //
          std::abs(translation[1]) > tolerance) {
        ACTS_ERROR("Volumes are not aligned: translation in x or y");
        throw std::invalid_argument(
            "Volumes are not aligned: translation in x or y");
      }
      ACTS_DEBUG(" -> Translation in x/y is ok!");

      n++;
    }

    ACTS_DEBUG("Sorting by volume z position");
    std::sort(volumeTuples.begin(), volumeTuples.end(),
              [](const auto& a, const auto& b) {
                return a.localTransform.translation()[eZ] <
                       b.localTransform.translation()[eZ];
              });

    ACTS_DEBUG("Checking for overlaps and attaching volumes in z");
    std::vector<VolumeTuple> gapVolumes;
    for (std::size_t i = 0; i < volumeTuples.size() - 1; i++) {
      std::size_t j = i + 1;
      auto& a = volumeTuples.at(i);
      auto& b = volumeTuples.at(j);

      auto gap = checkOverlapAndAttachInZ(a, b, strategy, logger);
      if (gap) {
        ACTS_DEBUG("Received gap volume, registering");
        auto* cylinderBounds =
            dynamic_cast<CylinderVolumeBounds*>(&gap->volumeBounds());
        if (cylinderBounds == nullptr) {
          throw std::invalid_argument{
              "CylinderVolumeStack constructed gap volume was not a cylinder "
              "volume"};
        }
        // @FIXME: This might not be safe, actually
        // Gap volume is in group coordinate system, transform to global
        Transform3 gapLocalTransform = gap->transform();
        gap->setTransform(groupTransform * gap->transform());
        gapVolumes.emplace_back(
            gap, cylinderBounds,
            std::make_shared<CylinderVolumeBounds>(*cylinderBounds),
            gapLocalTransform);
      }
    }

    ACTS_DEBUG("Appending " << gapVolumes.size()
                            << " gap volumes to the end of the volume vector");
    std::copy(gapVolumes.begin(), gapVolumes.end(),
              std::back_inserter(volumeTuples));

    ACTS_DEBUG("*** Volume configuration after z attachment:");
    printVolumeSequence(volumeTuples, logger);

    ActsScalar minR =
        std::min_element(volumeTuples.begin(), volumeTuples.end(),
                         [](const auto& a, const auto& b) {
                           return a.bounds->get(CylinderVolumeBounds::eMinR) <
                                  b.bounds->get(CylinderVolumeBounds::eMinR);
                         })
            ->bounds->get(CylinderVolumeBounds::eMinR);

    ActsScalar maxR =
        std::max_element(volumeTuples.begin(), volumeTuples.end(),
                         [](const auto& a, const auto& b) {
                           return a.bounds->get(CylinderVolumeBounds::eMaxR) <
                                  b.bounds->get(CylinderVolumeBounds::eMaxR);
                         })
            ->bounds->get(CylinderVolumeBounds::eMaxR);
    ACTS_DEBUG("Found: minR: " << minR << " maxR: " << maxR);

    for (auto& vt : volumeTuples) {
      vt.updatedBounds->set({
          {CylinderVolumeBounds::eMinR, minR},
          {CylinderVolumeBounds::eMaxR, maxR},
      });

      ACTS_DEBUG("Updated bounds for volume at z: "
                 << vt.localTransform.translation()[eZ]);
      ACTS_DEBUG(*vt.updatedBounds);

      // make a copy so we can't accidentally modify in-place
      auto copy = std::make_shared<CylinderVolumeBounds>(*vt.updatedBounds);
      vt.volume->assignVolumeBounds(std::move(vt.updatedBounds));
      vt.bounds = copy.get();
      vt.updatedBounds = std::move(copy);
    }

    ACTS_DEBUG("*** Volume configuration after r synchronization:");
    printVolumeSequence(volumeTuples, logger);

    // @TODO: Sort volume tuples again, then update the outer volume vector
    std::sort(volumeTuples.begin(), volumeTuples.end(),
              [](const auto& a, const auto& b) {
                return a.localTransform.translation()[eZ] <
                       b.localTransform.translation()[eZ];
              });

    // Update outer volume vector after sorting
    volumes.clear();
    for (const auto& vt : volumeTuples) {
      volumes.push_back(vt.volume);
    }

    ACTS_DEBUG("*** Volume configuration after final z sorting:");
    printVolumeSequence(volumeTuples, logger);

    ActsScalar minZ =
        volumeTuples.front().localTransform.translation()[eZ] -
        volumeTuples.front().bounds->get(CylinderVolumeBounds::eHalfLengthZ);
    ActsScalar maxZ =
        volumeTuples.back().localTransform.translation()[eZ] +
        volumeTuples.back().bounds->get(CylinderVolumeBounds::eHalfLengthZ);

    ActsScalar midZ = (minZ + maxZ) / 2.0;

    ACTS_DEBUG("Outer volume bounds are z: [ "
               << minZ << " <- " << midZ << " -> " << maxZ << " ], r: [ "
               << minR << " <-> " << maxR << " ]");

    Transform3 outerTransform = groupTransform * Translation3{0, 0, midZ};

    ACTS_DEBUG("Outer transform is:\n" << outerTransform.matrix());

    auto outerCylinderBounds =
        std::make_shared<CylinderVolumeBounds>(minR, maxR, (maxZ - minZ) / 2.0);
    ACTS_DEBUG("Outer bounds are:\n" << *outerCylinderBounds);

    return Volume(outerTransform, outerCylinderBounds);
  } else {
    throw std::invalid_argument(binningValueNames()[direction] +
                                " is not supported ");
  }
}

std::shared_ptr<Volume> CylinderVolumeStack::checkOverlapAndAttachInZ(
    VolumeTuple& a, VolumeTuple& b,
    CylinderVolumeStack::AttachmentStrategy strategy, const Logger& logger) {
  ActsScalar aZ = a.localTransform.translation()[eZ];
  ActsScalar aHlZ = a.bounds->get(CylinderVolumeBounds::eHalfLengthZ);
  ActsScalar aZMin = aZ - aHlZ;
  ActsScalar aZMax = aZ + aHlZ;

  ActsScalar bZ = b.localTransform.translation()[eZ];
  ActsScalar bHlZ = b.bounds->get(CylinderVolumeBounds::eHalfLengthZ);
  ActsScalar bZMin = bZ - bHlZ;
  ActsScalar bZMax = bZ + bHlZ;

  if (logger().doPrint(Acts::Logging::DEBUG)) {
    std::stringstream ss;
    ss << std::fixed;
    ss << std::setprecision(3);
    ss << std::setfill(' ');
    ACTS_DEBUG("Checking overlap between");
    int w = 9;
    ss << " - "
       << " z: [ " << std::setw(w) << aZMin << " <- " << std::setw(w) << aZ
       << " -> " << std::setw(w) << aZMax << " ]";
    ACTS_DEBUG(ss.str());

    ss.str("");
    ss << " - "
       << " z: [ " << std::setw(w) << bZMin << " <- " << std::setw(w) << bZ
       << " -> " << std::setw(w) << bZMax << " ]";
    ACTS_DEBUG(ss.str());
  }

  if (aZMax > bZMin) {
    ACTS_ERROR(" -> Overlap in z");
    throw std::invalid_argument("Volumes overlap in z");
  } else {
    ACTS_DEBUG(" -> No overlap");
  }

  std::shared_ptr<Volume> gap;

  constexpr auto tolerance = s_onSurfaceTolerance;
  if (std::abs(aZMax - bZMin) < tolerance) {
    ACTS_DEBUG("No gap between volumes, no attachment needed");
  } else {
    ACTS_DEBUG("Synchronizing bounds in z with strategy: " << strategy);

    switch (strategy) {
        // case AttachmentStrategy::Midpoint: {
        // ACTS_DEBUG(" -> Midpoint strategy");
        // // The midpoint of the two volumes
        // ActsScalar midpoint = (aZMax + bZMin) / 2.0;

        // ACTS_DEBUG(" - Midpoint: " << midpoint);

        // a.updatedBounds->set(CylinderVolumeBounds::eHalfLengthZ, halfLength);
        // a.updatedBounds->set(CylinderVolumeBounds::eZ, midpoint);
        // b.updatedBounds->set(CylinderVolumeBounds::eHalfLengthZ, halfLength);
        // b.updatedBounds->set(CylinderVolumeBounds::eZ, midpoint);
        // break;
        // }
      case AttachmentStrategy::Gap: {
        ACTS_DEBUG(" -> Gap strategy");
        ActsScalar gapHlZ = (bZMin - aZMax) / 2.0;
        ActsScalar gapMidZ = (bZMin + aZMax) / 2.0;

        Transform3 gapTransform{Translation3{0, 0, gapMidZ}};
        gap = std::make_shared<Volume>(
            gapTransform,
            std::make_shared<CylinderVolumeBounds>(
                a.bounds->get(CylinderVolumeBounds::eMinR),
                a.bounds->get(CylinderVolumeBounds::eMaxR), gapHlZ));

        break;
      }
      default:
        ACTS_ERROR("Attachment strategy " << strategy << " not implemented");
        std::stringstream ss;
        ss << strategy;
        throw std::invalid_argument("Attachment strategy " + ss.str() +
                                    " not implemented");
    }
  }

  return gap;
}

// std::shared_ptr<Volume> CylinderVolumeStack::createGapVolume(
// const VolumeTuple& a, const VolumeTuple& b, BinningValue direction,
// const Logger& logger) {
// if(direction
// // This calculation is duplicated, but the separation into a separate
// function
// // is worth it
// ActsScalar aZ = a.localTransform.translation()[eZ];
// ActsScalar aHlZ = a.bounds->get(CylinderVolumeBounds::eHalfLengthZ);
// ActsScalar aZMin = aZ - aHlZ;
// ActsScalar aZMax = aZ + aHlZ;

// ActsScalar bZ = b.localTransform.translation()[eZ];
// ActsScalar bHlZ = b.bounds->get(CylinderVolumeBounds::eHalfLengthZ);
// ActsScalar bZMin = bZ - bHlZ;
// ActsScalar bZMax = bZ + bHlZ;
// }

void CylinderVolumeStack::printVolumeSequence(
    const std::vector<VolumeTuple>& volumes, const Logger& logger) {
  if (!logger().doPrint(Acts::Logging::DEBUG)) {
    return;
  }
  for (const auto& vt : volumes) {
    std::stringstream ss;
    ss << std::fixed;
    ss << std::setprecision(3);
    ss << std::setfill(' ');
    ActsScalar z = vt.localTransform.translation()[eZ];
    ActsScalar hlZ = vt.updatedBounds->get(CylinderVolumeBounds::eHalfLengthZ);
    ActsScalar minZ = z - hlZ;
    ActsScalar maxZ = z + hlZ;

    ActsScalar minR = vt.updatedBounds->get(CylinderVolumeBounds::eMinR);
    ActsScalar maxR = vt.updatedBounds->get(CylinderVolumeBounds::eMaxR);

    int w = 9;
    ss << "z: [ " << std::setw(w) << minZ << " <- " << std::setw(w) << z
       << " -> " << std::setw(w) << maxZ << " ], r: [ " << std::setw(w) << minR
       << " <-> " << std::setw(w) << maxR << " ]";

    ACTS_DEBUG(ss.str());
  }
}

std::ostream& operator<<(std::ostream& os,
                         CylinderVolumeStack::AttachmentStrategy strategy) {
  switch (strategy) {
    case CylinderVolumeStack::AttachmentStrategy::First:
      os << "First";
      break;
    case CylinderVolumeStack::AttachmentStrategy::Second:
      os << "Second";
      break;
    case CylinderVolumeStack::AttachmentStrategy::Midpoint:
      os << "Midpoint";
      break;
    case CylinderVolumeStack::AttachmentStrategy::Gap:
      os << "Gap";
      break;
  }
  return os;
}

}  // namespace Acts
