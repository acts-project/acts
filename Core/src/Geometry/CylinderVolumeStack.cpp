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
  ACTS_DEBUG("Creating CylinderVolumeStack from "
             << volumes.size() << " volumes in direction "
             << binningValueNames()[direction]);
  if (volumes.empty()) {
    throw std::invalid_argument(
        "CylinderVolumeStack requires at least one volume");
  }

  if (direction != Acts::binZ && direction != Acts::binR) {
    throw std::invalid_argument(binningValueNames()[direction] +
                                " is not supported ");
  }

  // For alignment check, we have to pick one of the volumes as the base
  const auto groupTransform = volumes.front()->transform();
  ACTS_VERBOSE("Group transform is:\n" << groupTransform.matrix());

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
  printVolumeSequence(volumeTuples, logger, Acts::Logging::DEBUG);

  if (volumes.size() == 1) {
    ACTS_VERBOSE("Only one volume, returning")
    return Volume(*volumes.front());
  }

  if (direction == Acts::binZ) {
    ACTS_VERBOSE("Checking volume alignment");
    std::size_t n = 0;
    for (auto& vt : volumeTuples) {
      ACTS_VERBOSE("Checking volume #"
                   << n << " at z: " << vt.localTransform.translation()[eZ]);
      ACTS_VERBOSE("- Local transform is:\n" << vt.localTransform.matrix());

      // @TODO: What's a good tolerance here?
      constexpr auto tolerance = s_onSurfaceTolerance;

      // In the group coordinate system:

      // a) the volumes cannot rotate around x or y
      if (std::abs(vt.localTransform.rotation().col(eX)[eZ]) >= tolerance ||
          std::abs(vt.localTransform.rotation().col(eY)[eZ]) >= tolerance) {
        ACTS_ERROR("Volumes are not aligned: rotation is different");
        throw std::invalid_argument(
            "Volumes are not aligned: rotation is different");
      }

      ACTS_VERBOSE(" -> Rotation is ok!");

      // b) the volumes cannot have translation in x or y
      Vector2 translation = vt.localTransform.translation().head<2>();
      if (std::abs(translation[0]) > tolerance ||  //
          std::abs(translation[1]) > tolerance) {
        ACTS_ERROR("Volumes are not aligned: translation in x or y");
        throw std::invalid_argument(
            "Volumes are not aligned: translation in x or y");
      }
      ACTS_VERBOSE(" -> Translation in x/y is ok!");

      n++;
    }

    ACTS_VERBOSE("Sorting by volume z position");
    std::sort(volumeTuples.begin(), volumeTuples.end(),
              [](const auto& a, const auto& b) {
                return a.localTransform.translation()[eZ] <
                       b.localTransform.translation()[eZ];
              });

    ACTS_VERBOSE("Checking for overlaps and attaching volumes in z");
    std::vector<VolumeTuple> gapVolumes;
    for (std::size_t i = 0; i < volumeTuples.size() - 1; i++) {
      std::size_t j = i + 1;
      auto& a = volumeTuples.at(i);
      auto& b = volumeTuples.at(j);

      auto gap =
          checkOverlapAndAttachInZ(a, b, groupTransform, strategy, logger);
      if (gap) {
        ACTS_VERBOSE("Received gap volume, registering");
        auto* cylinderBounds =
            dynamic_cast<CylinderVolumeBounds*>(&gap->volumeBounds());
        if (cylinderBounds == nullptr) {
          throw std::invalid_argument{
              "CylinderVolumeStack constructed gap volume was not a cylinder "
              "volume"};
        }
        gapVolumes.emplace_back(
            gap, cylinderBounds,
            std::make_shared<CylinderVolumeBounds>(*cylinderBounds),
            groupTransform.inverse() * gap->transform());
      }
    }

    ACTS_VERBOSE("Appending "
                 << gapVolumes.size()
                 << " gap volumes to the end of the volume vector");
    std::copy(gapVolumes.begin(), gapVolumes.end(),
              std::back_inserter(volumeTuples));

    ACTS_VERBOSE("*** Volume configuration after z attachment:");
    printVolumeSequence(volumeTuples, logger, Acts::Logging::VERBOSE);

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
    ACTS_VERBOSE("Found: minR: " << minR << " maxR: " << maxR);

    for (auto& vt : volumeTuples) {
      vt.updatedBounds->set({
          {CylinderVolumeBounds::eMinR, minR},
          {CylinderVolumeBounds::eMaxR, maxR},
      });

      ACTS_VERBOSE("Updated bounds for volume at z: "
                   << vt.localTransform.translation()[eZ]);
      ACTS_VERBOSE(*vt.updatedBounds);

      // make a copy so we can't accidentally modify in-place
      auto copy = std::make_shared<CylinderVolumeBounds>(*vt.updatedBounds);
      vt.volume->assignVolumeBounds(std::move(vt.updatedBounds));
      vt.bounds = copy.get();
      vt.updatedBounds = std::move(copy);
    }

    ACTS_VERBOSE("*** Volume configuration after r synchronization:");
    printVolumeSequence(volumeTuples, logger, Acts::Logging::VERBOSE);

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
    printVolumeSequence(volumeTuples, logger, Acts::Logging::DEBUG);

    ActsScalar minZ =
        volumeTuples.front().localTransform.translation()[eZ] -
        volumeTuples.front().bounds->get(CylinderVolumeBounds::eHalfLengthZ);
    ActsScalar maxZ =
        volumeTuples.back().localTransform.translation()[eZ] +
        volumeTuples.back().bounds->get(CylinderVolumeBounds::eHalfLengthZ);

    ActsScalar midZ = (minZ + maxZ) / 2.0;
    ActsScalar hlZ = (maxZ - minZ) / 2.0;

    // ACTS_DEBUG("Outer volume bounds are z: [ "
    // << minZ << " <- " << midZ << " -> " << maxZ << " ], r: [ "
    // << minR << " <-> " << maxR << " ]");
    // ACTS_DEBUG(" -> Half length: " << hlZ);

    Transform3 outerTransform = groupTransform * Translation3{0, 0, midZ};

    auto outerCylinderBounds =
        std::make_shared<CylinderVolumeBounds>(minR, maxR, hlZ);
    ACTS_DEBUG("Outer bounds are:\n" << *outerCylinderBounds);
    ACTS_DEBUG("Outer transform is:\n" << outerTransform.matrix());

    return Volume(outerTransform, outerCylinderBounds);
  } else {
    throw std::invalid_argument(binningValueNames()[direction] +
                                " is not supported ");
  }
}

std::shared_ptr<Volume> CylinderVolumeStack::checkOverlapAndAttachInZ(
    VolumeTuple& a, VolumeTuple& b, const Transform3& groupTransform,
    CylinderVolumeStack::AttachmentStrategy strategy, const Logger& logger) {
  ActsScalar aZ = a.localTransform.translation()[eZ];
  ActsScalar aHlZ = a.updatedBounds->get(CylinderVolumeBounds::eHalfLengthZ);
  ActsScalar aZMin = aZ - aHlZ;
  ActsScalar aZMax = aZ + aHlZ;

  ActsScalar bZ = b.localTransform.translation()[eZ];
  ActsScalar bHlZ = b.updatedBounds->get(CylinderVolumeBounds::eHalfLengthZ);
  ActsScalar bZMin = bZ - bHlZ;
  ActsScalar bZMax = bZ + bHlZ;

  if (logger().doPrint(Acts::Logging::DEBUG)) {
    std::stringstream ss;
    ss << std::fixed;
    ss << std::setprecision(3);
    ss << std::setfill(' ');
    ACTS_VERBOSE("Checking overlap between");
    int w = 9;
    ss << " - "
       << " z: [ " << std::setw(w) << aZMin << " <- " << std::setw(w) << aZ
       << " -> " << std::setw(w) << aZMax << " ]";
    ACTS_VERBOSE(ss.str());

    ss.str("");
    ss << " - "
       << " z: [ " << std::setw(w) << bZMin << " <- " << std::setw(w) << bZ
       << " -> " << std::setw(w) << bZMax << " ]";
    ACTS_VERBOSE(ss.str());
  }

  if (aZMax > bZMin) {
    ACTS_ERROR(" -> Overlap in z");
    throw std::invalid_argument("Volumes overlap in z");
  } else {
    ACTS_VERBOSE(" -> No overlap");
  }

  std::shared_ptr<Volume> gap;

  constexpr auto tolerance = s_onSurfaceTolerance;
  if (std::abs(aZMax - bZMin) < tolerance) {
    ACTS_VERBOSE("No gap between volumes, no attachment needed");
  } else {
    ACTS_VERBOSE("Synchronizing bounds in z with strategy: " << strategy);

    switch (strategy) {
      case AttachmentStrategy::Midpoint: {
        ACTS_VERBOSE(" -> Strategy: Expand both volumes to midpoint");
        ActsScalar gapWidth = bZMin - aZMax;
        ACTS_VERBOSE("  - Gap width: " << gapWidth);

        ActsScalar aZMidNew = (aZMin + aZMax) / 2.0 + gapWidth / 4.0;
        ActsScalar aHlZNew = (aZMax - aZMin) / 2.0 + gapWidth / 4.0;
        ACTS_VERBOSE("  - New halflength for first volume: " << aHlZNew);
        ACTS_VERBOSE("  - New bounds for first volume: ["
                     << (aZMidNew - aHlZNew) << " <- " << aZMidNew << " -> "
                     << (aZMidNew + aHlZNew) << "]");

        assert(aZMin == aZMidNew - aHlZNew && "Volume shrunk");
        assert(aZMax <= aZMidNew + aHlZNew && "Volume shrunk");

        ActsScalar bZMidNew = (bZMin + bZMax) / 2.0 - gapWidth / 4.0;
        ActsScalar bHlZNew = (bZMax - bZMin) / 2.0 + gapWidth / 4.0;
        ACTS_VERBOSE("  - New halflength for second volume: " << bHlZNew);
        ACTS_VERBOSE("  - New bounds for second volume: ["
                     << (bZMidNew - bHlZNew) << " <- " << bZMidNew << " -> "
                     << (bZMidNew + bHlZNew) << "]");

        assert(bZMin >= bZMidNew - bHlZNew && "Volume shrunk");
        assert(bZMax == bZMidNew + bHlZNew && "Volume shrunk");

        a.localTransform = Translation3{0, 0, aZMidNew};
        a.volume->setTransform(groupTransform * a.localTransform);
        a.updatedBounds->set(CylinderVolumeBounds::eHalfLengthZ, aHlZNew);

        b.localTransform = Translation3{0, 0, bZMidNew};
        b.volume->setTransform(groupTransform * b.localTransform);
        b.updatedBounds->set(CylinderVolumeBounds::eHalfLengthZ, bHlZNew);

        break;
      }
      case AttachmentStrategy::First: {
        ACTS_VERBOSE(" -> Strategy: Expand first volume");
        ActsScalar gapWidth = bZMin - aZMax;
        ActsScalar aZMidNew = (aZMin + bZMin) / 2.0;
        ActsScalar aHlZNew = (bZMin - aZMin) / 2.0;
        ACTS_VERBOSE("  - Gap width: " << gapWidth);
        ACTS_VERBOSE("  - New halflength for first volume: " << aHlZNew);
        ACTS_VERBOSE("  - New bounds for first volume: ["
                     << (aZMidNew - aHlZNew) << " <- " << aZMidNew << " -> "
                     << (aZMidNew + aHlZNew) << "]");

        assert(aZMin == aZMidNew - aHlZNew && "Volume shrunk");
        assert(aZMax <= aZMidNew + aHlZNew && "Volume shrunk");

        a.localTransform = Translation3{0, 0, aZMidNew};
        a.volume->setTransform(groupTransform * a.localTransform);
        a.updatedBounds->set(CylinderVolumeBounds::eHalfLengthZ, aHlZNew);

        break;
      }
      case AttachmentStrategy::Second: {
        ACTS_VERBOSE(" -> Strategy: Expand second volume");
        ActsScalar gapWidth = bZMin - aZMax;
        ActsScalar bZMidNew = (aZMax + bZMax) / 2.0;
        ActsScalar bHlZNew = (bZMax - aZMax) / 2.0;
        ACTS_VERBOSE("  - Gap width: " << gapWidth);
        ACTS_VERBOSE("  - New halflength for second volume: " << bHlZNew);
        ACTS_VERBOSE("  - New bounds for second volume: ["
                     << (bZMidNew - bHlZNew) << " <- " << bZMidNew << " -> "
                     << (bZMidNew + bHlZNew) << "]");

        assert(bZMin >= bZMidNew - bHlZNew && "Volume shrunk");
        assert(bZMax == bZMidNew + bHlZNew && "Volume shrunk");

        b.localTransform = Translation3{0, 0, bZMidNew};
        b.volume->setTransform(groupTransform * b.localTransform);
        b.updatedBounds->set(CylinderVolumeBounds::eHalfLengthZ, bHlZNew);
        break;
      }
      case AttachmentStrategy::Gap: {
        ACTS_VERBOSE(" -> Strategy: Create a gap volume");
        ActsScalar gapHlZ = (bZMin - aZMax) / 2.0;
        ActsScalar gapMidZ = (bZMin + aZMax) / 2.0;

        ACTS_VERBOSE("  - Gap half length: " << gapHlZ << " at z: " << gapMidZ);

        Transform3 gapLocalTransform{Translation3{0, 0, gapMidZ}};
        gap = std::make_shared<Volume>(
            groupTransform * gapLocalTransform,
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

void CylinderVolumeStack::printVolumeSequence(
    const std::vector<VolumeTuple>& volumes, const Logger& logger,
    Acts::Logging::Level lvl) {
  if (!logger().doPrint(lvl)) {
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

    logger().log(lvl, ss.str());
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
