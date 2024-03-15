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
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <sstream>

namespace Acts {

struct CylinderVolumeStack::VolumeTuple {
  std::shared_ptr<Volume> volume;
  CylinderVolumeBounds* bounds{};
  std::shared_ptr<CylinderVolumeBounds> updatedBounds{};
  Transform3 localTransform;

  VolumeTuple(std::shared_ptr<Volume> volume_, const Transform3& groupTransform)
      : volume{std::move(volume_)},
        localTransform{groupTransform.inverse() * volume->transform()} {
    bounds = dynamic_cast<CylinderVolumeBounds*>(&volume->volumeBounds());
    assert(bounds != nullptr);
    updatedBounds = std::make_shared<CylinderVolumeBounds>(*bounds);
  }

  ActsScalar midZ() const { return localTransform.translation()[eZ]; }
  ActsScalar halfLengthZ() const {
    return updatedBounds->get(CylinderVolumeBounds::eHalfLengthZ);
  }
  ActsScalar minZ() const { return midZ() - halfLengthZ(); }
  ActsScalar maxZ() const { return midZ() + halfLengthZ(); }

  ActsScalar minR() const {
    return updatedBounds->get(CylinderVolumeBounds::eMinR);
  }
  ActsScalar maxR() const {
    return updatedBounds->get(CylinderVolumeBounds::eMaxR);
  }
  ActsScalar midR() const { return (minR() + maxR()) / 2.0; }

  void set(std::initializer_list<
           std::pair<CylinderVolumeBounds::BoundValues, ActsScalar>>
               keyValues) {
    updatedBounds->set(keyValues);
  }

  void setLocalTransform(const Transform3& transform,
                         const Transform3& groupTransform) {
    localTransform = transform;
    volume->setTransform(groupTransform * localTransform);
  }

  void commit() {
    // make a copy so we can't accidentally modify in-place
    auto copy = std::make_shared<CylinderVolumeBounds>(*updatedBounds);
    volume->assignVolumeBounds(std::move(updatedBounds));
    bounds = copy.get();
    updatedBounds = std::move(copy);
  }
};

CylinderVolumeStack::CylinderVolumeStack(
    std::vector<std::shared_ptr<Volume>>& volumes, BinningValue direction,
    AttachmentStrategy strategy, ResizeStrategy resizeStrategy,
    const Logger& logger)
    : Volume(createOuterVolume(volumes, direction, strategy, logger)),
      m_direction(direction),
      m_resizeStrategy(resizeStrategy),
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
  m_groupTransform = volumes.front()->transform();
  ACTS_VERBOSE("Initial group transform is:\n" << m_groupTransform.matrix());

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
    volumeTuples.emplace_back(volume, m_groupTransform);
  }

  ACTS_DEBUG("*** Initial volume configuration:");
  printVolumeSequence(volumeTuples, logger, Acts::Logging::DEBUG);

  if (volumes.size() == 1) {
    ACTS_VERBOSE("Only one volume, returning")
    return Volume(*volumes.front());
  }

  ACTS_VERBOSE("Checking volume alignment");
  checkVolumeAlignment(volumeTuples, logger);

  if (direction == Acts::binZ) {
    ACTS_VERBOSE("Sorting by volume z position");
    std::sort(volumeTuples.begin(), volumeTuples.end(),
              [](const auto& a, const auto& b) {
                return a.localTransform.translation()[eZ] <
                       b.localTransform.translation()[eZ];
              });

    ACTS_VERBOSE("Checking for overlaps and attaching volumes in z");
    std::vector<VolumeTuple> gapVolumes =
        checkOverlapAndAttachInZ(volumeTuples, strategy, logger);

    ACTS_VERBOSE("Appending "
                 << gapVolumes.size()
                 << " gap volumes to the end of the volume vector");
    std::copy(gapVolumes.begin(), gapVolumes.end(),
              std::back_inserter(volumeTuples));

    ACTS_VERBOSE("*** Volume configuration after z attachment:");
    printVolumeSequence(volumeTuples, logger, Acts::Logging::VERBOSE);

    ACTS_VERBOSE("Synchronizing bounds in r")
    const auto [minR, maxR] = synchronizeRBounds(volumeTuples, logger);

    for (auto& vt : volumeTuples) {
      ACTS_VERBOSE("Updated bounds for volume at z: "
                   << vt.localTransform.translation()[eZ]);
      ACTS_VERBOSE(*vt.updatedBounds);

      vt.commit();
    }

    ACTS_VERBOSE("*** Volume configuration after r synchronization:");
    printVolumeSequence(volumeTuples, logger, Acts::Logging::VERBOSE);

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

    ActsScalar minZ = volumeTuples.front().minZ();
    ActsScalar maxZ = volumeTuples.back().maxZ();

    ActsScalar midZ = (minZ + maxZ) / 2.0;
    ActsScalar hlZ = (maxZ - minZ) / 2.0;

    m_groupTransform = m_groupTransform * Translation3{0, 0, midZ};

    auto outerCylinderBounds =
        std::make_shared<CylinderVolumeBounds>(minR, maxR, hlZ);
    ACTS_DEBUG("Outer bounds are:\n" << *outerCylinderBounds);
    ACTS_DEBUG("Outer transform / new group transform is:\n"
               << m_groupTransform.matrix());

    return Volume(m_groupTransform, outerCylinderBounds);

  } else if (direction == Acts::binR) {
    ACTS_ERROR("Sorting by volume r middle point");
    std::sort(volumeTuples.begin(), volumeTuples.end(),
              [](const auto& a, const auto& b) { return a.midR() < b.midR(); });

    ACTS_VERBOSE("Checking for overlaps and attaching volumes in r");
    std::vector<VolumeTuple> gapVolumes =
        checkOverlapAndAttachInR(volumeTuples, strategy, logger);

    ACTS_VERBOSE("Appending "
                 << gapVolumes.size()
                 << " gap volumes to the end of the volume vector");
    std::copy(gapVolumes.begin(), gapVolumes.end(),
              std::back_inserter(volumeTuples));

    ACTS_VERBOSE("*** Volume configuration after r attachment:");
    printVolumeSequence(volumeTuples, logger, Acts::Logging::VERBOSE);

    ACTS_VERBOSE("Synchronizing bounds in z")
    const auto [minZ, maxZ] = synchronizeZBounds(volumeTuples, logger);

    for (auto& vt : volumeTuples) {
      ACTS_VERBOSE("Updated bounds for volume at r: " << vt.midR());
      ACTS_VERBOSE(*vt.updatedBounds);
      vt.commit();
    }

    ACTS_VERBOSE("*** Volume configuration after z synchronization:");
    printVolumeSequence(volumeTuples, logger, Acts::Logging::VERBOSE);

    std::sort(volumeTuples.begin(), volumeTuples.end(),
              [](const auto& a, const auto& b) { return a.midR() < b.midR(); });

    // Update outer volume vector after sorting
    volumes.clear();
    for (const auto& vt : volumeTuples) {
      volumes.push_back(vt.volume);
    }

    ACTS_DEBUG("*** Volume configuration after final r sorting:");
    printVolumeSequence(volumeTuples, logger, Acts::Logging::DEBUG);

    ActsScalar minR = volumeTuples.front().minR();
    ActsScalar maxR = volumeTuples.back().maxR();

    ActsScalar midZ = (minZ + maxZ) / 2.0;
    ActsScalar hlZ = (maxZ - minZ) / 2.0;

    Transform3 outerTransform = m_groupTransform * Translation3{0, 0, midZ};

    auto outerCylinderBounds =
        std::make_shared<CylinderVolumeBounds>(minR, maxR, hlZ);

    ACTS_DEBUG("Outer bounds are:\n" << *outerCylinderBounds);
    ACTS_DEBUG("Outer transform is:\n" << outerTransform.matrix());

    return Volume(outerTransform, outerCylinderBounds);

  } else {
    ACTS_ERROR("Binning in " << binningValueNames()[direction]
                             << " is not supported");
    throw std::invalid_argument(binningValueNames()[direction] +
                                " is not supported ");
  }
}

void CylinderVolumeStack::overlapPrint(
    const CylinderVolumeStack::VolumeTuple& a,
    const CylinderVolumeStack::VolumeTuple& b, const Logger& logger) {
  if (logger().doPrint(Acts::Logging::DEBUG)) {
    std::stringstream ss;
    ss << std::fixed;
    ss << std::setprecision(3);
    ss << std::setfill(' ');
    ACTS_VERBOSE("Checking overlap between");
    int w = 9;
    ss << " - "
       << " z: [ " << std::setw(w) << a.minZ() << " <- " << std::setw(w)
       << a.midZ() << " -> " << std::setw(w) << a.maxZ() << " ]";
    ACTS_VERBOSE(ss.str());

    ss.str("");
    ss << " - "
       << " z: [ " << std::setw(w) << b.minZ() << " <- " << std::setw(w)
       << b.midZ() << " -> " << std::setw(w) << b.maxZ() << " ]";
    ACTS_VERBOSE(ss.str());
  }
}

std::vector<CylinderVolumeStack::VolumeTuple>
CylinderVolumeStack::checkOverlapAndAttachInZ(
    std::vector<VolumeTuple>& volumes,
    CylinderVolumeStack::AttachmentStrategy strategy, const Logger& logger) {
  // Preconditions: volumes are sorted by z
  std::vector<VolumeTuple> gapVolumes;
  for (std::size_t i = 0; i < volumes.size() - 1; i++) {
    std::size_t j = i + 1;
    auto& a = volumes.at(i);
    auto& b = volumes.at(j);

    overlapPrint(a, b, logger);

    if (a.maxZ() > b.minZ()) {
      ACTS_ERROR(" -> Overlap in z");
      throw std::invalid_argument("Volumes overlap in z");
    } else {
      ACTS_VERBOSE(" -> No overlap");
    }

    constexpr auto tolerance = s_onSurfaceTolerance;
    if (std::abs(a.maxZ() - b.minZ()) < tolerance) {
      ACTS_VERBOSE("No gap between volumes, no attachment needed");
    } else {
      ActsScalar gapWidth = b.minZ() - a.maxZ();
      ACTS_VERBOSE("Gap width: " << gapWidth);

      ACTS_VERBOSE("Synchronizing bounds in z with strategy: " << strategy);
      switch (strategy) {
        case AttachmentStrategy::Midpoint: {
          ACTS_VERBOSE(" -> Strategy: Expand both volumes to midpoint");

          ActsScalar aZMidNew = (a.minZ() + a.maxZ()) / 2.0 + gapWidth / 4.0;
          ActsScalar aHlZNew = a.halfLengthZ() + gapWidth / 4.0;
          ACTS_VERBOSE("  - New halflength for first volume: " << aHlZNew);
          ACTS_VERBOSE("  - New bounds for first volume: ["
                       << (aZMidNew - aHlZNew) << " <- " << aZMidNew << " -> "
                       << (aZMidNew + aHlZNew) << "]");

          assert(a.minZ() == aZMidNew - aHlZNew && "Volume shrunk");
          assert(a.maxZ() <= aZMidNew + aHlZNew && "Volume shrunk");

          ActsScalar bZMidNew = (b.minZ() + b.maxZ()) / 2.0 - gapWidth / 4.0;
          ActsScalar bHlZNew = b.halfLengthZ() + gapWidth / 4.0;
          ACTS_VERBOSE("  - New halflength for second volume: " << bHlZNew);
          ACTS_VERBOSE("  - New bounds for second volume: ["
                       << (bZMidNew - bHlZNew) << " <- " << bZMidNew << " -> "
                       << (bZMidNew + bHlZNew) << "]");

          assert(b.minZ() >= bZMidNew - bHlZNew && "Volume shrunk");
          assert(b.maxZ() == bZMidNew + bHlZNew && "Volume shrunk");

          a.localTransform = Translation3{0, 0, aZMidNew};
          a.volume->setTransform(m_groupTransform * a.localTransform);
          a.updatedBounds->set(CylinderVolumeBounds::eHalfLengthZ, aHlZNew);

          b.localTransform = Translation3{0, 0, bZMidNew};
          b.volume->setTransform(m_groupTransform * b.localTransform);
          b.updatedBounds->set(CylinderVolumeBounds::eHalfLengthZ, bHlZNew);

          break;
        }
        case AttachmentStrategy::First: {
          ACTS_VERBOSE(" -> Strategy: Expand first volume");
          ActsScalar aZMidNew = (a.minZ() + b.minZ()) / 2.0;
          ActsScalar aHlZNew = (b.minZ() - a.minZ()) / 2.0;
          ACTS_VERBOSE("  - Gap width: " << gapWidth);
          ACTS_VERBOSE("  - New bounds for first volume: ["
                       << (aZMidNew - aHlZNew) << " <- " << aZMidNew << " -> "
                       << (aZMidNew + aHlZNew) << "]");

          assert(a.minZ() == aZMidNew - aHlZNew && "Volume shrunk");
          assert(a.maxZ() <= aZMidNew + aHlZNew && "Volume shrunk");

          a.localTransform = Translation3{0, 0, aZMidNew};
          a.volume->setTransform(m_groupTransform * a.localTransform);
          a.updatedBounds->set(CylinderVolumeBounds::eHalfLengthZ, aHlZNew);

          break;
        }
        case AttachmentStrategy::Second: {
          ACTS_VERBOSE(" -> Strategy: Expand second volume");
          ActsScalar bZMidNew = (a.maxZ() + b.maxZ()) / 2.0;
          ActsScalar bHlZNew = (b.maxZ() - a.maxZ()) / 2.0;
          ACTS_VERBOSE("  - New halflength for second volume: " << bHlZNew);
          ACTS_VERBOSE("  - New bounds for second volume: ["
                       << (bZMidNew - bHlZNew) << " <- " << bZMidNew << " -> "
                       << (bZMidNew + bHlZNew) << "]");

          assert(b.minZ() >= bZMidNew - bHlZNew && "Volume shrunk");
          assert(b.maxZ() == bZMidNew + bHlZNew && "Volume shrunk");

          b.localTransform = Translation3{0, 0, bZMidNew};
          b.volume->setTransform(m_groupTransform * b.localTransform);
          b.updatedBounds->set(CylinderVolumeBounds::eHalfLengthZ, bHlZNew);
          break;
        }
        case AttachmentStrategy::Gap: {
          ACTS_VERBOSE(" -> Strategy: Create a gap volume");
          ActsScalar gapHlZ = (b.minZ() - a.maxZ()) / 2.0;
          ActsScalar gapMidZ = (b.minZ() + a.maxZ()) / 2.0;

          ACTS_VERBOSE("  - Gap half length: " << gapHlZ
                                               << " at z: " << gapMidZ);

          Transform3 gapLocalTransform{Translation3{0, 0, gapMidZ}};
          Transform3 gapGlobalTransform = m_groupTransform * gapLocalTransform;
          auto gapBounds = std::make_shared<CylinderVolumeBounds>(
              a.minR(), b.maxR(), gapHlZ);
          gapVolumes.emplace_back(std::make_shared<Volume>(
                                      gapGlobalTransform, std::move(gapBounds)),
                                  m_groupTransform);

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
  }

  return gapVolumes;
}

std::vector<CylinderVolumeStack::VolumeTuple>
CylinderVolumeStack::checkOverlapAndAttachInR(
    std::vector<VolumeTuple>& volumes,
    CylinderVolumeStack::AttachmentStrategy strategy, const Logger& logger) {
  std::vector<VolumeTuple> gapVolumes;
  for (std::size_t i = 0; i < volumes.size() - 1; i++) {
    std::size_t j = i + 1;
    auto& a = volumes.at(i);
    auto& b = volumes.at(j);

    overlapPrint(a, b, logger);

    if (a.maxR() > b.minR()) {
      ACTS_ERROR(" -> Overlap in r");
      throw std::invalid_argument("Volumes overlap in r");
    } else {
      ACTS_VERBOSE(" -> No overlap");
    }

    constexpr auto tolerance = s_onSurfaceTolerance;
    if (std::abs(a.maxR() - b.minR()) < tolerance) {
      ACTS_VERBOSE("No gap between volumes, no attachment needed");
    } else {
      ActsScalar gapWidth = b.minR() - a.maxR();
      ACTS_VERBOSE("Gap width: " << gapWidth);

      ACTS_VERBOSE("Synchronizing bounds in r with strategy: " << strategy);
      switch (strategy) {
        case AttachmentStrategy::Midpoint: {
          ACTS_VERBOSE(" -> Strategy: Expand both volumes to midpoint");

          a.set({{CylinderVolumeBounds::eMaxR, a.maxR() + gapWidth / 2.0}});
          b.set({{CylinderVolumeBounds::eMinR, b.minR() - gapWidth / 2.0}});

          break;
        }
        case AttachmentStrategy::First: {
          ACTS_VERBOSE(" -> Strategy: Expand first volume");

          a.set({{CylinderVolumeBounds::eMaxR, b.minR()}});

          break;
        }
        case AttachmentStrategy::Second: {
          ACTS_VERBOSE(" -> Strategy: Expand second volume");

          b.set({{CylinderVolumeBounds::eMinR, a.maxR()}});

          break;
        }
        case AttachmentStrategy::Gap: {
          ACTS_VERBOSE(" -> Strategy: Create a gap volume");

          auto gapBounds = std::make_shared<CylinderVolumeBounds>(
              a.maxR(), b.minR(), a.halfLengthZ());
          auto gap = std::make_shared<Volume>(m_groupTransform, gapBounds);

          gapVolumes.emplace_back(std::move(gap), m_groupTransform);
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
  }

  return gapVolumes;
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

void CylinderVolumeStack::checkVolumeAlignment(
    const std::vector<VolumeTuple>& volumes, const Logger& logger) {
  std::size_t n = 0;
  for (auto& vt : volumes) {
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
}

std::pair<ActsScalar, ActsScalar> CylinderVolumeStack::synchronizeRBounds(
    std::vector<VolumeTuple>& volumes, const Logger& logger) {
  const ActsScalar minR =
      std::min_element(volumes.begin(), volumes.end(),
                       [](const auto& a, const auto& b) {
                         return a.bounds->get(CylinderVolumeBounds::eMinR) <
                                b.bounds->get(CylinderVolumeBounds::eMinR);
                       })
          ->bounds->get(CylinderVolumeBounds::eMinR);

  const ActsScalar maxR =
      std::max_element(volumes.begin(), volumes.end(),
                       [](const auto& a, const auto& b) {
                         return a.bounds->get(CylinderVolumeBounds::eMaxR) <
                                b.bounds->get(CylinderVolumeBounds::eMaxR);
                       })
          ->bounds->get(CylinderVolumeBounds::eMaxR);
  ACTS_VERBOSE("Found: minR: " << minR << " maxR: " << maxR);

  for (auto& vt : volumes) {
    vt.set({
        {CylinderVolumeBounds::eMinR, minR},
        {CylinderVolumeBounds::eMaxR, maxR},
    });
  }

  return {minR, maxR};
}

std::pair<ActsScalar, ActsScalar> CylinderVolumeStack::synchronizeZBounds(
    std::vector<VolumeTuple>& volumes, const Logger& logger) {
  const ActsScalar minZ = std::min_element(volumes.begin(), volumes.end(),
                                           [](const auto& a, const auto& b) {
                                             return a.minZ() < b.minZ();
                                           })
                              ->minZ();

  const ActsScalar maxZ = std::max_element(volumes.begin(), volumes.end(),
                                           [](const auto& a, const auto& b) {
                                             return a.maxZ() < b.maxZ();
                                           })
                              ->maxZ();
  const ActsScalar midZ = (minZ + maxZ) / 2.0;
  const ActsScalar hlZ = (maxZ - minZ) / 2.0;
  ACTS_DEBUG("Found overall z bounds: [ " << minZ << " <- " << midZ << " -> "
                                          << maxZ << " ]");
  const Transform3 transform{Translation3{0, 0, midZ}};

  for (auto& vt : volumes) {
    vt.set({{CylinderVolumeBounds::eHalfLengthZ, hlZ}});
    vt.setLocalTransform(transform, m_groupTransform);
  }

  return {minZ, maxZ};
}

void CylinderVolumeStack::assignVolumeBounds(
    std::shared_ptr<VolumeBounds> volbounds) {
  assignVolumeBounds(std::move(volbounds), Acts::getDummyLogger());
}

void CylinderVolumeStack::assignVolumeBounds(
    std::shared_ptr<VolumeBounds> volbounds, const Logger& logger) {
  ACTS_DEBUG(
      "Resizing CylinderVolumeStack with strategy: " << m_resizeStrategy);

  const auto* newBounds = dynamic_cast<CylinderVolumeBounds*>(volbounds.get());
  if (newBounds == nullptr) {
    ACTS_ERROR(
        "Tried to assign non-CylinderVolumeBounds to "
        "CylinderVolumeStack");
    throw std::invalid_argument(
        "CylinderVolumeStack requires CylinderVolumeBounds");
  }

  const auto* oldBounds =
      dynamic_cast<CylinderVolumeBounds*>(m_volumeBounds.get());
  const ActsScalar newMinR = newBounds->get(CylinderVolumeBounds::eMinR);
  const ActsScalar newMaxR = newBounds->get(CylinderVolumeBounds::eMaxR);
  const ActsScalar newHlZ = newBounds->get(CylinderVolumeBounds::eHalfLengthZ);

  const ActsScalar oldMinR = oldBounds->get(CylinderVolumeBounds::eMinR);
  const ActsScalar oldMaxR = oldBounds->get(CylinderVolumeBounds::eMaxR);
  const ActsScalar oldHlZ = oldBounds->get(CylinderVolumeBounds::eHalfLengthZ);

  ACTS_VERBOSE("Previous bounds are: z: [ " << -oldHlZ << " <- 0 -> " << oldHlZ
                                            << " ], r: [ " << oldMinR << " <-> "
                                            << oldMaxR << " ]");
  ACTS_VERBOSE("New bounds are: z:      [ " << -newHlZ << " <- 0 -> " << newHlZ
                                            << " ], r: [ " << newMinR << " <-> "
                                            << newMaxR << " ]");

  if (newHlZ < oldHlZ) {
    ACTS_ERROR("Shrinking the stack size in z is not supported");
    throw std::invalid_argument("Shrinking the stack is not supported");
  }

  if (newMinR > oldMinR) {
    ACTS_ERROR("Shrinking the stack size in r is not supported");
    throw std::invalid_argument("Shrinking the stack is not supported");
  }

  if (newMaxR < oldMaxR) {
    ACTS_ERROR("Shrinking the stack size in r is not supported");
    throw std::invalid_argument("Shrinking the stack is not supported");
  }

  if (*newBounds == *oldBounds) {
    ACTS_VERBOSE("Bounds are the same, no resize needed");
    return;
  }

  if (m_direction == BinningValue::binZ) {
    ACTS_VERBOSE("Stack direction is z");

    std::vector<VolumeTuple> volumes;
    volumes.reserve(m_volumes.size());
    std::transform(m_volumes.begin(), m_volumes.end(),
                   std::back_inserter(volumes), [this](const auto& volume) {
                     return VolumeTuple{volume, m_groupTransform};
                   });

    ACTS_VERBOSE("*** Initial volume configuration:");
    printVolumeSequence(volumes, logger, Acts::Logging::DEBUG);

    ACTS_VERBOSE("Resize all volumes to new r bounds");
    for (auto& volume : volumes) {
      volume.set({
          {CylinderVolumeBounds::eMinR, newMinR},
          {CylinderVolumeBounds::eMaxR, newMaxR},
      });
    }

    ACTS_VERBOSE("*** Volume configuration after r resizing:");
    printVolumeSequence(volumes, logger, Acts::Logging::DEBUG);

    if (newHlZ == oldHlZ) {
      ACTS_VERBOSE("Halflength z is the same, no z resize needed");
    } else {
      if (m_resizeStrategy == ResizeStrategy::Expand) {
        ACTS_VERBOSE("Expanding first and last volume to new z bounds");
        auto& first = volumes.front();
        ActsScalar newMinZFirst = -newHlZ;
        ActsScalar newMidZFirst = (newMinZFirst + first.maxZ()) / 2.0;
        ActsScalar newHlZFirst = (first.maxZ() - newMinZFirst) / 2.0;

        ACTS_VERBOSE(" -> first z: [ " << newMinZFirst << " <- " << newMidZFirst
                                       << " -> " << first.maxZ() << " ]");

        first.set({{CylinderVolumeBounds::eHalfLengthZ, newHlZFirst}});
        first.setLocalTransform(Transform3{Translation3{0, 0, newMidZFirst}},
                                m_groupTransform);

        auto& last = volumes.back();
        ActsScalar newMaxZLast = newHlZ;
        ActsScalar newMidZLast = (last.minZ() + newMaxZLast) / 2.0;
        ActsScalar newHlZLast = (newMaxZLast - last.minZ()) / 2.0;

        ACTS_VERBOSE(" -> last z: [ " << last.minZ() << " <- " << newMidZLast
                                      << " -> " << newMaxZLast << " ]");

        last.set({{CylinderVolumeBounds::eHalfLengthZ, newHlZLast}});
        last.setLocalTransform(Transform3{Translation3{0, 0, newMidZLast}},
                               m_groupTransform);
      } else if (m_resizeStrategy == ResizeStrategy::Gap) {
        ACTS_VERBOSE("Creating gap volumes to fill the new z bounds");

        ActsScalar gap1MinZ = -newHlZ;
        ActsScalar gap1MaxZ = volumes.front().minZ();
        ActsScalar gap1HlZ = (gap1MaxZ - gap1MinZ) / 2.0;
        ActsScalar gap1PZ = (gap1MaxZ + gap1MinZ) / 2.0;

        ACTS_VERBOSE(" -> gap1 z: [ " << gap1MinZ << " <- " << gap1PZ << " -> "
                                      << gap1MaxZ << " ]");

        auto gap1Bounds =
            std::make_shared<CylinderVolumeBounds>(newMinR, newMaxR, gap1HlZ);
        auto gap1Transform = m_groupTransform * Translation3{0, 0, gap1PZ};
        volumes.insert(volumes.begin(),
                       VolumeTuple{std::make_shared<Volume>(
                                       gap1Transform, std::move(gap1Bounds)),
                                   m_groupTransform});

        ActsScalar gap2MinZ = volumes.back().maxZ();
        ActsScalar gap2MaxZ = newHlZ;
        ActsScalar gap2HlZ = (gap2MaxZ - gap2MinZ) / 2.0;
        ActsScalar gap2PZ = (gap2MaxZ + gap2MinZ) / 2.0;

        ACTS_VERBOSE(" -> gap2 z: [ " << gap2MinZ << " <- " << gap2PZ << " -> "
                                      << gap2MaxZ << " ]");

        auto gap2Bounds =
            std::make_shared<CylinderVolumeBounds>(newMinR, newMaxR, gap2HlZ);
        auto gap2Transform = m_groupTransform * Translation3{0, 0, gap2PZ};
        volumes.emplace_back(
            std::make_shared<Volume>(gap2Transform, std::move(gap2Bounds)),
            m_groupTransform);
      }
    }

    ACTS_VERBOSE("*** Volume configuration after z resizing:");
    printVolumeSequence(volumes, logger, Acts::Logging::DEBUG);

    // Commit and update outer vector
    m_volumes.clear();
    for (auto& vt : volumes) {
      vt.commit();
      m_volumes.push_back(vt.volume);
    }
  }

  m_volumeBounds = std::move(volbounds);
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

std::ostream& operator<<(std::ostream& os,
                         CylinderVolumeStack::ResizeStrategy strategy) {
  switch (strategy) {
    case CylinderVolumeStack::ResizeStrategy::Expand:
      os << "Expand";
      break;
    case CylinderVolumeStack::ResizeStrategy::Gap:
      os << "Gap";
      break;
  }
  return os;
}

}  // namespace Acts
