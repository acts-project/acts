// This file is part of the Acts project.
//
// Copyright (C) 2022-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Utilities/Delegate.hpp"

#include <memory>
#include <vector>

namespace Acts {
namespace Experimental {

class DetectorVolume;
class Portal;

class IPortalGenerator {
 public:
  virtual ~IPortalGenerator() = default;

  virtual std::vector<std::shared_ptr<Portal>> generate(
      const Transform3& dTransform, const VolumeBounds& dBounds,
      const std::shared_ptr<DetectorVolume>& dVolume) const = 0;
};

/// The Portal genertor definition
using PortalGenerator =
    OwningDelegate<std::vector<std::shared_ptr<Portal>>(
                       const Transform3&, const VolumeBounds&,
                       const std::shared_ptr<DetectorVolume>&),
                   IPortalGenerator>;

struct DefaultPortalGenerator final : public IPortalGenerator {
  /// @brief Generator function for creation of portal surfaces
  ///
  /// @param dTransform a contextually resolved transform
  /// @param dBounds the detecor volume bounds
  /// @param dVolume the reference to the detector volume which generates this volume
  ///
  /// @return a vector of newly created portals with registered inside volume
  std::vector<std::shared_ptr<Portal>> generate(
      const Transform3& dTransform, const VolumeBounds& dBounds,
      const std::shared_ptr<DetectorVolume>& dVolume) const final;
};

struct PortalAndSubPortalGenerator final : public IPortalGenerator {
  DefaultPortalGenerator defaultGenerator;

  /// @brief Calls the portal generation and adds registration to sub portals
  ///
  /// This code is split off the PortalGenerator code in order to allow
  /// unit testing of the portal generation wihtout detector volume construction
  ///
  /// @param dTransform a contextually resolved transform
  /// @param dBounds the detecor volume bounds
  /// @param dVolume the reference to the detector volume which generates this volume
  ///
  /// @return a vector of newly created portals with registered inside volume
  std::vector<std::shared_ptr<Portal>> generate(
      const Transform3& dTransform, const VolumeBounds& dBounds,
      const std::shared_ptr<DetectorVolume>& dVolume) const final;
};

template <typename Derived, typename... Args>
inline static PortalGenerator makePortalGenerator(Args... args) {
  PortalGenerator delegate;
  delegate.template connect<&Derived::generate>(
      std::make_unique<Derived>(std::forward<Args>(args)...));
  return delegate;
}

}  // namespace Experimental
}  // namespace Acts
