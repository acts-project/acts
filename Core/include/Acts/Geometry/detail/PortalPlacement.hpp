// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.



#pragma once

#include "Acts/Geometry/DetectorElementBase.hpp"
#include "Acts/Surfaces/RegularSurface.hpp"

namespace Acts {
class VolumePlacementBase;

namespace detail {

class PortalPlacement : public DetectorElementBase {
 public:
  friend class Acts::VolumePlacementBase;
  const Transform3& localToGlobalTransform(
      const GeometryContext& gctx) const final;

  const Surface& surface() const final;

  Surface& surface() final;

  double thickness() const final;

  bool isSensitive() const final;

  std::size_t index() const;

  const Transform3& portalToVolumeCenter() const;

 protected:
  PortalPlacement(const std::size_t portalIdx, const Transform3& portalTrf,
                  VolumePlacementBase* parent,
                  std::shared_ptr<RegularSurface>&& surface);

 private:
  Transform3 m_interalTrf{Transform3::Identity()};
  std::shared_ptr<RegularSurface> m_surface{};
  VolumePlacementBase* m_parent{nullptr};

  std::size_t m_portalIdx{0ul};
};
}  // namespace detail
}  // namespace Acts
