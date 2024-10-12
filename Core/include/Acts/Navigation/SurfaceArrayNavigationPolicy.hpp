// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/SurfaceArrayCreator.hpp"
#include "Acts/Navigation/INavigationPolicy.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"

#pragma once

namespace Acts {

class SurfaceArrayNavigationPolicy : public INavigationPolicy {
 public:
  enum class LayerType { Cylinder, Disc, Plane };

  struct Config {
    SurfaceArrayCreator::Config surfaceArrayConfig{};
    LayerType layerType = LayerType::Cylinder;
  };

  explicit SurfaceArrayNavigationPolicy(const GeometryContext& gctx,
                                        const TrackingVolume& volume,
                                        const Logger& logger, Config config);

  void updateState(const NavigationArguments& args) const;

  void connect(NavigationDelegate& delegate) const override;

  friend std::ostream& operator<<(std::ostream& os,
                                  const LayerType& layerType) {
    switch (layerType) {
      case LayerType::Cylinder:
        os << "Cylinder";
        break;
      case LayerType::Disc:
        os << "Disc";
        break;
      case LayerType::Plane:
        os << "Plane";
        break;
    }
    return os;
  }

 private:
  std::unique_ptr<SurfaceArray> m_surfaceArray{};
};

}  // namespace Acts
