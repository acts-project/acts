// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// ILayerBuilder.h, Acts project
///////////////////////////////////////////////////////////////////

#pragma once
#include <memory>
#include <string>
#include <vector>

namespace Acts {

class Layer;
using LayerPtr    = std::shared_ptr<const Layer>;
using LayerVector = std::vector<LayerPtr>;

/// @class ILayerBuilder
///
/// Interface class for ILayerBuilders in a typical
/// | EC- | Central | EC+ |
/// detector setup.
///
class ILayerBuilder
{
public:
  /// Virtual destructor
  virtual ~ILayerBuilder() = default;
  /// LayerBuilder interface method
  /// @return  the layers at negative side
  virtual const LayerVector
  negativeLayers() const = 0;

  /// LayerBuilder interface method
  /// @return the layers at the central sector
  virtual const LayerVector
  centralLayers() const = 0;

  /// LayerBuilder interface method
  /// @return  the layers at positive side
  virtual const LayerVector
  positiveLayers() const = 0;

  /// Name identification
  /// @return the string based identification
  virtual const std::string&
  identification() const = 0;
};

}  // namespace
