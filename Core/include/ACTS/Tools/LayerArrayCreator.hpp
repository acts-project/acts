// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// LayerArrayCreator.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_TOOLS_LAYERARRAYCREATOR_H
#define ACTS_TOOLS_LAYERARRAYCREATOR_H 1

#ifndef ACTS_TOOLS_TAKESMALLERBIGGER
#define ACTS_TOOLS_TAKESMALLERBIGGER
#define takeSmaller(current, test) current = current < test ? current : test
#define takeBigger(current, test) current  = current > test ? current : test
#define takeSmallerBigger(cSmallest, cBiggest, test)                           \
  takeSmaller(cSmallest, test);                                                \
  takeBigger(cBiggest, test)
#endif

#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/Logger.hpp"
#include "ACTS/Tools/ILayerArrayCreator.hpp"
#include <algorithm>

namespace Acts {

class Surface;
class Layer;

/// @class LayerArrayCreator
///  The LayerArrayCreator is a simple Tool that helps to construct
///  LayerArrays from std::vector of Acts::CylinderLayer, Acts::DiscLayer,
/// Acts::PlaneLayer.
///
///  It fills the gaps automatically with Acts::NavigationLayer to be processed
/// easily in the
///  Navigation of the Extrapolation process.
///
/// @TODO Julia: make private tools private again after Gaudi update (bug in
/// Gaudi), marked with //b

class LayerArrayCreator : public ILayerArrayCreator
{
public:
  /// @strcut Config
  /// This is the nested configuration struct for the Layer Array creator
  struct Config
  {
    std::shared_ptr<Logger> logger;  ///< logging instance

    Config() : logger(getDefaultLogger("LayerArrayCreator", Logging::INFO)) {}
  };

  /// Constructor 
  /// @param cfg is the configuration struct that steers the behavoir
  LayerArrayCreator(const Config& cfg) : m_config(cfg) {}
  
  /// Destructor 
  virtual ~LayerArrayCreator() = default;

  /// LayerArrayCreator interface method
  /// @param layers are the layers to be moved into an array
  /// @min is the minimul value for binning
  /// @max is the maximum value for binning
  /// @btype is the binning type 
  /// @bvalue is the value in which the binning should be done 
  std::unique_ptr<const LayerArray>
  layerArray(const LayerVector& layers,
             double             min,
             double             max,
             BinningType        btype  = arbitrary,
             BinningValue       bvalue = binX) const override;

  /// Set configuration method 
  void
  setConfiguration(const Config& c)
  {
    m_config = c;
  }

  /// Get configuration method 
  Config
  getConfiguration() const
  {
    return m_config;
  }

private:
  /// Private access method to the logging instance
  const Logger&
  logger() const
  {
    return *m_config.logger;
  }

  /// logging instance
  Config m_config;
  
  /// Private helper method for creating a surface for
  /// the NavigationLayer
  Surface*
  createNavigationSurface(const Layer& layer,
                          BinningValue bvalue,
                          double       offset) const;
};

}  // end of namespace

#endif  // ACTS_GEOMETRYTOOLS_LAYERARRAYCREATOR_H
