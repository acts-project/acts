// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_DD4HEPPLUGIN_DD4HEPLAYERBUILDER_H
#define ACTS_DD4HEPPLUGIN_DD4HEPLAYERBUILDER_H 1

#include "ACTS/Surfaces/Surface.hpp"
#include "ACTS/Tools/ILayerBuilder.hpp"
#include "ACTS/Utilities/BinUtility.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/Logger.hpp"

class TGeoMatrix;

namespace DD4hep {
namespace Geometry {
  class DetElement;
}
}

namespace Acts {

class ILayerCreator;

/// @brief build layers of one cylinder-endcap setup from DD4hep input
///
/// This class is an implementation of the Acts::ILayerBuilder,
/// creating the central (layers of barrel), the negative and positive layers
/// (layers of endcaps) of one hierarchy (e.g. PixelDetector, StripDetector,...)
/// with input from DD4hep.

class DD4hepLayerBuilder : public ILayerBuilder
{
public:
  /// @struct Config
  /// nested configuration struct for steering of the layer builder
  struct Config
  {
    /// string based identification
    std::string configurationName = "undefined";
    /// layer creator which is internally used to build layers
    std::shared_ptr<const ILayerCreator> layerCreator = nullptr;
    /// the binning type of the contained surfaces in phi
    /// (equidistant/arbitrary)
    BinningType bTypePhi = equidistant;
    /// the binning type of the contained surfaces in r
    /// (equidistant/arbitrary)
    BinningType bTypeR = equidistant;
    /// the binning type of the contained surfaces in z
    /// (equidistant/arbitrary)
    BinningType bTypeZ = equidistant;
    /// the DD4hep::DetElements of the layers of the negative volume (negative
    /// endcap)
    /// @note if the current volume has no endcaps or no layers this parameter
    /// will not be set
    std::vector<DD4hep::Geometry::DetElement> negativeLayers;
    /// the DD4hep::DetElements of the layers of the central volume (barrel)
    /// @note if the current volume has no layers this parameter will not be set
    std::vector<DD4hep::Geometry::DetElement> centralLayers;
    /// the DD4hep::DetElements of the layers of the positive volume (positive
    /// endcap)
    /// @note if the current volume has no endcaps or no layers this parameter
    /// will not be set
    std::vector<DD4hep::Geometry::DetElement> positiveLayers;
    /// @param buildDigitizationModules Flag indicating if the
    /// Acts::DigitizationModule (needed for Acts geometric digitization) will
    /// be build for every single sensitive DD4hep DetElement translating
    /// directly the DD4hep Segmentation.
    /// @note For more information please see Acts::convertDD4hepDetector() &
    /// Acts::ActsExtension.
    bool buildDigitizationModules = false;
  };

  /// Constructor
  /// @param config is the configuration struct
  /// @param logger is the logging instance
  DD4hepLayerBuilder(const Acts::DD4hepLayerBuilder::Config& config,
                     std::unique_ptr<const Logger>           logger);
  /// Destructor
  ~DD4hepLayerBuilder();

  /// LayerBuilder interface method
  /// @return  the layers at negative side
  virtual const LayerVector
  negativeLayers() const final;

  /// LayerBuilder interface method
  /// @return the layers at the central sector
  virtual const LayerVector
  centralLayers() const final;

  /// LayerBuilder interface method
  /// @return  the layers at positive side
  virtual const LayerVector
  positiveLayers() const final;

  /// Name identification
  /// @return the string based identification of this configuration
  virtual const std::string&
  identification() const final;

  /// set the configuration object
  /// @param cfg is the configuration struct
  void
  setConfiguration(const Config& cfg);

  /// get the configuration object
  Config
  getConfiguration() const;

  /// set logging instance
  void
  setLogger(std::unique_ptr<const Logger> logger);

private:
  /// configruation object
  Config m_cfg;

  /// logging instance
  std::unique_ptr<const Logger> m_logger;

  /// Private access to the logger
  const Logger&
  logger() const
  {
    return *m_logger;
  }

  /// Private helper function collecting all sensitive detector elements of a
  /// layer
  /// @param detElement the DD4hep::DetElement of the layer
  /// @param surfaces the vector of surfaces which should be filled with the
  /// sensitive detector elements
  /// @param axes the orientation of the modules to the ACTS frame
  void
  collectSensitive(const DD4hep::Geometry::DetElement& detElement,
                   std::vector<const Acts::Surface*>&  surfaces,
                   const std::string&                  axes = "XYZ") const;

  // Private helper function to convert the TGeo transformation matrix into a
  // ACTS transformation matrix
  /// @param tGeoTrans TGeo transformation matrix which should be converted
  std::shared_ptr<const Acts::Transform3D>
  convertTransform(const TGeoMatrix* tGeoTrans) const;
};

inline const std::string&
DD4hepLayerBuilder::identification() const
{
  return m_cfg.configurationName;
}

inline DD4hepLayerBuilder::Config
DD4hepLayerBuilder::getConfiguration() const
{
  return m_cfg;
}

}  // end of namespace

#endif  // ACTS_DD4HEPPLUGIN_DD4HEPLAYERBUILDER_H
