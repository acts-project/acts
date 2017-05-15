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

#ifndef ACTS_EXAMPLES_GENERICLAYERBUILDER_H
#define ACTS_EXAMPLES_GENERICLAYERBUILDER_H 1

#include "ACTS/Layers/Layer.hpp"
#include "ACTS/Material/Material.hpp"
#include "ACTS/Material/MaterialProperties.hpp"
#include "ACTS/Tools/ILayerBuilder.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/Logger.hpp"

namespace Acts {

class ILayerCreator;

class Surface;
class DetecorElementBase;
typedef std::pair<const Surface*, Vector3D> SurfacePosition;

/// @class GenericLayerBuilder
///
/// The GenericLayerBuilder is able to build cylinder & disc layers from python
/// input.
/// This is ment for the simple detector examples.
///
class GenericLayerBuilder : public ILayerBuilder
{
public:
  /// @struct Config
  /// Nested configuration struct for the GenericLayerBuilder
  struct Config
  {
    /// the string based identification
    std::string layerIdentification = "";
    /// a single paramater for the approach surface envelope
    double approachSurfaceEnvelope = 0.5;
    /// central layer specification
    /// bin multipliers in rphi,z for finer module binning
    std::pair<int, int> centralLayerBinMultipliers;
    /// layer radii for the sensitive layers
    std::vector<double> centralLayerRadii;
    /// the (additional) layer envelope in R/Z
    std::vector<std::pair<double, double>> centralLayerEnvelopes;
    /// the material concentration: -1 inner, 0 central, 1 outer
    std::vector<int> centralLayerMaterialConcentration;
    /// the assigned material propertis @todo change to surface material
    std::vector<MaterialProperties> centralLayerMaterialProperties;
    /// teh binning schema: nPhi x nZ
    std::vector<std::pair<int, int>> centralModuleBinningSchema;
    /// the module center positions
    std::vector<std::vector<Vector3D>> centralModulePositions;
    /// the module tilt for this layer
    std::vector<double> centralModuleTiltPhi;
    /// the module bounds: local x
    std::vector<double> centralModuleHalfX;
    /// the module bounds: local y
    std::vector<double> centralModuleHalfY;
    /// the module bounds: local z -> thickness
    std::vector<double> centralModuleThickness;
    /// the central volume readout schema
    std::vector<size_t> centralModuleReadoutBinsX;
    /// the central volume readout schema
    std::vector<size_t> centralModuleReadoutBinsY;
    /// the central volume readout schema
    std::vector<int> centralModuleReadoutSide;
    /// the central volume readout schema
    std::vector<double> centralModuleLorentzAngle;
    /// the module material @todo change to surface material
    std::vector<Material> centralModuleMaterial;
    /// the module front side stereo (if exists)
    std::vector<double> centralModuleFrontsideStereo;
    /// the module back side stereo (if exists)
    std::vector<double> centralModuleBacksideStereo;
    /// the module gap between frontside and backside
    std::vector<double> centralModuleBacksideGap;

    /// the layers at p/e side
    /// bin multipliers in r,phi for finer module binning
    std::pair<int, int> posnegLayerBinMultipliers;
    /// layer positions in Z
    std::vector<double> posnegLayerPositionsZ;
    /// the
    std::vector<double> posnegLayerEnvelopeR;
    /// the material concentration: -1 inner, 0 central, 1 outer
    std::vector<int> posnegLayerMaterialConcentration;
    /// the material prooperties @todo change to surface material
    std::vector<MaterialProperties> posnegLayerMaterialProperties;
    /// the module center positions
    std::vector<std::vector<std::vector<Vector3D>>> posnegModulePositions;
    /// the phi binning
    std::vector<std::vector<size_t>> posnegModulePhiBins;
    /// the module bounds: min halfx
    std::vector<std::vector<double>> posnegModuleMinHalfX;
    /// the module bounds: max halfx
    std::vector<std::vector<double>> posnegModuleMaxHalfX;
    /// the module bounds: local y
    std::vector<std::vector<double>> posnegModuleHalfY;
    /// the module bounds: local z -> thickness
    std::vector<std::vector<double>> posnegModuleThickness;
    /// the central volume readout schema
    std::vector<std::vector<size_t>> posnegModuleReadoutBinsX;
    /// the central volume readout schema
    std::vector<std::vector<size_t>> posnegModuleReadoutBinsY;
    /// the central volume readout schema
    std::vector<std::vector<int>> posnegModuleReadoutSide;
    /// the central volume readout schema
    std::vector<std::vector<double>> posnegModuleLorentzAngle;
    /// the module material @todo change to surface material
    std::vector<std::vector<Material>> posnegModuleMaterial;
    /// the module front side stereo (if exists)
    std::vector<std::vector<double>> posnegModuleFrontsideStereo;
    /// the module back side stereo (if exists)
    std::vector<std::vector<double>> posnegModuleBacksideStereo;
    /// the module gap between frontside and backside
    std::vector<std::vector<double>> posnegModuleBacksideGap;

    /// helper tools: layer creator
    std::shared_ptr<const ILayerCreator> layerCreator = nullptr;
    /// helper tools: central passiva layer builder
    std::shared_ptr<const ILayerBuilder> centralPassiveLayerBuilder = nullptr;
    /// helper tools: p/n passive layer builder
    std::shared_ptr<const ILayerBuilder> posnegPassiveLayerBuilder = nullptr;
  };

  /// Constructor
  /// @param glbConfig is the configuration class
  GenericLayerBuilder(const Config&                 glbConfig,
                      std::unique_ptr<const Logger> logger
                      = getDefaultLogger("GenericLayerBuilder", Logging::INFO));

  /// Destructor
  ~GenericLayerBuilder();

  /// LayerBuilder interface method - returning the layers at negative side
  const LayerVector
  negativeLayers() const override;

  /// LayerBuilder interface method - returning the central layers
  const LayerVector
  centralLayers() const override;

  /// LayerBuilder interface method - returning the layers at negative side
  const LayerVector
  positiveLayers() const override;

  /// ILayerBuilder method
  const std::string&
  identification() const override
  {
    return m_cfg.layerIdentification;
  }

  /// set the configuration object
  void
  setConfiguration(const Config& glbConfig);

  /// get the configuration object
  Config
  getConfiguration() const;

  /// set logging instance
  void
  setLogger(std::unique_ptr<const Logger> logger);

private:
  void
  constructLayers();

  LayerVector m_nLayers;  ///< layers on negative side
  LayerVector m_cLayers;  ///< layers on central side
  LayerVector m_pLayers;  ///< layers on positive side

  std::vector<const DetectorElementBase*>
      m_centralModule;  ///< acts as detector store

  std::vector<const DetectorElementBase*>
      m_posnegModule;  ///< acts as detector store

  /// Configuration member
  Config m_cfg;

  /// Private access to the looging instance
  const Logger&
  logger() const
  {
    return *m_logger;
  }

  /// the loging instance
  std::unique_ptr<const Logger> m_logger;
};

inline const LayerVector
GenericLayerBuilder::positiveLayers() const
{
  return m_pLayers;
}

inline const LayerVector
GenericLayerBuilder::negativeLayers() const
{
  return m_nLayers;
}

inline const LayerVector
GenericLayerBuilder::centralLayers() const
{
  return m_cLayers;
}

inline GenericLayerBuilder::Config
GenericLayerBuilder::getConfiguration() const
{
  return m_cfg;
}
}  // end of namespace

#endif  // ACTS_EXAMPLES_GENERICLAYERBUILDER_H
