// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// PassiveLayerBuilder.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_GEOMETRYTOOLS_PASSIVELAYERBUILDER_H
#define ACTS_GEOMETRYTOOLS_PASSIVELAYERBUILDER_H 1

#include "ACTS/Layers/Layer.hpp"
#include "ACTS/Material/Material.hpp"
#include "ACTS/Tools/ILayerBuilder.hpp"
#include "ACTS/Utilities/Logger.hpp"

namespace Acts {

///  @class PassiveLayerBuilder
///
///   The PassiveLayerBuilder is able to build cylinder & disc layers with given
///  detector
///
///   @todo Julia: make private tools private again after Gaudi update (bug in
///  Gaudi), marked with //b

class PassiveLayerBuilder : public ILayerBuilder
{
public:
  /// @struct Config
  /// Configuration struct for the passive layer builder
  /// This nested struct is used to configure the layer building
  struct Config
  {
    /// string based identification
    std::string layerIdentification;

    std::vector<double>   centralLayerRadii;        ///< central layer specs
    std::vector<double>   centralLayerHalflengthZ;  ///< central layer specs
    std::vector<double>   centralLayerThickness;    ///< central layer specs
    std::vector<Material> centralLayerMaterial;     ///< central layer specs

    // the layers at p/e side
    std::vector<double>   posnegLayerPositionZ;  ///< p/n layer specs
    std::vector<double>   posnegLayerRmin;       ///< p/n layer specs
    std::vector<double>   posnegLayerRmax;       ///< p/n layer specs
    std::vector<double>   posnegLayerThickness;  ///< p/n layer specs
    std::vector<Material> posnegLayerMaterial;   ///< p/n layer specs
  };

  /// Constructor
  /// @param plConfig is the ocnfiguration struct that steers behavior
  /// @param logger logging instance
  PassiveLayerBuilder(const Config&           plConfig,
                      std::unique_ptr<Logger> logger
                      = getDefaultLogger("PassiveLayerBuilder", Logging::INFO));

  /// Destructor
  virtual ~PassiveLayerBuilder() = default;

  /// LayerBuilder interface method
  /// @return  the layers at negative side
  const LayerVector
  negativeLayers() const override;

  /// LayerBuilder interface method
  /// @return the layers at the central sector
  const LayerVector
  centralLayers() const override;

  /// LayerBuilder interface method
  /// @return  the layers at positive side
  const LayerVector
  positiveLayers() const override;

  /// Name identification
  /// @return the string based identification
  const std::string&
  identification() const override
  {
    return m_cfg.layerIdentification;
  }

  /// Set configuration method
  /// @param plConfig is a configuration struct
  /// it overwrites the current configuration
  void
  setConfiguration(const Config& plConfig);

  /// Get configuration method
  Config
  getConfiguration() const;

  /// set logging instance
  void
  setLogger(std::unique_ptr<Logger> logger);

protected:
  Config m_cfg;  //!< configuration

private:
  const Logger&
  logger() const
  {
    return *m_logger;
  }

  /// logging instance
  std::unique_ptr<Logger> m_logger;

  bool
  constructLayers() const;

  mutable LayerVector m_nLayers;  ///< layers on negative side
  mutable LayerVector m_cLayers;  ///< layers on central side
  mutable LayerVector m_pLayers;  ///< layers on positive side

  mutable bool m_constructionFlag;  ///< indicator if the layer construction has
                                    /// been done already
};

inline PassiveLayerBuilder::Config
PassiveLayerBuilder::getConfiguration() const
{
  return m_cfg;
}

inline const LayerVector
PassiveLayerBuilder::positiveLayers() const
{
  if (not m_constructionFlag) constructLayers();
  return m_pLayers;
}

inline const LayerVector
PassiveLayerBuilder::negativeLayers() const
{
  if (not m_constructionFlag) constructLayers();
  return m_nLayers;
}

inline const LayerVector
PassiveLayerBuilder::centralLayers() const
{
  if (not m_constructionFlag) constructLayers();
  return m_cLayers;
}

}  // end of namespace

#endif  // ACTS_TOOLS_PASSIVELAYERBUILDER_H
