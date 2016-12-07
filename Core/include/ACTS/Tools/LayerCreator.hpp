// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// LayerCreator.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_GEOMETRYTOOLS_LAYERCREATOR_H
#define ACTS_GEOMETRYTOOLS_LAYERCREATOR_H 1

#include "ACTS/Tools/ILayerCreator.hpp"
#include "ACTS/Tools/ISurfaceArrayCreator.hpp"
#include "ACTS/Utilities/Logger.hpp"

#ifndef ACTS_LAYERCREATOR_TAKESMALLERBIGGER
#define ACTS_LAYERCREATOR_TAKESMALLERBIGGER
#define takeSmaller(current, test) current = current < test ? current : test
#define takeBigger(current, test) current  = current > test ? current : test
#define takeSmallerBigger(cSmallest, cBiggest, test)                           \
  takeSmaller(cSmallest, test);                                                \
  takeBigger(cBiggest, test)
#endif

namespace Acts {

class ISurfaceArrayCreator;
/// @class LayerCreator
///
/// The LayerCreator is able to build cylinde,r disc layers or plane layers from
/// detector elements
///
class LayerCreator : public ILayerCreator
{
public:
  ///  @struct Config
  ///  Configuration for the LayerCreator
  ///  This is the nexted configuration struct for the LayerCreator class
  struct Config
  {
    /// surface array helper
    std::shared_ptr<ISurfaceArrayCreator> surfaceArrayCreator = nullptr;
    /// cylinder module z tolerance : it counts at same z, if ...
    double cylinderZtolerance;
    /// cylinder module phi tolerance : it counts at same phi, if ...
    double cylinderPhiTolerance;

    // standard constructor
    Config() : cylinderZtolerance(10.), cylinderPhiTolerance(0.1) {}
  };

  /// Constructor
  ///
  /// @param lcConfig is the configuration object
  /// @param logger logging instance
  LayerCreator(const Config&           lcConfig,
               std::unique_ptr<Logger> logger
               = getDefaultLogger("LayerCreator", Logging::INFO));

  /// Destructor
  ~LayerCreator() = default;

  /// ILayerCreator interface method - returning a cylindrical layer
  ///
  /// @param surfaces is the vector of pointers to sensitive surfaces
  /// represented by this layer
  /// @pre the pointers to the sensitive surfaces in the surfaces vectors all
  /// need to be valid, since no check is performed
  /// @param envelopeR is the additional envelope applied in R
  /// @param envelopeZ is the additional envelope applied in z
  /// @param binsPhi is number of bins the sensitive surfaces are ordered in phi
  /// @param binsZ is number of bins the sensitive surfaces are ordered in Z
  /// @param transform is the (optional) transform of the layer
  ///
  /// @return shared pointer to a newly created layer
  LayerPtr
  cylinderLayer(const std::vector<const Surface*>&  surfaces,
                double                              envelopeR,
                double                              envelopeZ,
                size_t                              binsPhi,
                size_t                              binsZ,
                std::shared_ptr<Transform3D>        transform = nullptr,
                std::unique_ptr<ApproachDescriptor> ad
                = nullptr) const override;

  /// ILayerCreator interface method - returning a cylindrical layer
  ///
  /// @param surfaces is the vector of pointers to sensitive surfaces
  /// represented by this layer
  /// @pre the pointers to the sensitive surfaces in the surfaces vectors all
  /// need to be valid, since no check is performed
  /// @param layerRmin is the inner radius of the layer
  /// @param layerRmax is the outer radius of the layer
  /// @param layerHalfZ is the half length in z of the layer
  /// @param bTypePhi binning type in phi (equidistant/arbitrary)
  /// @param bTypeZ binning type in z (equidistant/arbitrary)
  /// @param transform is the (optional) transform of the layer
  ///
  /// @return shared pointer to a newly created layer
  LayerPtr
  cylinderLayer(const std::vector<const Surface*>&  surfaces,
                double                              layerRmin,
                double                              layerRmax,
                double                              layerHalfZ,
                BinningType                         bTypePhi,
                BinningType                         bTypeZ,
                std::shared_ptr<Transform3D>        transform = nullptr,
                std::unique_ptr<ApproachDescriptor> ad
                = nullptr) const override;

  /// ILayerCreator interface method - returning a cylindrical layer
  ///
  /// @param surfaces is the vector of pointers to sensitive surfaces
  /// represented by this layer
  /// @pre the pointers to the sensitive surfaces in the surfaces vectors all
  /// need to be valid, since no check is performed
  /// @param envelopeMinR is the additional envelope applied in R at Rmin
  /// @param envelopeMaxR is the additional envelope applied in R in Rmax
  /// @param envelopeZ is the additional envelope applied in z
  /// @param binsR is number of bins the sensitive surfaces are ordered in R
  /// @param binsPhi is number of bins the sensitive surfaces are ordered in Phi
  /// @param transform is the (optional) transform of the layer
  ///
  /// @return shared pointer to a newly created layer
  LayerPtr
  discLayer(const std::vector<const Surface*>&  surfaces,
            double                              envelopeMinR,
            double                              envelopeMaxR,
            double                              envelopeZ,
            size_t                              binsR,
            size_t                              binsPhi,
            std::shared_ptr<Transform3D>        transform = nullptr,
            std::unique_ptr<ApproachDescriptor> ad = nullptr) const override;

  /// ILayerCreator interface method - returning a cylindrical layer
  ///
  /// @param surfaces is the vector of pointers to sensitive surfaces
  /// represented by this layer
  /// @pre the pointers to the sensitive surfaces in the surfaces vectors all
  /// need to be valid, since no check is performed
  /// @param layerRmin is the inner radius of the layer
  /// @param layerRmax is the outer radius of the layer
  /// @param layerZmin is the minimum in z of the layer
  /// @param layerZmax is the maximum in z of the layer
  /// @param bTypeR binning type in r (equidistant/arbitrary)
  /// @param bTypePhi binning type in phi (equidistant/arbitrary)
  /// @param transform is the (optional) transform of the layer
  ///
  /// @return shared pointer to a newly created layer
  LayerPtr
  discLayer(const std::vector<const Surface*>&  surfaces,
            double                              layerZmin,
            double                              layerZmax,
            double                              layerRmin,
            double                              layerRmax,
            BinningType                         bTypeR,
            BinningType                         bTypePhi,
            std::shared_ptr<Transform3D>        transform = nullptr,
            std::unique_ptr<ApproachDescriptor> ad = nullptr) const override;

  /// ILayerCreator interface method - returning a cylindrical layer
  ///
  /// @param surfaces is the vector of pointers to sensitive surfaces
  /// represented by this layer
  /// @pre the pointers to the sensitive surfaces in the surfaces vectors all
  /// need to be valid, since no check is performed
  /// @param envelopeXY is the additional envelope applied in XY
  /// @param envelopeZ is the additional envelope applied in Z
  /// @param binsX is number of bins the sensitive surfaces are ordered in X
  /// @param binsY is number of bins the sensitive surfaces are ordered in Y
  /// @param transform is the (optional) transform of the layer
  ///
  /// @return shared pointer to a newly created layer
  LayerPtr
  planeLayer(const std::vector<const Surface*>&  surfaces,
             double                              envelopeXY,
             double                              envelopeZ,
             size_t                              binsX,
             size_t                              binsY,
             std::shared_ptr<Transform3D>        transform = nullptr,
             std::unique_ptr<ApproachDescriptor> ad = nullptr) const override;

  /// Set the configuration object
  /// @param lcConfig is the configuration struct
  void
  setConfiguration(const Config& lcConfig);

  /// Access th configuration object
  Config
  getConfiguration() const;

  /// set logging instance
  void
  setLogger(std::unique_ptr<Logger> logger);

private:
  /// Method to get the global extends in space for the module
  /// @todo shift to vertices of surfaces
  ///
  /// @param sf is the surface to be examinated
  /// @param minR minimal R extend
  /// @param maxR maximal R extend
  /// @param minPhi minimal phi extend
  /// @param maxPhi maximal phi extend
  /// @param minZ minimal z extend
  /// @param maxZ maximal z extend
  void
  moduleExtend(const Surface& sf,
               double&        minR,
               double&        maxR,
               double&        minPhi,
               double&        maxPhi,
               double&        minZ,
               double&        maxZ) const;

  /// Calculates the closest radial distance of a line
  ///
  /// @param pos1 is the first position on the line
  /// @param pos2 is the second position on the line
  ///
  /// @return is the closest distance
  double
  radialDistance(const Vector3D& pos1, const Vector3D& pos2) const;

  void
  binningParameters(const std::vector<const Acts::Surface*>& surfaces,
                    Acts::BinningValue                       bValue,
                    double&                                  minR,
                    double&                                  maxR,
                    double&                                  minPhi,
                    double&                                  maxPhi,
                    double&                                  minZ,
                    double&                                  maxZ,
                    size_t                                   nBins) const;

  /// configuration object
  Config m_cfg;

  /// Private acces method to the logger
  const Logger&
  logger() const
  {
    return *m_logger;
  }

  /// logging instance
  std::unique_ptr<Logger> m_logger;
};

inline LayerCreator::Config
LayerCreator::getConfiguration() const
{
  return m_cfg;
}

}  // end of namespace

#endif  // ACTS_GEOMETRYTOOLS_LAYERCREATOR_H
