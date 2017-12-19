// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// SurfaceArrayCreator.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_TOOLS_SURFACERARRAYCREATOR_H
#define ACTS_TOOLS_SURFACERARRAYCREATOR_H 1

#include <boost/none.hpp>
#include <boost/optional.hpp>
#include "ACTS/Layers/ProtoLayer.hpp"
#include "ACTS/Surfaces/PlanarBounds.hpp"
#include "ACTS/Surfaces/Surface.hpp"
#include "ACTS/Surfaces/SurfaceArray.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/Logger.hpp"
#include "ACTS/Utilities/Units.hpp"

namespace Acts {

namespace Test {
  struct SurfaceArrayCreatorFixture;
}

// typedef BinnedArray<const Surface*> SurfaceArray_old;
using SurfaceMatcher
    = std::function<bool(BinningValue, const Surface*, const Surface*)>;

typedef std::vector<const Surface*> SurfaceVector;
typedef std::vector<SurfaceVector>  SurfaceMatrix;
typedef std::vector<SurfaceMatrix>  SurfaceGrid_old;

typedef std::vector<Vector3D> V3Vector;
typedef std::vector<V3Vector> V3Matrix;

/// @class SurfaceArrayCreator
///
/// It is designed create sub surface arrays to be ordered on Surfaces
///
/// @todo write more documentation on how this is done
class SurfaceArrayCreator
{
public:
  friend Acts::Test::SurfaceArrayCreatorFixture;

  struct ProtoAxis
  {
    BinningType         bType;
    BinningValue        bValue;
    size_t              nBins;
    double              min;
    double              max;
    std::vector<double> binEdges;
  };

  // Configuration struct
  struct Config
  {
    SurfaceMatcher surfaceMatcher = SurfaceArrayCreator::isSurfaceEquivalent;
  };

  /// Constructor with default config
  ///
  /// @param logger logging instance
  SurfaceArrayCreator(std::unique_ptr<const Logger> logger
                      = getDefaultLogger("SurfaceArrayCreator", Logging::INFO))
    : m_cfg(Config()), m_logger(std::move(logger))
  {
  }
  /// Constructor with explicit config
  ///
  /// @param cfg Explicit config struct
  /// @param logger logging instance
  SurfaceArrayCreator(const Config&                 cfg,
                      std::unique_ptr<const Logger> logger
                      = getDefaultLogger("SurfaceArrayCreator", Logging::INFO))
    : m_cfg(cfg), m_logger(std::move(logger))
  {
  }

  /// Destructor
  virtual ~SurfaceArrayCreator() = default;

  /// SurfaceArrayCreator interface method
  ///
  /// - create an array in a cylinder, binned in phi, z when extremas and
  /// bin numbers are known
  /// @param surfaces is the vector of pointers to sensitive surfaces
  /// to be ordered on the cylinder
  /// @pre the pointers to the sensitive surfaces in the surfaces vectors all
  /// need to be valid, since no check is performed
  /// @param protoLayer The proto layer containing the layer size
  /// @param binsPhi is the number of bins in phi for the surfaces
  /// @param binsZ is the number of bin in Z for the surfaces
  /// @param transform is the (optional) additional transform applied
  ///
  /// @return a unique pointer to a new SurfaceArray
  std::unique_ptr<SurfaceArray>
  surfaceArrayOnCylinder(const std::vector<const Surface*>& surfaces,
                         size_t                             binsPhi,
                         size_t                             binsZ,
                         boost::optional<ProtoLayer> protoLayer = boost::none,
                         std::shared_ptr<const Transform3D> transform
                         = nullptr) const;

  /// SurfaceArrayCreator interface method
  ///
  /// - create an array in a cylinder, binned in phi, z when extremas and bin
  /// numbers are unknown - this method goes through the surfaces and finds out
  /// the needed information
  /// @param surfaces is the vector of pointers to sensitive surfaces
  /// to be ordered on the cylinder
  /// @pre the pointers to the sensitive surfaces in the surfaces vectors all
  /// need to be valid, since no check is performed
  /// @param protoLayer The proto layer containing the layer size
  /// @param bTypePhi the binning type in phi direction (equidistant/aribtrary)
  /// @param bTypeZ the binning type in z direction (equidistant/aribtrary)
  /// @param transform is the (optional) additional transform applied
  ///
  /// @return a unique pointer a new SurfaceArray
  std::unique_ptr<Acts::SurfaceArray>
  surfaceArrayOnCylinder(const std::vector<const Surface*>& surfaces,
                         BinningType                 bTypePhi   = equidistant,
                         BinningType                 bTypeZ     = equidistant,
                         boost::optional<ProtoLayer> protoLayer = boost::none,
                         std::shared_ptr<const Transform3D> transform
                         = nullptr) const;

  /// SurfaceArrayCreator interface method
  /// - create an array on a disc, binned in r, phi when extremas and
  /// bin numbers are known
  ///
  /// @param surfaces is the vector of pointers to sensitive surfaces
  /// to be ordered on the disc
  /// @pre the pointers to the sensitive surfaces in the surfaces vectors all
  /// need to be valid, since no check is performed
  /// @param protoLayer The proto layer containing the layer size
  /// @param binsPhi is the number of bins in phi for the surfaces
  /// @param binsR is the number of bin in R for the surfaces
  /// @param transform is the (optional) additional transform applied
  ///
  /// @return a unique pointer a new SurfaceArray
  std::unique_ptr<SurfaceArray>
  surfaceArrayOnDisc(const std::vector<const Surface*>& surfaces,
                     size_t                             binsR,
                     size_t                             binsPhi,
                     boost::optional<ProtoLayer> protoLayer = boost::none,
                     std::shared_ptr<const Transform3D> transform
                     = nullptr) const;

  /// SurfaceArrayCreator interface method
  ///
  /// - create an array in a cylinder, binned in phi, r when extremas and bin
  /// numbers are unknown - this method goes through the surfaces and finds out
  /// the needed information
  /// @param surfaces is the vector of pointers to sensitive surfaces
  /// to be ordered on the disc
  /// @pre the pointers to the sensitive surfaces in the surfaces vectors all
  /// need to be valid, since no check is performed
  /// @param protoLayer The proto layer containing the layer size
  /// @param bTypeR the binning type in r direction (equidistant/aribtrary)
  /// @param bTypePhi the binning type in phi direction (equidistant/aribtrary)
  /// @param transform is the (optional) additional transform applied
  ///
  /// @return a unique pointer a new SurfaceArray
  std::unique_ptr<Acts::SurfaceArray>
  surfaceArrayOnDisc(const std::vector<const Surface*>& surfaces,
                     BinningType                        bTypeR,
                     BinningType                        bTypePhi,
                     boost::optional<ProtoLayer> protoLayer = boost::none,
                     std::shared_ptr<const Transform3D> transform
                     = nullptr) const;

  /// SurfaceArrayCreator interface method
  /// - create an array on a plane
  ///
  /// @param surfaces is the vector of pointers to sensitive surfaces
  /// to be ordered on the plane
  /// @pre the pointers to the sensitive surfaces in the surfaces vectors all
  /// need to be valid, since no check is performed
  /// @param halflengthX is the half length in X
  /// @param halflengthY is the half length in Y
  /// @param binsX is the number of bins in X
  /// @param binsY is the number of bins in Y
  /// @param transform is the (optional) additional transform applied
  ///
  /// @return a unique pointer a new SurfaceArray
  std::unique_ptr<SurfaceArray>
  surfaceArrayOnPlane(const std::vector<const Surface*>& surfaces,
                      double                             halflengthX,
                      double                             halflengthY,
                      size_t                             binsX,
                      size_t                             binsY,
                      std::shared_ptr<const Transform3D> transform
                      = nullptr) const;

  static bool
  isSurfaceEquivalent(BinningValue bValue, const Surface* a, const Surface* b)
  {

    if (bValue == Acts::binPhi)
      return (std::abs(a->center().phi() - b->center().phi())
              < M_PI * 1 / 180.);

    if (bValue == Acts::binZ)
      return (std::abs(a->center().z() - b->center().z()) < Acts::units::_um);

    if (bValue == Acts::binR)
      return (std::abs(a->center().perp() - b->center().perp())
              < Acts::units::_um);

    return false;
  }

  /// Set logging instance
  /// @param logger is the logging instance to be set
  void
  setLogger(std::unique_ptr<const Logger> logger)
  {
    m_logger = std::move(logger);
  }

private:
  /// configuration object
  Config m_cfg;

  /// Private access to logger
  const Logger&
  logger() const
  {
    return *m_logger;
  }

  size_t
  determineBinCount(const std::vector<const Surface*>& surfaces,
                    BinningValue                       bValue) const;

  /// SurfaceArrayCreator internal method
  /// Creates a variable @c ProtoAxis from a vector of (unsorted) surfaces with
  /// PlanarBounds
  /// It loops through the surfaces and finds out the needed information
  /// First the surfaces are sorted in the binning direction and the so called
  /// "key" surfaces (surfaces with different positions in the binning
  /// direction) are extracted. The boundary value between two surfaces is the
  /// mean value of the two center position in the binning direction. The first
  /// and the last boundaries are calculated from the vertices of the first and
  /// last surface.
  /// @note currently implemented for phi, r and z bining
  /// @todo implement for x,y binning
  /// @param surfaces are the sensitive surfaces to be
  /// @param bValue the BinningValue in which direction should be binned
  /// (currently possible: binPhi, binR, binZ)
  /// @param transform is the (optional) additional transform applied
  /// @return Instance of @c ProtoAxis containing determined properties
  /// @note This only creates the @c ProtoAxis, this needs to be turned
  ///       into an actual @c Axis object to be used
  ProtoAxis
  createVariableAxis(const std::vector<const Surface*>& surfaces,
                     BinningValue                       bValue,
                     Transform3D&                       transform) const;

  /// SurfaceArrayCreator internal method
  /// Creates a equidistant @c ProtoAxis when the extremas and the bin number
  /// are
  /// It loops through the surfaces and finds out the needed information
  /// First the surfaces are sorted in the binning direction and the so called
  /// "key" surfaces (surfaces with different positions in the binning
  /// direction) are extracted. The number of key surfaces euqals the number of
  /// bins. Afterwards the minimum and maximum are calculated by
  /// subtracting/adding half of a bin size to the center position (in the
  /// binning direction) to the first/last surface.
  /// @note currently implemented for phi, r and z bining
  /// @todo implement for x,y binning
  /// @param surfaces are the sensitive surfaces to be
  /// @param bValue the BinningValue in which direction should be binned
  /// (currently possible: binPhi, binR, binZ)
  /// @param protoLayer Instance of @c ProtoLayer holding generic layer info
  /// @param transform is the (optional) additional transform applied
  /// @return Instance of @c ProtoAxis containing determined properties
  /// @note This only creates the @c ProtoAxis, this needs to be turned
  ///       into an actual @c Axis object to be used
  ProtoAxis
  createEquidistantAxis(const std::vector<const Surface*>& surfaces,
                        BinningValue                       bValue,
                        ProtoLayer                         protoLayer,
                        Transform3D&                       transform,
                        size_t                             nBins = 0) const;

  /// SurfaceArrayCreator internal method
  /// @brief Creates a SurfaceGridLookup instance within an any
  /// This is essentially a factory which absorbs some if/else logic
  /// that is required by the templating.
  /// @tparam wrapA AxisWrapping of axis A
  /// @tparam wrapB AxisWrapping of axis B
  /// @tparam F1 type-deducted value of g2l lambda
  /// @tparam F2 type-deducted value of l2g lambda
  /// @param globalToLocal transform callable
  /// @param localToGlobal transform callable
  /// @param pAxisA ProtoAxis object for axis A
  /// @param pAxisB ProtoAxis object for axis B
  template <detail::AxisWrapping wrapA,
            detail::AxisWrapping wrapB,
            typename F1,
            typename F2>
  SurfaceArray::AnySurfaceGridLookup_t
  makeSurfaceGridLookup2D(F1        globalToLocal,
                          F2        localToGlobal,
                          ProtoAxis pAxisA,
                          ProtoAxis pAxisB) const
  {

    // this becomes completely unreadable otherwise
    // clang-format off
    if (pAxisA.bType == equidistant && pAxisB.bType == equidistant) {

      detail::Axis<detail::AxisType::Equidistant, wrapA> axisA(pAxisA.min, pAxisA.max, pAxisA.nBins);
      detail::Axis<detail::AxisType::Equidistant, wrapB> axisB(pAxisB.min, pAxisB.max, pAxisB.nBins);
      return SurfaceArray::SurfaceGridLookup<decltype(axisA), decltype(axisB)>(globalToLocal,
                                                                                 localToGlobal,
                                                                                 std::make_tuple(axisA, axisB));

    } else if (pAxisA.bType == equidistant && pAxisB.bType == arbitrary) {

      detail::Axis<detail::AxisType::Equidistant, wrapA> axisA(pAxisA.min, pAxisA.max, pAxisA.nBins);
      detail::Axis<detail::AxisType::Variable, wrapB> axisB(pAxisB.binEdges);
      return SurfaceArray::SurfaceGridLookup<decltype(axisA), decltype(axisB)>(globalToLocal,
                                                                                 localToGlobal,
                                                                                 std::make_tuple(axisA, axisB));
      
    } else if (pAxisA.bType == arbitrary && pAxisB.bType == equidistant) {

      detail::Axis<detail::AxisType::Variable, wrapA> axisA(pAxisA.binEdges);
      detail::Axis<detail::AxisType::Equidistant, wrapB> axisB(pAxisB.min, pAxisB.max, pAxisB.nBins);
      return SurfaceArray::SurfaceGridLookup<decltype(axisA), decltype(axisB)>(globalToLocal,
                                                                                 localToGlobal,
                                                                                 std::make_tuple(axisA, axisB));

    } else /*if (pAxisA.bType == arbitrary && pAxisB.bType == arbitrary)*/ {

      detail::Axis<detail::AxisType::Variable, wrapA> axisA(pAxisA.binEdges);
      detail::Axis<detail::AxisType::Variable, wrapB> axisB(pAxisB.binEdges);
      return SurfaceArray::SurfaceGridLookup<decltype(axisA), decltype(axisB)>(globalToLocal,
                                                                                 localToGlobal,
                                                                                 std::make_tuple(axisA, axisB));
    }
    // clang-format on
  }

  /// logging instance
  std::unique_ptr<const Logger> m_logger;

  /// Private helper method to complete the binning
  ///
  ///
  ///  given a grid point o
  ///    |  0  |  1 |  2  |  3 |  4  |
  ///    ------------------------------
  ///  0 |  x  |    |     |    |  x  |
  ///  1 |     |    |  o  |    |     |
  ///  2 |  x  |    |     |    |  x  |
  ///
  /// This is being called when you chose to use more bins thans surfaces
  /// I.e. to put a finer granularity binning onto your surface
  /// Neighbour bins are then filled to contain pointers as well
  /// This method delegates to SurfaceGridLookup itself.
  /// @param sl The @c SurfaceGridLookup
  /// @param surfaces the surfaces
  template <class T>
  void
  completeBinning(T& sl, const std::vector<const Surface*>& surfaces) const
  {
    ACTS_VERBOSE("Complete binning by filling closest neighbour surfaces into "
                 "empty bins.");

    size_t binCompleted = sl.completeBinning(surfaces);

    ACTS_VERBOSE("       filled  : " << binCompleted
                                     << " (includes under/overflow)");
  }

  /// Private helper method to transform the  vertices of surface bounds into
  /// global coordinates
  /// @param surface the surface associated with the given vertices
  /// @param locVertices a vector of the vertices in local coordinates
  /// @return a vector of the vertices in global coordinates
  std::vector<Acts::Vector3D>
  makeGlobalVertices(const Acts::Surface&               surface,
                     const std::vector<Acts::Vector2D>& locVertices) const;
};

}  // namespace Acts

#endif  // ACTS_TOOLS_SURFACERARRAYCREATOR_H
