#include <memory>

#include "ACTS/Seeding/SeedSPGrid.hpp"
#include "ACTS/Utilities/detail/Axis.hpp"


std::unique_ptr<Acts::Seeding::SPGrid> Acts::Seeding::SPGridCreator::createGrid(std::shared_ptr<Acts::Seeding::Config> config)
{
  //calculate circle intersections of helix and max detector radius
  float minHelixRadius = config->minPt/(300.*config->bFieldInZ); // in mm
  float maxR2 = config->rMax*config->rMax; 
  float xOuter = maxR2/(2*minHelixRadius);
  float yOuter = sqrt(maxR2-xOuter*xOuter);
  float outerAngle = atan(xOuter/yOuter);
  // intersection of helix and max detector radius minus maximum R distance per seed
  float innerCircleR2 = (config->rMax - config->deltaRMax) * (config->rMax - config->deltaRMax);
  float xInner = innerCircleR2/(2*minHelixRadius);
  float yInner = std::sqrt(innerCircleR2-xInner*xInner);
  float innerAngle = std::atan(xInner/yInner);

  // divide 2pi by angle delta to get number of phi-bins
  // size is always 2pi even for regions of interest
  int phiBins = std::ceil(2*M_PI/(outerAngle-innerAngle));
  Acts::detail::Axis<detail::AxisType::Equidistant, detail::AxisBoundaryType::Closed> phiAxis(-M_PI, M_PI, phiBins);

  // TODO: can probably be optimized using smaller z bins 
  // and returning (multiple) neighbors only in one z-direction for forward seeds
  // FIXME: zBinSize must include scattering
  float zBinSize = config->cotThetaMax * config->deltaRMax;
  int zBins = std::ceil((config->zMax - config->zMin)/zBinSize);
  detail::Axis<detail::AxisType::Equidistant, detail::AxisBoundaryType::Bound> zAxis(config->zMin, config->zMax, zBins);

  Acts::Seeding::SPGrid grid(std::make_tuple(phiAxis, zAxis));
  return std::make_unique<Acts::Seeding::SPGrid>(grid);
}

