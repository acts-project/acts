#pragma once

#include <array>
#include <cstdint>

namespace Acts {

// class of an space point for the purpose of covariance scaling
enum class CosmicSpacePointClass : std::uint8_t {
  BarrelPixel = 0,
  EndcapPixel = 1,
  BarrelStrip = 2,
  EndcapStrip = 3,
};

// cuts aplied by CosmicTripletSeedFinder when forming top-middle-bottom triplets from cosmic ray space points
struct CosmicSeedCuts {
  // order: top.y() > middle.y() > bottom.y()
  bool requireTopMiddleBottom = true;
  // with useSubdetectorAngularCuts is the maximum kink for strip SPs
  // and pixel SPs are cut by pixelMinCosDeltaPhi instead
  double maxKinkAngle = 0.05;
  // vertical separation for each segment [mm], for any middle SP (pixel or strip):
  // |y_top - y_middle| and |y_middle - y_bottom| in [minDeltaY, maxDeltaY]
  // these are mindRadius/maxdRadius of SiSpacePointsSeedMaker_Cosmic (m_drmin = 5, m_drmax = 1000),
  // where dR is Δy
  double minDeltaY = 5.0;
  double maxDeltaY = 1000.0;
  // minimal alignment with the y direction: |dy| / |dr|, where dr = sqrt(dx^2 + dy^2 + dz^2)
  double minCosThetaY = 0.3;
  // maximal distance of the middle SP to the line top-bottom [mm]
  // with usePtDependentKinkResidual it is the floor of the tolerance of the
  // component out of the transverse plane (see passCosmicCuts)
  // negative = no cut
  double maxResidual = 5.0;
  // zMin/zMax [mm] of the space points to be considered for seeding
  double zMin = -3000.0;
  double zMax = 3000.0;
  // aceptance in |eta| < maxAbsEta
  double maxAbsEta = 2.3;
  // momentum cut: pT > minPt [MeV] (requires magneticFieldBz != 0)
  double minPt = 300.0;
  bool usePtCut = false;
  double magneticFieldBz = 2.0;  // T
  double multipleScatteringCof = 60.3;
  // if is true, rejects triplets whose longitudinal discriminant dT is positive,
  // i.e. incompatible with a straight (or curved) trajectory within the errors
  bool useLongitudinalCompatibility = false;
  bool usePtDependentKinkResidual = false;
  // maximum pT for scattering: if the estimated pT is above this value,
  // it is replaced by this value in the scattering term (same as in Athena)
  double maxPtForScattering = 5000.0;

  // active by default, as in Athena, but can be disabled for testing
  bool useCovarianceScale = true;
  std::array<double, 4> covarianceScale = {16., 16., 100., 16.};

  bool requireSameSubdetector = false;

  // if is true, the top-middle and middle-bottom SPs cannot be in the same surface (module)
  // this requires passing the surface of each SP to createSeeds()
  bool rejectSameSurface = true;

  // if is true, the cuts on the angular separation of the segments are applied depending
  // on the subdetector type (pixel or strip)
  bool useSubdetectorAngularCuts = false;
  double pixelMinCosDeltaPhi = 0.96;

  // if is true, the positions of strip space points in the seeds are recalculated with the direction of their best seed
  // this requires passing a copy of the space points with the same indices to createSeeds()
  bool useStripCalibration = false;
  double stripCalibrationTolerance = 1.1;

  // quality assigned to the seeds, used for sorting and filtering in the CKF
  bool useSeedQuality = true;
  float seedQuality = 0.0f;
};

}
