#include <memory>

#include "ACTS/Examples/GenericDetectorExample/BuildGenericDetector.hpp"
#include "ACTS/Detector/TrackingGeometry.hpp"

int main()
{
  std::unique_ptr<const Acts::TrackingGeometry> geo = Acts::trackingGeometry();
  return 0;
}
