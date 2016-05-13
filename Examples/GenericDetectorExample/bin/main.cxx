#include <memory>

#include "ACTS/Examples/GenericDetectorExample/BuildGenericDetector.h"
#include "ACTS/Detector/TrackingGeometry.h"

int main()
{
  std::unique_ptr<const Acts::TrackingGeometry> geo = Acts::trackingGeometry();
  return 0;
}
