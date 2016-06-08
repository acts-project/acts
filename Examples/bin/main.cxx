#include <memory>

#include "ACTS/Examples/BuildGenericDetector.hpp"
#include "ACTS/Detector/TrackingGeometry.hpp"
#include "ACTS/Utilities/Logger.hpp"

int main()
{
  std::unique_ptr<const Acts::TrackingGeometry> geo = Acts::trackingGeometry(Acts::Logging::DEBUG);
  return 0;
}
