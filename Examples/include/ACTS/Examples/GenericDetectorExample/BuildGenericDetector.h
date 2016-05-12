#ifndef ACTS_BUILDGENERICDETECTOR_H
#define ACTS_BUILDGENERICDETECTOR_H 1

// STL include(s)
#include <memory>

// ACTS include(s)
namespace Acts
{
  class TrackingGeometry;

  std::unique_ptr<const Acts::TrackingGeometry> trackingGeometry();
}

#endif // ACTS_BUILDGENERICDETECTOR_H
