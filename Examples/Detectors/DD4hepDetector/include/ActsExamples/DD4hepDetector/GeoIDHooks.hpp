#pragma once

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Surfaces/Surface.hpp"

/// Geometry identifiers hooks to be used with DD4Hep detectors
/// to had some extra identifier to sensitives surfaces.
namespace det {
namespace geoIDHook {

/// Use the extra identifier in the endcap to separate the two row of modules
/// @param identifier geometry identifier
/// @param surface coresponding surface
Acts::GeometryIdentifier stripEndcapODD(Acts::GeometryIdentifier identifier,
                                        const Acts::Surface& surface);
}  // namespace geoIDHook
}  // namespace det
