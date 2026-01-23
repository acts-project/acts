


#include "Acts/Geometry/detail/PortalPlacement.hpp"
#include "Acts/Geometry/VolumePlacementBase.hpp"

namespace Acts::detail {
    PortalPlacement::PortalPlacement(
    const std::size_t portalIdx, const Transform3& portalTrf, 
    VolumePlacementBase* parent,
    std::shared_ptr<RegularSurface>&& surface)
    : m_interalTrf{portalTrf},
      m_surface{std::move(surface)},
      m_parent{parent},
      m_portalIdx{portalIdx} {
  assert(m_surface != nullptr);
  m_surface->assignDetectorElement(*this);
}

const Transform3& PortalPlacement::localToGlobalTransform(
    const GeometryContext& gctx) const {
  if (!m_parent->portalTransformCached(m_portalIdx)) {
    m_parent->cachePortalTransform(
        m_portalIdx, m_parent->localToGlobalTransform(gctx) * m_interalTrf);
  }
  return m_parent->portalLocalToGlobal(gctx, m_portalIdx);
}

const Surface& PortalPlacement::surface() const {
  return *m_surface;
}

Surface& PortalPlacement::surface() {
  return *m_surface;
}

double PortalPlacement::thickness() const {
  return 0.;
}

bool PortalPlacement::isSensitive() const {
  return false;
}

std::size_t PortalPlacement::index() const {
  return m_portalIdx;
}

const Transform3& PortalPlacement::portalToVolumeCenter()
    const {
  return m_interalTrf;
}
}