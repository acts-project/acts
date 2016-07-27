#include "ACTS/Plugins/MaterialPlugins/MaterialTrackRecord.hpp"

Acts::MaterialTrackRecord::MaterialTrackRecord(
    const MaterialStep::Position& startPos,
    double                        eta,
    double                        phi,
    std::vector<MaterialStep>     materialSteps)
  : m_startPosition(startPos)
  , m_eta(eta)
  , m_phi(phi)
  , m_materialSteps(materialSteps)
{
}

Acts::MaterialTrackRecord::MaterialTrackRecord(
    const MaterialTrackRecord& mtrecord)
  : m_startPosition(mtrecord.m_startPosition)
  , m_eta(mtrecord.m_eta)
  , m_phi(mtrecord.m_phi)
  , m_materialSteps(mtrecord.m_materialSteps)
{
}

const Acts::MaterialTrackRecord*
Acts::MaterialTrackRecord::clone() const
{
  return (new MaterialTrackRecord(*this));
}

Acts::MaterialTrackRecord&
Acts::MaterialTrackRecord::operator=(const MaterialTrackRecord& mtrecord)
{
  if (this != &mtrecord) {
    m_startPosition = mtrecord.m_startPosition;
    m_eta           = mtrecord.m_eta;
    m_phi           = mtrecord.m_phi;
    m_materialSteps = mtrecord.m_materialSteps;
  }
  return (*this);
}