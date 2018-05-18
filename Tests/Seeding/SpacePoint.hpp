#pragma once
struct SpacePoint{
  float m_x;
  float m_y;
  float m_z;
  float m_r;
  int surface;
  float covr;
  float covz;
  float x() const {
    return m_x;
  }
  float y() const {
    return m_y;
  }
  float z() const {
    return m_z;
  }
  float r() const {
    return m_r;
  }
};
