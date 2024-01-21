//
// Created by Alon Levi on 21/01/2024.
//

namespace Acts {
// Make Mock Detector
class MockDetector {
    public:
        MockTrack(const Vector3 &mom, const Vector3 &pos) : m_mom(mom), m_pos(pos) {
            // nop
        }

        Vector3 momentum() const { return m_mom; }

        Vector3 position() const { return m_pos; }

    private:
        Vector3 m_mom;
        Vector3 m_pos;
    };


// Make Mock Tracking Geometry

} // namespace acts