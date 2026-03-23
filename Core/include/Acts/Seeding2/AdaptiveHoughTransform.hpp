#include <vector>
#include <concepts>

namespace Acts {

// Helper class describing one section of the accumulator space
class AccumulatorSection {
 public:
  AccumulatorSection() = default;

  AccumulatorSection(float xw, float yw, float xBegin, float yBegin,
                     int div = 0, const std::vector<unsigned> &indices = {},
                     const std::vector<float> &history = {});

  /// @brief keep indices and update parameters of the box
  /// This method is useful when changing direction of the search
  void updateDimensions(float xw, float yw, float xBegin, float yBegin);

  /// @brief keep indices and update parameters of the box by scalling
  /// @param xs - scale in x direction, if bigger than 1 the size increases
  /// @param ys - scale in y direction
  /// The box is recentred
  void expand(float xs, float ys);

  inline unsigned int count() const {
    return static_cast<unsigned>(m_indices.size());
  }
  inline const std::vector<unsigned> &indices() const { return m_indices; }
  inline std::vector<unsigned> &indices() { return m_indices; }

  // create section that is bottom part this this one
  // +------+
  // |      |
  // +------+
  // |     <|-- this part
  // +------+
  AccumulatorSection bottom(float yFraction = 0.5) const;
  // see @bottom
  AccumulatorSection top(float yFraction = 0.5) const;
  // @see @bottom
  AccumulatorSection left(float xFraction = 0.5) const;
  // @see @bottom
  AccumulatorSection right(float xFraction = 0.5) const;

  // create section that is bottom left corner of this this one
  // by default the section is divided into 4 quadrants,
  // if parameters are provided the quadrants size can be adjusted
  // +---+---+
  // |   |   |
  // +---+---+
  // |   |  <|-- this part
  // +---+---+
  AccumulatorSection bottomRight(float xFraction = 0.5,
                                 float yFraction = 0.5) const;
  AccumulatorSection bottomLeft(float xFraction = 0.5,
                                float yFraction = 0.5) const;

  AccumulatorSection topLeft(float xFraction = 0.5,
                             float yFraction = 0.5) const;
  AccumulatorSection topRight(float xFraction = 0.5,
                              float yFraction = 0.5) const;

  /// @brief true if the line defined by given parameters passes the section
  /// @param function is callable used to check crossing at the edges
  template <typename F>
  inline bool isLineInside(F &&function) const &
    requires std::invocable<F, float>
  {
    const float yB = function(m_xBegin);
    const float yE = function(m_xBegin + m_xSize);
    return (yE > yB) ? yB < m_yBegin + m_ySize && yE > m_yBegin
                     : yB > m_yBegin && yE < m_yBegin + m_ySize;
  }

  /// @brief check if the lines cross inside the section
  /// @brief line1 - functional form of line 1
  /// @brief line1 - functional form of line 2
  /// @warning note that this function is assuming that these are lines and the derivative is positive.
  /// It may be incorrect assumption for rapidly changing function or large
  /// sections
  /// @return true if the lines cross in the section
  template <typename F1, typename F2>
  inline bool isCrossingInside(F1 &&line1, F2 &&line2) const &
    requires std::invocable<F1, float> and std::invocable<F2, float>
  {
    // this microalgorithm idea is illustrated below
    // section left section right
    // example with crossing
    //                                       |            +2
    // line 1 crossing left section edge     +1          _|
    // left edge mid point                   |_           |
    //                                       |            +1
    // line 2crossing left section           +2           |
    //
    // example with no crossing
    //                                       |            +1
    // line 1 crossing left section edge     +1          _|
    // left edge mid point                   |_           |
    //                                       |            +2
    // line 2crossing left section           +2           |
    //


    // if for any of the two lines the condition
    // (line1_left_y-middle_on_the_left_y)*(line1_right_y-middle_on_the_right_y)
    // < 0 means that there is crossing

    const float line1_left_y = line1(xBegin());
    const float line1_right_y = line1(xBegin() + xSize());
    const float line2_left_y = line2(xBegin());
    const float line2_right_y = line2(xBegin() + xSize());
    const float left_mid = 0.5f * (line1_left_y + line2_left_y);
    const float right_mid = 0.5f * (line1_right_y + line2_right_y);
    return (line1_left_y - left_mid) * (line1_right_y - right_mid) < 0 
      and left_mid ;
  }

  // sizes
  float xSize() const { return m_xSize; }
  float ySize() const { return m_ySize; }
  float xBegin() const { return m_xBegin; }
  float yBegin() const { return m_yBegin; }
  unsigned divisionLevel() const { return m_divisionLevel; }

  /// store additional (arbitrary) info in indexed array
  /// @param index - identifier
  /// @param value - value to store
  void setHistory(unsigned index, float value) {
    m_history.resize(index + 1);
    m_history.at(index) = value;
  }
  /// @brief retrieve history info
  /// @param index - item index
  /// @return value stored by @see setHistory
  float history(unsigned index) const { return m_history[index]; }

 private:
  float m_xSize = 0;
  float m_ySize = 0;
  float m_xBegin = 0;
  float m_yBegin = 0;
  unsigned m_divisionLevel =
      0;  // number of times the starting section was already divided
  std::vector<unsigned>
      m_indices;  // indices of measurements contributing to this section
  std::vector<float> m_history;  // additional record where an arbitrary
                                 // information can be stored
};
}