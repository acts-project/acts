#include <cmath>
#include <concepts>
#include <vector>

namespace Acts {

// Helper class describing one section of the accumulator space
class HoughAccumulatorSection {
 public:
  HoughAccumulatorSection() = default;

  HoughAccumulatorSection(float xw, float yw, float xBegin, float yBegin,
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

  /// create section that is bottom part this this one
  /// @arg copyIndices - copies indices from the parent
  /// +------+
  /// |      |
  /// +------+
  /// |     <|-- this part
  /// +------+
  HoughAccumulatorSection bottom(bool copyIndices = false) const;
  /// see @bottom
  HoughAccumulatorSection top(bool copyIndices = false) const;
  /// @see @bottom
  HoughAccumulatorSection left(bool copyIndices = false) const;
  /// @see @bottom
  HoughAccumulatorSection right(bool copyIndices = false) const;

  /// create section that is bottom left corner of this this one
  /// by default the section is divided into 4 quadrants,
  /// if parameters are provided the quadrants size can be adjusted
  /// +---+---+
  /// |   |   |
  /// +---+---+
  /// |   |  <|-- this part
  /// +---+---+
  HoughAccumulatorSection bottomRight(bool copyIndices = false) const;
  HoughAccumulatorSection bottomLeft(bool copyIndices = false) const;

  HoughAccumulatorSection topLeft(bool copyIndices = false) const;
  HoughAccumulatorSection topRight(bool copyIndices = false) const;

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
  template <typename F>
  inline bool isCrossingInside(F &&line1, F &&line2) const &
    requires std::invocable<F, float>;

  /// sizes access
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

template <typename F>
inline bool HoughAccumulatorSection::isCrossingInside(F &&line1, F &&line2) const &
  requires std::invocable<F, float>
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
  // The above covers most of the cases. 
  // Additional precautions are made when both lines cross
  // left & right (x) bounds outside of vertical (y) bounds.


  const float xL = xBegin();
  const float xR = xBegin() + xSize();

  // Evaluate both lines at section boundaries
  const float y1L = line1(xL);
  const float y1R = line1(xR);
  const float y2L = line2(xL);
  const float y2R = line2(xR);

  // --- Step 1: Quick rejection based on ordering ---
  // If one line is consistently above the other → no crossing
  const float dL = y1L - y2L;
  const float dR = y1R - y2R;

  if (dL * dR > 0) {
    return false;
  }

  // --- Step 2: Check if either line fully spans inside vertically ---
  auto is_in = [this](float y) {
    return yBegin() < y && y < yBegin() + ySize();
  };

  if (is_in(y1L) && is_in(y1R)) {
    return true;
  }

  if (is_in(y2L) && is_in(y2R)) {
    return true;
  }

  // --- Step 3: Approximate crossing position ---
  // Linear interpolation assuming near-linear behavior
  const float abs_dL = std::abs(dL);
  const float abs_dR = std::abs(dR);

  // Avoid division by zero (parallel & overlapping edge case)
  if (abs_dL + abs_dR == 0.0f) {
    return false;
  }

  const float t = abs_dL / (abs_dL + abs_dR);  // fraction from left
  const float xCross = xL + t * (xR - xL);

  // --- Step 4: Check if crossing lies inside vertical bounds ---
  const float yCross1 = line1(xCross);
  const float yCross2 = line2(xCross);
  const float yCross = 0.5f * (yCross1 + yCross2);  // intersection approx

  return is_in(yCross);

  // const float line1_left_y = line1(xBegin());
  // const float line1_right_y = line1(xBegin() + xSize());
  // const float line2_left_y = line2(xBegin());
  // const float line2_right_y = line2(xBegin() + xSize());
  // const float left_mid = 0.5f * (line1_left_y + line2_left_y);
  // const float right_mid = 0.5f * (line1_right_y + line2_right_y);
  // // they do not cross for sure
  // if ((line1_left_y - left_mid) * (line1_right_y - right_mid) > 0) {
  //   return false;
  // }
  // //
  // auto is_in = [this](float v) {
  //   return yBegin() < v and v < yBegin() + ySize();
  // };

  // if (is_in(line1_left_y) and is_in(line1_right_y)) {
  //   return true;
  // }

  // if (is_in(line2_left_y) and is_in(line2_right_y)) {
  //   return true;
  // }

  // // check where they cross assuming thy are lines
  // // assume that the functions are close to linear
  // const float d_y_left = std::abs(line1_left_y - line2_left_y);
  // const float d_y_right = std::abs(line1_right_y - line2_right_y);
  // const float x_cross = xSize() * d_y_right / (d_y_left + d_y_right);
  // // if any of the lines at point xBegin()+x_cross is in box we consider
  // the
  // // crossing to be close

  // if (is_in(line1(xBegin() + x_cross))) {
  //   return true;
  // }
  // if (is_in(line2(xBegin() + x_cross))) {
  //   return true;
  // }
  // return false;
}

template <typename M, typename Options, typename Functor>
void exploreParametersSpace(std::vector<HoughAccumulatorSection> &sectionsStack,
                            const std::vector<M> &measurements,
                            const Options &opt,
                            Functor lineFunctor,
                            std::vector<HoughAccumulatorSection> &results) {
    
    using Decision = typename Options::Decision;
    
    while (!sectionsStack.empty()) {
        //ACTS_VERBOSE("Stack size " << sectionsStack.size());
        
        HoughAccumulatorSection thisSection = std::move(sectionsStack.back());
        sectionsStack.pop_back();

        Decision whatNext = opt.decisionFunctor(thisSection, measurements);

        if (whatNext == Decision::Discard) {
            continue;
        } else if (whatNext == Decision::Accept) {
            results.push_back(std::move(thisSection));
        } else {
            bool splitX = thisSection.xSize() > opt.xMinBinSize;
            bool splitY = thisSection.ySize() > opt.yMinBinSize;

            if (splitX && splitY) {
                // Split into 4 sections
                HoughAccumulatorSection bL = thisSection.bottomLeft();
                HoughAccumulatorSection tL = thisSection.topLeft();
                HoughAccumulatorSection bR = thisSection.bottomRight();
                HoughAccumulatorSection tR = thisSection.topRight();

                if (whatNext == Decision::DrillAndExpand) {
                    bL.expand(opt.expandX, opt.expandY);
                    tL.expand(opt.expandX, opt.expandY);
                    bR.expand(opt.expandX, opt.expandY);
                    tR.expand(opt.expandX, opt.expandY);
                }
                
                // checking if lines are crossing 4 new sections 
                for (unsigned idx : thisSection.indices()) {
                    const auto &m = measurements[idx];
                    auto line = [&](float x){return lineFunctor(m,x);};

                    if (bL.isLineInside(line)) bL.indices().push_back(idx);
                    if (tL.isLineInside(line)) tL.indices().push_back(idx);
                    if (bR.isLineInside(line)) bR.indices().push_back(idx);
                    if (tR.isLineInside(line)) tR.indices().push_back(idx);
                }
                sectionsStack.push_back(std::move(bL));
                sectionsStack.push_back(std::move(tL));
                sectionsStack.push_back(std::move(bR));
                sectionsStack.push_back(std::move(tR));

            } else if (!splitX && splitY) {
                // Split into 2 horizontal sections
                HoughAccumulatorSection b = thisSection.bottom();
                HoughAccumulatorSection t = thisSection.top();

                if (whatNext == Decision::DrillAndExpand) {
                    b.expand(opt.expandX, opt.expandY);
                    t.expand(opt.expandX, opt.expandY);
                }

                for (unsigned idx : thisSection.indices()) {
                    const auto &m = measurements[idx];
                    auto line = [&](float x){return lineFunctor(m,x);};

                    if (b.isLineInside(line)) b.indices().push_back(idx);
                    if (t.isLineInside(line)) t.indices().push_back(idx);
                }
                sectionsStack.push_back(std::move(b));
                sectionsStack.push_back(std::move(t));

            } else if (splitX && !splitY) {
                // Split into 2 vertical sections
                HoughAccumulatorSection l = thisSection.left();
                HoughAccumulatorSection r = thisSection.right();

                if (whatNext == Decision::DrillAndExpand) {
                    l.expand(opt.expandX, opt.expandY);
                    r.expand(opt.expandX, opt.expandY);
                }

                for (unsigned idx : thisSection.indices()) {
                    const auto &m = measurements[idx];
                    auto line = [&](float x){return lineFunctor(m,x);};

                    if (l.isLineInside(line)) l.indices().push_back(idx);
                    if (r.isLineInside(line)) r.indices().push_back(idx);
                }
                sectionsStack.push_back(std::move(l));
                sectionsStack.push_back(std::move(r));
            }
        }
    }
}

template <typename SectionType, typename MeasurementType, typename Functor>
bool passIntersectionsCheck(
    const SectionType &section,
    const std::vector<MeasurementType> &measurements,
    Functor lineFunctor, 
    const unsigned threshold) {
    
    const unsigned count = section.count();
    const float xLeft = section.xBegin();
    const float xRight = xLeft + section.xSize();
    const auto& indices = section.indices();

    // Small Buffer Optimization
    constexpr unsigned kMaxStackLines = 64; 
    
    if (count <= kMaxStackLines) {
        float yLeft[kMaxStackLines];
        float yRight[kMaxStackLines];        
        
        for (unsigned i = 0; i < count; ++i) {
            const auto &m = measurements[indices[i]];
            yLeft[i] = lineFunctor(m, xLeft);
            yRight[i] = lineFunctor(m, xRight);
        }
        
        unsigned inside = 0;
        for (unsigned i = 0; i < count; ++i) {
            for (unsigned j = i + 1; j < count; ++j) {
                if ((yLeft[i] - yLeft[j]) * (yRight[i] - yRight[j]) < 0.0f) {
                    inside++;
                    if (inside >= threshold) return true; // Early exit
                }
            }
        }
        return false;
    } 
    else {
        std::vector<float> yLeft(count);
        std::vector<float> yRight(count);
        for (unsigned i = 0; i < count; ++i) {
            const auto &m = measurements[indices[i]];
            yLeft[i] = lineFunctor(m, xLeft);
            yRight[i] = lineFunctor(m, xRight);
        }
        unsigned inside = 0;
        for (unsigned i = 0; i < count; ++i) {
            for (unsigned j = i + 1; j < count; ++j) {
                if ((yLeft[i] - yLeft[j]) * (yRight[i] - yRight[j]) < 0.0f) {
                    inside++;
                    if (inside >= threshold) return true; // Early exit
                }
            }
        }
        return inside >= threshold;
    }
}

}  // namespace Acts