
#include "Acts/Seeding/GbtsLutParser.hpp"
#include <fstream>

namespace Acts::Experimental{

   void  GbtsLutParser::parseLutFile(std::string& lutInputFile, bool useML){

    if (useML) {
      if (lutInputFile.empty()) {
        throw std::runtime_error("Cannot find ML predictor LUT file");
      } else {
        m_mlLUT.reserve(100);
        std::ifstream ifs(std::string(lutInputFile).c_str());

        if (!ifs.is_open()) {
          throw std::runtime_error("Failed to open LUT file");
        }

        float cl_width{}, min1{}, max1{}, min2{}, max2{};

        while (ifs >> cl_width >> min1 >> max1 >> min2 >> max2) {
          std::array<float, 5> lut_line = {cl_width, min1, max1, min2, max2};
          m_mlLUT.emplace_back(lut_line);
        }
        if (!ifs.eof()) {
          // ended if parse error present, not clean EOF

          throw std::runtime_error(
              "Stopped reading LUT file due to parse error");
        }

        ifs.close();
      }
    }
   }
}