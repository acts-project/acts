#ifndef ACTS_LOGFILE_ERASURE_FIXTURE_HPP
#define ACTS_LOGFILE_ERASURE_FIXTURE_HPP 1

#include <fstream>
#include <string>

namespace Acts {

namespace Test {

  struct logfile_erasure_fixture
  {
    logfile_erasure_fixture(const std::string& name)
    {
      std::ofstream f(name);
      f.close();
    }
  };

}  // namespace Test

}  // namespace Acts
#endif  // ACTS_LOGFILE_ERASURE_FIXTURE_HPP
