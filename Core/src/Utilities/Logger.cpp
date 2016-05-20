#include "ACTS/Utilities/Logger.hpp"

namespace Acts
{
  std::unique_ptr<Logger> getDefaultLogger(const std::string& name,const Logging::Level& lvl)
  {
    using namespace Logging;
    auto output= std::make_unique<LevelOutputDecorator>(std::make_unique<NamedOutputDecorator>(std::make_unique<TimedOutputDecorator>(std::make_unique<DefaultOutputPolicy>(stdout)),name));
    auto print= std::make_unique<DefaultPrintPolicy>(lvl);
    return std::make_unique<Logger>(std::move(output),std::move(print));
  }
}  // end of namespace Acts
