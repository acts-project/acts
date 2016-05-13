#ifndef ACTS_LOGGER_H
#define ACTS_LOGGER_H 1

// STL include(s)
#include <type_traits>
#include <sstream>
#include <cstdio>
#include <ios>
#include <iomanip>
#include <ctime>

namespace Acts
{
  template<typename Stream>
  class Logger
  {
  public:
    virtual ~Logger() = default;

    virtual bool printVerbose() const = 0;
    virtual bool printDebug() const = 0;
    virtual bool printInfo() const = 0;
    virtual bool printWarning() const = 0;
    virtual bool printError() const = 0;
    virtual bool printFatal() const = 0;

    virtual Stream verbose() const = 0;
    virtual Stream debug() const = 0;
    virtual Stream info() const = 0;
    virtual Stream warning() const = 0;
    virtual Stream error() const = 0;
    virtual Stream fatal() const = 0;

    typedef typename std::remove_reference<Stream>::type& StreamRef;
    typedef StreamRef (&EndMsgFunc)(StreamRef);
    virtual EndMsgFunc endmsg() const = 0;
  };

  enum Level {VERBOSE = 0, DEBUG, INFO, WARNING, ERROR, FATAL};

  class LoggerImplementation
  {
  public:
    LoggerImplementation(const Level& lvl):
      m_stream()
    {
      m_stream << std::left << std::setw(10) << now() << std::left << std::setw(10) << toString(lvl);
    }

    LoggerImplementation(const LoggerImplementation& copy):
      m_stream()
    {
      m_stream << copy.m_stream.str();
    }

    ~LoggerImplementation()
    {
      fprintf(stdout, "%s", m_stream.str().c_str());
      fflush(stdout);
    }

    template<typename T>
    LoggerImplementation& operator<<(T&& input)
    {
      m_stream << std::forward<T>(input);

      return *this;
    }

    template<typename T>
    LoggerImplementation& operator<<(T& (*f)(T&))
    {
      f(m_stream);

      return *this;
    }

    LoggerImplementation& operator<<(LoggerImplementation& (*f)(LoggerImplementation&))
    {
      f(*this);

      return *this;
    }

    friend LoggerImplementation& EndMsg(LoggerImplementation& s);

  private:
    std::string toString(const Level& lvl) const
    {
      static const char* const buffer[] = {"VERBOSE", "DEBUG", "INFO", "WARNING", "ERROR", "FATAL"};
      return buffer[lvl];
    }

    std::string now() const
    {
      char buffer[20];
      time_t t;
      std::time(&t);
      std::strftime(buffer, sizeof(buffer), "%X", localtime(&t));
      return buffer;
    }

    mutable std::ostringstream m_stream;
  };

  LoggerImplementation& EndMsg(LoggerImplementation& s)
  {
    s.m_stream << std::endl;
    return s;
  }

  class DefaultLogger: public Logger<LoggerImplementation>
  {
  public:
    DefaultLogger(const Level& lvl):
      m_level(lvl)
    {}

    bool printVerbose() const override {return m_level <= VERBOSE;}
    bool printDebug() const override {return m_level <= DEBUG;}
    bool printInfo() const override {return m_level <= INFO;}
    bool printWarning() const override {return m_level <= WARNING;}
    bool printError() const override {return m_level <= ERROR;}
    bool printFatal() const override {return m_level <= FATAL;}

    LoggerImplementation verbose() const override {return LoggerImplementation(VERBOSE);}
    LoggerImplementation debug() const override   {return LoggerImplementation(DEBUG);}
    LoggerImplementation info() const override    {return LoggerImplementation(INFO);}
    LoggerImplementation warning() const override {return LoggerImplementation(WARNING);}
    LoggerImplementation error() const override   {return LoggerImplementation(ERROR);}
    LoggerImplementation fatal() const override   {return LoggerImplementation(FATAL);}

    typedef LoggerImplementation& (&EndMsgFunc)(LoggerImplementation&);
    EndMsgFunc endmsg() const override
    {return EndMsg;}

  private:
    Level        m_level;
  };
}  // end of namespace Acts

#endif // ACTS_LOGGER_H
