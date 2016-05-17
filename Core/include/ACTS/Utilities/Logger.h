#ifndef ACTS_LOGGER_H
#define ACTS_LOGGER_H 1

// STL include(s)
#include <sstream>
#include <cstdio>
#include <ios>
#include <iomanip>
#include <ctime>
#include <string>
#include <functional>
#include <ostream>

namespace Acts
{
  namespace Logging
  {
    enum Level {VERBOSE = 0, DEBUG, INFO, WARNING, ERROR, FATAL};

    class OutputPolicy
    {
    public:
      virtual ~OutputPolicy() = default;

      virtual void flush(const Level& lvl,const std::ostringstream& input) = 0;
    };

    class PrintPolicy
    {
    public:
      virtual ~PrintPolicy() = default;

      virtual bool doPrint(const Level& lvl) const = 0;
    };

    class OutStream
    {
      typedef std::function<void(const std::ostringstream&)> OutputFunc;
    public:
      OutStream(OutputFunc output):
        m_stream(),
        m_outputFunctor(output)
    {}

      OutStream(const OutStream& copy):
        m_stream(),
        m_outputFunctor(copy.m_outputFunctor)
      {
        m_stream << copy.m_stream.str();
      }

      ~OutStream()
      {
        m_outputFunctor(m_stream);
      }

      template<typename T>
      OutStream& operator<<(T&& input)
      {
        m_stream << std::forward<T>(input);
        return *this;
      }

      template<typename T>
      OutStream& operator<<(T& (*f)(T&))
      {
        f(m_stream);
        return *this;
      }

    private:
      std::ostringstream m_stream;
      OutputFunc  m_outputFunctor;
    };

    class DefaultPrintPolicy: public PrintPolicy
    {
    public:
      DefaultPrintPolicy(const Level& lvl):
        m_level(lvl)
    {}

      bool doPrint(const Level& lvl) const override {return m_level <= lvl;}

    private:
      Level m_level;
    };

    class OutputDecorator: public OutputPolicy
    {
    public:
      OutputDecorator(std::unique_ptr<OutputPolicy> wrappee):
        m_wrappee(std::move(wrappee))
    {}

      void flush(const Level& lvl,const std::ostringstream& input) override
          {
        m_wrappee->flush(lvl,input);
          }

    private:
      std::unique_ptr<OutputPolicy> m_wrappee;
    };

    class NamedOutputDecorator final : public OutputDecorator
    {
    public:
      NamedOutputDecorator(std::unique_ptr<OutputPolicy> wrappee,const std::string& name,unsigned int maxWidth=15):
        OutputDecorator(std::move(wrappee)),
        m_name(name),
        m_maxWidth(maxWidth)
    {}

      void flush(const Level& lvl,const std::ostringstream& input) override
          {
        std::ostringstream os;
        os << std::left << std::setw(m_maxWidth) << m_name.substr(0,m_maxWidth - 3) << input.str();
        OutputDecorator::flush(lvl,os);
          }

    private:
      std::string m_name;
      unsigned int m_maxWidth;
    };

    class TimedOutputDecorator final : public OutputDecorator
    {
    public:
      TimedOutputDecorator(std::unique_ptr<OutputPolicy> wrappee,const std::string& format = "%X"):
        OutputDecorator(std::move(wrappee)),
        m_format(format)
    {}

      void flush(const Level& lvl,const std::ostringstream& input) override
          {
        std::ostringstream os;
        os << std::left << std::setw(12) << now() << input.str();
        OutputDecorator::flush(lvl,os);
          }

    private:
      std::string now() const
      {
        char buffer[20];
        time_t t;
        std::time(&t);
        std::strftime(buffer, sizeof(buffer), m_format.c_str(), localtime(&t));
        return buffer;
      }

      std::string m_format;
    };

    class ThreadOutputDecorator final : public OutputDecorator
    {
    public:
      ThreadOutputDecorator(std::unique_ptr<OutputPolicy> wrappee):
        OutputDecorator(std::move(wrappee))
    {}

      void flush(const Level& lvl,const std::ostringstream& input) override
          {
        std::ostringstream os;
        os << std::left << std::setw(20) << std::this_thread::get_id() << input.str();
        OutputDecorator::flush(lvl,os);
          }
    };

    class LevelOutputDecorator final : public OutputDecorator
    {
    public:
      LevelOutputDecorator(std::unique_ptr<OutputPolicy> wrappee):
        OutputDecorator(std::move(wrappee))
    {}

      void flush(const Level& lvl,const std::ostringstream& input) override
          {
        std::ostringstream os;
        os << std::left << std::setw(10) << toString(lvl) << input.str();
        OutputDecorator::flush(lvl,os);
          }

    private:
      std::string toString(const Level& lvl) const
      {
        static const char* const buffer[] = {"VERBOSE", "DEBUG", "INFO", "WARNING", "ERROR", "FATAL"};
        return buffer[lvl];
      }
    };

    class DefaultOutputPolicy: public OutputPolicy
    {
    public:
      DefaultOutputPolicy(std::FILE* out = stdout):
        m_out(out)
    {}

      void flush(const Level& lvl,const std::ostringstream& input) final
          {
        fprintf(m_out, "%s\n", input.str().c_str());
        fflush(m_out);
          }

    private:
      std::FILE* m_out;
    };
  }  // end of namespace Logging

  class Logger
  {
  public:
    template<typename Output,typename Print>
    Logger(std::unique_ptr<Output> pOutput,std::unique_ptr<Print> pPrint):
    m_outputPolicy(std::move(pOutput)),
    m_printPolicy(std::move(pPrint))
    {}

    bool print(const Logging::Level& lvl) const {return m_printPolicy->doPrint(lvl);}

    Logging::OutStream log(const Logging::Level& lvl) const
    {
      return Logging::OutStream(std::bind(&Logging::OutputPolicy::flush,m_outputPolicy.get(),lvl,std::placeholders::_1));
    }

  private:
    std::unique_ptr<Logging::OutputPolicy> m_outputPolicy;
    std::unique_ptr<Logging::PrintPolicy> m_printPolicy;
  };
}  // end of namespace Acts

#endif // ACTS_LOGGER_H
