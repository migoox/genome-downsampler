#ifndef LOG_HPP
#define LOG_HPP

#include <sstream>
#include <string>

#define LOG_WITH_LEVEL(level) logging::Log().get(level)
#define SET_LOG_LEVEL(level) logging::Log::ReportingLevel = level

namespace logging {
enum LogLevel { ERROR, INFO, DEBUG };

class Log {
   public:
    Log() = default;
    virtual ~Log();
    std::ostringstream& get(LogLevel level = INFO);
    static inline LogLevel ReportingLevel = INFO;

   private:
    Log(const Log&);
    Log& operator=(const Log&);
    LogLevel messageLevel_;
    std::ostringstream os_;
    std::string logLevelPrefixes_[3] = {
        "Error",
        "Info",
        "Debug",
    };
};
}  // namespace logging

#endif
