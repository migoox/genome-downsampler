#pragma once

#include <sstream>
#include <string>

#define LOG(level) logging::Log().get(level)
#define SET_LOG_LEVEL(level) logging::Log::ReportingLevel = level

namespace logging {
enum LogLevel { kError, kInfo, kDebug };

class Log {
   public:
    Log() = default;
    virtual ~Log();
    std::ostringstream& get(LogLevel level = kInfo);
    static inline LogLevel ReportingLevel = kInfo;

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
