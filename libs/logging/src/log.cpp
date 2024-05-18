#include "../include/logging/log.hpp"

#include <iostream>
#include <ostream>

logging::Log::~Log() {
    if (messageLevel_ <= ReportingLevel) {
        os_ << std::endl;
        std::cerr << os_.str();
   }
}

std::ostringstream& logging::Log::Get(LogLevel level) {
    os_ << logLevelPrefixes_[level] << ": ";
    messageLevel_ = level;
    return os_;
}
