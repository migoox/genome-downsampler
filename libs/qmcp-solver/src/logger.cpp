#include "../include/qmcp-solver/logger.hpp"

#include <iostream>

void qmcp::Logger::Info(const std::string& message) {
    std::cout << kInfoPrefix << additional_prefix_ << kDivider << message
              << std::endl;
}

void qmcp::Logger::Debug(const std::string& message) {
    if (!verbose_mode_) return;

    std::cout << kDebugPrefix << additional_prefix_ << kDivider << message
              << std::endl;
}
