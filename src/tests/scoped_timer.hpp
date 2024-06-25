#ifndef SCOPED_TIMER
#define SCOPED_TIMER
#include <chrono>

#include "logging/log.hpp"
class ScopedTimer {
   public:
    ScopedTimer() : start_(std::chrono::high_resolution_clock::now()) {}
    ~ScopedTimer() {
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start_;
        LOG_WITH_LEVEL(logging::INFO) << "Took: " << duration.count() << " seconds.";
    }

   private:
    std::chrono::time_point<std::chrono::high_resolution_clock> start_;
};

#endif