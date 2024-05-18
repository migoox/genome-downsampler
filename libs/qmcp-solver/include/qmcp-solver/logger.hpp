#include <string>

namespace qmcp {
class Logger {
   public:
    static constexpr const char* kDebugPrefix = "[Debug]";
    static constexpr const char* kInfoPrefix = "[Info]";
    static constexpr const char* kDivider = ": ";

    explicit Logger(bool verbose_mode,
                    const std::string& additional_prefix = "")
        : verbose_mode_(verbose_mode), additional_prefix_(additional_prefix) {}
    void Info(const std::string& message);
    void Debug(const std::string& message);

   private:
    bool verbose_mode_;
    std::string additional_prefix_;
};
}  // namespace qmcp
