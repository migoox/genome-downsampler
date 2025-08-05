#ifndef HELPERS_HPP
#define HELPERS_HPP

#include <map>
#include <string>
#include <vector>

namespace helpers {

template <typename T>
std::vector<std::string> get_names_from_map(const std::map<std::string, T>& input_map) {
    std::vector<std::string> names;
    names.reserve(input_map.size());

    for (const auto& mapping : input_map) {
        names.push_back(mapping.first);
    }

    return names;
}

}  // namespace helpers

#endif
