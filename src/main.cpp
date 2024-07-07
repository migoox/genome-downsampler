#include <cstdlib>

#include "app.hpp"

int main(int argc, char** argv) {
    App app;

    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError& e) {
        return app.exit(e);
    }

    return EXIT_SUCCESS;
}
