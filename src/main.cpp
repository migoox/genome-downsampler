#include <cstdlib>

#include "app.hpp"

int main(int argc, char** argv) {
    App app;

    try {
        app.Parse(argc, argv);
    } catch (const CLI::ParseError& e) {
        return app.Exit(e);
    }

    app.Solve();

    return EXIT_SUCCESS;
}
