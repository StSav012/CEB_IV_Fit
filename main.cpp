#include <cmath>
#include <chrono>
#include <iostream>

#include "constants.h"

#include "cfoo.h"

int main() {
    const std::chrono::time_point<std::chrono::steady_clock> beginning = std::chrono::steady_clock::now();
    try {
        CFoo foo(0);
        foo.loadExperimentData(fname);
        foo.CEB_2eq_parallel_lite();

        auto [Irex, Vrex] = foo.resample();

        foo.SeqFit(1, Irex);
    } catch (std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }
    std::clog << '\n' << "Finished." << std::endl;
    std::clog << "Spent " << std::chrono::duration<double>(std::chrono::steady_clock::now() - beginning) << std::endl;

    // std::getchar();
    return 0;
}
