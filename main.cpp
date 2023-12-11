#include <cmath>
#include <ctime>
#include <iostream>
#include <limits>

#include "cfoo.h"

int main() {
    //brent.bracket(0.0, dPbgPB, foo);
    //GetDBStartPoint(par, "DB_300mK.txt");
    const time_t begin = time(nullptr);
    CFoo foo(0);

    foo.SeqFit(3);
    std::clog << '\n' << "Finished." << std::endl;
    const time_t end = time(nullptr);
    const double time_spent = difftime(end, begin);
    std::clog << '\n' << "spent " << time_spent << " seconds" << std::endl;

    std::getchar();
    return 0;
}
