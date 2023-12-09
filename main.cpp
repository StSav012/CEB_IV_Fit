#include <cmath>
#include <ctime>
#include <iostream>
#include <limits>

#include "cfoo.h"

int main()
{
    //brent.bracket(0.0, dPbgPB, foo);
    //GetDBStartPoint(par, "DB_300mK.txt");
    time_t begin = time(NULL);
    CFoo foo(0);

    foo.SeqFit(3);
    std::clog << '\n' << "Finished." << std::endl;
    time_t end = time(NULL);
    double time_spent = difftime(end, begin);
    std::clog << '\n' << "spent " << time_spent << " seconds" << std::endl;

    std::getchar();
    return 0;
}
