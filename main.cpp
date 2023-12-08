#include "classfuncs.h"

int main()
{
    //brent.bracket(0.0, dPbgPB, foo);
    //GetDBStartPoint(par, "DB_300mK.txt");
    clock_t begin = clock();
    CFoo foo(0);

    foo.SeqFit(3);
    printf("\nFinished.");
    clock_t end = clock();
    double time_spent = (double) (end - begin) / CLOCKS_PER_SEC;
    printf("\n%f", time_spent);

    double tmp = 0;
    scanf("%lf", &tmp);
    return 0;
}
