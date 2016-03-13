#include "gapmanager.h"

GapManager::GapManager(unsigned gaps, unsigned intervals, double overlap):
    gaps(gaps),
    intervals(intervals),
    overlap(overlap)
{
    double o=overlap,m=1-o;
    double a=(intervals*m)/(gaps*m+o);
    double b=(intervals*o)/(gaps*m+o);
    for (unsigned i=0; i<gaps; i++){
        froms.push_back(i*a);
    }
    for (unsigned i=0; i<gaps-1; i++){
        tos.push_back(i*a+a+b);
    }
    tos.push_back(intervals-1); //последний гап обязательно последний интервал
}

double GapManager::from(unsigned gap)
{
    return froms[gap];
}

double GapManager::to(unsigned gap)
{
    return tos[gap];
}

double GapManager::inRange(unsigned interval, unsigned gap)
{
    return interval>=froms[gap] && interval<=tos[gap];
}
