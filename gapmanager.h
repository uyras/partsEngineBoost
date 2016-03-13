#ifndef GAPMANAGER_H
#define GAPMANAGER_H

#include <vector>

using namespace std;
class GapManager
{
public:
    GapManager(unsigned gaps, unsigned intervals, double overlap);
    double from(unsigned gap);
    double to(unsigned gap);
    double inRange(unsigned interval, unsigned gap);

    inline unsigned Gaps(){ return gaps; }
    inline unsigned Intervals(){return intervals;}
    inline double Overlap(){return overlap;}

private:
    unsigned gaps, intervals;
    double overlap;
    vector<unsigned> froms,tos;

};

#endif // GAPMANAGER_H
