#include "gapmanager.h"

GapManager::GapManager():
    gaps(0),
    intervals(0)
{

}

GapManager::GapManager(unsigned gaps, unsigned intervals, double overlap):
    gaps(0),
    intervals(0)
{
    setUniform(gaps,intervals,overlap);
}

void GapManager::setUniform(unsigned gaps, unsigned intervals, double overlap)
{
    this->gaps = gaps;
    this->intervals=intervals;
    froms.resize(gaps);
    tos.resize(gaps);
    double o=overlap,m=1-o;
    double a=(intervals*m)/(gaps*m+o);
    double b=(intervals*o)/(gaps*m+o);
    for (unsigned i=0; i<gaps; i++){
        froms[i]=i*a;
    }
    for (unsigned i=0; i<gaps-1; i++){
        tos[i]=i*a+a+b;
    }
    tos[gaps-1]=intervals-1; //последний гап обязательно последний интервал
}

void GapManager::setLogarithmic(unsigned gaps, unsigned intervals, double alfa, double beta)
{
    this->gaps = gaps;
    this->intervals=intervals;
    froms.resize(gaps);
    tos.resize(gaps);
    double a=2., k=a-1.;

    froms[gaps-1]=(1.-k/a)*intervals;
    tos[gaps-1]=intervals;

    for (unsigned i=1; i<gaps; i++){
        double e=1.-1./std::pow(a,i);
        double sigm=k/std::pow(a,i+1);
        froms[gaps-i-1]=(1.-(e+sigm))*intervals;
        tos[gaps-i-1]=(1.-(e-sigm))*intervals;
    }

    froms[0]=0; //первый гап обязательно первый интервал
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
