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

void GapManager::setLogarithmic(const unsigned gaps, const unsigned intervals, const double overlap, const double reductionSpeed)
{
    this->gaps = gaps;
    this->intervals=intervals;
    froms.resize(gaps);
    tos.resize(gaps);

    double kPow=0;
    const double ooverlap=1.-overlap;
    for (unsigned i=0; i<gaps; i++)
        kPow=reductionSpeed*kPow+1.;

    double l=(double)intervals/(ooverlap*kPow+overlap*pow(reductionSpeed,gaps-1));

    for (unsigned i=0; i<gaps; i++){
        kPow=0;
        for (unsigned j=0; j<i; j++){
            kPow=reductionSpeed*kPow+1.;
        }
        tos[gaps-1-i]=(double)(intervals-1)-(l*kPow*ooverlap);
        froms[gaps-1-i]=(double)(intervals-1)-(l*kPow*ooverlap+l*pow(reductionSpeed,i));
    }

    froms[0]=0;
    tos[gaps-1]=intervals-1;
}

void GapManager::setLogarithmicBothSides(const unsigned gaps, const unsigned intervals, const double overlap, const double reductionSpeed)
{
    this->gaps = gaps;
    this->intervals=intervals;
    froms.resize(gaps);
    tos.resize(gaps);

    unsigned halfGaps=floor(gaps/2);
    unsigned halfIntervals=floor(intervals/2);

    double kPow=0;
    const double ooverlap=1.-overlap;

    for (unsigned i=0; i<halfGaps; i++)
        kPow=reductionSpeed*kPow+1.;

    double l=(double)halfIntervals/(ooverlap*kPow+overlap*pow(reductionSpeed,halfGaps-1));

    //первая половина гапов - увеличение интервала
    for (unsigned i=0; i<halfGaps; i++){
        kPow=0;
        for (unsigned j=0; j<i; j++){
            kPow=reductionSpeed*kPow+1.;
        }
        tos[halfGaps-1-i]=(double)(halfIntervals-1)-(l*kPow*ooverlap);
        froms[halfGaps-1-i]=(double)(halfIntervals-1)-(l*kPow*ooverlap+l*pow(reductionSpeed,i));
    }

    unsigned middlePoint=halfGaps;
    if (gaps%2)
        middlePoint=halfGaps+1;

    //вторая половина - уменьшение гапа
    for (unsigned i=0; i<halfGaps; i++){
        kPow=0;
        for (unsigned j=0; j<i; j++){
            kPow=reductionSpeed*kPow+1.;
        }
        froms[middlePoint+i]=halfIntervals+(l*kPow*ooverlap);
        tos[middlePoint+i]=halfIntervals+(l*kPow*ooverlap+l*pow(reductionSpeed,i));
    }

    if (gaps%2){ //если четное число гапов, задаем среднее
        froms[halfGaps]=froms[halfGaps+1]-(tos[halfGaps+1]-froms[halfGaps+1])*overlap;
        tos[halfGaps]=tos[halfGaps-1]+(tos[halfGaps-1]-froms[halfGaps-1])*overlap;
    } else {//если нечетное число гапов, расширяем средние
        tos[halfGaps-1]+=(tos[halfGaps-1]-froms[halfGaps-1])*overlap;
        froms[halfGaps]-=(tos[halfGaps]-froms[halfGaps])*overlap;
    }

    froms[0]=0;
    tos[gaps-1]=intervals-1;
}

void GapManager::setLogarithmic2(unsigned gaps, unsigned intervals, double alfa, double beta)
{
    this->gaps = gaps;
    this->intervals=intervals;
    froms.resize(gaps);
    tos.resize(gaps);
    double a=alfa, k=beta;

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

void GapManager::setLinear(unsigned gaps, unsigned intervals, double overlap)
{
    this->gaps = gaps;
    this->intervals=intervals;

    froms.resize(gaps);
    tos.resize(gaps);

    double m=(double)gaps;

    double mingap=(double)intervals/m;

    for (unsigned i=0; i<gaps; i++){
        double j=i*intervals;
        froms[i]=j/(m*(1.+overlap));
        tos[i]=j/m + mingap;
    }
    froms[0]=0;
    tos[gaps-1]=intervals-1;
}

double GapManager::inRange(unsigned interval, unsigned gap)
{
    return interval>=froms[gap] && interval<=tos[gap];
}

string GapManager::toWolframString()
{
    ostringstream s;
    bool first=true;
    for (unsigned i=0;i<Gaps();i++){
        if (!first)
            s<<",";
        s<<"{{"<<this->from(i)<<","<<i<<"},";
        s<<"{"<<this->to(i)<<","<<i<<"}}";
        first=false;
    }
    return s.str();
}
