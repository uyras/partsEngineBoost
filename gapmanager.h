#ifndef GAPMANAGER_H
#define GAPMANAGER_H

#include <vector>
#include <cmath>

using namespace std;
class GapManager
{
public:
    GapManager();
    GapManager(unsigned gaps, unsigned intervals, double overlap);

    void setUniform(unsigned gaps, unsigned intervals, double overlap);
    /**
     * @brief setLogarithmic Устанавливает логарифмически увеличивающуюся шкалу
     * @param gaps Число энергетических зон шкалы
     * @param intervals Число интервалов шкалы
     * @param alfa Скоростьсдвига интервала относительно соседнего
     * @param beta Скорость увеличения ширины интервала
     */
    void setLogarithmic(unsigned gaps, unsigned intervals, double alfa=2., double beta=1.);

    void setLinear(unsigned gaps, unsigned intervals, double overlap=0.8);

    double from(unsigned gap);
    double to(unsigned gap);
    double inRange(unsigned interval, unsigned gap);

    inline unsigned Gaps(){ return gaps; }
    inline unsigned Intervals(){return intervals;}

private:
    unsigned gaps, intervals;
    vector<unsigned> froms,tos;

};

#endif // GAPMANAGER_H
