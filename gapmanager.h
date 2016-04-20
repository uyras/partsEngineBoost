#ifndef GAPMANAGER_H
#define GAPMANAGER_H

#include <vector>
#include <cmath>
#include <sstream>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>

using namespace std;
using namespace boost::serialization;
class GapManager
{
public:
    GapManager();
    GapManager(unsigned gaps, unsigned intervals, double overlap);

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & gaps;
        ar & intervals;
        ar & froms;
        ar & tos;
    }

    void setUniform(unsigned gaps, unsigned intervals, double overlap);

    void setLogarithmic(const unsigned gaps, const unsigned intervals, const double overlap=0.6, const double reductionSpeed=0.8);
    void setLogarithmicBothSides(const unsigned gaps, const unsigned intervals, const double overlap=0.6, const double reductionSpeed=0.8);

    /**
     * @brief setLogarithmic2 Устанавливает логарифмически увеличивающуюся шкалу
     * @param gaps Число энергетических зон шкалы
     * @param intervals Число интервалов шкалы
     * @param alfa Скоростьсдвига интервала относительно соседнего
     * @param beta Скорость увеличения ширины интервала
     */
    void setLogarithmic2(unsigned gaps, unsigned intervals, double alfa=2., double beta=1.);

    void setLinear(unsigned gaps, unsigned intervals, double overlap=0.8);

    inline unsigned& from(unsigned gap){ return froms[gap]; }
    inline unsigned& to(unsigned gap){ return tos[gap]; }
    double inRange(unsigned interval, unsigned gap);

    inline unsigned Gaps(){ return gaps; }
    inline unsigned Intervals(){return intervals;}

    std::string toWolframString();

private:
    unsigned gaps, intervals;
    vector<unsigned> froms,tos;

};

#endif // GAPMANAGER_H
