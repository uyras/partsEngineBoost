#ifndef WANGLANDAUMPI_H
#define WANGLANDAUMPI_H

#include <vector>
#include <PartArray.h>
#include <cmath>
#include <QString>
#include <QDebug>
#include "partarrayboost.h"
#include "random.h"
#include <boost/mpi.hpp>

using namespace std;
using namespace boost::mpi;

class WangLandauMPI
{
public:
    WangLandauMPI(PartArray *system, unsigned int intervals, unsigned int gaps=10, double overlap=0.75, double accuracy=0.75);
    ~WangLandauMPI();

    vector<double> dos();
    void walk(unsigned stepsPerWalk);

    inline unsigned int getIntervalNumber(const double Energy);
    inline double getEnergyByInterval(const unsigned int interval);
    bool inRange(const double E);

    bool isFlat();//критерий плоскости гистограммы

    /**
     * @brief processWalk Выполняется только когда h стала плоской
     */
    void processWalk();

    inline bool finished();

    inline double getG(double e);

    void updateGH(double E);

    void makeNormalInitState();

    double calcAverageH();

    bool recieveSystem(); //Получить систему из любого узла, неблокирующая операция
    void sendSystem(); //Отправить систему случайному блуждателю своего окна

private:
    void averageHistogramms(); //усреднить гистограмму между блуждателями своего окна, блокирующая
    void averageMaster(); //хост усреднения
    void averageSlave(int host); // ведомые процессы усреднения

    bool allFinished(); //Возвращает true если все процессы завершили работу.
    void sygnaliseFinish(); //сообщить всем узлам об окончании
    int finishedProcesses; //число финишировавших

    void makeAnalyseFile();
    void saveToFile();

    inline void setValues(vector<double> &_h, double _v);

    PartArray *sys; //экземпляр вычисляемой системы
    double eMin,eMax, //максимальная и минимальная энергии системы
        dE, fMin, f;
    double from,to; //энергетические интервалы текущего блуждателя
    unsigned int nFrom, nTo; //численные значения интервалов

    double eInit; // энергия для быстрого вычисления энергии
    unsigned int gaps, //число энергетических пробелов
        walkersByGap, //число блуждателей на интервал
        gapNumber, //номер интервала, в котором работает текущий блуждатель
        intervals; //число интервалов в плотности состояний
    double overlap,//степень перекрытия энергетических окон
        accuracy; //величина погрешности для степени плоскости гистограммы

    double average; //подсчитывает среднее число для h
    unsigned hCount; //количество ненулевых элементов h, нужно для подсчета среднего

    vector<int>
        neightbourWalkers, //валкеры из соседнего окна, для которых возможен обмен
        sameWalkers; //валкеры из того же энергетического окна

    vector<double> h,g;//g - логарифм плотности состояний (энтропия), h - вспомогательная гистограмма, которая должна быть плоской

    environment env;
    communicator world;
    int root;

    const int
        tag_configSwap = 1,//
        tag_finish = 2,//
        tag_averageQuery=3,//
        tag_averageHistogramm=4,//
        tag_averagedHistogramm=5;//
};

#endif // WANGLANDAUMPI_H
