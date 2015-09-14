#ifndef WANGLANDAUMPI_H
#define WANGLANDAUMPI_H

#include <vector>
#include <PartArray.h>
#include <cmath>
#include <ctime>
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
    inline bool inRange(const double _E);
    inline bool inRange(const int _E);

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


    bool checkFlat();//критерий плоскости гистограммы
    bool allFlatted(); //Возвращает true, если все блуждатели в окне плоские
    void sygnaliseFlat(); //сообщить всем узлам о достижении плоскости
    int flatedProcesses;
    bool thisFlatted;


    void saveToFile();

    void checkStop();
    void callStop();

    inline void setValues(vector<double> &_h, double _v);

    std::string dump();

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
    int size; //xbckj задействованных блуждателей
    int root;

    const int
        tag_swapEnergy, //отправка энергий
        tag_swapEnergyPack, //отправка энергии, и двух точек гистограммы
        tag_swapConfig, //отправка конфигурации системы
        tag_swapFalse, //ошибка отправки (вне энергетической зоны или малая вероятность обмена)
        tag_finish,//
        tag_averageQuery,//
        tag_averageHistogramm,//
        tag_averagedHistogramm, //
        tag_complete_swap,     //отправляется из хоста всем узлам, когда процесс обмена завершился
        tag_stopsignal,
        tag_flatSignal; //сигнал о том, что узел плоский
};

#endif // WANGLANDAUMPI_H
