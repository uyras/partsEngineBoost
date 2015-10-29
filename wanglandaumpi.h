#ifndef WANGLANDAUMPI_H
#define WANGLANDAUMPI_H

#include <vector>
#include <PartArray.h>
#include <cmath>
#include <ctime>
#include <QString>
#include <QDebug>
#include <QtDebug>
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
    void testDos(); //Функция для тестирования разных функций в динамике
    void walk(unsigned stepsPerWalk);

    inline unsigned int getIntervalNumber(const double Energy);
    inline double getEnergyByInterval(const unsigned int interval);
    inline bool inRange(const double _E);
    inline bool inRange(const unsigned int _E);

    /**
     * @brief processWalk Выполняется только когда h стала плоской
     */
    void processWalk();

    inline bool finished();

    inline double getG(double e);

    void updateGH(double E=0.0);

    /**
     * @brief makeNormalInitState
     * генерирует состояние системы, входящее в интервал энергий (this->from, this->to),
     * начиная из любого состояния (не обязательно GS).
     * Просто блуждает по пространству состояний.
     */
    void makeNormalInitState();

    /**
     * @brief makeNormalInitStateFromGS
     * генерирует состояние системы, входящее в интервал (this->from, this->to),
     * корректно работает только начиная из GS.
     * Блуждает по пространству состояний, постоянно повышая энергию.
     * @param revert если false, то метод пытается понизить энергию. Если true, то пытается повысить.
     */
    void makeNormalInitStateFromGS(bool revert = false);

    /**
     * @brief makeNormalInitStateBothSides
     * генерирует состояние системы, входящее в интервал (this->from, this->to),
     * при работе определяет, энергетический интервал находится ближе к концу или к началу энергий.
     * Если к концу, то приводит систему в максимальное состояние и понижает энергию.
     * Если к началу, приводит систему в минимальное состояние и повышает энергию.
     * Для работы должны быть оптимизированы setToGS и setToMaximalState системы.
     */
    void makeNormalInitStateBothSides();

    double calcAverageH();

    bool recieveSystem(); //Получить систему из любого узла, неблокирующая операция
    void sendSystem(int pair=-1); //Отправить систему случайному блуждателю своего окна

    void setMinMaxEnergy(double eMin, double eMax);

    /**
     * @brief save сохранить гистограммы в файл
     * @param filename Имя файла для сохранения. По умолчанию сохраняет в формате g_<number_of_parts>_<intervals>.dat.
     */
    void save(std::string filename="");

private:
    void averageHistogramms(); //усреднить гистограмму между блуждателями своего окна, блокирующая
    void averageMaster(); //хост усреднения
    void averageSlave(int host); // ведомые процессы усреднения

    bool allFinished(); //Возвращает true если все процессы завершили работу.
    void sygnaliseFinish(); //сообщить всем узлам об окончании
    unsigned int finishedProcesses; //число финишировавших


    bool checkFlat();//критерий плоскости гистограммы
    bool allFlatted(); //Возвращает true, если все блуждатели в окне плоские
    void sygnaliseFlat(); //сообщить всем узлам о достижении плоскости
    unsigned int flatedProcesses;
    bool thisFlatted;

    void checkStop();
    void callStop();

    inline void setValues(vector<double> &_h, double _v);

    /**
     * @brief getFromTo Получить индивидуальные границы волкера
     * @param gap номер интервала, для которого получать границы
     * @param from Сюда записывается граница "от"
     * @param to Сюда записывается граница "до"
     */
    void getFromTo(double gap, double & from, double & to);

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
    int size; //число задействованных блуждателей
    int rank;
    int root;

    const int
        tag_swapEnergy, //отправка энергий
        tag_swapEnergyPack, //отправка энергии, и двух точек гистограммы
        tag_swapConfig, //отправка конфигурации системы
        tag_finish,//
        tag_averageHistogramm,//
        tag_averagedHistogramm, //
        tag_complete_swap,     //отправляется из хоста всем узлам, когда процесс обмена завершился
        tag_stopsignal,
        tag_flatSignal, //сигнал о том, что узел плоский
        tag_saveGistogramm; //cигнал сохранения гистограммы
};

#endif // WANGLANDAUMPI_H
