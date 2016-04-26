#ifndef WANGLANDAUMPI_H
#define WANGLANDAUMPI_H

#include <vector>
#include <PartArray.h>
#include <cmath>
#include <ctime>
#include <QString>
#include <QDebug>
#include <QtDebug>
#include <iomanip>
#include "partarrayboost.h"
#include "random.h"
#include <boost/mpi.hpp>
#include <boost/serialization/vector.hpp>
#include "dos2.h"
#include "gapmanager.h"

using namespace std;
using namespace boost::mpi;

class WangLandauMPI
{
public:
    WangLandauMPI(PartArray *system, unsigned int intervals, double accuracy=0.75, double fmin=1.00001);
    ~WangLandauMPI();

    void run(unsigned stepCount=10000);
    void testDos(); //Функция для тестирования разных функций в динамике
    void walk(unsigned stepsPerWalk);

    void checkParams();
    void init();

    /**
     * @brief processWalk Выполняется только когда h стала плоской
     */
    void processWalk();

    inline bool finished();

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

    bool recieveSystem(int pair); //Получить систему из узла pair
    void sendSystem(int pair); //Отправить систему узлу pair

    void save();

    /**
     * @brief save сохранить гистограммы в файл
     * @param filename Имя файла для сохранения. По умолчанию сохраняет в формате g_<number_of_parts>_<intervals>.dat.
     */
    void save2(std::string filename="");

    /**
     * @brief balanceGaps Балансирует энергетические интервалы в GapManager классе с учетом сложности энергетического ландшафта
     * @param mcSteps Число шагов за одну WL-итерацию
     *
     * Процесс балансировки проходит следующим образом:
     * 1. Запускается mcSteps шагов метрополиса
     * 2. Проверяются на плоскость гистограммы H кадого гапа
     * 3. Если не достигли плоскости, идем в шаг 1. Иначе увеличиваем завершившийся гап на 1%.
     */
    void balanceGaps(unsigned mcSteps=10000);

    /**
     * @brief balanceGaps2 Второй алгоритм балансировки
     * @param mcSteps
     *
     * Алгоритм выполняется только на первом узле (в дальнейшем можно добавить выполнение на всех и дальнейшее усреднение)
     * 1. Выполняется mcSteps WL шагов над системой и формируется g(E) гистограмма
     * 2. По пути от минимума к максимуму проверяется плоскость g(E) на возрастающем интервале (не менее минимального)
     * 3. Если критерий плоскости нарушен, значит в предыдущей точке гистограмма завершается.
     * 4. Следующая начинается с точки с учетом перекрытия
     */
    void balanceGaps2(unsigned mcSteps=100000);


    int size; ///число задействованных блуждателей
    int rank; ///число блуждателей всего

private:
    void averageHistogramms(); //усреднить гистограмму между блуждателями своего окна, блокирующая
    void averageMaster(); //хост усреднения
    void averageSlave(); // ведомые процессы усреднения

    bool allFinished(bool showMessage=true); //Возвращает true если все процессы завершили работу.
    void sygnaliseFinish(); //сообщить всем узлам об окончании
    unsigned int finishedProcesses; //число финишировавших
    bool finishSignalSended; //флаг, отправлен ли сигнал финиша


    bool checkFlat();//критерий плоскости гистограммы
    bool allFlatted(); //Возвращает true, если все блуждатели в окне плоские
    void sygnaliseFlat(); //сообщить всем узлам о достижении плоскости
    unsigned int flatedProcesses;
    bool thisFlatted;

    void checkStop();
    void callStop();

    inline void resetH();

    std::string dump();

    PartArray *sys; //экземпляр вычисляемой системы
    double fMin, f;

    unsigned int walkersByGap, //число блуждателей на интервал
        gapNumber, //номер интервала, в котором работает текущий блуждатель
        intervals; //число интервалов в плотности состояний
    double accuracy; //величина погрешности для степени плоскости гистограммы

    double average; //подсчитывает среднее число для h
    unsigned hCount; //количество ненулевых элементов h, нужно для подсчета среднего
public:
    GapManager gaps; //энергетические интервалы

private:
    vector<int> neightbourWalkers; //валкеры из соседнего окна, для которых возможен обмен
    communicator sameWalkers; //валкеры из того же энергетического окна

    Dos2<double> g;//g - логарифм плотности состояний (энтропия), h - вспомогательная гистограмма, которая должна быть плоской
    Dos2<unsigned> h;

    environment env;
    communicator world;
    int root;
    bool inited;

    enum {
        tag_swapEnergy, //отправка энергий
        tag_swapEnergyPack, //отправка энергии, и двух точек гистограммы
        tag_swapConfig, //отправка конфигурации системы
        tag_finish,//
        tag_averageHistogramm,//
        tag_averagedHistogramm, //
        tag_complete_swap,     //отправляется из хоста всем узлам, когда процесс обмена завершился
        tag_stopsignal,
        tag_flatSignal, //сигнал о том, что узел плоский
        tag_saveGistogramm //cигнал сохранения гистограммы
    };
};

#endif // WANGLANDAUMPI_H
