#ifndef PARTARRAYMPI_H
#define PARTARRAYMPI_H

#include <boost/mpi.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <partarrayboost.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include "PartArray.h"
#include "config.h"
#include "Vect.h"
#include "statemachinegmp.h"


class PartArrayMPI :
        public PartArray
{
public:
    PartArrayMPI();

    /**
    * Создает пустой массив частиц с подложкой размером X,Y,Z (в нанометрах)
    * @param x
    * @param y
    * @param z
    */
    PartArrayMPI(double x, double y, double z);

    /**
    * Создает подложку размером x,y,z и заполняет её случайным набором частиц до заданной плотности
    */
    PartArrayMPI(double x, double y, double z, double density);

    /**
    * Создает подложку размером x,y,z и заполняет её указанным количеством частиц
    */
    PartArrayMPI(double x, double y, double z, int count);

    PartArrayMPI(char* file);

    /**
      * Сериализация класса
    **/
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        (void)version; //Заглушка для компиллятора
        ar & boost::serialization::base_object<PartArray>(*this);
    }

    // Бросает частицы на плоскость в многопоточном режиме, заполняет плоскость до определенного насыщения
    // Может дополнительно добрасывать частицы на уже набросанную плоскость
    void dropRandomMPI(double maxDestiny, int x=10., int y=10.);

    // удаляет пересекающиеся элементы из швов после dropRandomMPI
    void filterInterMPI();

    //проверяет, принадлежит ли сектор потоку. startFrom указывает сколько потоков отведено под root нужды.
    //Обычно это один нулевой поток, который не участвует в распределении секторовы
    bool isMySector(int sector, int startFrom=1);

    bool setToGroundState(int thread=0); //переводим систему в GS методом полного перебора. Пространство состояний разбиваем динамически

    void getMinMaxEnergy(double & eMin, double & eMax); //получает минимум и максимум энергии, записывает их в eMin и eMax


    /*
    //проверяем сколько частиц стоят по полю, а сколько против поля
    void checkFM2(char* file, double c, double realC);
    void calcEnergyMPI1(); //считает энергию в многопоточном режиме
    */

    Vect sector; //размер одного сектора при разделении плоскости на подзадачи (в MPI)

private:
    //преобразует массив ссылок на частицы в массив частиц
    vector<Part> transformToParts();
    void transformFromParts(vector<Part> &temp);
    void _construct();
    boost::mpi::environment env;
    boost::mpi::communicator world;
};

#endif // PARTARRAY_H
