#ifndef PARTARRAYMPI_H
#define PARTARRAYMPI_H

#include <boost/mpi.hpp>
#include <boost/serialization/access.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <partarrayboost.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <QString>

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
      * Сериализация класса
    **/
    friend class boost::serialization::access;
    template<class Archive>
    void save(Archive & ar, const unsigned int version) const
    {
        (void)version;
        int sysSize = this->size();
        ar & sysSize;
        double ir = this->interactionRange();
        ar & ir;
        for (int i=0; i<this->size(); i++)
            ar & (*this).parts[i];
        ar & this->eMin;
        ar & this->eMax;
        ar & this->eInit;
        ar & this->eTemp;
        ar & this->minstate;
        ar & this->maxstate;
        ar & this->state;
        string t = this->type().toStdString();
        ar & t;
    }

    template<class Archive>
    void load(Archive & ar, const unsigned int version)
    {
        (void)version;
        int arraySize;
        ar & arraySize;
        double intRange=0;
        ar & intRange;
        this->setInteractionRange(intRange);
        Part* temp;
        for (int i=0; i<arraySize; i++){
            ar & temp;
            this->insert(temp);
        }
        ar & this->eMin;
        ar & this->eMax;
        ar & this->eInit;
        ar & this->eTemp;
        ar & this->minstate;
        ar & this->maxstate;
        ar & this->state;
        string t;
        ar & t;
        this->_type = QString::fromStdString(t);
    }

    BOOST_SERIALIZATION_SPLIT_MEMBER()

    // Бросает частицы на плоскость в многопоточном режиме, заполняет плоскость до определенного насыщения
    // Может дополнительно добрасывать частицы на уже набросанную плоскость
    //void dropRandomMPI(double maxDestiny, int x=10., int y=10.);

    // удаляет пересекающиеся элементы из швов после dropRandomMPI
    //void filterInterMPI();

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

BOOST_CLASS_VERSION(PartArrayMPI, 1)

#endif // PARTARRAY_H
