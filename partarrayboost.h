#ifndef PARTARRAYBOOST_H
#define PARTARRAYBOOST_H

#define UNUSED(x) (void)x; //Макрос для того, чтобы компилятор не ругался на неиспользуемые параметры

#include <boost/mpi.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/split_free.hpp>
#include <string>
#include "Vect.h"
#include "Part.h"
#include "StateMachine.h"
#include "StateMachineFree.h"
#include "PartArray.h"
#include "dos2.h"

namespace boost {
namespace serialization {

template<class Archive> void serialize(Archive & ar, Vect & v, const unsigned int version)
{
    UNUSED(version)

    ar & v.x;
    ar & v.y;
    ar & v.z;
}

template<class Archive> void serialize(Archive & ar, Part & p, const unsigned int version)
{
    UNUSED(version)

    ar & p.pos;
    ar & p.m;
    ar & p.h;
    ar & p.e;
    ar & p.r;
    ar & p._vol;
    ar & p.sector;
    ar & p.state;
    ar & p.h1;
    ar & p.w1;
}

template<class Archive> void load(Archive & ar, StateMachineFree & v, const unsigned int version)
{
    UNUSED(version)

    string s;
    ar & s;
    v.fromString(s);
}

template<class Archive> void save(Archive & ar,const StateMachineFree & v, const unsigned int version)
{
    UNUSED(version)

    string s = v.toString();
    ar & s;
}


template<class Archive> void load(Archive & ar, StateMachine & v, const unsigned int version)
{
    UNUSED(version)

    string s;
    ar & s;
    v.fromString(s);
}

template<class Archive> void save(Archive & ar,const StateMachine & v, const unsigned int version)
{
    UNUSED(version)
    string s = v.toString();
    ar & s;
}

template<class Archive, class T> void load(Archive & ar, Dos2<T> & v, const unsigned int version)
{
    UNUSED(version)
    double min,max;
    unsigned intervals;
    ar & min;
    ar & max;
    ar & intervals;
    v.resize(min,max,intervals);
    T temp;
    for (unsigned i=0; i<intervals; i++){
        ar & temp;
        v.at(i) = temp;
    }
}

template<class Archive, class T> void save(Archive & ar,const Dos2<T> & v, const unsigned int version)
{
    UNUSED(version)
    double min = v.Min(),max = v.Max();
    unsigned intervals = v.Intervals();
    ar & min;
    ar & max;
    ar & intervals;
    for (unsigned i=0; i<v.Intervals(); i++){
        ar & v.at(i);
    }
}

template<class Archive, class T>                         \
inline void serialize(                          \
        Archive & ar,                               \
        Dos2<T> & t,                                      \
        const unsigned int file_version             \
){                                              \
        split_free(ar, t, file_version);            \
}

} // namespace serialization
} // namespace boost

BOOST_SERIALIZATION_SPLIT_FREE(StateMachine)
BOOST_SERIALIZATION_SPLIT_FREE(StateMachineFree)
BOOST_IS_MPI_DATATYPE(Vect)

#endif // PARTARRAYBOOST_H
