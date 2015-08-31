#ifndef STATEMACHINEGMP_H
#define STATEMACHINEGMP_H

#include <boost/multiprecision/cpp_int.hpp>
#include <vector>
#include "StateMachineFree.h"

using namespace boost::multiprecision;

class StateMachineGmp: public StateMachineFree
{
public:
    StateMachineGmp();
    StateMachineGmp(const unsigned long int size);
    StateMachineGmp(const cpp_int & b);
    StateMachineGmp(const StateMachine &state);
    StateMachineGmp(const StateMachine *state);

    void add(const cpp_int & b);
    StateMachineGmp &operator += (const cpp_int &b);
};

#endif // STATEMACHINEGMP_H

