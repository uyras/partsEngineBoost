#include "statemachinegmp.h"

void StateMachineGmp::add(const cpp_int & b)
{
    std::vector<char>::iterator iter = this->_state.begin();

    unsigned i=0;
    unsigned short int m=0;
    while (iter!=this->_state.end()){
        m=(int)bit_test(b,i)+(int)(*iter)+m;
        switch (m){
        case 1:
            if ((*iter)!=1)
                (*iter)=!(*iter);
            m=0;
            break;
        case 2:
            if ((*iter)==1)
                (*iter)=!(*iter);
            m=1;
            break;
        case 3:
            if ((*iter)!=1)
                (*iter)=!(*iter);
            m=1;
            break;
        }
        iter++; i++;
    }
    return;
}


StateMachineGmp::StateMachineGmp():StateMachineFree()
{

}

StateMachineGmp::StateMachineGmp(const unsigned long size):StateMachineFree(size)
{

}

StateMachineGmp::StateMachineGmp(const cpp_int &b)
{
    this->add(b);
}

StateMachineGmp::StateMachineGmp(const StateMachine &state):StateMachineFree(state)
{

}

StateMachineGmp::StateMachineGmp(const StateMachine *state):StateMachineFree(*state)
{

}

StateMachineGmp &StateMachineGmp::operator +=(const boost::multiprecision::cpp_int &b)
{
    this->add(b);
    return *this;
}
