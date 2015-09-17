#include <QCoreApplication>
#include <ctime>

#include "config.h"
#include "PartArray.h"
#include "squarespinicearray.h"
#include "honeycombspinicearray.h"

#include "StateMachine.h"
#include "StateMachineFree.h"
#include "wanglandaumpi.h"

#include "random.h"

#include <boost/mpi.hpp>


using namespace std;
using namespace boost::mpi;

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);
    environment env;
    communicator comm;


    Random::Instance(time(NULL)+comm.rank());
    config::Instance()->m = 1;


    PartArray *sys;

    sys = new PartArray("honeycomb_30_circle_1.sys");
    //sys->PartArray::setToGroundState();


    comm.barrier();
    if (comm.rank()==0) qDebug()<<"init Wang Landau Parallel";
    WangLandauMPI w(sys,1000,3,0.8);
    w.setMinMaxEnergy(-78.1466,112.283);
    if (comm.rank()==0) qDebug()<<"start Wang Landau DOS";
    w.dos();



    delete sys;

    cout<<comm.rank()<<" finish"<<endl;
    return 0;
    //return a.exec();
}
