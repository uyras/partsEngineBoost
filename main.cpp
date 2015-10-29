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

#include <QSettings>
//переделать сохранение дампов и файла

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

    const int
            width = 50,
            height = 50,
            lattice = 1;
    SquareSpinIceArray *ssi = new SquareSpinIceArray();
    ssi->dropSpinIce(width, height, lattice);
    sys = ssi;

    qDebug()<<sys->count();

    const int
            intervals = 1000,
            gaps = 20;

    const double
            overlap=0.8,
            accuracy=0.1;


    WangLandauMPI w(sys, intervals, gaps, overlap, accuracy);

    if (comm.rank()==0) qInfo()<<"start Wang Landau DOS";
    w.testDos();

    if (comm.rank()==0) qInfo()<<"WL DOS finished, saving data";
    w.save();

    delete sys;

    qInfo()<<"finish";

    return 0;
}
