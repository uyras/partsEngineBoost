#include "PartArrayMPI.h"


PartArrayMPI::PartArrayMPI() : PartArray(){}

/*
void PartArrayMPI::checkFM2(char* file, double c, double realC){

    if (config::Instance()->rank==0)
        std::cout<<"start checking ferromanetizm with destiny="<<c<<endl;

    this->sendParticlesBcast(0);

    std::vector<Part>::iterator iter1, iter2;
    const int total = this->parts.size() * (this->parts.size()-1); //итоговое количество взаимодействий
    double cosFi, angleErr = cos((90.-0.)*M_PI / 180.);

    iter1 = this->parts.begin();
    long int parall = 0, aparall = 0, perpend = 0;//количество параллельных, антипараллельных и перпендикулярных частиц
    Vect h; ///сюда пишется рассчитанное поле взаимодействия

    while (iter1!= this->parts.end()) {

        //MPI добавочка, рассчитывает только частицы из определенного (своего) квадрата
        if (this->isMySector(iter1->sector)) {

            iter2 = this->parts.begin();

            //считаем количество сонаправленных и обратнонаправленных частиц
            while (iter2 != this->parts.end()){
                if (iter1!=iter2){
                    h = iter1->interact(&(*iter2)); //считаем поле взаимодействия двух частиц
                    cosFi = h.scalar(iter1->m)/( h.length() * iter1->m.length() );

                    if (cosFi > angleErr) parall++; else
                        if (cosFi < -1*angleErr) aparall++; else
                            perpend++;

                }
                iter2++;
            }

            //this->calcInteraction(&(*iter1));
        }

        iter1++;
    }

    double par=0,apar=0,perp=0;

    long int *a = new long int[config::Instance()->size];

    //считаем количество параллельных
    MPI_Gather(&parall,1,MPI_LONG,a,1,MPI_LONG,0,MPI_COMM_WORLD);
    if (config::Instance()->rank==0){
        for (int i=1; i<config::Instance()->size; i++){
            par += (double)a[i] / (double)total;
        }
    }

    //считаем количество антипараллельных
    MPI_Gather(&aparall,1,MPI_LONG,a,1,MPI_LONG,0,MPI_COMM_WORLD);
    if (config::Instance()->rank==0){
        for (int i=1; i<config::Instance()->size; i++){
            apar += (double)a[i] / (double)total;
        }
    }

    //считаем количество перпендикулярных
    MPI_Gather(&perpend,1,MPI_LONG,a,1,MPI_LONG,0,MPI_COMM_WORLD);
    if (config::Instance()->rank==0){
        for (int i=1; i<config::Instance()->size; i++){
            perp += (double)a[i] / (double)total;
        }
    }

    if (config::Instance()->rank==0){
        std::ofstream f(file, ios_base::app);
        f
            << par << "\t"
            << apar << "\t"
            << perp << "\t"
            << c << "\t"
            << realC <<endl;
        f.close();
        std::cout<<"end checking ferromanetizm with destiny="<<c<<endl;
    }

}


void PartArrayMPI::calcEnergyMPI1(){
    this->E1 = 0;

    if (config::Instance()->rank!=0){
        std::vector < Part >::iterator iterator2;
        iterator2 = this->parts.begin();
        while (iterator2 != this->parts.end()) {
            if (this->isMySector(iterator2->sector)){
                //std::cout<<"calcing energy for x="<<iterator2->pos.x<<" y="<<iterator2->pos.y<<" sector="<<iterator2->sector<<endl;
                this->E1 += this->calcEnergy1(&*iterator2);
            }
            iterator2++;
        }
        this->E1 *= 0.5; //@todo уточнить почему надо делить на два

        MPI_Send(&(this->E1),1,MPI_DOUBLE,0,96,MPI_COMM_WORLD);
    } else {
        std::cout<<"calculating energy start"<<endl;
        MPI_Status status;
        double buf;
        for (int i=1;i<config::Instance()->size;i++) {
            MPI_Recv(&buf,1,MPI_DOUBLE,MPI_ANY_SOURCE,96,MPI_COMM_WORLD,&status);
            this->E1+=buf;
        }
        std::cout<<"calculating energy complete"<<endl;
    }
}
*/

bool PartArrayMPI::isMySector(int sector, int startFrom){
    if (sector % (config::Instance()->size-startFrom) == config::Instance()->rank-startFrom) return true;
    else return false;
}

/*
void PartArrayMPI::dropRandomMPI(double maxDestiny, int x, int y){

    if (config::Instance()->rank==0) std::cout<<"dropping particles with destiny="<<maxDestiny<<" start"<<endl;

    this->sector.setXYZ(x,y,0);

    int squareNum = 0;
    const int squareCount = (int)(this->size.x * this->size.y) / (x*y);

    double surfVol = x * y;

    Part* temp; //временная частица

    int partCount = 0; //количество сброшеных частиц
    double destiny; //плотность заполнения образца

    std::vector < Part* >::iterator iterator1; // итератор для обхода массива частиц
    bool dropError; //флаг ошибки сброса частицы. true, если 2 частицы пересеклись
    int dropErrorCount=0; //количество неудачных попыток наброса подряд
    //int redropCount=0;
    if (config::Instance()->rank>0) {
        while(squareNum<squareCount){//главный цикл, обходит квадраты
            if (squareNum % (config::Instance()->size-1) != config::Instance()->rank-1) { squareNum++; continue; } //проверяем, чтобы поток взял только свой квадрат

            //std::cout<<"thread "<<config::rank<<" take "<<squareNum<<" sector of "<<squareCount<<endl;

            do {
                //далее наброс частиц
                temp = new Part(); //создаем частицу в памяти

                //пытаемся поместить её на массив (задаём ей координаты)
                dropErrorCount=0;
                do {

                    //генерим координаты, цифра 100 - сколько знаков учитывать после запятой в random
                    temp->pos.x = rand() % (x * 100 );
                    temp->pos.y = rand() % (y * 100 );
                    temp->pos.z = 0;

                    temp->pos.x /= 100;
                    temp->pos.y /= 100;
                    temp->pos.z /= 100;

                    //сдвигаем точку в нужный квадрат
                    //@todo эту часть можно упростить рассчитывая сдвиг для каждого кадрата отдельно
                    temp->pos.x += ( squareNum % (int)(this->size.x/x) ) * x;
                    temp->pos.y += (int)( squareNum / ((int)this->size.x/x) ) * y;


                    //проверяем чтобы сгенная точка не пересекалась ни с какой другой в рассматриваемом квадрате (это значит что площади их сфер не пересекались)
                    iterator1 = this->parts.begin();
                    dropError = false;
                    while (iterator1 != this->parts.end()) {
                        if (temp->pos.radius((*iterator1)->pos).length() <= config::Instance()->partR*2){
                            dropError = true;	break;
                        }
                        ++iterator1;
                    }
                    dropErrorCount++;
                } while (dropError && dropErrorCount<50000); //@TODO переделать, количество попыток не должно задаваться численно

                if (dropErrorCount>=50000) break; //если частицы уже не кидаются, брэйк


                //генерируем вектор M для частицы
                double longitude=(double)rand()/(double)config::Instance()->rand_max*2*M_PI;
                double lattitude=0;
                switch (config::Instance()->dimensions()){
                case 2: lattitude=0; break;
                case 3:
                    lattitude=(double)rand()/(double)config::Instance()->rand_max*2-1; // если частица 2-х мерная то угол отклонения должен быть 0
                    break;
                }

                temp->m.x = config::Instance()->m * cos(longitude)*sqrt(1-lattitude*lattitude);
                temp->m.y = config::Instance()->m * sin(longitude)*sqrt(1-lattitude*lattitude);
                temp->m.z = config::Instance()->m * lattitude;

                //отмесаем к какому квадрату принадлежит частица
                temp->sector = squareNum;

                //добавляем частицу на экземпляр
                this->parts.push_back(temp);
                partCount++;


                //считаем плотность заполнения экземпляра
                destiny = (float)(config::Instance()->vol * partCount) / surfVol;
            } while (destiny < maxDestiny);


            //передаем массив в главный поток и чистим память для следующего прохода
            world.send(0,1,this);
            this->parts.clear();
            partCount=0;
            squareNum++;
        }
    }
    else {
        squareNum=0;
        while(squareNum<squareCount){
            world.recv(boost::mpi::any_source,1,*this);
            squareNum++;
            //std::cout<<"recieve parts from "<<squareNum<<" of "<<squareCount<<" sector with destiny="<<maxDestiny<<endl;
        }
    }


    if (config::Instance()->rank==0) std::cout<<"dropping particles with destiny="<<maxDestiny<<" complete"<<endl;

}
*/
/*
void PartArrayMPI::filterInterMPI(){
    //границы точки вектора
    const float
            seamMinX = (float) ( (float)config::Instance()->partR*2 / this->sector.x ),
            seanMaxX = 1 - seamMinX,
            seamMinY = (float) ( (float)config::Instance()->partR*2 / this->sector.y ),
            seamMaxY = 1 - seamMinY;
    const int
            sectorXCount = (int) ( this->size.x / this->sector.x ),
            sectorYCount = (int) ( this->size.y / this->sector.y );
    (void)sectorYCount;

    vector<Part*>::iterator iter, iter2;

    //Шаг 1. Убираем те элементы, которые выходят за край подложки
    std::cout<<"start filtering particles: step 1, delete border collapces"<<endl;
    iter = this->parts.begin();
    while (iter != this->parts.end()){
        if (
                (*iter)->pos.x < config::Instance()->partR ||
                (*iter)->pos.x > this->size.x - config::Instance()->partR ||
                (*iter)->pos.y < config::Instance()->partR ||
                (*iter)->pos.y > this->size.y - config::Instance()->partR
                ) {
            iter = this->parts.erase(iter);
        } else
            iter++;
    }


    //Шаг 2. убираем пересекающиеся частицы.
    std::cout<<"start filtering particles: step 2, remove crossed particles"<<endl;
    bool deleted; //флаг, удалена ли частица (нужен где-то глубоко внутри)
    iter = this->parts.begin();
    while (iter != this->parts.end()){

        deleted = false;
        //сначала смотрим лежит ли она на краю (необходимо только верхние правые частицы одного квадрата сравнивать с нижними левыми другого
        if (
                //iter->pos.x / this->sector.x - std::floor( iter->pos.x / this->sector.x ) < seamMinX || //выше линии X
                ( (*iter)->pos.x / this->sector.x ) - std::floor( (*iter)->pos.x / this->sector.x ) > seanMaxX || //ниже линии X
                //iter->pos.y / this->sector.y - std::floor( iter->pos.y / this->sector.y ) < seamMinY || //правее линии Y
                ( (*iter)->pos.y / this->sector.y ) - std::floor( (*iter)->pos.y / this->sector.y ) > seamMaxY //левее линии Y
                ) {
            //если лежит, то проверяем, лежит ли её оппонент на краю и не в том же секторе, и не дальше нужного сектора

            //сначала задаем оппонентов
            iter2 = this->parts.begin();
            while (iter2!=this->parts.end()){
                //проверяем частицу
                if (
                        (*iter2)->sector != (*iter)->sector && //должны лежать не в одном и том же секторе
                        ( //должны находиться в соседних секторах
                          (*iter2)->sector + sectorXCount + 1 == (*iter)->sector ||
                          //iter2->sector + sectorXCount == iter->sector ||
                          //iter2->sector + sectorXCount - 1 == iter->sector ||
                          //iter2->sector + 1 == iter->sector ||
                          (*iter2)->sector - 1 == (*iter)->sector || //сравниваем только с верхним левым, верхним, верхним правым, средним правым и нижним правым квадратами
                          (*iter2)->sector - sectorXCount + 1 == (*iter)->sector ||
                          (*iter2)->sector - sectorXCount == (*iter)->sector ||
                          (*iter2)->sector - sectorXCount - 1 == (*iter)->sector
                          ) &&
                        ( //если оппонент лежит на швах
                          (*iter2)->pos.x / this->sector.x - std::floor( (*iter2)->pos.x / this->sector.x ) < seamMinX || //выше линии X
                          //iter2->pos.x / this->sector.x - std::floor( iter2->pos.x / this->sector.x ) > seanMaxX || //ниже линии X
                          (*iter2)->pos.y / this->sector.y - std::floor( (*iter2)->pos.y / this->sector.y ) < seamMinY //|| //правее линии Y
                          //iter2->pos.y / this->sector.y - std::floor( iter2->pos.y / this->sector.y ) > seamMaxY //левее линии Y
                          ) &&
                        (*iter2) != (*iter) && //не должны быть одной и той же точкой (само собой, на всякий случай)
                        (*iter2)->pos.radius((*iter)->pos).length()<=config::Instance()->partR*2 // и наконец, если их окружности пересекаются...
                        ) {
                    //std::cout<<"i1 x="<<iter->pos.x<<" y="<<iter->pos.y<<"; with x="<<iter2->pos.x<<" y="<<iter2->pos.y<<endl;
                    iter = this->parts.erase(iter); //удаляем наш элемент
                    deleted = true;
                    break;
                }
                iter2++;
            }
        }

        if (!deleted) iter++; //если частица была удалена, итератор сам переключится на следующую частицу

    }
    std::cout << "complete filtering particles" << endl;
}
*/

bool PartArrayMPI::setToGroundState(int thread){
    boost::mpi::broadcast(world,*this,thread);//передали набросанную систему всем для работы

    StateMachineFree minstate;
    //const unsigned long long int maxstate =  pow(2.,this->count()-1);

    double eTemp, minE=9999999;

    this->state += (config::Instance()->rank-1);

    //считаем минимум для каждого потока
    do {
        eTemp = this->E();
        if (eTemp<minE) { minE=eTemp; minstate = this->state; }

        for (int j=0;j<config::Instance()->size-1;j++){
            if (!this->state.halfNext())
                break;
        }
    } while (true);

    //передаем минимумы из всех потоков
    //unsigned long long int *stateBuf = new unsigned long long int[config::Instance()->size];
    //double *eBuf = new double[config::Instance()->size];

    /* MPI_Gather(&minstate,1,MPI_UNSIGNED_LONG_LONG,stateBuf,1,MPI_UNSIGNED_LONG_LONG,thread,MPI_COMM_WORLD);
    MPI_Gather(&minE,1,MPI_DOUBLE,eBuf,1,MPI_DOUBLE,thread,MPI_COMM_WORLD);

    if (config::Instance()->rank==0){
        for (int i=1; i<config::Instance()->size; i++){
            if (eBuf[i]<minE) {
                minE=eBuf[i];
                minstate=stateBuf[i];
            }
        }
        (&this->state) = (&minstate);
    }*/
    boost::mpi::broadcast(world,*this,thread);
}

void PartArrayMPI::getMinMaxEnergy(double &eMin, double &eMax)
{
    //делим пространство состояний на подпространства
    boost::multiprecision::cpp_int dState=1, //дельта каждого состояния
            dTotal=1; //Всего состояний
    dTotal<<=this->count()-1;
    dState=dTotal/(world.size()-1);

    double eeMin=DBL_MAX,eeMax=DBL_MIN,eTemp=0;

    this->state.reset();
    StateMachineGmp tempState = this->state;
    if (world.rank()!=0){
        tempState+=(dState*(world.rank()-1)); //потоки обрабатывают свою часть работы, главный остатки
    } else {
        tempState+=(dState*(world.size()-1));
        dState = dTotal - (dState*(world.size()-1));
    }
    this->state=tempState;
    while(dState!=0){
        eTemp = this->E();
        if (eeMin>eTemp){
            eeMin = eTemp;
        }
        if (eeMax<eTemp){
            eeMax = eTemp;
        }
        dState--; this->state.next();
    }
    boost::mpi::all_reduce(world,eeMin,eMin,boost::mpi::minimum<double>());
    boost::mpi::all_reduce(world,eeMax,eMax,boost::mpi::maximum<double>());
    this->state.reset();
}

vector<Part> PartArrayMPI::transformToParts()
{
    vector<Part> temp;
    vector<Part*>::iterator iter = this->parts.begin();
    while (iter!=this->parts.end()){
        temp.push_back(**iter);
        iter++;
    }
    return temp;
}

void PartArrayMPI::transformFromParts(vector<Part>& temp)
{
    this->parts.clear();
    this->parts.resize(temp.size());
    vector<Part>::iterator iter = temp.begin();
    while (iter!=temp.end()){
        this->parts.push_back(&(*iter));
        iter++;
    }
}

void PartArrayMPI::_construct()
{

}
