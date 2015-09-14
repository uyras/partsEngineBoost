#include "wanglandaumpi.h"

WangLandauMPI::WangLandauMPI(PartArray *system, unsigned int intervals, unsigned int gaps, double overlap, double accuracy):
    sys(system->copy()),
    gaps(gaps),
    intervals(intervals),
    overlap(overlap),
    accuracy(accuracy),
    average(0.0),
    hCount(0),
    root(0),

    tag_swapEnergy(1),
    tag_swapEnergyPack(2),
    tag_swapConfig(3),
    tag_swapFalse(4),
    tag_finish(5),
    tag_averageQuery(6),
    tag_averageHistogramm(7),
    tag_averagedHistogramm(8),
    tag_complete_swap(9),
    tag_stopsignal(10),
    tag_flatSignal(11)
{
    //@todo Выдавать ошибку если валкеров на гап хотя бы меньше чем три
    //@todo Выдавать ошибку если валкеров на гап нечерное число (например 10 валкеров и 3 гапа

    qDebug()<<"(init) start init";
    walkersByGap = floor(world.size()/gaps);
    size = walkersByGap*gaps;
    gapNumber = floor(world.rank()/walkersByGap);

    //Определяем соседние узлы следующего (соседнего) окна. Для последнего окна соседей не будет (принимает от предыдущих)
    if (gapNumber!=gaps-1) {
        for (unsigned int i=0;i<walkersByGap;i++)
            neightbourWalkers.push_back((gapNumber+1)*walkersByGap+i);
    }

    //определяем узлы своего окна
    for (unsigned int i=0;i<walkersByGap;i++)
        sameWalkers.push_back(gapNumber*walkersByGap+i);

    //считаем минимум и максимум системы
    qDebug()<< "(init) calculating maximal state";
    sys->setToMaximalState();
    this->eMax = sys->calcEnergy1();

    qDebug()<<"(init) calculating ground state";
    sys->setToGroundState();
    this->eMin = sys->calcEnergy1();

    //слегка расширяем границы, дабы нормально работало сравнение double
    double delta=(eMax-eMin)*0.001/(double)(intervals-1);
    eMax+=delta; eMin-=delta;

    this->dE = (eMax-eMin)/(intervals-1);
    this->f = exp(1);
    this->fMin=1.00001;

    qDebug()<<"(init) init starting energy";
    sys->state->reset();
    this->eInit = sys->calcEnergy1FastIncrementalFirst();

    //Добавляем интервалы в гистограмму
    qDebug()<<"(init) init histogramms";
    for (unsigned i=0;i<intervals;i++){
        g.push_back(0);
        h.push_back(0);
    }

    //устанавливаем границы волкера
    double x = (this->eMax-this->eMin)/((double)gaps*(1.-overlap)+overlap); //ширина одного интервала
    this->from = this->eMin + (double)gapNumber * x * (1.-overlap);
    this->to = this->from+x;
    this->nFrom = getIntervalNumber(from);
    this->nTo = getIntervalNumber(to);


    qDebug()<<"(init) averaging init state";
    this->makeNormalInitState();

    finishedProcesses = 0;
    flatedProcesses = 0; thisFlatted = false;
}

WangLandauMPI::~WangLandauMPI()
{

}

vector<double> WangLandauMPI::dos()
{
    qDebug()<<"prepare to start DOS";
    if (world.rank()==0) qDebug()<<"start DOS";
    world.barrier();

    bool continueFlag=true;
    //выполняем несколько шагов WL на каждом walker'е. Каждый walker имеет свой критерий плоскости
    while (continueFlag){
        checkStop();

        qDebug()<<world.rank()<<"start step";
        this->walk(10000);

        if (this->checkFlat() && !this->thisFlatted) //проверяем на плоскость
            this->sygnaliseFlat();

        checkStop();

        //qDebug()<<"swap configurations";
        //проверяем нет ли заявок на обмен конфигами
        //если есть заявка, принимаем ее, отправляем и принимаем конфиг
        qDebug()<<"recieve system";
        this->recieveSystem();
        //если нет заявок на обмен, отправляем заявку и ждем обмена
        qDebug()<<"send system";
        this->sendSystem();
        //qDebug()<<"swap complete";
        checkStop();


        if (this->allFlatted()){ //если все гистограммы в окне плоские,
            this->averageHistogramms(); //усредняем ее
            this->processWalk(); //обнуляем H, усредняем f
        }
        checkStop();

        //если ВСЕ блуждатели завершили работу, сохраняем результат и выходим
        if (this->allFinished()){
            qDebug()<<"All precesses finished!";
            this->saveToFile();

            this->recieveSystem(); //получаем все сообщения перед выходом

            return this->g;
        }

        ofstream f(QString("dump_%1.txt").arg(world.rank()).toStdString().c_str());
        f<<this->dump();
        f.close();

    }
    return this->g;
}

void WangLandauMPI::walk(unsigned stepsPerWalk)
{
    //qDebug()<<"MC step";
    double eOld = sys->calcEnergy1FastIncremental(eInit);
    double eNew;

    //повторяем алгоритм сколько-то шагов
    for (unsigned i=1;i<=stepsPerWalk;i++){

        int partNum = sys->state->randomize();

        eNew = sys->calcEnergy1FastIncremental(eInit);

        double randnum = Random::Instance()->nextDouble();
        if (randnum==0.0)
            randnum=10e-20;
        randnum = std::log(randnum);

        if (inRange(eNew) && randnum <= this->getG(eOld)-this->getG(eNew) ) {
            eOld = eNew;
        } else {
            sys->parts[partNum]->rotate(); //откатываем состояние
        }

        this->updateGH(eOld);
    }
}

unsigned int WangLandauMPI::getIntervalNumber(const double Energy)
{
    return floor((Energy-this->eMin)/this->dE);
}

double WangLandauMPI::getEnergyByInterval(const unsigned int interval)
{
    return dE*(double)interval+eMin;
}

bool WangLandauMPI::inRange(const double _E)
{
    unsigned _ir = getIntervalNumber(_E);
    return _ir<=nTo && _ir>=nFrom;
}

bool WangLandauMPI::inRange(const int _E)
{
    return _E<=nTo && _E>=nFrom;
}

bool WangLandauMPI::checkFlat()
{
    qDebug()<<"check the flattiness of h";
    vector<double>::iterator iter;

    if (average==0.0){
        average = this->calcAverageH();
    }

    iter = h.begin()+nFrom;
    while (iter!=h.begin()+this->nTo+1){ //плоскость гистограммы только в своем интервале
        if ((*iter)!=0. && fabs(*iter-average)/average > (1.0 - accuracy)) //критерий плоскости
            return false;
        iter++;
    }
    return true;
}

bool WangLandauMPI::allFlatted()
{
    while (boost::optional<status> s = this->world.iprobe(any_source,tag_flatSignal)){
        this->world.recv((*s).source(),tag_flatSignal); //получаем запрос, чтобы не валялся в буфере
        flatedProcesses++;
    }

    qDebug()<<"check flatted:"<<flatedProcesses;
    return flatedProcesses==walkersByGap;
}

void WangLandauMPI::sygnaliseFlat()
{
    for (int i=0;i<sameWalkers.size();i++){
        if (sameWalkers[i]!=world.rank())
            world.isend(sameWalkers[i],tag_flatSignal);
    }
    this->flatedProcesses++;
    thisFlatted=true;
    qDebug()<<"Send about flat";
}

void WangLandauMPI::processWalk()
{
    bool beforeFinished = this->finished();
    setValues(h,0); //обнуляем гистограмму
    f=sqrt(f);
    average = 0;
    hCount = 0;
    this->flatedProcesses = 0;
    this->thisFlatted = false;
    qDebug()<<"modify f="<<f<<" on "<<world.rank()<<" walker";
    qDebug()<<"process walk parameters, modify f as "<<f;
    if (!beforeFinished && this->finished()) //если процесс финишировал в процессе изменения параметра f, сигнализируем об этом всем
        this->sygnaliseFinish();
}

bool WangLandauMPI::finished()
{
    return f<=fMin;
}

double WangLandauMPI::getG(double e)
{
    return this->g[this->getIntervalNumber(e)];
}

void WangLandauMPI::updateGH(double E)
{
    if (E==0.){
        E = this->sys->calcEnergy1FastIncremental(this->eInit);
    }

    g[this->getIntervalNumber(E)]+=log(f);

    bool increased = (h[this->getIntervalNumber(E)]+=1) == 1;

    if (inRange(E)){//считаем среднее только в разрешенном интервале
        if (increased){ //прибавляем h и одновременно считаем среднее значение
            //случай если изменилось число ненулевых элементов
            hCount++;
            average = (average * (hCount-1) + 1) / hCount;
        } else {
            average += (1./(double)hCount);
        }
    }
}

void WangLandauMPI::makeNormalInitState()
{
    unsigned long int i=0;
    double eTemp=0;
    while (!inRange(eTemp=sys->calcEnergy1FastIncremental(eInit))){
        this->sys->state->randomize();
        i++;
    }
    qDebug()<<": normalize init state takes "<<i<<" steps";
    this->updateGH(eTemp);
}

double WangLandauMPI::calcAverageH()
{
    double avg=0;
    //считаем среднее значение
    int step=0;
    vector<double>::iterator iter = h.begin()+this->getIntervalNumber(this->from);
    while (iter!=h.begin()+this->getIntervalNumber(this->to)+1){ //плоскость гистограммы только в своем интервале
        if (*iter != 0){
            avg = (avg*step+(*iter))/(step+1);
            step++;
        }
        iter++;
    }
    return avg;
}

bool WangLandauMPI::recieveSystem()
{
    bool recieved=false;
    while (boost::optional<status> s = this->world.iprobe(any_source,tag_swapEnergy)){ //получаем все входящие сообщения
        //qDebug()<<"Recieve system: signal found from"<<(*s).source();
        int pair = (*s).source();
        int ex,ey;
        double gjex,gjey;

        qDebug()<<QString("(recv) %2 recieve pair energy from %1").arg(pair).arg(clock());
        world.recv(pair,tag_swapEnergy,ex);

        if (!this->inRange(ex)){
            qDebug()<<QString("(recv) %2 send that energy not in range %1").arg(pair).arg(clock());
            world.send(pair,tag_swapFalse);
        } else {
            qDebug()<<QString("(recv) %2 send energypack to %1").arg(pair).arg(clock());
            string newState = this->sys->state->toString();
            ey = this->getIntervalNumber(this->sys->calcEnergy1FastIncremental(this->eInit));
            gjex = this->g[ex];
            gjey = this->g[ey];
            packed_oarchive pack(world);
            pack<<ey;
            pack<<gjex;
            pack<<gjey;
            pack<<newState;
            world.send(pair,tag_swapEnergyPack,pack);


            qDebug()<<QString("(recv) %2 waiting for new state from %1").arg(pair).arg(clock());
            status s1 = this->world.probe(pair,any_tag);

            if (s1.tag()==tag_swapConfig){
                string newState;
                world.recv(pair,tag_swapConfig,newState);
                sys->state->fromString(newState);
                qDebug()<<QString("(recv) %2 new state applied from %1").arg(pair).arg(clock());
                recieved = true;
            } else if (s1.tag()==tag_swapFalse) {
                world.recv(pair,tag_swapFalse);
                qDebug()<<QString("(recv) %2 new state rejected by %1").arg(pair).arg(clock());
            }
        }
    }
    return recieved;
}

void WangLandauMPI::sendSystem()
{
    if (neightbourWalkers.size()>0) {
        int pair = neightbourWalkers[Random::Instance()->next(neightbourWalkers.size())];
        qDebug()<<QString("(send) %2 offer to exchange with %1").arg(pair).arg(clock());

        int ex,ey;
        double giex,giey,gjex,gjey;

        ex = this->getIntervalNumber(this->sys->calcEnergy1FastIncremental(this->eInit));
        giex = this->g[ex];

        //отправлем свою энергию и ждем пока блуждатель ответит
        world.send(pair,tag_swapEnergy,ex);

        //получаем ответ
        qDebug()<<QString("(send) %2 waiting for energypack from %1").arg(pair).arg(clock());
        status s = world.probe(pair);
        if (s.tag()==tag_swapFalse){ //если энергия не подходит паре, завершаем процесс обмена
            qDebug()<<QString("(send) %2 exchange rejected by %1").arg(pair).arg(clock());
            world.recv(s.source(),tag_swapFalse);
            return;
        }
        else if (s.tag()==tag_swapEnergyPack){
            //получаем пакет с новыми энергиями
            string newState;
            packed_iarchive pack(world);
            world.recv(pair,tag_swapEnergyPack,pack);
            pack>>ey;
            pack>>gjex;
            pack>>gjey;
            pack>>newState;
            giey = this->g[ey];

            qDebug()<<QString("(send) %2 calculating probability for exchange with %1").arg(pair).arg(clock());
            double p = (giex+gjey) - (giey+gjex);
            double randnum = Random::Instance()->nextDouble();
            if (0.0==randnum) //логарифма нуля не существует
                randnum=10e-20;
            randnum = std::log(randnum);
            if (inRange(ey) && randnum <= p){ //если вероятность получилась, отправляем свой конфиг
                qDebug()<<QString("(send) %2 probability confirmed, send config %1").arg(pair).arg(clock());
                world.send(pair,tag_swapConfig,sys->state->toString());
                sys->state->fromString(newState);
            } else {
                qDebug()<<QString("(send) %2 probability rejected, send signal %1").arg(pair).arg(clock());
                world.send(pair,tag_swapFalse);
            }
        }
    }
}

void WangLandauMPI::averageHistogramms()
{
    qDebug()<<"average histogramms";
    //проверить, достиг ли кто-то еще этой точки. Если нет, стать "хостом".
    if (boost::optional<status> s = this->world.iprobe(any_source,tag_averageQuery)){
        int n=0;
        this->world.recv((*s).source(),tag_averageQuery,n); //получаем запрос, чтобы не валялся в буфере
        this->averageSlave((*s).source());
    } else {
        this->averageMaster();
    }
}

void WangLandauMPI::averageMaster()
{
    qDebug()<<"average as master";
    //отправляем всем сигнал что стал хостом
    for (unsigned int i=0;i<sameWalkers.size();i++){
        if (sameWalkers[i]!=this->world.rank()){
            world.isend(sameWalkers[i],tag_averageQuery,0);
        }
    }

    vector< vector<double> > allG; //массив где будут храниться все гистограммы
    allG.push_back(this->g); //добавляем свою гистограмму

    //получаем все остальные гистограммы
    for (unsigned int i=0;i<sameWalkers.size()-1;i++){
        vector<double> tempG;
        this->world.recv(any_source,tag_averageHistogramm,tempG);
        allG.push_back(tempG);
    }

    //усредняем все гистограммы на примере своей
    double average=0;
    for (unsigned i=0;i<this->intervals;i++){
        average=0;
        for (unsigned j=0;j<allG.size();j++){
            average = (average * (double)j + allG[j][i]) / (double)(j+1);
        }

        this->g[i] = average;
    }

    //отплавляем всем гистограммы
    for (unsigned int i=0;i<sameWalkers.size();i++){
        if (sameWalkers[i]!=this->world.rank()){
            world.isend(sameWalkers[i],tag_averagedHistogramm,this->g);
        }
    }
}

void WangLandauMPI::averageSlave(int host)
{
    qDebug()<<"average as slave";
    world.send(host,tag_averageHistogramm,this->g);
    world.recv(host,tag_averagedHistogramm,this->g);
}

bool WangLandauMPI::allFinished()
{
    while (boost::optional<status> s = this->world.iprobe(any_source,tag_finish)){
        this->world.recv((*s).source(),tag_finish); //получаем запрос, чтобы не валялся в буфере
        finishedProcesses++;
    }

    qDebug()<<"check finished:"<<finishedProcesses;
    return finishedProcesses==size;
}

void WangLandauMPI::sygnaliseFinish()
{
    for (int i=0;i<size;i++){
        if (i!=world.rank())
            world.isend(i,tag_finish);
    }
    this->finishedProcesses++;
    qDebug()<<"Send about finish";
}

void WangLandauMPI::saveToFile()
{
    ofstream f(QString("g_%1.txt").arg(world.rank()).toStdString().c_str());
    vector<double>::iterator iter=this->g.begin();
    int i=0;
    while (iter!=this->g.end()){
        f<<i<<"\t"<<this->getEnergyByInterval(i)<<"\t"<<*iter<<endl;
        iter++; i++;
    }
    f.close();
}

void WangLandauMPI::checkStop()
{
    if (boost::optional<status> s = this->world.iprobe(any_source,tag_stopsignal)){
        this->world.recv((*s).source(),tag_stopsignal); //получаем запрос, чтобы не валялся в буфере
        qDebug()<<"Process stopped by"<<(*s).source();
        std::exit(0);
    }
}

void WangLandauMPI::callStop()
{
    for (int i=0;i<world.size();i++){
        if (i!=world.rank())
            world.isend(i,tag_stopsignal);
    }
    qDebug()<<"Send stop signal to all";
    std::exit(0);
}

void WangLandauMPI::setValues(vector<double> &_h, double _v)
{
    vector<double>::iterator _i = _h.begin();
    while (_i!=_h.end()){
        (*_i)=_v;
        _i++;
    }
}

string WangLandauMPI::dump()
{
    stringstream ss;
    ss<<"rank="<<world.rank()<<endl;
    ss<<"w.size="<<size<<endl;
    ss<<"size="<<size<<endl;
    ss<<"eMin="<<this->eMin<<endl;
    ss<<"dE="<<this->dE<<endl;
    ss<<"fMin="<<this->fMin<<endl;
    ss<<"f="<<this->f<<endl;

    ss<<"from="<<this->from<<endl;
    ss<<"to="<<this->to<<endl;

    ss<<"nFrom="<<this->nFrom<<endl;
    ss<<"nTo="<<this->nTo<<endl;

    ss<<"eInit="<<this->eInit<<endl;
    ss<<"gaps="<<this->gaps<<endl;
    ss<<"walkersByGap="<<this->walkersByGap<<endl;
    ss<<"gapNumber="<<this->gapNumber<<endl;
    ss<<"intervals="<<this->intervals<<endl;
    ss<<"overlap="<<this->overlap<<endl;
    ss<<"accuracy="<<this->accuracy<<endl;

    ss<<"finishedProcesses="<<this->finishedProcesses<<endl;

    ss<<"average="<<this->average<<endl;
    ss<<"hCount="<<this->hCount<<endl;

    ss<<"neightbourWalkers=";
    for (unsigned int i=0;i<neightbourWalkers.size();i++){
        ss<<neightbourWalkers[i]<<",";
    }
    ss<<endl;
    ss<<"sameWalkers=";
    for (unsigned int i=0;i<sameWalkers.size();i++){
        ss<<sameWalkers[i]<<",";
    }
    ss<<endl;

    ss<<"system size="<<this->sys->count()<<endl;
    ss<<"system state="<<this->sys->state->toString()<<endl;
    return ss.str();
}

