#include "wanglandaumpi.h"

WangLandauMPI::WangLandauMPI(PartArray *system, unsigned int intervals, unsigned int gaps, double overlap, double accuracy):
    sys(system->copy()),
    gaps(gaps),
    intervals(intervals),
    overlap(overlap),
    accuracy(accuracy),
    average(0.0),
    hCount(0),
    root(0)
{
    //@todo Выдавать ошибку если валкеров на гап хотя бы меньше чем три
    //@todo Выдавать ошибку если валкеров на гап нечерное число (например 10 валкеров и 3 гапа

    walkersByGap = floor(world.size()/gaps);
    size = walkersByGap*gaps;
    gapNumber = floor(world.rank()/walkersByGap);

    //Определяем соседние узлы следующего (соседнего) окна
    if (gapNumber!=gaps-1) {
        for (unsigned int i=0;i<walkersByGap;i++)
            neightbourWalkers.push_back((gapNumber+1)*walkersByGap+i);
    } else {
        for (unsigned int i=0;i<walkersByGap;i++)
            neightbourWalkers.push_back(i);
    }

    //определяем узлы своего окна
    for (unsigned int i=0;i<walkersByGap;i++)
        sameWalkers.push_back(gapNumber*walkersByGap+i);

    //считаем минимум и максимум системы
    sys->setToMaximalState();
    this->eMax = sys->calcEnergy1();

    sys->setToGroundState();
    this->eMin = sys->calcEnergy1();

    //слегка расширяем границы, дабы нормально работало сравнение double
    double delta=(eMax-eMin)*0.001/(double)(intervals-1);
    eMax+=delta; eMin-=delta;

    this->dE = (eMax-eMin)/(intervals-1);
    this->f = exp(1);
    this->fMin=1.00001;

    sys->state->reset();
    this->eInit = sys->calcEnergy1FastIncrementalFirst();

    //Добавляем интервалы в гистограмму
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

    this->makeNormalInitState();

    finishedProcesses = 0;
}

WangLandauMPI::~WangLandauMPI()
{

}

vector<double> WangLandauMPI::dos()
{
    if (world.rank()==0) qDebug()<<"start DOS";
    world.barrier();

    bool continueFlag=true;
    //выполняем несколько шагов WL на каждом walker'е. Каждый walker имеет свой критерий плоскости
    while (continueFlag){

        //analys
        makeAnalyseFile();

        //qDebug()<<"start step";
        this->walk(10000);

        //qDebug()<<"swap configurations";
        //проверяем нет ли заявок на обмен конфигами
        //если есть заявка, принимаем ее, отправляем и принимаем конфиг
        if (!this->recieveSystem()){
            //если нет заявок на обмен, отправляем заявку и ждем обмена
            this->sendSystem();
        }


        if (this->isFlat()){ //если гистограмма плоская,
            this->averageHistogramms(); //усредняем ее
            this->processWalk(); //обнуляем H, усредняем f
        }

        //если ВСЕ блуждатели завершили работу, сохраняем результат и выходим
        if (this->allFinished()){
            this->saveToFile();
            return this->g;
        }

    }
    return this->g;
}

void WangLandauMPI::walk(unsigned stepsPerWalk)
{
    qDebug()<<"MC step";
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

bool WangLandauMPI::inRange(const double E)
{
    unsigned i = getIntervalNumber(E);
    return i<=nTo && i>=nFrom;
}

bool WangLandauMPI::isFlat()
{
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

void WangLandauMPI::processWalk()
{
    bool beforeFinished = this->finished();
    setValues(h,0); //обнуляем гистограмму
    f=sqrt(f);
    average = 0;
    hCount = 0;
    qDebug()<<"modify f="<<f<<" on "<<world.rank()<<" walker";
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
    if (boost::optional<status> s = this->world.iprobe(any_source,tag_configSwap)){
        //qDebug()<<"Recieve system: signal found from"<<(*s).source();
        string state;
        world.recv((*s).source(),tag_configSwap,state);
        world.send((*s).source(),tag_configSwap,sys->state->toString());
        sys->state->fromString(state);
        qDebug()<<"Recieve system: swap complete with"<<(*s).source();
        return true;
    } else {
        //qDebug()<<"Recieve system: signal not found";
        return false;
    }
}

void WangLandauMPI::sendSystem()
{
    int rnum=Random::Instance()->next(neightbourWalkers.size());
    //cout<<"random gives "<<(rnum = Random::Instance()->next(neightbourWalkers.size()))<<endl;
    int pair = neightbourWalkers[rnum];
    qDebug()<<"Send system: exchange with"<<pair;
    world.send(pair,tag_configSwap,sys->state->toString());
    string state;
    world.recv(pair,tag_configSwap,state);
    sys->state->fromString(state);
    //qDebug()<<"Send system: exchange complete with"<<pair;
}

void WangLandauMPI::averageHistogramms()
{
    qDebug()<<"Start averaging";
    //проверить, достиг ли кто-то еще этой точки. Если нет, стать "хостом".
    if (boost::optional<status> s = this->world.iprobe(any_source,tag_averageQuery)){
        int n=0;
        this->world.recv((*s).source(),tag_averageQuery,n); //получаем запрос, чтобы не валялся в буфере
        this->averageSlave((*s).source());
    } else {
        this->averageMaster();
    }
    qDebug()<<"Finish averaging";
}

void WangLandauMPI::averageMaster()
{
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
    world.send(host,tag_averageHistogramm,this->g);
    world.recv(host,tag_averagedHistogramm,this->g);
}

bool WangLandauMPI::allFinished()
{
    while (boost::optional<status> s = this->world.iprobe(any_source,tag_finish)){
        this->world.recv((*s).source(),tag_averageQuery); //получаем запрос, чтобы не валялся в буфере
        finishedProcesses++;
    }
    return finishedProcesses==size;
}

void WangLandauMPI::sygnaliseFinish()
{
    for (int i=0;i<size;i++){
        world.isend(i,tag_finish);
    }
}

void WangLandauMPI::makeAnalyseFile()
{
    ofstream f(QString("analys_%1.txt").arg(world.rank()).toStdString().c_str());
    f<<
            this->gapNumber<<"\t"<<
            this->from<<"\t"<<
            this->to<<"\t"<<
        this->f<<"\t"<<endl;
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

    ss<<"average="<<this->average<<endl;
    ss<<"hCount="<<this->hCount<<endl;

    ss<<"neightbourWalkers=";
    for (int i=0;i<neightbourWalkers.size();i++){
        ss<<neightbourWalkers[i]<<",";
    }
    ss<<endl;
    ss<<"sameWalkers=";
    for (int i=0;i<sameWalkers.size();i++){
        ss<<sameWalkers[i]<<",";
    }
    ss<<endl;

    ss<<"system size="<<this->sys->count()<<endl;
    ss<<"system state="<<this->sys->state->toString()<<endl;
    return ss.str();
}

