#include "wanglandaumpi.h"

WangLandauMPI::WangLandauMPI(PartArray *system, unsigned int intervals, unsigned int gaps, double overlap, double accuracy):
    finishedProcesses(0),
    flatedProcesses(0),
    thisFlatted(false),
    sys(new PartArray(*system)),
    eMin(0),
    eMax(0),
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
    tag_finish(5),
    tag_averageHistogramm(7),
    tag_averagedHistogramm(8),
    tag_complete_swap(9),
    tag_stopsignal(10),
    tag_flatSignal(11),
    tag_saveGistogramm(12)
{
    qDebug()<<"(init) start init";
    size = world.size();
    rank = world.rank();

    //@todo Выдавать ошибку если валкеров на гап хотя бы меньше чем три

    if (size%gaps!=0) {
        qFatal("Error! Process number is not divisible to gap number.");
    }

    if (system->count()==0){
        qFatal("Error! The system is empty.");
    }

    walkersByGap = floor(size/gaps);

    if (walkersByGap<3){
        qFatal("Too small number of processes, minimal is %d", gaps*3);
    }

    gapNumber = floor(rank/walkersByGap);

    //Определяем соседние узлы следующего (соседнего) окна. Для последнего окна соседей не будет (принимает от предыдущих)
    if (gapNumber!=gaps-1) {
        for (unsigned int i=0;i<walkersByGap;i++)
            neightbourWalkers.push_back((gapNumber+1)*walkersByGap+i);
    }

    //определяем узлы своего окна
    for (unsigned int i=0;i<walkersByGap;i++)
        sameWalkers.push_back(gapNumber*walkersByGap+i);


    this->f = exp(1);
    this->fMin=1.00001;

    //Добавляем интервалы в гистограмму
    for (unsigned i=0;i<intervals;i++){
        g.push_back(0);
        h.push_back(0);
    }
}

WangLandauMPI::~WangLandauMPI()
{

}

vector<double> WangLandauMPI::dos()
{
    if (rank==0) qInfo()<<"prepare to start DOS";

    //считаем минимум и максимум системы
    if (sys->Minstate().size()==0 || sys->Maxstate().size()==0)
        qFatal("Min or max state is unknown. DOS calculation is impossible.");

    this->setMinMaxEnergy(sys->EMin(),sys->EMax());

    if (rank==0) qInfo()<<"(init) init starting energy";

    //устанавливаем границы волкера
    this->getFromTo(gapNumber,this->from,this->to);
    this->nFrom = getIntervalNumber(this->from);
    this->nTo = getIntervalNumber(this->to);

    if (rank==0) qInfo()<<"(init) make normal init state";

    this->makeNormalInitStateBothSides();

    world.barrier();
    if (rank==0) qInfo()<<"init state completes;";
    if (rank==0) qInfo()<<"start DOS";

    bool continueFlag=true;
    //выполняем несколько шагов WL на каждом walker'е. Каждый walker имеет свой критерий плоскости
    while (continueFlag){
        checkStop();

        //qDebug()<<rank<<"start step";
        this->walk(10000);

        if (this->checkFlat() && !this->thisFlatted) //проверяем на плоскость
            this->sygnaliseFlat();

        checkStop();

        //проверяем нет ли заявок на обмен конфигами
        //если есть заявка, принимаем ее, отправляем и принимаем конфиг
        //qDebug()<<"recieve system";
        this->recieveSystem();
        //если нет заявок на обмен, отправляем заявку и ждем обмена
        //qDebug()<<"send system";
        this->sendSystem();
        //qDebug()<<"swap complete";
        checkStop();


        if (this->allFlatted()){ //если все гистограммы в окне плоские,
            this->averageHistogramms(); //усредняем ее
            this->processWalk(); //обнуляем H, усредняем f
        }
        checkStop();

        ofstream f(QString("dump_%1.txt").arg(rank).toStdString().c_str());
        f<<this->dump();
        f.close();

        //если ВСЕ блуждатели завершили работу, сохраняем результат и выходим
        if (this->allFinished()){
            if (rank==0) qInfo("WL completed");

            //qDebug()<<"Sleeping 1 seconds before quit";
            sleep(1);
            this->recieveSystem(); //получаем все сообщения перед выходом

            return this->g;
        }

    }
    return this->g;
}

void WangLandauMPI::testDos()
{

    this->setMinMaxEnergy(sys->EMin(),sys->EMax());

    sys->state.hardReset();

    //устанавливаем границы волкера
    this->getFromTo(gapNumber,this->from,this->to);
    this->nFrom = getIntervalNumber(this->from);
    this->nTo = getIntervalNumber(this->to);

    qDebug("Normalize states");

}

void WangLandauMPI::walk(unsigned stepsPerWalk)
{
    //qDebug()<<"MC step";
    double eOld = sys->E();
    double eNew;

    //повторяем алгоритм сколько-то шагов
    for (unsigned i=1;i<=stepsPerWalk;i++){

        int partNum = sys->state.randomize();

        eNew = sys->E();

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
    return dE*(double)interval+this->eMin;
}

bool WangLandauMPI::inRange(const double _E)
{
    unsigned _ir = getIntervalNumber(_E);
    return _ir<=nTo && _ir>=nFrom;
}

bool WangLandauMPI::inRange(const unsigned int _E)
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

    qDebug("check flatted: %d",flatedProcesses);
    return flatedProcesses==walkersByGap;
}

void WangLandauMPI::sygnaliseFlat()
{
    for (unsigned i=0;i<sameWalkers.size();i++){
        if (sameWalkers[i]!=rank)
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
    qDebug()<<"modify f="<<f<<" on "<<rank<<" walker";
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
        E = this->sys->E();
    }

    if (!this->inRange(E)){
        qFatal("Error! Trying to add e=%f out of range (%f(%d), %f(%d)). Real e=%f.",E,from,nFrom,to,nTo,sys->E());
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
    while (!inRange(eTemp=sys->E())){
        this->sys->state.randomize();
        i++;
    }
    qDebug()<<"normalize init state takes "<<i<<" steps";
    this->updateGH(eTemp);
}

void WangLandauMPI::makeNormalInitStateFromGS(bool revert)
{
    unsigned long int i=0;
    int rotated = 0;
    double eTemp, eTempPrev;
    eTemp = eTempPrev = sys->E();


    if (!inRange(eTemp)){ //защита, если система уже в интервале
        do {
            eTempPrev = eTemp;
            rotated = this->sys->state.randomize();
            eTemp=sys->E();
            if (revert){//если новое состояние ниже, отменяем переворот
                if (eTemp>eTempPrev){
                    this->sys->parts[rotated]->rotate();
                    eTemp = eTempPrev;
                }
            } else {
                if (eTemp<eTempPrev){
                    this->sys->parts[rotated]->rotate();
                    eTemp = eTempPrev;
                }
            }
            i++;
        } while (!inRange(eTemp));
    }

    qDebug()<<"normalize init state takes "<<i<<" steps";
    this->updateGH(eTemp);
}

void WangLandauMPI::makeNormalInitStateBothSides()
{
    if (gapNumber<=floor(gaps/2.)){
        qDebug("calc from GS");
        sys->setState(sys->Minstate());

        this->makeNormalInitStateFromGS(false);
    } else {
        qDebug("calc from maximal");
        sys->setState(sys->Maxstate());

        this->makeNormalInitStateFromGS(true);
    }
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
        unsigned ex,ey;
        double gjex,gjey;

        qDebug()<<QString("(recv) recieve pair energy from %1").arg(pair);
        world.recv(pair,tag_swapEnergy,ex);
        qDebug("(recv) recieved from %d E=%f(%d)",pair,getEnergyByInterval(ex),ex);

        string newState = this->sys->state.toString();
        bool rejected = !this->inRange(ex);
        ey = this->getIntervalNumber(this->sys->E());
        gjex = this->g[ex];
        gjey = this->g[ey];
        packed_oarchive pack(world);
        pack<<rejected;
        pack<<ey;
        pack<<gjex;
        pack<<gjey;
        pack<<newState;

        if (rejected){
            qDebug()<<QString("(recv) send that energy not in range %1").arg(pair);
            world.send(pair,tag_swapEnergyPack,pack);
            return recieved;
        } else {
            qDebug()<<QString("(recv) send energypack to %1").arg(pair);
            world.send(pair,tag_swapEnergyPack,pack);
        }

        qDebug()<<QString("(recv) waiting for new state from %1").arg(pair);


        world.recv(pair,tag_swapConfig,newState);
        if (newState.compare("false")==0){
            qDebug()<<QString("(recv) new state rejected by %1").arg(pair);
        } else {
            sys->state.fromString(newState);
            qDebug()<<
                       QString("(recv) new state (%2) applied from %1")
                       .arg(pair)
                       .arg(sys->state.toString().c_str());
            this->updateGH();
            recieved = true;
        }
    }
    return recieved;
}

void WangLandauMPI::sendSystem(int pair)
{
    if (neightbourWalkers.size()>0) {
        if (pair==(-1))
            pair = neightbourWalkers[Random::Instance()->next(neightbourWalkers.size())];

        unsigned ex,ey;
        double giex,giey,gjex,gjey;

        {
            double eTemp = this->sys->E();
            ex = this->getIntervalNumber(eTemp);
            giex = this->g[ex];
            qDebug("(send) offer to exchange with %d, E=%f (%d)", pair, eTemp, ex);
        }

        //отправлем свою энергию и ждем пока блуждатель ответит
        world.send(pair,tag_swapEnergy,ex);

        //получаем ответ
        qDebug()<<QString("(send) waiting for energypack from %1").arg(pair);
        packed_iarchive pack(world);
        world.recv(pair,tag_swapEnergyPack,pack);

        bool rejected; string newState;
        pack>>rejected;
        pack>>ey;
        pack>>gjex;
        pack>>gjey;
        pack>>newState;
        giey = this->g[ey];

        if (rejected){ //если энергия не подходит паре, завершаем процесс обмена
            qDebug()<<QString("(send) exchange rejected by %1").arg(pair);
            return;
        } else {
            qDebug()<<QString("(send) calculating probability for exchange with %1").arg(pair);
            double p = (giex+gjey) - (giey+gjex);
            double randnum = Random::Instance()->nextDouble();
            if (0.0==randnum) //логарифма нуля не существует
                randnum=10e-20;
            randnum = std::log(randnum);
            if (inRange(ey) && randnum <= p){ //если вероятность получилась, отправляем свой конфиг
                qDebug()<<QString("(send) probability confirmed, send config %3 to %1").arg(pair).arg(sys->state.toString().c_str());
                world.send(pair,tag_swapConfig,sys->state.toString());
                sys->state.fromString(newState);
                this->updateGH();
            } else {
                qDebug()<<QString("(send) probability rejected, send signal %1").arg(pair);
                string s = "false";
                world.send(pair,tag_swapConfig,s);
            }
        }
    }
}

void WangLandauMPI::setMinMaxEnergy(double eMin, double eMax)
{
    //Выставляем границы
    this->eMin = eMin;
    this->eMax = eMax;

    //слегка расширяем границы, дабы нормально работало сравнение double
    double delta=(this->eMax-this->eMin)*0.01/(double)(intervals-1);
    this->eMax+=delta; this->eMin-=delta;

    this->dE = (this->eMax-this->eMin)/(intervals-1);
}

void WangLandauMPI::averageHistogramms()
{
    qDebug()<<"average histogramms";
    //Если узел первый в окне, стать хостом
    if (rank==this->sameWalkers[0]){
        this->averageMaster();
    } else {
        this->averageSlave(this->sameWalkers[0]);
    }
}

void WangLandauMPI::averageMaster()
{
    qDebug()<<"average as master";

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
        if (sameWalkers[i]!=this->rank){
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
        int buff=0;
        this->world.recv((*s).source(),tag_finish,buff); //получаем запрос, чтобы не валялся в буфере
        finishedProcesses++;
    }

    qDebug("(%d) check finished: %d",rank,finishedProcesses);
    return finishedProcesses==(unsigned)size;
}

void WangLandauMPI::sygnaliseFinish()
{
    for (int i=0;i<size;i++){
        //не отправлять своей группе сигнал о завершении
        bool sendFlag = true;
        for (unsigned j=0;j<this->sameWalkers.size();j++){
            if (i==sameWalkers[j])
                sendFlag=false;
        }
        if (sendFlag)
            world.isend(i,tag_finish,0);
    }
    this->finishedProcesses+=sameWalkers.size();
    qDebug()<<"Send about finish";
}

void WangLandauMPI::save(string filename)
{
    if (rank!=0){
        world.send(0,tag_saveGistogramm,this->g);
    } else {
        QString fname = QString::fromStdString(filename);
        if (fname=="")
            fname=QString("g_%1_%2.dat").arg(sys->count()).arg(intervals);
        ofstream f(fname.toStdString().c_str());

        f<<"eMin="<<this->eMin<<endl;
        f<<"eMax="<<this->eMax<<endl;
        f<<"dE="<<this->dE<<endl;
        f<<"intervals="<<this->intervals<<endl;
        f<<"gaps="<<this->gaps<<endl;
        f<<"walkers="<<size<<endl;
        f<<"walkersByGap="<<this->walkersByGap<<endl;
        f<<"from=";
        for (unsigned i=0;i<gaps;i++){
            double from,to;
            this->getFromTo(i,from,to);
            f<<from<<",";
        }
        f<<endl;
        f<<"to=";
        for (unsigned i=0;i<gaps;i++){
            double from,to;
            this->getFromTo(i,from,to);
            f<<to<<",";
        }
        f<<endl;
        f<<"nfrom=";
        for (unsigned i=0;i<gaps;i++){
            double from,to;
            this->getFromTo(i,from,to);
            f<<getIntervalNumber(from)<<",";
        }
        f<<endl;
        f<<"nto=";
        for (unsigned i=0;i<gaps;i++){
            double from,to;
            this->getFromTo(i,from,to);
            f<<getIntervalNumber(to)<<",";
        }
        f<<endl;

        for (int i=0;i<size;i++){
            vector<double> gg;
            int sourceNode=0;
            if (i==0){
                gg = this->g;
                sourceNode=0;
            } else {
                status s = world.recv(any_source,tag_saveGistogramm,gg);
                sourceNode = s.source();
            }
            f<<"-----"<<endl;
            f<<sourceNode<<endl;
            vector<double>::iterator iter=gg.begin();
            int j=0;
            while (iter!=gg.end()){
                f<<j<<"\t"<<this->getEnergyByInterval(j)<<"\t"<<*iter<<endl;
                iter++; j++;
            }
        }
        f.close();
    }
    world.barrier();
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
    for (int i=0;i<size;i++){
        if (i!=rank)
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

void WangLandauMPI::getFromTo(double gap, double &from, double &to)
{
    double x = (this->eMax-this->eMin)/((double)gaps*(1.-overlap)+overlap); //ширина одного интервала
    from = this->eMin + gap * x * (1.-overlap);
    to = from+x;
}

string WangLandauMPI::dump()
{
    stringstream ss;
    ss<<"rank="<<rank<<endl;
    ss<<"size="<<size<<endl;
    ss<<"eMin="<<this->eMin<<endl;
    ss<<"eMax="<<this->eMax<<endl;
    ss<<"dE="<<this->dE<<endl;
    ss<<"fMin="<<this->fMin<<endl;
    ss<<"f="<<this->f<<endl;

    ss<<"from="<<this->from<<endl;
    ss<<"to="<<this->to<<endl;

    ss<<"nFrom="<<this->nFrom<<endl;
    ss<<"nTo="<<this->nTo<<endl;

    ss<<"gaps="<<this->gaps<<endl;
    ss<<"walkersByGap="<<this->walkersByGap<<endl;
    ss<<"gapNumber="<<this->gapNumber<<endl;
    ss<<"intervals="<<this->intervals<<endl;
    ss<<"overlap="<<this->overlap<<endl;
    ss<<"accuracy="<<this->accuracy<<endl;

    ss<<"finishedProcesses="<<this->finishedProcesses<<endl;
    ss<<"flfattedProcesses="<<this->flatedProcesses<<endl;

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
    ss<<"system state="<<this->sys->state.toString()<<endl;
    return ss.str();
}

