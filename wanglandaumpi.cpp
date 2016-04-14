#include "wanglandaumpi.h"

WangLandauMPI::WangLandauMPI(
        PartArray *system,
        unsigned int intervals,
        unsigned int gapCount,
        double overlap,
        double accuracy,
        double fmin
        ):
    finishedProcesses(0),
    flatedProcesses(0),
    thisFlatted(false),
    sys(new PartArray(*system)),
    fMin(fmin),
    intervals(intervals),
    overlap(overlap),
    accuracy(accuracy),
    average(0.0),
    hCount(0),
    gaps(gapCount,intervals,overlap),
    g(0,0,0),
    h(0,0,0),
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

    if (size%gapCount!=0) {
        qFatal("Error! Process number is not divisible to gap number.");
    }

    if (system->count()==0){
        qFatal("Error! The system is empty.");
    }

    walkersByGap = floor(size/gapCount);

    if (walkersByGap<3){
        qFatal("Too small number of processes, minimal is %d", gapCount*3);
    }

    gapNumber = floor(rank/walkersByGap);

    //Определяем соседние узлы следующего (соседнего) окна. Для последнего окна соседей не будет (принимает от предыдущих)
    if (gapNumber!=gapCount-1) {
        for (unsigned int i=0;i<walkersByGap;i++)
            neightbourWalkers.push_back((gapNumber+1)*walkersByGap+i);
    }

    //определяем узлы своего окна
    for (unsigned int i=0;i<walkersByGap;i++)
        sameWalkers.push_back(gapNumber*walkersByGap+i);

    //считаем минимум и максимум системы
    if (sys->Minstate().size()==0 || sys->Maxstate().size()==0)
        qFatal("Min or max state is unknown. DOS calculation is impossible.");

    //инициируем DOS
    h.resize(sys->E(sys->Minstate()), sys->E(sys->Maxstate()), intervals);
    g.resize(sys->E(sys->Minstate()), sys->E(sys->Maxstate()), intervals);
}

WangLandauMPI::~WangLandauMPI()
{
    delete sys;
}

void WangLandauMPI::run(unsigned stepCount)
{

    if (!gaps.inRange(g.num(sys->E()),gapNumber))
        qFatal("Init state is not in range");

    this->f = std::exp(1);
    this->updateGH(sys->E());

    world.barrier();
    if (rank==0) cout<<"start DOS with "<<stepCount<<" steps per single walk"<<endl;

    bool continueFlag=true;
    //выполняем несколько шагов WL на каждом walker'е. Каждый walker имеет свой критерий плоскости
    while (continueFlag){

        checkStop();

        //qDebug()<<rank<<"start step";
        this->walk(stepCount);

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

            return;
        }

    }
    return;
}

void WangLandauMPI::testDos()
{

    sys->state.hardReset();

    qDebug("Normalize states");

}

void WangLandauMPI::walk(unsigned stepsPerWalk)
{
    //qDebug()<<"MC step";
    double eOld = sys->E();
    double eNew;

    unsigned accepted=0, reverted=0;

    //повторяем алгоритм сколько-то шагов
    for (unsigned i=1;i<=stepsPerWalk;i++){

        int partNum = sys->state.randomize();

        eNew = sys->E();

        double randnum = Random::Instance()->nextDouble();
        if (randnum==0.0)
            randnum=10e-20;
        randnum = std::log(randnum);

        if (gaps.inRange(g.num(eNew),gapNumber) && randnum <= this->g[eOld] - this->g[eNew] ) {
            eOld = eNew;
            accepted++;
        } else {
            sys->parts[partNum]->rotate(true); //откатываем состояние
            reverted++;
        }

        this->updateGH(eOld);
    }

    qDebug()<<"accepted"<<accepted<<"rejected"<<reverted;
}

bool WangLandauMPI::checkFlat()
{
    qDebug()<<"check the flattiness of h";

    if (average==0.0){
        average = this->calcAverageH();
    }

    for (unsigned i=gaps.from(gapNumber); i<=gaps.to(gapNumber); i++){//плоскость гистограммы только в своем интервале
        if (h.at(i)!=0. && fabs(h.at(i)-average)/average > (1.0 - accuracy)) //критерий плоскости
            return false;
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
    resetH(); //обнуляем гистограмму
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

void WangLandauMPI::updateGH(double E)
{
    if (E==0.){
        E = this->sys->E();
    }

    if (!gaps.inRange(g.num(E),gapNumber)){
        qFatal(
                    "Error! Trying to add e=%f out of range (%f(%d), %f(%d)). Real e=%f.",
                    E,
                    g.val(gaps.from(gapNumber)),
                    (int)gaps.from(gapNumber),
                    g.val(gaps.to(gapNumber)),
                    (int)gaps.to(gapNumber),
                    sys->E());
    }

    g[E]+=log(f);

    bool increased = (h[E]+=1) == 1;

    if (gaps.inRange(h.num(E),gapNumber)){//считаем среднее только в разрешенном интервале
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
    while (!gaps.inRange(g.num(eTemp=sys->E()),gapNumber)){
        this->sys->state.randomize();
        i++;
    }
    qDebug()<<"normalize init state takes "<<i<<" steps";
}

void WangLandauMPI::makeNormalInitStateFromGS(bool revert)
{
    unsigned long int i=0;
    int rotated = 0;
    double eTemp, eTempPrev;
    eTemp = eTempPrev = sys->E();


    if (!gaps.inRange(g.num(eTemp),gapNumber)){ //защита, если система уже в интервале
        do {
            eTempPrev = eTemp;
            rotated = this->sys->state.randomize();
            eTemp=sys->E();
            if (revert){//идем сверху вниз по энергии: если новое состояние выше, отменяем переворот
                if (eTemp>eTempPrev){
                    this->sys->parts[rotated]->rotate();
                    eTemp = eTempPrev;
                }
            } else {//идем снизу вверх по энергии: если новое состояние ниже, отменяем переворот
                if (eTemp<eTempPrev){
                    this->sys->parts[rotated]->rotate();
                    eTemp = eTempPrev;
                }
            }
            //если случайно проскочили энергию, меняем направление работы
            if (revert==true) {
                if ((g.num(eTemp)) < (gaps.from(gapNumber))){
                    revert=!revert;
                }
            } else {
                if ((g.num(eTemp)) > (gaps.to(gapNumber))){
                    revert=!revert;
                }
            }

            i++;
            if (i==1000000)
                qInfo()<<"init state makes too long on gap "<<gapNumber<<", E="<<eTemp<<
                         " ("<<g.val(gaps.from(gapNumber))<<";"<<g.val(gaps.to(gapNumber)+1)<<")";
        } while (!gaps.inRange(g.num(eTemp),gapNumber));
    }

    qDebug()<<"normalize init state takes "<<i<<" steps, E="<<eTemp<<
              " ("<<g.val(gaps.from(gapNumber))<<";"<<g.val(gaps.to(gapNumber)+1)<<")";
}

void WangLandauMPI::makeNormalInitStateBothSides()
{
    if (gapNumber<=floor(gaps.Gaps()/2)){
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
    for (unsigned i=gaps.from(gapNumber); i<=gaps.to(gapNumber); i++ ){//плоскость гистограммы только в своем интервале
        if (h.at(i) != 0){
            avg = (avg*step+h.at(i))/(step+1);
            step++;
        }
    }
    return avg;
}

bool WangLandauMPI::recieveSystem()
{
    bool recieved=false;
    while (boost::optional<status> s = this->world.iprobe(any_source,tag_swapEnergy)){ //получаем все входящие сообщения
        //qDebug()<<"Recieve system: signal found from"<<(*s).source();
        int pair = (*s).source();
        double ex,ey;
        double gjex,gjey;

        qDebug()<<QString("(recv) recieve pair energy from %1").arg(pair);
        world.recv(pair,tag_swapEnergy,ex);
        qDebug("(recv) recieved from %d E=%f",pair,ex);

        string newState = this->sys->state.toString();
        bool rejected = !this->gaps.inRange(g.num(ex),gapNumber);
        ey = this->sys->E();
        gjex = this->g[ex];
        gjey = this->g[ey];
        packed_oarchive pack(world);
        pack<<rejected;
        pack<<ey;
        pack<<gjex;
        pack<<gjey;
        pack<<newState;


        world.send(pair,tag_swapEnergyPack,pack);
        if (rejected){
            qDebug()<<QString("(recv) send that energy not in range %1").arg(pair);
            return recieved;
        } else {
            qDebug()<<QString("(recv) send energypack to %1").arg(pair);
        }

        qDebug()<<QString("(recv) waiting for new state from %1").arg(pair);

        //ждем ответ, параллельно проверяем завершились ли остальные процессы
        boost::mpi::request r = world.irecv(pair,tag_swapConfig,newState);
        while (!r.test()){
            if (this->allFinished(false))
                return false;
        }

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

        double ex,ey;
        double giex,giey,gjex,gjey;

        {
            ex = this->sys->E();
            giex = this->g[ex];
            qDebug("(send) offer to exchange with %d, E=%f", pair, ex);
        }

        //отправлем свою энергию и ждем пока блуждатель ответит
        world.send(pair,tag_swapEnergy,ex);

        //получаем ответ
        qDebug()<<QString("(send) waiting for energypack from %1").arg(pair);
        packed_iarchive pack(world);
        //ждем ответ, параллельно проверяем завершились ли остальные процессы
        boost::mpi::request r = world.irecv(pair,tag_swapEnergyPack,pack);
        while (!r.test()){
            if (this->allFinished(false))
                return;
        }

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
            if (gaps.inRange(g.num(ey),gapNumber) && randnum <= p){ //если вероятность получилась, отправляем свой конфиг
                qDebug()<<QString("(send) probability confirmed, send config %2 to %1").arg(pair).arg(sys->state.toString().c_str());
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

    vector< Dos2<double> > allG; //массив где будут храниться все гистограммы
    allG.push_back(this->g); //добавляем свою гистограмму

    //получаем все остальные гистограммы
    for (unsigned int i=0;i<sameWalkers.size()-1;i++){
        Dos2<double> tempG(0,0,0);
        this->world.recv(any_source,tag_averageHistogramm,tempG);
        allG.push_back(tempG);
    }

    //усредняем все гистограммы на примере своей
    double average=0;
    for (unsigned i=0;i<this->intervals;i++){
        average=0;
        for (unsigned j=0;j<allG.size();j++){
            average = (average * (double)j + allG[j].at(i)) / (double)(j+1);
        }

        this->g.at(i) = average;
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

bool WangLandauMPI::allFinished(bool showMessage)
{
    while (boost::optional<status> s = this->world.iprobe(any_source,tag_finish)){
        int buff=0;
        this->world.recv((*s).source(),tag_finish,buff); //получаем запрос, чтобы не валялся в буфере
        finishedProcesses++;
    }

    if (showMessage)
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
            fname=QString("g_%1_%2_%3_%4_%5.dat")
                    .arg(sys->count())
                    .arg(intervals)
                    .arg(gaps.Gaps())
                    .arg(overlap)
                    .arg(accuracy);
        ofstream f(fname.toStdString().c_str());

        f<<"eMin="<<sys->EMin()<<endl;
        f<<"eMax="<<sys->EMax()<<endl;
        f<<"intervals="<<this->intervals<<endl;
        f<<"gaps="<<this->gaps.Gaps()<<endl;
        f<<"walkers="<<size<<endl;
        f<<"walkersByGap="<<this->walkersByGap<<endl;
        f<<"from=";
        for (unsigned i=0;i<gaps.Gaps();i++){
            f<<h.val(gaps.from(i))<<",";
        }
        f<<endl;
        f<<"to=";
        for (unsigned i=0;i<gaps.Gaps();i++){
            f<<h.val(gaps.to(i))<<",";
        }
        f<<endl;
        f<<"nfrom=";
        for (unsigned i=0;i<gaps.Gaps();i++){
            f<<gaps.from(i)<<",";
        }
        f<<endl;
        f<<"nto=";
        for (unsigned i=0;i<gaps.Gaps();i++){
            f<<gaps.to(i)<<",";
        }
        f<<endl;

        for (int i=0;i<size;i++){
            Dos2<double> gg(0,0,0);
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
            for (unsigned i=0;i<gg.Intervals();i++){
                f<<i<<"\t"<<gg.val(i)<<"\t"<<gg.at(i)<<endl;
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

void WangLandauMPI::resetH()
{
    for (unsigned i=0; i<h.Intervals(); i++){
        h.at(i)=0;
    }
}

string WangLandauMPI::dump()
{
    stringstream ss;
    ss<<std::setprecision (15);
    ss<<"rank="<<rank<<endl;
    ss<<"size="<<size<<endl;
    ss<<"eMin="<<sys->EMin()<<endl;
    ss<<"eMax="<<sys->EMax()<<endl;
    ss<<"fMin=" <<this->fMin<<endl;
    ss<<"f=" <<this->f<<endl;

    ss<<"from="<<g.val(gaps.from(gapNumber))<<endl;
    ss<<"to="<<g.val(gaps.to(gapNumber))<<endl;

    ss<<"nFrom="<<gaps.from(gapNumber)<<endl;
    ss<<"nTo="<<gaps.to(gapNumber)<<endl;

    ss<<"gaps="<<this->gaps.Gaps()<<endl;
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

