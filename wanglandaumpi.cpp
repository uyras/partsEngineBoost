#include "wanglandaumpi.h"

WangLandauMPI::WangLandauMPI(
        PartArray *system,
        unsigned int intervals,
        double accuracy,
        double fmin
        ):
    showMessages(false),
    finishedProcesses(0),
    flatedProcesses(0),
    thisFlatted(false),
    sys(system),
    fMin(fmin),
    f(0),
    intervals(intervals),
    accuracy(accuracy),
    average(0.0),
    hCount(0),
    g(0,0,0),
    h(0,0,0),
    root(0),
    inited(false)
{
    qDebug()<<"(init) start init";
    size = world.size();
    rank = world.rank();

    //инициируем DOS
    h.resize(sys->E(sys->Minstate()), sys->E(sys->Maxstate()), intervals);
    g.resize(sys->E(sys->Minstate()), sys->E(sys->Maxstate()), intervals);

    f = std::exp(1);
}

WangLandauMPI::~WangLandauMPI()
{

}

void WangLandauMPI::run(unsigned stepCount)
{
    this->checkParams();
    this->init();
    if (walkersByGap<2){
        qFatal("Too small number of processes, minimal is %d", gaps.Gaps()*2);
    }

    this->updateGH(sys->E());

    world.barrier();
    if (rank==0) cout<<"start DOS with "<<stepCount<<" steps per single walk"<<endl;

    bool continueFlag=true;
    //выполняем несколько шагов WL на каждом walker'е. Каждый walker имеет свой критерий плоскости
    while (continueFlag){

        msg("walk");

        //qDebug()<<rank<<"start step";
        this->walk(stepCount);

        msg("checkFlat");
        if (this->checkFlat() && !this->thisFlatted){ //проверяем на плоскость
            this->sygnaliseFlat();
            msg("sayFlat");
        }

        //проверяем нет ли заявок на обмен конфигами
        //если есть заявка, принимаем ее, отправляем и принимаем конфиг
        //qDebug()<<"recieve system";
        //начинаем процедуру обмена конфигурациями
        msg("barrier before sendrecv");
        world.barrier();
        int pair=0;
        for (int i=0; i<world.size(); i++){
            if (world.rank()==i){
                pair=neightbourWalkers[Random::Instance()->next(neightbourWalkers.size())];
            }
            broadcast(world,pair,i);
            if (world.rank()==i){
                this->sendSystem(pair);
            }
            if (world.rank()==pair){
                this->recieveSystem(i);
            }
        }
        world.barrier();
        qDebug()<<"swap complete";

        msg("sheckAllFlatted");
        if (this->allFlatted()){ //если все гистограммы в окне плоские,
            msg("average G");
            this->averageHistogramms(); //усредняем ее
            this->processWalk(); //обнуляем H, усредняем f
        }

        msg("save dump");
        ofstream f(QString("dump_%1.txt").arg(rank).toStdString().c_str());
        f<<this->dump();
        f.close();

        msg("checkAllFinished");
        //если ВСЕ блуждатели завершили работу, сохраняем результат и выходим
        if (this->allFinished()){
            if (rank==0) qInfo("WL completed");
            world.barrier();
            return;
        }


        msg("check am i finished");
        if (!finishSignalSended && this->finished()){
            msg("say about finish");
            this->sygnaliseFinish();
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

void WangLandauMPI::checkParams()
{
    if (gaps.Gaps()==0 || gaps.Intervals() ==0)
        qFatal("gaps is not inited. Please, one of WangLandauMPI::gaps::set***() function");

    if (size%gaps.Gaps()!=0) {
        qFatal("Error! Process number is not divisible to gap number.");
    }

    if (!gaps.inRange(g.num(sys->E()),gapNumber))
        qFatal("Init state is not in range");

    if (sys->count()==0){
        qFatal("Error! The system is empty.");
    }

    if (sys->Minstate().size()==0 || sys->Maxstate().size()==0)
        qFatal("Min or max state is unknown. DOS calculation is impossible.");

}

void WangLandauMPI::init()
{
    walkersByGap = floor(size/gaps.Gaps());
    gapNumber = floor(rank/walkersByGap);

    //Определяем соседние узлы следующего (соседнего) окна. Для последнего окна соседей не будет (принимает от предыдущих)
    neightbourWalkers.clear();
    if (gapNumber!=0) {
        for (unsigned int i=0;i<walkersByGap;i++)
            neightbourWalkers.push_back((gapNumber-1)*walkersByGap+i);
    }

    if (gapNumber!=gaps.Gaps()-1) {
        for (unsigned int i=0;i<walkersByGap;i++)
            neightbourWalkers.push_back((gapNumber+1)*walkersByGap+i);
    }


    //определяем узлы своего окна
    sameWalkers = world.split(gapNumber);

    //обнуляем G
    for (unsigned i=0; i<g.Intervals(); i++){
        g[i]=0;
    }
    this->resetH();

    this->f = std::exp(1);
    inited=true;
    finishSignalSended=false;
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
    while (boost::optional<status> s = this->sameWalkers.iprobe(any_source,tag_flatSignal)){
        char buff=false;
        this->sameWalkers.recv((*s).source(),tag_flatSignal,buff); //получаем запрос, чтобы не валялся в буфере
        flatedProcesses++;
    }

    qDebug("check flatted: %d",flatedProcesses);
    return flatedProcesses==walkersByGap;
}

void WangLandauMPI::sygnaliseFlat()
{
    for (int i=0;i<sameWalkers.size();i++){
        if (i!=sameWalkers.rank())
            sameWalkers.isend(i,tag_flatSignal,true);
    }
    this->flatedProcesses++;
    thisFlatted=true;
    qDebug()<<"Send about flat";
}

void WangLandauMPI::processWalk()
{
    resetH(); //обнуляем гистограмму
    f=sqrt(f);
    this->flatedProcesses = 0;
    this->thisFlatted = false;
    qDebug()<<"modify f="<<f<<" on "<<rank<<" walker";
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
    this->init();
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
    this->init();

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
    this->init();
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

bool WangLandauMPI::recieveSystem(int pair)
{
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
        return false;
    } else {
        qDebug()<<QString("(recv) send energypack to %1").arg(pair);
    }

    qDebug()<<QString("(recv) waiting for new state from %1").arg(pair);

    //ждем ответ
    world.recv(pair,tag_swapConfig,newState);

    if (newState.compare("false")==0){
        qDebug()<<QString("(recv) new state rejected by %1").arg(pair);
        return false;
    } else {
        sys->state.fromString(newState);
        qDebug()<<
                   QString("(recv) new state applied from %1")
                   .arg(pair);
        this->updateGH();
        return true;
    }
}

void WangLandauMPI::sendSystem(int pair)
{
    if (neightbourWalkers.size()>0) {

        double ex,ey;
        double giex,giey,gjex,gjey;

        ex = this->sys->E();
        giex = this->g[ex];
        qDebug("(send) offer to exchange with %d, E=%f", pair, ex);

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
                randnum=10e-15;
            randnum = std::log(randnum);
            if (gaps.inRange(g.num(ey),gapNumber) && randnum <= p){ //если вероятность получилась, отправляем свой конфиг
                qDebug()<<QString("(send) probability confirmed, send config to %1").arg(pair);
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

void WangLandauMPI::save()
{
    QString fname = QString("g_%1_%2_%3_%4.dat_%5")
            .arg(sys->count())
            .arg(intervals)
            .arg(gaps.Gaps())
            .arg(accuracy)
            .arg(rank);
    ofstream f(fname.toStdString().c_str());

    if (rank==0){
        f<<"intervals="<<this->intervals<<endl;
        f<<"gaps="<<this->gaps.Gaps()<<endl;
        f<<"walkers="<<size<<endl;
        f<<"walkersByGap="<<this->walkersByGap<<endl;
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
        f<<"energies:"<<endl;
        for (unsigned i=0;i<g.Intervals();i++){
            f<<i<<"\t"<<g.val(i)<<"\t"<<g.val(i+1)<<endl;
        }
    }

    f<<"-----"<<endl;
    f<<rank<<endl;
    for (unsigned i=0;i<g.Intervals();i++){
        if (g.at(i)>0)
            f<<i<<"\t"<<g.val(i)<<"\t"<<g.at(i)<<endl;
    }

    f.close();

    world.barrier();
}

void WangLandauMPI::averageHistogramms()
{
    qDebug()<<"average histogramms";
    //Если узел первый в окне, стать хостом
    if (this->sameWalkers.rank()==0){
        this->averageMaster();
    } else {
        this->averageSlave();
    }
}

void WangLandauMPI::averageMaster()
{
    qDebug()<<"average as master";

    vector< Dos2<double> > allG; //массив где будут храниться все гистограммы
    allG.push_back(this->g); //добавляем свою гистограмму

    //получаем все остальные гистограммы
    for (int i=1;i<sameWalkers.size();i++){
        Dos2<double> tempG(0,0,0);
        this->sameWalkers.recv(any_source,tag_averageHistogramm,tempG);
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
    for (int i=1;i<sameWalkers.size();i++){
        sameWalkers.send(i,tag_averagedHistogramm,this->g);
    }
}

void WangLandauMPI::averageSlave()
{
    qDebug()<<"average as slave";
    sameWalkers.send(0,tag_averageHistogramm,this->g);
    sameWalkers.recv(0,tag_averagedHistogramm,this->g);
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
    qDebug()<<"Send about finish";
    for (int i=0;i<size;i++){
        if (i!=rank)
            world.isend(i,tag_finish,0);
    }
    this->finishedProcesses++;
    this->finishSignalSended=true;
}

void WangLandauMPI::save2(string filename)
{
    if (rank!=0){
        world.send(0,tag_saveGistogramm,this->g);
    } else {
        QString fname = QString::fromStdString(filename);
        if (fname=="")
            fname=QString("g_%1_%2_%3_%4.dat")
                    .arg(sys->count())
                    .arg(intervals)
                    .arg(gaps.Gaps())
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

void WangLandauMPI::balanceGaps(unsigned mcSteps)
{

    StateMachineFree initState=sys->state;
    const double minSize=(gaps.Intervals()/gaps.Gaps())*0.5;
    const double flatCriterion=0.2;
    bool finished=false;
    double scalePercent=0.01;

    while (!finished){
        //очищаем систему
        this->makeNormalInitStateBothSides();
//        this->checkParams();
        this->init();
        //if (walkersByGap<2){
        //    qFatal("Too small number of processes, minimal is %d", gaps.Gaps()*2);
        //}
        this->updateGH(sys->E());

        //проходим МК-шаги
        this->walk(mcSteps);

        //собираем информацию о плоскости гистограмм
        vector<double> flattiness;
        double disp=0, m=0, m2=0;
        unsigned iCount=0;
        for (unsigned i=gaps.from(gapNumber); i<=gaps.to(gapNumber); i++){//считаем матожидание
            if (h.at(i)!=0.){
                m=((double)m*(double)iCount+(double)h.at(i))/(double)(iCount+1);
                m2=((double)m2*(double)iCount+(double)(h.at(i)*h.at(i)))/(double)(iCount+1);
                iCount++;
            }
        }
        disp=m2-m*m;

        gather(world,disp,flattiness,0);

        if (rank==0) {
            vector<double> flatAvg(gaps.Gaps(),0);
            //усредняем значения блуждателей
            for (unsigned i=0; i<gaps.Gaps(); i++){
                for (unsigned j=0; j<walkersByGap; j++){
                    flatAvg[i]=(flatAvg[i]*j+flattiness[i*gaps.Gaps()+j])/(double)(j+1);
                }
            }


            //выбираем максимальное
            unsigned flatNum=0; double maxFlat=0;
            for (unsigned i=0; i<gaps.Gaps(); i++){
                if (flatAvg[i]>maxFlat){
                    flatNum=i;
                    maxFlat=flattiness[i];
                }
            }

            for (unsigned i=0;i<gaps.Gaps();i++){
                if (
                        gaps.from(i)>(gaps.Intervals()-1) ||
                        gaps.from(i)<0 ||
                        gaps.to(i)>(gaps.Intervals()-1) ||
                        gaps.to(i)<0
                        ){
                    int j=0;
                }
            }

            //если размер позволяет уменьшить гап, уменьшаем с сохранением уровня перекрытия
            if (gaps.to(flatNum)-gaps.from(flatNum) > minSize){

                if (flatNum!=(size-1)){ //если несбалансированный гап не последний, уменьшаем справа
                    gaps.to(flatNum)-=(double)(gaps.to(flatNum)-gaps.from(flatNum))*scalePercent;
                    gaps.from(flatNum+1)-=(double)(gaps.from(flatNum+1)-gaps.from(flatNum))*scalePercent;
                }
                if (flatNum!=0) { //если несбалансированный гап не первый, уменьшаем слева
                    gaps.to(flatNum-1) +=
                            ((double)(gaps.to(flatNum-1)-gaps.from(flatNum-1))*
                            (double)(gaps.from(flatNum)-gaps.to(flatNum))*
                            scalePercent) /
                            (double)(gaps.from(flatNum-1)-gaps.to(flatNum));
                    gaps.from(flatNum)+=(double)(gaps.to(flatNum)-gaps.from(flatNum))*scalePercent;
                }
            } else { //если уменьшаемый гап достиг минимального разера, завершаем баланс
                finished=true;
            }

            for (unsigned i=0;i<gaps.Gaps();i++){
                if (
                        gaps.from(i)>(gaps.Intervals()-1) ||
                        gaps.from(i)<0 ||
                        gaps.to(i)>(gaps.Intervals()-1) ||
                        gaps.to(i)<0
                        ){
                    int j=0;
                }
            }

            //определяем, завершилась ли балансировка
            double balanceAvg=0;
            for (unsigned i=0; i<flatAvg.size(); i++){
                balanceAvg=(balanceAvg*i+flatAvg[i])/(double)(i+1);
            }
            bool allBalanced=true;
            for (unsigned i=0; i<flatAvg.size(); i++){
                if ((flatAvg[i]/balanceAvg)>(1.+flatCriterion) || (flatAvg[i]/balanceAvg)<(1.-flatCriterion)){
                    allBalanced=false;
                    break;
                }
            }
            if (allBalanced)
                finished=true;
        }

        //отправляем всем флаг finished
        broadcast(world,finished,0);

        //отправляем всем новое распределение гапов
        broadcast(world,gaps,0);
    }
    sys->state = initState;
}

void WangLandauMPI::balanceGaps2(unsigned mcSteps)
{

}

void WangLandauMPI::resetH()
{
    for (unsigned i=0; i<h.Intervals(); i++){
        h.at(i)=0;
    }
    average = 0;
    hCount = 0;
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

    ss<<"system size="<<this->sys->count()<<endl;
    //ss<<"system state="<<this->sys->state.toString()<<endl;
    return ss.str();
}

