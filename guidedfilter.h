#ifndef GUIDED_FILTER_H
#define GUIDED_FILTER_H
#include <iostream>
#include <stdint.h>
#include <thread>
#include <mutex>
#include <string>
#include <vector>
#include <unordered_map>
#include <condition_variable>

template <typename T> class ThreadArg;
template <typename T> class WorkerThread;
template <typename T> void runThread(WorkerThread<T> *t);

template <typename T> void splitC(ThreadArg<T> *arg, int sI, int eI);
template <typename T> void mergeC(ThreadArg<T> *arg, int sI, int eI);
template <typename T> void blurIH(ThreadArg<T> *arg, int sI, int eI);
template <typename T> void blurFH(ThreadArg<T> *arg, int sI, int eI);
template <typename T> void blurFV(ThreadArg<T> *arg, int sI, int eI);
template <typename T> void mulMatI(ThreadArg<T> *arg, int sI, int eI);
template <typename T> void mulMatIF(ThreadArg<T> *arg, int sI, int eI);
template <typename T> void calcVar(ThreadArg<T> *arg, int sI, int eI);
template <typename T> void calcInvVar(ThreadArg<T> *arg, int sI, int eI);
template <typename T> void calcCovDet(ThreadArg<T> *arg, int sI, int eI);
template <typename T> void normVar(ThreadArg<T> *arg, int sI, int eI);
template <typename T> void calcA(ThreadArg<T> *arg, int sI, int eI);
template <typename T> void calcB(ThreadArg<T> *arg, int sI, int eI);
template <typename T> void calcResult(ThreadArg<T> *arg, int sI, int eI);
template <typename T> void calcResultColor(ThreadArg<T> *arg, int sI, int eI);

template <typename T> class ThreadArg
{
public:
    std::mutex mutex;
    std::string funcName;
    T *inI1, *inI2, *result;
    int width, height, channels, r, taskRem, numThreads, taskSize;
    float *inF11, *inF12, *inF21, *inF22, *inF31, *inF32, *output, *inter;
    float eps;
    ThreadArg():eps(0)
    {}
};

template <typename T> class WorkerThread
{
public:
    int tid;
    std::thread *thread;
    std::mutex mutex;
    std::condition_variable cond;
    ThreadArg<T> *arg;
    bool ready, processed, stop;

    WorkerThread(ThreadArg<T> *arg, int tid)
        : arg(arg), tid(tid), ready(false), processed(false), stop(false)
    {}

    ~WorkerThread()
    {
        {
            std::lock_guard<std::mutex> lock(mutex);
            ready=true;
            processed=true;
            stop=true;
        }
        cond.notify_one();
        thread->join();
        delete thread;
    }

    void startThread()
    {
        thread = new std::thread(runThread<T>, this);
    }
};

template <typename T> class ThreadPool
{
public:
    ThreadArg<T> *arg;
    int numThreads, taskSize;
    std::vector<WorkerThread<T>*> threads;
    ThreadPool(int numThreads, int taskSize): numThreads(numThreads), taskSize(taskSize), arg(new ThreadArg<T>())
    {
        for (int i=0; i<numThreads; i++)
        {
            threads.push_back(new WorkerThread<T>(arg, i));
            threads[i]->startThread();
        }
    }
    ~ThreadPool()
    {
        delete arg;
        for (int i=0; i<numThreads; i++)
            delete threads[i];
    }
    void launchThreads()
    {
        arg->taskSize = taskSize;
        arg->taskRem = taskSize;
        arg->numThreads = numThreads;
        for (int i=0; i<numThreads; i++)
        {
            {
                std::lock_guard<std::mutex> lock(threads[i]->mutex);
                threads[i]->ready=true;
                threads[i]->processed=false;
            }
            threads[i]->cond.notify_one();
        }

        for (int i=0; i<numThreads; i++)
        {
            auto *t=threads[i];
            std::unique_lock<std::mutex> lock(t->mutex);
            t->cond.wait(lock, [t]{return t->processed;});
        }
    }
};

template <typename T> class FuncList
{
public:
    FuncList()
    {
        funcMap["splitC"]=splitC, funcMap["mergeC"]=mergeC, funcMap["mulMatI"]=mulMatI,
        funcMap["mulMatIF"]=mulMatIF, funcMap["calcVar"]=calcVar, funcMap["calcInvVar"]=calcInvVar,
        funcMap["calcCovDet"]=calcCovDet, funcMap["normVar"]=normVar, funcMap["calcA"]=calcA,
        funcMap["calcB"]=calcB, funcMap["calcResult"]=calcResult, funcMap["calcResultColor"]=calcResultColor,
        funcMap["blurIH"]=blurIH, funcMap["blurFH"]=blurFH, funcMap["blurFV"]=blurFV;
    }
    std::unordered_map<std::string, void (*)(ThreadArg<T>*,int,int)> funcMap;
};

template <typename T> class AuxMono
{
public:
    T *pc;
    float *mean_I, *var_I, *tempArr, *mean_p, *cov_Ip, *a, *b, *mean_a, *mean_b, *inter;
    int size, width, height, nchannelsP, numThreads;
    ThreadPool<T> *pool;

    AuxMono(int width, int height, int nchannelsP, int numThreads)
        : width(width), height(height), size(width*height), nchannelsP(nchannelsP), numThreads(numThreads)
    {
        pc = new T[size*nchannelsP], mean_I = new float[size], var_I = new float[size],
        tempArr = new float[size], mean_p = new float[size], cov_Ip = new float[size],
        a = new float[size], b = new float[size], mean_a = new float[size],
        mean_b = new float[size], inter = new float[size];
        pool = new ThreadPool<T>(numThreads, height);
    }
    ~AuxMono()
    {
        delete[] mean_I, var_I, tempArr, mean_p, cov_Ip, a, b, mean_a, mean_b, pc, inter;
        delete pool;
    }
};

template <typename T> class AuxColor
{
public:
    T *pc, *tempArrC;
    float *covDet, *tempArr, *mean_p, *mean_I_arr, *var_I_arr, *inv_arr,
        *cov_Ip_arr, *a_arr, *b, *inter;
    int size, width, height, nchannelsP, nchannelsI, numThreads;
    ThreadPool<T> *pool;

    AuxColor(int width, int height, int nchannelsP, int nchannelsI, int numThreads)
        : width(width), height(height), size(width*height), nchannelsP(nchannelsP),
         nchannelsI(nchannelsI), numThreads(numThreads)
    {
        covDet = new float[size], tempArr = new float[size],
        mean_p = new float[size], mean_I_arr = new float[size*nchannelsI],
        var_I_arr = new float[size*nchannelsI*2], inv_arr = new float[size*nchannelsI*2],
        cov_Ip_arr = new float[size*nchannelsI], a_arr = new float[size*nchannelsI],
        b = new float[size], pc = new T[size*nchannelsP],
        tempArrC = new T[size*nchannelsI], inter = new float[size];
        pool = new ThreadPool<T>(numThreads, height);
    }
    ~AuxColor()
    {
        delete[] mean_I_arr, var_I_arr, covDet, tempArr, mean_p, inv_arr, cov_Ip_arr,
                 a_arr, b, pc, tempArrC, inter;
        delete pool;
    }
};

template <typename T> void guidedFilterMono(T *origI, T *p, int r, double eps, T *result, AuxMono<T> *am);
template <typename T> void guidedFilterColor(T *origI, T *p, int r, double eps, T *result, AuxColor<T> *ac);

#endif
