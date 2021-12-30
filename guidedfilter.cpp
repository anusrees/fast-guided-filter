#include "guidedfilter.h"

#define MIN2(a,b) ((a)<(b)?(a):(b))
#define MAX2(a,b) ((a)>(b)?(a):(b))

template <typename T> void runThread(WorkerThread<T> *t)
{
    int size, sI, eI;
    FuncList<T>* list=new FuncList<T>();
    bool stop=false;
    do
    {
        {
            std::unique_lock<std::mutex> lock(t->mutex);
            t->cond.wait(lock, [t]{return t->ready;});
            t->ready=false;
            stop=t->stop;
            while (!t->processed)
            {
                {
                    std::unique_lock<std::mutex> arglock(t->arg->mutex);
                    if (t->arg->taskRem == 0)
                        t->processed = true;
                    else
                    {
                        if (t->arg->taskRem < t->arg->numThreads)
                            size=t->arg->taskRem;
                        else
                            size=t->arg->taskRem/t->arg->numThreads;
                        sI=t->arg->taskSize-t->arg->taskRem;
                        eI=sI+size;
                        t->arg->taskRem-=size;
                    }
                }
                if (sI < eI)
                    list->funcMap[t->arg->funcName](t->arg, sI, eI);
            }
        }
        t->cond.notify_one();
    }
    while (!stop);
}

template <typename T> void splitC(ThreadArg<T> *arg, int sI, int eI)
{
    auto *input=arg->inI1, *output=arg->result;
    int width=arg->width, height=arg->height, channels=arg->channels;
    int index, size=width*height;
    for (int c=0; c<channels; c++)
    {
        for (int i=sI; i<eI; i++)
        {
            for (int j=0; j<width; j++)
            {
                index=i*width+j;
                output[c*size+index]=input[index*channels+c];
            }
        }
    }
}

template <typename T> void mergeC(ThreadArg<T> *arg, int sI, int eI)
{
    auto *input=arg->inI1, *output=arg->result;
    int width=arg->width, height=arg->height, channels=arg->channels;
    int index, size=width*height;
    for (int i=sI; i<eI; i++)
    {
        for (int j=0; j<width; j++)
        {
            index=i*width+j;
            for (int c=0; c<channels; c++)
            {
                output[index*channels+c]=input[c*size+index];
            }
        }
    }
}

template <typename T> void blurIH(ThreadArg<T> *arg, int sI, int eI)
{
    auto *input=arg->inI1;
    float *output=arg->inter;
    int r=arg->r, width=arg->width, height=arg->height;
    int hr=r/2, index=0;
    float sumVal=0, denom=r;
    for (int i=sI; i<eI; i++)
    {
        sumVal=0;

        for (int j=-hr; j<hr; j++)
        {
            index=i*width+MAX2(j,0);
            sumVal+=input[index];
        }

        for (int j=0; j<width; j++)
        {
            index=i*width+MIN2(MAX2((j+hr),0),(width-1));
            sumVal+=input[index];
            output[j*height+i]=sumVal/denom;
            index=i*width+MIN2(MAX2((j-hr),0),(width-1));
            sumVal-=input[index];
        }
    }
}

template <typename T> void blurFH(ThreadArg<T> *arg, int sI, int eI)
{
    float *input=arg->inF11, *output=arg->inter;
    int r=arg->r, width=arg->width, height=arg->height;
    int hr=r/2, index;
    float sumVal=0, denom=r;
    for (int i=sI; i<eI; i++)
    {
        sumVal=0;

        for (int j=-hr; j<hr; j++)
        {
            index=i*width+MAX2(j,0);
            sumVal+=input[index];
        }

        for (int j=0; j<width; j++)
        {
            index=i*width+MIN2(MAX2((j+hr),0),(width-1));
            sumVal+=input[index];
            output[j*height+i]=sumVal/denom;
            index=i*width+MIN2(MAX2((j-hr),0),(width-1));
            sumVal-=input[index];
        }
    }
}

template <typename T> void blurFV(ThreadArg<T> *arg, int sI, int eI)
{
    float *input=arg->inter, *output=arg->output;
    int r=arg->r, width=arg->width, height=arg->height;
    int hr=r/2, index;
    float sumVal=0, denom=r;
    for (int i=0; i<width; i++)
    {
        sumVal=0;

        for (int j=-hr; j<hr; j++)
        {
            index=i*height+MAX2((sI+j),0);
            sumVal+=input[index];
        }

        for (int j=sI; j<eI; j++)
        {
            index=i*height+MIN2(MAX2((j+hr),0),(height-1));
            sumVal+=input[index];
            output[j*width+i]=sumVal/denom;
            index=i*height+MIN2(MAX2((j-hr),0),(height-1));
            sumVal-=input[index];
        }
    }
}

template <typename T> void mulMatI(ThreadArg<T> *arg, int sI, int eI)
{
    auto *in1=arg->inI1, *in2=arg->inI2;
    float *output=arg->output;
    int width=arg->width, height=arg->height, index;
    for (int i=sI; i<eI; i++)
    {
        for (int j=0; j<width; j++)
        {
            index=j+i*width;
            output[index]=static_cast<float>(in1[index])*in2[index];
        }
    }
}

template <typename T> void mulMatIF(ThreadArg<T> *arg, int sI, int eI)
{
    auto *in1=arg->inI1;
    float *in2=arg->inF11, *output=arg->output;
    int width=arg->width, height=arg->height;
    int index;
    for (int i=sI; i<eI; i++)
    {
        for (int j=0; j<width; j++)
        {
            index=j+i*width;
            output[index]=in1[index]*in2[index];
        }
    }
}

template <typename T> void calcVar(ThreadArg<T> *arg, int sI, int eI)
{
    float *in1=arg->inF11, *in2=arg->inF12, *output=arg->output;
    int width=arg->width, height=arg->height;
    float eps=arg->eps;
    int index;
    for (int i=sI; i<eI; i++)
    {
        for (int j=0; j<width; j++)
        {
            index=j+i*width;
            output[index]=output[index]-in1[index]*in2[index]+eps;
        }
    }
}

template <typename T> void calcInvVar(ThreadArg<T> *arg, int sI, int eI)
{
    float *in11=arg->inF11, *in12=arg->inF12, *in21=arg->inF21, *in22=arg->inF22,
            *output=arg->output;
    int width=arg->width, height=arg->height;
    int index;
    for (int i=sI; i<eI; i++)
    {
        for (int j=0; j<width; j++)
        {
            index=j+i*width;
            output[index]=in11[index]*in12[index] - in21[index]*in22[index];
        }
    }
}

template <typename T> void calcCovDet(ThreadArg<T> *arg, int sI, int eI)
{
    float *in11=arg->inF11, *in12=arg->inF12, *in21=arg->inF21, *in22=arg->inF22;
    float *in31=arg->inF31, *in32=arg->inF32, *output=arg->output;
    int width=arg->width, height=arg->height;
    int index;
    for (int i=sI; i<eI; i++)
    {
        for (int j=0; j<width; j++)
        {
            index=j+i*width;
            output[index]=in11[index]*in12[index]+in21[index]*in22[index]+in31[index]*in32[index];
        }
    }    
}

template <typename T> void normVar(ThreadArg<T> *arg, int sI, int eI)
{
    float *dst=arg->output, *src=arg->inF11;
    int width=arg->width, height=arg->height;
    int index;
    for (int i=sI; i<eI; i++)
    {
        for (int j=0; j<width; j++)
        {
            index=j+i*width;
            dst[index]/=src[index];
        }
    }
}

template <typename T> void calcA(ThreadArg<T> *arg, int sI, int eI)
{
    float *in1=arg->inF11, *in2=arg->inF12, *output=arg->output;
    float eps=arg->eps;
    int width=arg->width, height=arg->height;
    int index;
    for (int i=sI; i<eI; i++)
    {
        for (int j=0; j<width; j++)
        {
            index=j+i*width;
            output[index]=in1[index]/(in2[index] + eps);
        }
    }
}

template <typename T> void calcB(ThreadArg<T> *arg, int sI, int eI)
{
    float *in1=arg->inF11, *in2=arg->inF12, *in3=arg->inF21, *output=arg->output;
    int width=arg->width, height=arg->height;
    int index;
    for (int i=sI; i<eI; i++)
    {
        for (int j=0; j<width; j++)
        {
            index=j+i*width;
            if (in3==NULL)
                output[index]=in1[index] - in2[index];
            else
                output[index]=in1[index] - in2[index]*in3[index];
        }
    }
}

template <typename T> void calcResult(ThreadArg<T> *arg, int sI, int eI)
{
    float *in1=arg->inF11, *in3=arg->inF12;
    auto *in2=arg->inI1, *output=arg->result;
    int width=arg->width, height=arg->height;
    int index;
    for (int i=sI; i<eI; i++)
    {
        for (int j=0; j<width; j++)
        {
            index=j+i*width;
            output[index]=MAX2(MIN2((in1[index]*in2[index] + in3[index]), 255), 0);
        }
    }
}

template <typename T> void calcResultColor(ThreadArg<T> *arg, int sI, int eI)
{
    float *in1=arg->inF11, *in2=arg->inF12, *in3=arg->inF21, *in4=arg->inF22;
    auto *output=arg->result;
    int width=arg->width, height=arg->height;
    int index;
    for (int i=sI; i<eI; i++)
    {
        for (int j=0; j<width; j++)
        {
            index=j+i*width;
            output[index]=MAX2(MIN2((in1[index]+in2[index]+in3[index]+in4[index]), 255), 0);
        }
    }
}

template <typename C, typename T> void threadProcess(ThreadArg<T> *arg, std::string str, C *aux)
{
    arg->funcName = str;
    aux->pool->launchThreads();
}

template <typename T> void guidedFilterMono(T *origI, T *p, int r, double eps, T *result, AuxMono<T> *am)
{
    int width=am->width, height=am->height, size=am->size, nchannelsP=am->nchannelsP;
    T *pc = am->pc;
    float *mean_I = am->mean_I, *var_I = am->var_I, *tempArr = am->tempArr;
    float *mean_p = am->mean_p, *cov_Ip = am->cov_Ip, *a = am->a, *b = am->b;
    float *mean_a = am->mean_a, *mean_b = am->mean_b, *inter = am->inter;
    ThreadArg<T> *arg=am->pool->arg;

    r = 2 * r + 1;
    arg->width = width;
    arg->height = height;
    arg->inter = inter;
    arg->r = r;

    arg->inI1 = origI;
    arg->output = mean_I;
    threadProcess<AuxMono<T>, T>(arg, "blurIH", am);
    threadProcess<AuxMono<T>, T>(arg, "blurFV", am);

    arg->inI1=origI;
    arg->inI2=origI;
    arg->output=tempArr;
    threadProcess<AuxMono<T>, T>(arg, "mulMatI", am);

    arg->inF11=tempArr;
    arg->output=var_I;
    threadProcess<AuxMono<T>, T>(arg, "blurFH", am);
    threadProcess<AuxMono<T>, T>(arg, "blurFV", am);

    arg->inF11=mean_I;
    arg->inF12=mean_I;
    arg->output=var_I;
    arg->eps=0;
    threadProcess<AuxMono<T>, T>(arg, "calcVar", am);

    if (nchannelsP == 1)
    {
        arg->inI1 = p;
        arg->output = mean_p;
        threadProcess<AuxMono<T>, T>(arg, "blurIH", am);
        threadProcess<AuxMono<T>, T>(arg, "blurFV", am);

        arg->inI1=origI;
        arg->inI2=p;
        arg->output=tempArr;
        threadProcess<AuxMono<T>, T>(arg, "mulMatI", am);
 
        arg->inF11=tempArr;
        arg->output=cov_Ip;
        threadProcess<AuxMono<T>, T>(arg, "blurFH", am);
        threadProcess<AuxMono<T>, T>(arg, "blurFV", am);

        arg->inF11=mean_I;
        arg->inF12=mean_p;
        arg->output=cov_Ip;
        arg->eps=0;
        threadProcess<AuxMono<T>, T>(arg, "calcVar", am);

        arg->inF11=cov_Ip;
        arg->inF12=var_I;
        arg->output=a;
        arg->eps=eps;
        threadProcess<AuxMono<T>, T>(arg, "calcA", am);

        arg->inF11=mean_p;
        arg->inF12=a;
        arg->inF21=mean_I;
        arg->output=b;
        threadProcess<AuxMono<T>, T>(arg, "calcB", am);

        arg->inF11=a;
        arg->output=mean_a;
        threadProcess<AuxMono<T>, T>(arg, "blurFH", am);
        threadProcess<AuxMono<T>, T>(arg, "blurFV", am);

        arg->inF11=b;
        arg->output=mean_b;
        threadProcess<AuxMono<T>, T>(arg, "blurFH", am);
        threadProcess<AuxMono<T>, T>(arg, "blurFV", am);

        arg->inF11=mean_a;
        arg->inF12=mean_b;
        arg->inI1=origI;
        arg->result=result;
        threadProcess<AuxMono<T>, T>(arg, "calcResult", am);
    }
    else
    {
        arg->inI1=p;
        arg->result=pc;
        arg->channels=nchannelsP;
        threadProcess<AuxMono<T>, T>(arg, "splitC", am);

        for (int i = 0; i < nchannelsP; ++i)
        {
            arg->inI1 = pc + i*size;
            arg->output = mean_p;
            threadProcess<AuxMono<T>, T>(arg, "blurIH", am);
            threadProcess<AuxMono<T>, T>(arg, "blurFV", am);

            arg->inI1=origI;
            arg->inI2=pc + i*size;
            arg->output=tempArr;
            threadProcess<AuxMono<T>, T>(arg, "mulMatI", am);

            arg->inF11=tempArr;
            arg->output=cov_Ip;
            threadProcess<AuxMono<T>, T>(arg, "blurFH", am);
            threadProcess<AuxMono<T>, T>(arg, "blurFV", am);

            arg->inF11=mean_I;
            arg->inF12=mean_p;
            arg->output=cov_Ip;
            arg->eps=0;
            threadProcess<AuxMono<T>, T>(arg, "calcVar", am);

            arg->inF11=cov_Ip;
            arg->inF12=var_I;
            arg->output=a;
            arg->eps=eps;
            threadProcess<AuxMono<T>, T>(arg, "calcA", am);

            arg->inF11=mean_p;
            arg->inF12=a;
            arg->inF21=mean_I;
            arg->output=b;
            threadProcess<AuxMono<T>, T>(arg, "calcB", am);

            arg->inF11=a;
            arg->output=mean_a;
            threadProcess<AuxMono<T>, T>(arg, "blurFH", am);
            threadProcess<AuxMono<T>, T>(arg, "blurFV", am);

            arg->inF11=b;
            arg->output=mean_b;
            threadProcess<AuxMono<T>, T>(arg, "blurFH", am);
            threadProcess<AuxMono<T>, T>(arg, "blurFV", am);

            arg->inF11=mean_a;
            arg->inF12=mean_b;
            arg->inI1=origI;
            arg->result=pc + i*size;
            threadProcess<AuxMono<T>, T>(arg, "calcResult", am);
        }

        arg->inI1=pc;
        arg->result=result;
        arg->channels=nchannelsP;
        threadProcess<AuxMono<T>, T>(arg, "mergeC", am);
    }
}

template <typename T> void guidedFilterColor(T *origI, T *p, int r, double eps, T *result, AuxColor<T> *ac)
{
    int nchannelsP=ac->nchannelsP, size=ac->size;
    int width=ac->width, height=ac->height, nchannelsI=ac->nchannelsI;
    float *covDet = ac->covDet, *tempArr = ac->tempArr, *mean_p = ac->mean_p;
    float *mean_I_arr = ac->mean_I_arr, *var_I_arr = ac->var_I_arr;
    float *inv_arr = ac->inv_arr, *cov_Ip_arr=ac->cov_Ip_arr;
    float *a_arr=ac->a_arr, *b=ac->b, *inter=ac->inter;
    T *pc=ac->pc, *tempArrC = ac->tempArrC;
    std::vector<T*> Ichannels;
    std::vector<float*> mean_I, var_I, inv, cov_Ip, a;
    ThreadArg<T> *arg = ac->pool->arg;

    r = 2 * r + 1;
    arg->width = width;
    arg->height = height;
    arg->inter = inter;
    arg->r = r;

    arg->inI1=origI;
    arg->result=tempArrC;
    arg->channels=nchannelsI;
    threadProcess<AuxColor<T>, T>(arg, "splitC", ac);

    for (int i=0; i<nchannelsI; i++)
    {
        Ichannels.push_back(tempArrC+i*size);
        mean_I.push_back(mean_I_arr+i*size);
        var_I.push_back(var_I_arr+2*i*size);
        var_I.push_back(var_I_arr+(2*i+1)*size);
        inv.push_back(inv_arr+2*i*size);
        inv.push_back(inv_arr+(2*i+1)*size);
        cov_Ip.push_back(cov_Ip_arr+i*size);
        a.push_back(a_arr+i*size);
    }

    for (int i=0; i<nchannelsI; i++)
    {
        arg->inI1 = Ichannels[i];
        arg->output = mean_I[i];
        threadProcess<AuxColor<T>, T>(arg, "blurIH", ac);
        threadProcess<AuxColor<T>, T>(arg, "blurFV", ac);
    }

    arg->inI1=Ichannels[0];
    arg->inI2=Ichannels[0];
    arg->output=tempArr;
    threadProcess<AuxColor<T>, T>(arg, "mulMatI", ac);

    arg->inF11=tempArr;
    arg->output=var_I[0];
    threadProcess<AuxColor<T>, T>(arg, "blurFH", ac);
    threadProcess<AuxColor<T>, T>(arg, "blurFV", ac);

    arg->inF11=mean_I[0];
    arg->inF12=mean_I[0];
    arg->output=var_I[0];
    arg->eps=eps;
    threadProcess<AuxColor<T>, T>(arg, "calcVar", ac);

    arg->inI1=Ichannels[0];
    arg->inI2=Ichannels[1];
    arg->output=tempArr;
    threadProcess<AuxColor<T>, T>(arg, "mulMatI", ac);

    arg->inF11=tempArr;
    arg->output=var_I[1];
    threadProcess<AuxColor<T>, T>(arg, "blurFH", ac);
    threadProcess<AuxColor<T>, T>(arg, "blurFV", ac);

    arg->inF11=mean_I[0];
    arg->inF12=mean_I[1];
    arg->output=var_I[1];
    arg->eps=0;
    threadProcess<AuxColor<T>, T>(arg, "calcVar", ac);

    arg->inI1=Ichannels[0];
    arg->inI2=Ichannels[2];
    arg->output=tempArr;
    threadProcess<AuxColor<T>, T>(arg, "mulMatI", ac);

    arg->inF11=tempArr;
    arg->output=var_I[2];
    threadProcess<AuxColor<T>, T>(arg, "blurFH", ac);
    threadProcess<AuxColor<T>, T>(arg, "blurFV", ac);

    arg->inF11=mean_I[0];
    arg->inF12=mean_I[2];
    arg->output=var_I[2];
    arg->eps=0;
    threadProcess<AuxColor<T>, T>(arg, "calcVar", ac);

    arg->inI1=Ichannels[1];
    arg->inI2=Ichannels[1];
    arg->output=tempArr;
    threadProcess<AuxColor<T>, T>(arg, "mulMatI", ac);

    arg->inF11=tempArr;
    arg->output=var_I[3];
    threadProcess<AuxColor<T>, T>(arg, "blurFH", ac);
    threadProcess<AuxColor<T>, T>(arg, "blurFV", ac);

    arg->inF11=mean_I[1];
    arg->inF12=mean_I[1];
    arg->output=var_I[3];
    arg->eps=eps;
    threadProcess<AuxColor<T>, T>(arg, "calcVar", ac);

    arg->inI1=Ichannels[1];
    arg->inI2=Ichannels[2];
    arg->output=tempArr;
    threadProcess<AuxColor<T>, T>(arg, "mulMatI", ac);

    arg->inF11=tempArr;
    arg->output=var_I[4];
    threadProcess<AuxColor<T>, T>(arg, "blurFH", ac);
    threadProcess<AuxColor<T>, T>(arg, "blurFV", ac);

    arg->inF11=mean_I[1];
    arg->inF12=mean_I[2];
    arg->output=var_I[4];
    arg->eps=0;
    threadProcess<AuxColor<T>, T>(arg, "calcVar", ac);

    arg->inI1=Ichannels[2];
    arg->inI2=Ichannels[2];
    arg->output=tempArr;
    threadProcess<AuxColor<T>, T>(arg, "mulMatI", ac);

    arg->inF11=tempArr;
    arg->output=var_I[5];
    threadProcess<AuxColor<T>, T>(arg, "blurFH", ac);
    threadProcess<AuxColor<T>, T>(arg, "blurFV", ac);

    arg->inF11=mean_I[2];
    arg->inF12=mean_I[2];
    arg->output=var_I[5];
    arg->eps=eps;
    threadProcess<AuxColor<T>, T>(arg, "calcVar", ac);

    arg->inF11=var_I[3];
    arg->inF12=var_I[5];
    arg->inF21=var_I[4];
    arg->inF22=var_I[4];
    arg->output=inv[0];
    threadProcess<AuxColor<T>, T>(arg, "calcInvVar", ac);

    arg->inF11=var_I[4];
    arg->inF12=var_I[2];
    arg->inF21=var_I[1];
    arg->inF22=var_I[5];
    arg->output=inv[1];
    threadProcess<AuxColor<T>, T>(arg, "calcInvVar", ac);

    arg->inF11=var_I[1];
    arg->inF12=var_I[4];
    arg->inF21=var_I[3];
    arg->inF22=var_I[2];
    arg->output=inv[2];
    threadProcess<AuxColor<T>, T>(arg, "calcInvVar", ac);

    arg->inF11=var_I[0];
    arg->inF12=var_I[5];
    arg->inF21=var_I[2];
    arg->inF22=var_I[2];
    arg->output=inv[3];
    threadProcess<AuxColor<T>, T>(arg, "calcInvVar", ac);

    arg->inF11=var_I[2];
    arg->inF12=var_I[1];
    arg->inF21=var_I[0];
    arg->inF22=var_I[4];
    arg->output=inv[4];
    threadProcess<AuxColor<T>, T>(arg, "calcInvVar", ac);

    arg->inF11=var_I[0];
    arg->inF12=var_I[3];
    arg->inF21=var_I[1];
    arg->inF22=var_I[1];
    arg->output=inv[5];
    threadProcess<AuxColor<T>, T>(arg, "calcInvVar", ac);

    arg->inF11=inv[0];
    arg->inF12=var_I[0];
    arg->inF21=inv[1];
    arg->inF22=var_I[1];
    arg->inF31=inv[2];
    arg->inF32=var_I[2];
    arg->output=covDet;
    threadProcess<AuxColor<T>, T>(arg, "calcCovDet", ac);

    for (int i=0; i<nchannelsI*2; i++)
    {
        arg->inF11=covDet;
        arg->output=inv[i];
        threadProcess<AuxColor<T>, T>(arg, "normVar", ac);
    }

    if (nchannelsP == 1)
    {
        arg->inI1 = p;
        arg->output = mean_p;
        threadProcess<AuxColor<T>, T>(arg, "blurIH", ac);
        threadProcess<AuxColor<T>, T>(arg, "blurFV", ac);

        arg->inI2=p;
        for (int i=0; i<nchannelsI; i++)
        {
            arg->inI1=Ichannels[i];
            arg->output=tempArr;
            threadProcess<AuxColor<T>, T>(arg, "mulMatI", ac);

            arg->inF11=tempArr;
            arg->output=cov_Ip[i];
            threadProcess<AuxColor<T>, T>(arg, "blurFH", ac);
            threadProcess<AuxColor<T>, T>(arg, "blurFV", ac);
        }

        arg->eps=0;
        for (int i=0; i<nchannelsI; i++)
        {
            arg->inF11=mean_I[i];
            arg->inF12=mean_p;
            arg->output=cov_Ip[i];
            threadProcess<AuxColor<T>, T>(arg, "calcVar", ac);
        }

        arg->inF11=inv[0];
        arg->inF12=cov_Ip[0];
        arg->inF21=inv[1];
        arg->inF22=cov_Ip[1];
        arg->inF31=inv[2];
        arg->inF32=cov_Ip[2];
        arg->output=a[0];
        threadProcess<AuxColor<T>, T>(arg, "calcCovDet", ac);

        arg->inF11=inv[1];
        arg->inF12=cov_Ip[0];
        arg->inF21=inv[3];
        arg->inF22=cov_Ip[1];
        arg->inF31=inv[4];
        arg->inF32=cov_Ip[2];
        arg->output=a[1];
        threadProcess<AuxColor<T>, T>(arg, "calcCovDet", ac);

        arg->inF11=inv[2];
        arg->inF12=cov_Ip[0];
        arg->inF21=inv[4];
        arg->inF22=cov_Ip[1];
        arg->inF31=inv[5];
        arg->inF32=cov_Ip[2];
        arg->output=a[2];
        threadProcess<AuxColor<T>, T>(arg, "calcCovDet", ac);

        arg->inF11=a[0];
        arg->inF12=mean_I[0];
        arg->inF21=a[1];
        arg->inF22=mean_I[1];
        arg->inF31=a[2];
        arg->inF32=mean_I[2];
        arg->output=tempArr;
        threadProcess<AuxColor<T>, T>(arg, "calcCovDet", ac);

        arg->inF11=mean_p;
        arg->inF12=tempArr;
        arg->inF21=nullptr;
        arg->output=b;
        threadProcess<AuxColor<T>, T>(arg, "calcB", ac);

        for (int i=0; i<nchannelsI; i++)
        {
            arg->inF11=a[i];
            arg->output=tempArr;
            threadProcess<AuxColor<T>, T>(arg, "blurFH", ac);
            threadProcess<AuxColor<T>, T>(arg, "blurFV", ac);

            arg->inI1=Ichannels[i];
            arg->inF11=tempArr;
            arg->output=a[i];
            threadProcess<AuxColor<T>, T>(arg, "mulMatIF", ac);
        }

        arg->inF11=b;
        arg->output=tempArr;
        threadProcess<AuxColor<T>, T>(arg, "blurFH", ac);
        threadProcess<AuxColor<T>, T>(arg, "blurFV", ac);

        arg->inF11=a[0];
        arg->inF12=a[1];
        arg->inF21=a[2];
        arg->inF22=tempArr;
        arg->result=result;
        threadProcess<AuxColor<T>, T>(arg, "calcResultColor", ac);
    }
    else
    {
        arg->inI1=p;
        arg->result=pc;
        arg->channels=nchannelsP;
        threadProcess<AuxColor<T>, T>(arg, "splitC", ac);

        for (int i=0; i<nchannelsP; i++)
        {
            arg->inI1 = pc+i*size;
            arg->output = mean_p;
            threadProcess<AuxColor<T>, T>(arg, "blurIH", ac);
            threadProcess<AuxColor<T>, T>(arg, "blurFV", ac);

            for (int j=0; j<nchannelsI; j++)
            {
                arg->inI2=Ichannels[j];
                arg->output=tempArr;
                threadProcess<AuxColor<T>, T>(arg, "mulMatI", ac);

                arg->inF11=tempArr;
                arg->output=cov_Ip[j];
                threadProcess<AuxColor<T>, T>(arg, "blurFH", ac);
                threadProcess<AuxColor<T>, T>(arg, "blurFV", ac);
            }

            arg->eps=0;
            for (int j=0; j<nchannelsI; j++)
            {
                arg->inF11=mean_I[j];
                arg->inF12=mean_p;
                arg->output=cov_Ip[j];
                threadProcess<AuxColor<T>, T>(arg, "calcVar", ac);
            }

            arg->inF11=inv[0];
            arg->inF12=cov_Ip[0];
            arg->inF21=inv[1];
            arg->inF22=cov_Ip[1];
            arg->inF31=inv[2];
            arg->inF32=cov_Ip[2];
            arg->output=a[0];
            threadProcess<AuxColor<T>, T>(arg, "calcCovDet", ac);

            arg->inF11=inv[1];
            arg->inF12=cov_Ip[0];
            arg->inF21=inv[3];
            arg->inF22=cov_Ip[1];
            arg->inF31=inv[4];
            arg->inF32=cov_Ip[2];
            arg->output=a[1];
            threadProcess<AuxColor<T>, T>(arg, "calcCovDet", ac);

            arg->inF11=inv[2];
            arg->inF12=cov_Ip[0];
            arg->inF21=inv[4];
            arg->inF22=cov_Ip[1];
            arg->inF31=inv[5];
            arg->inF32=cov_Ip[2];
            arg->output=a[2];
            threadProcess<AuxColor<T>, T>(arg, "calcCovDet", ac);

            arg->inF11=a[0];
            arg->inF12=mean_I[0];
            arg->inF21=a[1];
            arg->inF22=mean_I[1];
            arg->inF31=a[2];
            arg->inF32=mean_I[2];
            arg->output=tempArr;
            threadProcess<AuxColor<T>, T>(arg, "calcCovDet", ac);

            arg->inF11=mean_p;
            arg->inF12=tempArr;
            arg->inF21=nullptr;
            arg->output=b;
            threadProcess<AuxColor<T>, T>(arg, "calcB", ac);

            for (int j=0; j<nchannelsI; j++)
            {
                arg->inF11=a[j];
                arg->output=tempArr;
                threadProcess<AuxColor<T>, T>(arg, "blurFH", ac);
                threadProcess<AuxColor<T>, T>(arg, "blurFV", ac);

                arg->inI1=Ichannels[j];
                arg->inF11=tempArr;
                arg->output=a[j];
                threadProcess<AuxColor<T>, T>(arg, "mulMatIF", ac);
            }

            arg->inF11=b;
            arg->output=tempArr;
            threadProcess<AuxColor<T>, T>(arg, "blurFH", ac);
            threadProcess<AuxColor<T>, T>(arg, "blurFV", ac);

            arg->inF11=a[0];
            arg->inF12=a[1];
            arg->inF21=a[2];
            arg->inF22=tempArr;
            arg->result=pc+i*size;
            threadProcess<AuxColor<T>, T>(arg, "calcResultColor", ac);
        }

        arg->inI1=pc;
        arg->result=result;
        arg->channels=nchannelsP;
        threadProcess<AuxColor<T>, T>(arg, "mergeC", ac);
    }
}

template void runThread(WorkerThread<uint8_t> *t);
template void runThread(WorkerThread<uint16_t> *t);
template void runThread(WorkerThread<uint32_t> *t);
template void runThread(WorkerThread<uint64_t> *t);
template void runThread(WorkerThread<int8_t> *t);
template void runThread(WorkerThread<int16_t> *t);
template void runThread(WorkerThread<int32_t> *t);
template void runThread(WorkerThread<int64_t> *t);
template void runThread(WorkerThread<float> *t);
template void runThread(WorkerThread<double> *t);

template void guidedFilterMono(uint8_t *origI, uint8_t *p, int r, double eps, uint8_t *result,
 AuxMono<uint8_t> *aux);
template void guidedFilterMono(uint16_t *origI, uint16_t *p, int r, double eps, uint16_t *result,
 AuxMono<uint16_t> *aux);
template void guidedFilterMono(uint32_t *origI, uint32_t *p, int r, double eps, uint32_t *result,
 AuxMono<uint32_t> *aux);
template void guidedFilterMono(uint64_t *origI, uint64_t *p, int r, double eps, uint64_t *result,
 AuxMono<uint64_t> *aux);
template void guidedFilterMono(int8_t *origI, int8_t *p, int r, double eps, int8_t *result,
 AuxMono<int8_t> *aux);
template void guidedFilterMono(int16_t *origI, int16_t *p, int r, double eps, int16_t *result,
 AuxMono<int16_t> *aux);
template void guidedFilterMono(int32_t *origI, int32_t *p, int r, double eps, int32_t *result,
 AuxMono<int32_t> *aux);
template void guidedFilterMono(int64_t *origI, int64_t *p, int r, double eps, int64_t *result,
 AuxMono<int64_t> *aux);
template void guidedFilterMono(float *origI, float *p, int r, double eps, float *result,
 AuxMono<float> *aux);
template void guidedFilterMono(double *origI, double *p, int r, double eps, double *result,
 AuxMono<double> *aux);

template void guidedFilterColor(uint8_t *origI, uint8_t *p, int r, double eps, uint8_t *result,
 AuxColor<uint8_t> *aux);
template void guidedFilterColor(uint16_t *origI, uint16_t *p, int r, double eps, uint16_t *result,
 AuxColor<uint16_t> *aux);
template void guidedFilterColor(uint32_t *origI, uint32_t *p, int r, double eps, uint32_t *result,
 AuxColor<uint32_t> *aux);
template void guidedFilterColor(uint64_t *origI, uint64_t *p, int r, double eps, uint64_t *result,
 AuxColor<uint64_t> *aux);
template void guidedFilterColor(int8_t *origI, int8_t *p, int r, double eps, int8_t *result,
 AuxColor<int8_t> *aux);
template void guidedFilterColor(int16_t *origI, int16_t *p, int r, double eps, int16_t *result,
 AuxColor<int16_t> *aux);
template void guidedFilterColor(int32_t *origI, int32_t *p, int r, double eps, int32_t *result,
 AuxColor<int32_t> *aux);
template void guidedFilterColor(int64_t *origI, int64_t *p, int r, double eps, int64_t *result,
 AuxColor<int64_t> *aux);
template void guidedFilterColor(float *origI, float *p, int r, double eps, float *result,
 AuxColor<float> *aux);
template void guidedFilterColor(double *origI, double *p, int r, double eps, double *result,
 AuxColor<double> *aux);
