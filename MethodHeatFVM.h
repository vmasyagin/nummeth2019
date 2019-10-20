#pragma once

#include "method.h"


class MethodHeatFVM : public Method
{
public:
    virtual void convertToParam(int, Point, Param&);
    virtual void init();
    virtual void run();

    ~MethodHeatFVM();
protected:

    double *T;
    double *intT;

    double TMAX;
    double TAU;
};



