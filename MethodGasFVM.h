#pragma once

#include "method.h"


class MethodGasFVM : public Method
{
public:
    virtual void convertToParam(int, Point, Param&);
    virtual void init();
    virtual void run();

    ~MethodGasFVM();
protected:

    void initValues();
    void flux(Param pl, Param pr, Vector n, double &fr, double &fu, double &fv, double &fe);
    void bnd(Edge *e, Param p1, Param &p2);

    double *ro, *ru, *rv, *re;
    double *intRO, *intRU, *intRV, *intRE;

    double TMAX;
    double TAU;

    const double GAM = 1.4;
};



