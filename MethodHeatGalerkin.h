#pragma once

#include "method.h"

class MethodHeatGalerkin : public Method {
public:
    virtual void convertToParam(int, Point, Param &);

    virtual void init();

    virtual void run();

    ~MethodHeatGalerkin();

protected:

    void calcMassMatrix();
    void calcGaussParam();
    double getT(int iCell, Point pt);
    double getQx(int iCell, Point pt);
    double getQy(int iCell, Point pt);
    double getF(int i, int iCell, Point pt);
    double getDfDx(int i, int iCell, Point pt);
    double getDfDy(int i, int iCell, Point pt);
    void calcGradients();
    void calcDiffusionVol();
    void calcDiffusionSurf();
    void calcNewValues();
    void bnd(Edge *e, Point pt, Param p1, Param &p2);
    void flux(Param pl, Param pr, Vector n, double &fT, double &fqx, double &fqy);
    double **T;
    double **qx;
    double **qy;

    double **intT;
    double **intQx;
    double **intQy;

    double ***A;
    double ***invA;

    double **cellGW;
    Point **cellGP;
    double **edgeGW;
    Point **edgeGP;
    double *cellJ;
    double *edgeJ;


    double TMAX;
    double TAU;
    int BASE_FUNC_COUNT;
    int GP_CELL_COUNT;
    int GP_EDGE_COUNT;
    double C11;
};
