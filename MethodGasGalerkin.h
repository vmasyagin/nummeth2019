#pragma once

#include "method.h"

const double	LIMITER_ALFA	= 1.5;

class MethodGasGalerkin : public Method
{
public:
    virtual void convertToParam(int, Point, Param&);
    virtual void init();
    virtual void run();

    ~MethodGasGalerkin();
protected:
    void calcMassMatrix();
    void calcGaussParam();
    double getRO(int iCell, Point pt, DATA_LAYER layer);
    double getRO(int iCell, double x, double y, DATA_LAYER layer);
    double getRU(int iCell, Point pt, DATA_LAYER layer);
    double getRU(int iCell, double x, double y, DATA_LAYER layer);
    double getRV(int iCell, Point pt, DATA_LAYER layer);
    double getRV(int iCell, double x, double y, DATA_LAYER layer);
    double getRE(int iCell, Point pt, DATA_LAYER layer);
    double getRE(int iCell, double x, double y, DATA_LAYER layer);
    double getF(int i, int iCell, Point pt);
    double getF(int i, int iCell, double x, double y);
    double getDFDX(int i, int iCell, Point pt);
    double getDFDY(int i, int iCell, Point pt);
    void bnd(Edge *e, Point pt, Param p1, Param &p2, Vector n);
    void calcConvectionVol();
    void calcConvectionSurf();
    void calcNewValues();
    void calcLimiterCockburn();
    void calcLimiter_II();
    void initValues();
    void initLimiterParameters();
    void copyToOld();
    void getMatrLR(double** L, double** R,  double u, double v, double c);
    int getEdgeByCells(int c1, int c2);
    void flux(Param pL, Param pR, Vector n, double& fr, double& fu, double& fv, double& fe);
    void choiseDirection(int& nn1, int& nn2, double& a1, double& a2, int n0, int n1, int n2, int n3, Point pm, int mm);

    double **ro, **ru, **rv, **re;
    double **roOld, **ruOld, **rvOld, **reOld;
    double **intRO, **intRU, **intRV, **intRE;

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

    /*!
	параметры для лимитера
    */
    double	***	limAlfa;
    int		***	limNeigh;
    double	**	deltaU1;
    double  **	deltaU2;
    double	*	deltaS1;
    double  *	deltaS2;
    double	**	matrL;
    double	**	matrR;
    Vector	**	limLm;
    Vector	**	limLmN;
    Point	**	limPm;

    const double GAM = 1.4;
};


