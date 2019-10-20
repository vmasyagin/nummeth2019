#include "MethodGasGalerkin.h"
#include <cstdio>
#include <cstring>
#include <cmath>

#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))

inline double MINMOD(double a, double b)
{
    if (a*b>0.0)
    {
        if (fabs(a)<fabs(b))
            return a;
        else
            return b;
    }
    else {
        return 0;
    }
}

inline double MINMOD_B(double a, double b)
{
    return MINMOD(a, b);
}

static void inverseMatrix(double** a_src, double **am, int N)
{
    int	*	mask;
    double	fmaxval;
    int		maxind;
    int		tmpi;
    double	tmp;

    double	**a;

    mask = new int[N];
    a = new double*[N];
    for (int i = 0; i < N; i++)
    {
        a[i] = new double[N];
        for (int j = 0; j < N; j++)
        {
            a[i][j] = a_src[i][j];
        }
    }

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if (i == j)
            {
                am[i][j] = 1.0;
            } else {
                am[i][j] = 0.0;
            }
        }
    }
    for (int i = 0; i < N; i++)
    {
        mask[i] = i;
    }
    for (int i = 0; i < N; i++)
    {
        maxind = i;
        fmaxval = fabs(a[i][i]);
        for (int ni = i+1; ni < N; ni++)
        {
            if (fabs(fmaxval) <= fabs(a[ni][i]))
            {
                fmaxval	= fabs(a[ni][i]);
                maxind	= ni;
            }
        }
        fmaxval = a[maxind][i];
        if (fmaxval == 0)
        {
            printf("ERROR! Determinant of mass matrix is zero...\n");
            return;
        }
        if (i != maxind)
        {
            for (int nj = 0; nj < N; nj++)
            {
                tmp				= a[i][nj];
                a[i][nj]		= a[maxind][nj];
                a[maxind][nj]	= tmp;

                tmp				= am[i][nj];
                am[i][nj]		= am[maxind][nj];
                am[maxind][nj]	= tmp;
            }
            tmpi = mask[i];
            mask[i] = mask[maxind];
            mask[maxind] = tmpi;
        }
        double aii = a[i][i];
        for (int j = 0; j < N; j++)
        {
            a[i][j]	=  a[i][j] / aii;
            am[i][j]	= am[i][j] / aii;
        }
        for (int ni = 0; ni < N; ni++)
        {
            if (ni != i)
            {
                double fconst = a[ni][i];
                for (int nj = 0; nj < N; nj++)
                {
                    a[ni][nj]	=  a[ni][nj] - fconst *  a[i][nj];
                    am[ni][nj]	= am[ni][nj] - fconst * am[i][nj];
                }
            }
        }
    }
    for (int i = 0; i < N; i++)
    {
        delete[] a[i];
    }
    delete[] a;
    delete[] mask;
    return;
}

void MethodGasGalerkin::convertToParam(int i, Point pt, Param& p)
{
    p.r = getRO(i, pt, NEW);
    p.u = getRU(i, pt, NEW)/p.r;
    p.v = getRV(i, pt, NEW)/p.r;
    p.e = getRE(i, pt, NEW)/p.r-p.U2()/2.;
    p.p = p.r*p.e*(GAM-1.0);
    p.T = 0.0;
    p.M = p.U2()/sqrt(GAM*p.p/p.r);
}

void MethodGasGalerkin::init()
{
    mesh = new Mesh();
    mesh->initFromFiles((char*)"step.1");
    BASE_FUNC_COUNT = 3;
    GP_CELL_COUNT = 3;
    GP_EDGE_COUNT = 2;

    cellGP = new Point*[mesh->cCount];
    cellGW = new double*[mesh->cCount];
    cellJ = new double[mesh->cCount];
    edgeGP = new Point*[mesh->eCount];
    edgeGW = new double*[mesh->eCount];
    edgeJ = new double[mesh->eCount];

    ro = new double*[mesh->cCount];
    ru = new double*[mesh->cCount];
    rv = new double*[mesh->cCount];
    re = new double*[mesh->cCount];
    roOld = new double*[mesh->cCount];
    ruOld = new double*[mesh->cCount];
    rvOld = new double*[mesh->cCount];
    reOld = new double*[mesh->cCount];
    intRO = new double*[mesh->cCount];
    intRU = new double*[mesh->cCount];
    intRV = new double*[mesh->cCount];
    intRE = new double*[mesh->cCount];
    for(int i = 0; i < mesh->cCount; i++){
        ro[i] = new double[BASE_FUNC_COUNT];
        ru[i] = new double[BASE_FUNC_COUNT];
        rv[i] = new double[BASE_FUNC_COUNT];
        re[i] = new double[BASE_FUNC_COUNT];

        roOld[i] = new double[BASE_FUNC_COUNT];
        ruOld[i] = new double[BASE_FUNC_COUNT];
        rvOld[i] = new double[BASE_FUNC_COUNT];
        reOld[i] = new double[BASE_FUNC_COUNT];

        intRO[i] = new double[BASE_FUNC_COUNT];
        intRU[i] = new double[BASE_FUNC_COUNT];
        intRV[i] = new double[BASE_FUNC_COUNT];
        intRE[i] = new double[BASE_FUNC_COUNT];

        cellGP[i] = new Point[GP_CELL_COUNT];
        cellGW[i] = new double[GP_CELL_COUNT];
    }

    initValues();
    saveVTK(0);

    TMAX    = 4.0;
    TAU     = 1.e-4;

    for(int i = 0; i < mesh->eCount; i++){
        edgeGP[i] = new Point[GP_EDGE_COUNT];
        edgeGW[i] = new double[GP_EDGE_COUNT];
    }

    calcGaussParam();
    A = new double**[mesh->cCount];
    invA = new double**[mesh->cCount];
    for (int i = 0; i < mesh->cCount; i++){
        A[i] = new double*[BASE_FUNC_COUNT];
        invA[i] = new double*[BASE_FUNC_COUNT];
        for(int j = 0; j < BASE_FUNC_COUNT; j++){
            A[i][j] = new double[BASE_FUNC_COUNT];
            invA[i][j] = new double[BASE_FUNC_COUNT];
        }
    }
    calcMassMatrix();

    initLimiterParameters();
};

void MethodGasGalerkin::calcMassMatrix() {
    for (int i = 0; i < mesh->cCount; i++) {
        double **mA = A[i];
        double **invA_ = invA[i];
        for (int i = 0; i < BASE_FUNC_COUNT; i++) {
            for (int j = 0; j < BASE_FUNC_COUNT; j++) {
                mA[i][j] = 0.0;
                for (int iGP = 0; iGP < GP_CELL_COUNT; iGP++) {
                    mA[i][j] += cellGW[i][iGP] * getF(i, i, cellGP[i][iGP])
                                * getF(j, i, cellGP[i][iGP]);
                }
                mA[i][j] *= cellJ[i];
            }
        }
        inverseMatrix(mA, invA_, BASE_FUNC_COUNT);
    }
}

void MethodGasGalerkin::calcGaussParam() {
    // для ячеек
    if (GP_CELL_COUNT == 4) {
        for (int i = 0; i < mesh->cCount; i++) {

            double a = 1.0 / 3.0;
            double b = 1.0 / 5.0;
            double c = 3.0 / 5.0;
            double x1 = mesh->nodes[mesh->cells[i].nodesInd[0]].x;
            double y1 = mesh->nodes[mesh->cells[i].nodesInd[0]].y;
            double x2 = mesh->nodes[mesh->cells[i].nodesInd[1]].x;
            double y2 = mesh->nodes[mesh->cells[i].nodesInd[1]].y;
            double x3 = mesh->nodes[mesh->cells[i].nodesInd[2]].x;
            double y3 = mesh->nodes[mesh->cells[i].nodesInd[2]].y;
            double a1 = x1 - x3;
            double a2 = y1 - y3;
            double b1 = x2 - x3;
            double b2 = y2 - y3;
            double c1 = x3;
            double c2 = y3;

            cellGW[i][0] = -27.0 / 48.0;
            cellGP[i][0].x = a1 * a + b1 * a + c1;
            cellGP[i][0].y = a2 * a + b2 * a + c2;

            cellGW[i][1] = 25.0 / 48.0;
            cellGP[i][1].x = a1 * c + b1 * b + c1;
            cellGP[i][1].y = a2 * c + b2 * b + c2;

            cellGW[i][2] = 25.0 / 48.0;
            cellGP[i][2].x = a1 * b + b1 * c + c1;
            cellGP[i][2].y = a2 * b + b2 * c + c2;

            cellGW[i][3] = 25.0 / 48.0;
            cellGP[i][3].x = a1 * b + b1 * b + c1;
            cellGP[i][3].y = a2 * b + b2 * b + c2;

            cellJ[i] = 0.5 * fabs(a1 * b2 - a2 * b1);
        }
    } else if (GP_CELL_COUNT == 3) {
        for (int i = 0; i < mesh->cCount; i++) {
            double a = 1.0 / 6.0;
            double b = 2.0 / 3.0;
            double x1 = mesh->nodes[mesh->cells[i].nodesInd[0]].x;
            double y1 = mesh->nodes[mesh->cells[i].nodesInd[0]].y;
            double x2 = mesh->nodes[mesh->cells[i].nodesInd[1]].x;
            double y2 = mesh->nodes[mesh->cells[i].nodesInd[1]].y;
            double x3 = mesh->nodes[mesh->cells[i].nodesInd[2]].x;
            double y3 = mesh->nodes[mesh->cells[i].nodesInd[2]].y;
            double a1 = x1 - x3;
            double a2 = y1 - y3;
            double b1 = x2 - x3;
            double b2 = y2 - y3;
            double c1 = x3;
            double c2 = y3;

            cellGW[i][0] = 1.0 / 6.0;
            cellGP[i][0].x = a1 * a + b1 * a + c1;
            cellGP[i][0].y = a2 * a + b2 * a + c2;

            cellGW[i][1] = 1.0 / 6.0;
            cellGP[i][1].x = a1 * a + b1 * b + c1;
            cellGP[i][1].y = a2 * a + b2 * b + c2;

            cellGW[i][2] = 1.0 / 6.0;
            cellGP[i][2].x = a1 * b + b1 * a + c1;
            cellGP[i][2].y = a2 * b + b2 * a + c2;

            cellJ[i] = fabs(a1 * b2 - a2 * b1);

        }
    }

    // для ребер
    if (GP_EDGE_COUNT == 3) {
        for (int i = 0; i < mesh->eCount; i++) {
            double gp1 = -3.0 / 5.0;
            double gp2 = 0.0;
            double gp3 = 3.0 / 5.0;
            double x1 = mesh->nodes[mesh->edges[i].n1].x;
            double y1 = mesh->nodes[mesh->edges[i].n1].y;
            double x2 = mesh->nodes[mesh->edges[i].n2].x;
            double y2 = mesh->nodes[mesh->edges[i].n2].y;

            edgeGW[i][0] = 5.0 / 9.0;
            edgeGP[i][0].x = ((x1 + x2) + gp1 * (x2 - x1)) / 2.0;
            edgeGP[i][0].y = ((y1 + y2) + gp1 * (y2 - y1)) / 2.0;

            edgeGW[i][1] = 8.0 / 9.0;
            edgeGP[i][1].x = (x1 + x2) / 2.0;
            edgeGP[i][1].y = (y1 + y2) / 2.0;

            edgeGW[i][2] = 5.0 / 9.0;
            edgeGP[i][2].x = ((x1 + x2) + gp3 * (x2 - x1)) / 2.0;
            edgeGP[i][2].y = ((y1 + y2) + gp3 * (y2 - y1)) / 2.0;

            edgeJ[i] = sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1)) * 0.5;

        }
    } else if (GP_EDGE_COUNT == 2) {
        for (int i = 0; i < mesh->eCount; i++) {
            double _sqrt3 = 1.0 / sqrt(3.0);
            double x1 = mesh->nodes[mesh->edges[i].n1].x;
            double y1 = mesh->nodes[mesh->edges[i].n1].y;
            double x2 = mesh->nodes[mesh->edges[i].n2].x;
            double y2 = mesh->nodes[mesh->edges[i].n2].y;

            edgeGW[i][0] = 1.0;
            edgeGP[i][0].x = ((x1 + x2) - _sqrt3 * (x2 - x1)) * 0.5;
            edgeGP[i][0].y = ((y1 + y2) - _sqrt3 * (y2 - y1)) * 0.5;

            edgeGW[i][1] = 1.0;
            edgeGP[i][1].x = ((x1 + x2) + _sqrt3 * (x2 - x1)) * 0.5;
            edgeGP[i][1].y = ((y1 + y2) + _sqrt3 * (y2 - y1)) * 0.5;

            edgeJ[i] = sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1)) * 0.5;

        }
    }
}


void MethodGasGalerkin::initValues()
{
    for (int i = 0; i < mesh->cCount; i++) {
        ro[i][0] = 1.;
        ru[i][0] = 3.;
        rv[i][0] = 0.0;
        re[i][0] = (1./1.4)/(GAM-1.0)+0.5*(ru[i][0]*ru[i][0])/ro[i][0];
        for (int j = 1; j < BASE_FUNC_COUNT; j++){
            ro[i][j] = 0.0;
            ru[i][j] = 0.0;
            rv[i][j] = 0.0;
            re[i][j] = 0.0;
        }
    }
}


void MethodGasGalerkin::run()
{
    double t = 0.0;
    int step = 0;

    copyToOld();
    calcLimiterCockburn();
    calcLimiter_II();

    while (t < TMAX) {
        t += TAU;
        step++;
        // обнуляем правые части системы
        for (int i = 0; i < mesh->cCount; i++) {
            memset(intRO[i], 0, sizeof(double) * BASE_FUNC_COUNT);
            memset(intRU[i], 0, sizeof(double) * BASE_FUNC_COUNT);
            memset(intRV[i], 0, sizeof(double) * BASE_FUNC_COUNT);
            memset(intRE[i], 0, sizeof(double) * BASE_FUNC_COUNT);
        }

        calcConvectionVol();
        calcConvectionSurf();

        calcNewValues();
        copyToOld();

        calcLimiterCockburn();
        calcLimiter_II();

        copyToOld();

        if (step % 1000 == 0)
        {
            saveVTK(step);
            printf("Calculation results for step %d are saved.\n", step);
        }
    }
}


void MethodGasGalerkin::calcConvectionVol() {
    double gpRO,  gpRU,  gpRV,  gpRE;
    double gpFRx, gpFUx, gpFVx, gpFEx;
    double gpFRy, gpFUy, gpFVy, gpFEy;
    double tmpU, tmpV, tmpP, tmpE, tmpe;
    double intROtmp, intRUtmp, intRVtmp, intREtmp;
    for (int iCell = 0; iCell < mesh->cCount; iCell++)
    {
        for (int iF = 0; iF < BASE_FUNC_COUNT; iF++) {
            intROtmp = 0.0;
            intRUtmp = 0.0;
            intRVtmp = 0.0;
            intREtmp = 0.0;
            for (int iGP = 0; iGP < GP_CELL_COUNT; iGP++)
            {
                gpRO = getRO(iCell, cellGP[iCell][iGP], NEW);
                gpRU = getRU(iCell, cellGP[iCell][iGP], NEW);
                gpRV = getRV(iCell, cellGP[iCell][iGP], NEW);
                gpRE = getRE(iCell, cellGP[iCell][iGP], NEW);

                tmpU = gpRU/gpRO;
                tmpV = gpRV/gpRO;
                tmpE = gpRE/gpRO;
                tmpe = tmpE-(tmpU*tmpU+tmpV*tmpV)*0.5;
                tmpP = tmpe*gpRO*(GAM-1);

                gpFRx = gpRU;
                gpFUx = gpRU*tmpU+tmpP;
                gpFVx = gpRV*tmpU;
                gpFEx = (gpRE+tmpP)*tmpU;

                gpFRy = gpRV;
                gpFUy = gpRU*tmpV;
                gpFVy = gpRV*tmpV+tmpP;
                gpFEy = (gpRE+tmpP)*tmpV;

                double bfDFDx = getDFDX(iF, iCell, cellGP[iCell][iGP]);
                double bfDFDy = getDFDY(iF, iCell, cellGP[iCell][iGP]);

                intROtmp += cellGW[iCell][iGP]*( gpFRx*bfDFDx + gpFRy*bfDFDy );
                intRUtmp += cellGW[iCell][iGP]*( gpFUx*bfDFDx + gpFUy*bfDFDy );
                intRVtmp += cellGW[iCell][iGP]*( gpFVx*bfDFDx + gpFVy*bfDFDy );
                intREtmp += cellGW[iCell][iGP]*( gpFEx*bfDFDx + gpFEy*bfDFDy );
            }

            if (fabs(intROtmp) <= EPS) intROtmp = 0.0;
            if (fabs(intRUtmp) <= EPS) intRUtmp = 0.0;
            if (fabs(intRVtmp) <= EPS) intRVtmp = 0.0;
            if (fabs(intREtmp) <= EPS) intREtmp = 0.0;

            intRO[iCell][iF] += intROtmp*cellJ[iCell];
            intRU[iCell][iF] += intRUtmp*cellJ[iCell];
            intRV[iCell][iF] += intRVtmp*cellJ[iCell];
            intRE[iCell][iF] += intREtmp*cellJ[iCell];
        }
    }

}

void MethodGasGalerkin::calcConvectionSurf() {
    double fr,  fu,  fv,  fe;	// потоки на границе ячейки
    int c1, c2;
    for (int i = 0; i < mesh->eCount; i++) {
        Edge &e = mesh->edges[i];
        c1 = e.c1;
        c2 = e.c2;
        Vector n = e.n;

        for (int iF = 0;
             iF < BASE_FUNC_COUNT; iF++) // каждый случай при скалярном умножении уравнений на baseF( iF, ... )
        {
            double tmpintRO1 = 0.0, tmpintRU1 = 0.0, tmpintRV1 = 0.0, tmpintRE1 = 0.0; // для ячейки С1
            double tmpintRO2 = 0.0, tmpintRU2 = 0.0, tmpintRV2 = 0.0, tmpintRE2 = 0.0; // для ячейки С2
            for (int iGP = 0; iGP < GP_EDGE_COUNT; iGP++) {
                Point &pt = edgeGP[i][iGP];

                Param p1, p2;
                convertToParam(c1, pt, p1);

                if (c2 >= 0) {
                    convertToParam(c2, pt, p2);
                } else {
                    bnd(&e, pt, p1, p2, n);
                }

                flux(p1, p2, n, fr, fu, fv, fe);

                double cGP1 = edgeGW[i][iGP] * getF(iF, c1, edgeGP[i][iGP].x, edgeGP[i][iGP].y);
                tmpintRO1 += fr * cGP1;
                tmpintRU1 += fu * cGP1;
                tmpintRV1 += fv * cGP1;
                tmpintRE1 += fe * cGP1;

                if (c2 >= 0) {
                    double cGP2 = edgeGW[i][iGP] * getF(iF, c2, edgeGP[i][iGP].x, edgeGP[i][iGP].y);
                    tmpintRO2 += fr * cGP2;
                    tmpintRU2 += fu * cGP2;
                    tmpintRV2 += fv * cGP2;
                    tmpintRE2 += fe * cGP2;
                }
            }

            if (fabs(tmpintRO1) <= EPS) tmpintRO1 = 0.0;
            if (fabs(tmpintRU1) <= EPS) tmpintRU1 = 0.0;
            if (fabs(tmpintRV1) <= EPS) tmpintRV1 = 0.0;
            if (fabs(tmpintRE1) <= EPS) tmpintRE1 = 0.0;

            intRO[c1][iF] -= tmpintRO1 * edgeJ[i];
            intRU[c1][iF] -= tmpintRU1 * edgeJ[i];
            intRV[c1][iF] -= tmpintRV1 * edgeJ[i];
            intRE[c1][iF] -= tmpintRE1 * edgeJ[i];

            if (c2 >= 0) {
                if (fabs(tmpintRO2) <= EPS) tmpintRO2 = 0.0;
                if (fabs(tmpintRU2) <= EPS) tmpintRU2 = 0.0;
                if (fabs(tmpintRV2) <= EPS) tmpintRV2 = 0.0;
                if (fabs(tmpintRE2) <= EPS) tmpintRE2 = 0.0;

                intRO[c2][iF] += tmpintRO2 * edgeJ[i];
                intRU[c2][iF] += tmpintRU2 * edgeJ[i];
                intRV[c2][iF] += tmpintRV2 * edgeJ[i];
                intRE[c2][iF] += tmpintRE2 * edgeJ[i];
            }
        }
    }
}

double MethodGasGalerkin::getF(int i, int iCell, Point pt)
{
    if (BASE_FUNC_COUNT == 3)
    {
        switch (i) {
            case 0:
                return 1.0;
            case 1:
                return (pt.x - mesh->cells[iCell].c.x) / mesh->cells[iCell].HX;
            case 2:
                return (pt.y - mesh->cells[iCell].c.y) / mesh->cells[iCell].HY;
        }

    }
    else if (BASE_FUNC_COUNT == 6)
    {
        switch(i) {
            case 0:
                return 1.0;
            case 1:
                return (pt.x - mesh->cells[iCell].c.x) / mesh->cells[iCell].HX;
            case 2:
                return (pt.y - mesh->cells[iCell].c.y) / mesh->cells[iCell].HY;
            case 3:
                return (pt.x - mesh->cells[iCell].c.x) * (pt.x - mesh->cells[iCell].c.x) /
                       (mesh->cells[iCell].HX * mesh->cells[iCell].HX);
            case 4:
                return (pt.y - mesh->cells[iCell].c.y) * (pt.y - mesh->cells[iCell].c.y) /
                       (mesh->cells[iCell].HY * mesh->cells[iCell].HY);
            case 5:
                return (pt.x - mesh->cells[iCell].c.x) * (pt.y - mesh->cells[iCell].c.y) /
                       (mesh->cells[iCell].HX * mesh->cells[iCell].HY);
        }
    }
    else {
        printf("ERROR: wrong basic functions count!\n");
        exit(1);
    }
}

double MethodGasGalerkin::getF(int i, int iCell, double x, double y)
{
    if (BASE_FUNC_COUNT == 3)
    {
        switch (i) {
            case 0:
                return 1.0;
            case 1:
                return (x - mesh->cells[iCell].c.x) / mesh->cells[iCell].HX;
            case 2:
                return (y - mesh->cells[iCell].c.y) / mesh->cells[iCell].HY;
        }

    }
    else if (BASE_FUNC_COUNT == 6)
    {
        switch(i) {
            case 0:
                return 1.0;
            case 1:
                return (x - mesh->cells[iCell].c.x) / mesh->cells[iCell].HX;
            case 2:
                return (y - mesh->cells[iCell].c.y) / mesh->cells[iCell].HY;
            case 3:
                return (x - mesh->cells[iCell].c.x) * (x - mesh->cells[iCell].c.x) /
                       (mesh->cells[iCell].HX * mesh->cells[iCell].HX);
            case 4:
                return (y - mesh->cells[iCell].c.y) * (y - mesh->cells[iCell].c.y) /
                       (mesh->cells[iCell].HY * mesh->cells[iCell].HY);
            case 5:
                return (x - mesh->cells[iCell].c.x) * (y - mesh->cells[iCell].c.y) /
                       (mesh->cells[iCell].HX * mesh->cells[iCell].HY);
        }
    }
    else {
        printf("ERROR: wrong basic functions count!\n");
        exit(1);
    }
}

double MethodGasGalerkin::getDFDX(int i, int iCell, Point pt)
{
    if (BASE_FUNC_COUNT == 3)
    {
        switch (i) {
            case 0:
                return 0.0;
            case 1:
                return 1.0 / mesh->cells[iCell].HX;
            case 2:
                return 0.0;
        }

    }
    else if (BASE_FUNC_COUNT == 6)
    {
        switch(i) {
            case 0:
                return 0.0;
            case 1:
                return 1.0 / mesh->cells[iCell].HX;
            case 2:
                return 0.0;
            case 3:
                return 2.0*(pt.x - mesh->cells[iCell].c.x)/(mesh->cells[iCell].HX*mesh->cells[iCell].HX);
            case 4:
                return 0.0;
            case 5:
                return (pt.y - mesh->cells[iCell].c.y)/(mesh->cells[iCell].HX*mesh->cells[iCell].HY);
        }
    }
    else {
        printf("ERROR: wrong basic functions count!\n");
        exit(1);
    }
}

double MethodGasGalerkin::getDFDY(int i, int iCell, Point pt)
{
    if (BASE_FUNC_COUNT == 3)
    {
        switch (i) {
            case 0:
                return 0.0;
            case 1:
                return 0.0;
            case 2:
                return 1.0/mesh->cells[iCell].HY;
        }

    }
    else if (BASE_FUNC_COUNT == 6)
    {
        switch(i) {
            case 0:
                return 0.0;
            case 1:
                return 1.0 / mesh->cells[iCell].HX;
            case 2:
                return 0.0;
            case 3:
                return 0.0;
            case 4:
                return 2.0*(pt.y - mesh->cells[iCell].c.y)/(mesh->cells[iCell].HY*mesh->cells[iCell].HY);
            case 5:
                return (pt.x - mesh->cells[iCell].c.x)/(mesh->cells[iCell].HX*mesh->cells[iCell].HY);
        }
    }
    else {
        printf("ERROR: wrong basic functions count!\n");
        exit(1);
    }
}

double MethodGasGalerkin::getRO(int iCell, Point pt, DATA_LAYER layer) {
    double RO, *fieldRO;
    switch (layer)
    {
        case OLD:
            fieldRO = roOld[iCell];
            break;
        case NEW:
        default:
            fieldRO = ro[iCell];
            break;
    }
    RO = fieldRO[0];
    for (int j = 1; j < BASE_FUNC_COUNT; j++)
    {
        double bF = getF(j, iCell, pt);
        RO += fieldRO[j]*bF;
    }
    return RO;
}

double MethodGasGalerkin::getRO(int iCell, double x, double y, DATA_LAYER layer) {
    double RO, *fieldRO;
    switch (layer)
    {
        case OLD:
            fieldRO = roOld[iCell];
            break;
        case NEW:
        default:
            fieldRO = ro[iCell];
            break;
    }
    RO = fieldRO[0];
    for (int j = 1; j < BASE_FUNC_COUNT; j++)
    {
        double bF = getF(j, iCell, x, y);
        RO += fieldRO[j]*bF;
    }
    return RO;
}

double MethodGasGalerkin::getRU(int iCell, Point pt, DATA_LAYER layer) {
    double RU, *fieldRU;
    switch (layer)
    {
        case OLD:
            fieldRU = ruOld[iCell];
            break;
        case NEW:
        default:
            fieldRU = ru[iCell];
            break;
    }
    RU = fieldRU[0];
    for (int j = 1; j < BASE_FUNC_COUNT; j++)
    {
        double bF = getF(j, iCell, pt);
        RU += fieldRU[j]*bF;
    }
    return RU;
}

double MethodGasGalerkin::getRU(int iCell, double x, double y, DATA_LAYER layer) {
    double RU, *fieldRU;
    switch (layer)
    {
        case OLD:
            fieldRU = ruOld[iCell];
            break;
        case NEW:
        default:
            fieldRU = ru[iCell];
            break;
    }
    RU = fieldRU[0];
    for (int j = 1; j < BASE_FUNC_COUNT; j++)
    {
        double bF = getF(j, iCell, x, y);
        RU += fieldRU[j]*bF;
    }
    return RU;
}

double MethodGasGalerkin::getRV(int iCell, Point pt, DATA_LAYER layer) {
    double RV, *fieldRV;
    switch (layer)
    {
        case OLD:
            fieldRV = rvOld[iCell];
            break;
        case NEW:
        default:
            fieldRV = rv[iCell];
            break;
    }
    RV = fieldRV[0];
    for (int j = 1; j < BASE_FUNC_COUNT; j++)
    {
        double bF = getF(j, iCell, pt);
        RV += fieldRV[j]*bF;
    }
    return RV;
}

double MethodGasGalerkin::getRV(int iCell, double x, double y, DATA_LAYER layer) {
    double RV, *fieldRV;
    switch (layer)
    {
        case OLD:
            fieldRV = rvOld[iCell];
            break;
        case NEW:
        default:
            fieldRV = rv[iCell];
            break;
    }
    RV = fieldRV[0];
    for (int j = 1; j < BASE_FUNC_COUNT; j++)
    {
        double bF = getF(j, iCell, x, y);
        RV += fieldRV[j]*bF;
    }
    return RV;
}

double MethodGasGalerkin::getRE(int iCell, Point pt, DATA_LAYER layer) {
    double RE, *fieldRE;
    switch (layer)
    {
        case OLD:
            fieldRE = reOld[iCell];
            break;
        case NEW:
        default:
            fieldRE = re[iCell];
            break;
    }
    RE = fieldRE[0];
    for (int j = 1; j < BASE_FUNC_COUNT; j++)
    {
        double bF = getF(j, iCell, pt);
        RE += fieldRE[j]*bF;
    }
    return RE;
}

double MethodGasGalerkin::getRE(int iCell, double x, double y, DATA_LAYER layer) {
    double RE, *fieldRE;
    switch (layer)
    {
        case OLD:
            fieldRE = reOld[iCell];
            break;
        case NEW:
        default:
            fieldRE = re[iCell];
            break;
    }
    RE = fieldRE[0];
    for (int j = 1; j < BASE_FUNC_COUNT; j++)
    {
        double bF = getF(j, iCell, x, y);
        RE += fieldRE[j]*bF;
    }
    return RE;
}

MethodGasGalerkin::~MethodGasGalerkin()
{
    for (int i = 0; i < mesh->cCount; i++) {
        delete[] ro[i]; delete[] ru[i]; delete[] rv[i]; delete[] re[i];
        delete[] roOld[i]; delete[] ruOld[i]; delete[] rvOld[i]; delete[] reOld[i];
        delete[] intRO[i]; delete[] intRU[i]; delete[] intRV[i]; delete[] intRE[i];
        delete[] cellGP[i], cellGW[i];

        delete[] limAlfa[i]; delete[] limNeigh[i]; delete[] limLm[i]; delete[] limLmN[i];
        delete[] limPm[i]; delete[] matrL[i]; delete[] matrR[i];

        for(int j = 0; j < BASE_FUNC_COUNT; j++){
            for(int k = 0; k < BASE_FUNC_COUNT; k++){
                delete[] A[j][k];
                delete[] invA[j][k];
            }
            delete[] A[j];
            delete[] invA[j];
        }
    }

    for (int i = 0; i < mesh->eCount; i++) {
        delete[] edgeGP[i], edgeGW[i];
    }

    delete[] ro, ru, rv, re;
    delete[] roOld, ruOld, rvOld, reOld;
    delete[] intRO, intRU, intRV, intRE;
    delete[] cellGP, cellGW, cellJ;
    delete[] edgeGP, edgeGW, edgeJ;
    delete[] A, invA;

    delete[] limAlfa, limNeigh, limLm, limLmN, limPm;
    delete[] deltaS1, deltaS2, deltaU1, deltaU2, matrL, matrR;
    delete mesh;
}

void MethodGasGalerkin::calcLimiterCockburn() {
    for (int iCell = 0; iCell < mesh->cCount; iCell++) // цикл по ячейкам
    {
        int n0 = mesh->cells[iCell].neigh[0];
        int n1 = mesh->cells[iCell].neigh[1];
        int n2 = mesh->cells[iCell].neigh[2];
        if ( (n0 < 0) || (n1 < 0) || (n2 < 0) )
        {
            ro[iCell][1] = 0.0;
            ru[iCell][1] = 0.0;
            rv[iCell][1] = 0.0;
            re[iCell][1] = 0.0;

            ro[iCell][2] = 0.0;
            ru[iCell][2] = 0.0;
            rv[iCell][2] = 0.0;
            re[iCell][3] = 0.0;

            continue;
        }

        double ROc, RUc, RVc, REc;
        ROc = getRO(iCell, mesh->cells[iCell].c, OLD);
        RUc = getRU(iCell, mesh->cells[iCell].c, OLD);
        RVc = getRV(iCell, mesh->cells[iCell].c, OLD);
        REc = getRE(iCell, mesh->cells[iCell].c, OLD);

        for (int m = 0; m < 3; m++)
        {
            double ROm, RUm, RVm, REm;
            double RO1, RU1, RV1, RE1;
            double RO2, RU2, RV2, RE2;

            ROm = getRO(iCell, limPm[iCell][m], OLD);
            RUm = getRU(iCell, limPm[iCell][m], OLD);
            RVm = getRV(iCell, limPm[iCell][m], OLD);
            REm = getRE(iCell, limPm[iCell][m], OLD);

            RO1 = getRO(limNeigh[iCell][m][0], mesh->cells[limNeigh[iCell][m][0]].c, OLD);
            RU1 = getRU(limNeigh[iCell][m][0], mesh->cells[limNeigh[iCell][m][0]].c, OLD);
            RV1 = getRV(limNeigh[iCell][m][0], mesh->cells[limNeigh[iCell][m][0]].c, OLD);
            RE1 = getRE(limNeigh[iCell][m][0], mesh->cells[limNeigh[iCell][m][0]].c, OLD);

            RO2 = getRO(limNeigh[iCell][m][1], mesh->cells[limNeigh[iCell][m][1]].c, OLD);
            RU2 = getRU(limNeigh[iCell][m][1], mesh->cells[limNeigh[iCell][m][1]].c, OLD);
            RV2 = getRV(limNeigh[iCell][m][1], mesh->cells[limNeigh[iCell][m][1]].c, OLD);
            RE2 = getRE(limNeigh[iCell][m][1], mesh->cells[limNeigh[iCell][m][1]].c, OLD);

            deltaU1[m][0] = ROm-ROc;
            deltaU1[m][1] = RUm-RUc;
            deltaU1[m][2] = RVm-RVc;
            deltaU1[m][3] = REm-REc;

            deltaU2[m][0] = limAlfa[iCell][m][0]*(RO1-ROc)+limAlfa[iCell][m][1]*(RO2-ROc);
            deltaU2[m][1] = limAlfa[iCell][m][0]*(RU1-RUc)+limAlfa[iCell][m][1]*(RU2-RUc);
            deltaU2[m][2] = limAlfa[iCell][m][0]*(RV1-RVc)+limAlfa[iCell][m][1]*(RV2-RVc);
            deltaU2[m][3] = limAlfa[iCell][m][0]*(RE1-REc)+limAlfa[iCell][m][1]*(RE2-REc);

            double un, ut, c, utmp, vtmp;
            double nx = limLmN[iCell][m].x;
            double ny = limLmN[iCell][m].y;
            utmp = RUc/ROc;
            vtmp = RVc/ROc;
            un = utmp*nx+vtmp*ny;
            ut = utmp*ny-vtmp*nx;
            c  = sqrt(GAM*(GAM - 1.0)*(REc/ROc-(utmp*utmp+vtmp*vtmp)/2.0));

            utmp = deltaU1[m][1]*nx+deltaU1[m][2]*ny;
            vtmp = deltaU1[m][1]*ny-deltaU1[m][2]*nx;
            deltaU1[m][1] = utmp;
            deltaU1[m][2] = vtmp;

            utmp = deltaU2[m][1]*nx+deltaU2[m][2]*ny;
            vtmp = deltaU2[m][1]*ny-deltaU2[m][2]*nx;
            deltaU2[m][1] = utmp;
            deltaU2[m][2] = vtmp;

            getMatrLR(matrL, matrR, un,ut,c);
            memset(deltaS1, 0, 4*sizeof(double));
            memset(deltaS2, 0, 4*sizeof(double));
            for (int k = 0; k < 4; k++)
            {
                deltaS1[0] += matrL[0][k]*deltaU1[m][k];
                deltaS1[1] += matrL[1][k]*deltaU1[m][k];
                deltaS1[2] += matrL[2][k]*deltaU1[m][k];
                deltaS1[3] += matrL[3][k]*deltaU1[m][k];

                deltaS2[0] += matrL[0][k]*deltaU2[m][k];
                deltaS2[1] += matrL[1][k]*deltaU2[m][k];
                deltaS2[2] += matrL[2][k]*deltaU2[m][k];
                deltaS2[3] += matrL[3][k]*deltaU2[m][k];
            }

            deltaS1[0] = MINMOD_B( deltaS1[0], LIMITER_ALFA*deltaS2[0] );
            deltaS1[1] = MINMOD_B( deltaS1[1], LIMITER_ALFA*deltaS2[1] );
            deltaS1[2] = MINMOD_B( deltaS1[2], LIMITER_ALFA*deltaS2[2] );
            deltaS1[3] = MINMOD_B( deltaS1[3], LIMITER_ALFA*deltaS2[3] );

            memset(deltaU1[m], 0, 4*sizeof(double));
            for (int k = 0; k < 4; k++)
            {
                deltaU1[m][0] += matrR[0][k]*deltaS1[k];
                deltaU1[m][1] += matrR[1][k]*deltaS1[k];
                deltaU1[m][2] += matrR[2][k]*deltaS1[k];
                deltaU1[m][3] += matrR[3][k]*deltaS1[k];
            }

            utmp = deltaU1[m][1]*nx+deltaU1[m][2]*ny;
            vtmp = deltaU1[m][1]*ny-deltaU1[m][2]*nx;
            deltaU1[m][1] = utmp;
            deltaU1[m][2] = vtmp;
        }

        double pos, neg;
        for (int k = 0; k < 4; k++)
        {

            double & d1n = deltaU1[0][k];
            double & d2n = deltaU1[1][k];
            double & d3n = deltaU1[2][k];

            if (fabs(deltaU1[0][k]+deltaU1[1][k]+deltaU1[2][k]) > EPS)
            {
                pos = 0; neg = 0;
                for (int m = 0; m < 3; m++)
                {
                    pos += MAX(0.0,  deltaU1[m][k]);
                    neg += MAX(0.0, -deltaU1[m][k]);
                }
                double thetaP, thetaM;
                if (fabs(pos)<EPS)
                {
                    thetaP = 1.0;
                    thetaM = 0.0;
                }
                else if (fabs(neg)<EPS)
                {
                    thetaP = 0.0;
                    thetaM = 1.0;
                }
                else
                {
                    thetaP = MIN(1.0, neg/pos);
                    thetaM = MIN(1.0, pos/neg);
                }
                deltaU1[0][k] = thetaP * MAX( 0.0,  deltaU1[0][k] ) - thetaM * MAX( 0.0, -deltaU1[0][k] );
                deltaU1[1][k] = thetaP * MAX( 0.0,  deltaU1[1][k] ) - thetaM * MAX( 0.0, -deltaU1[1][k] );
                deltaU1[2][k] = thetaP * MAX( 0.0,  deltaU1[2][k] ) - thetaM * MAX( 0.0, -deltaU1[2][k] );
            }

            double& xb0 = mesh->cells[iCell].c.x;
            double& yb0 = mesh->cells[iCell].c.y;
            double& xm1 = limPm[iCell][0].x;
            double& ym1 = limPm[iCell][0].y;
            double& xm2 = limPm[iCell][1].x;
            double& ym2 = limPm[iCell][1].y;
            double& xm3 = limPm[iCell][2].x;
            double& ym3 = limPm[iCell][2].y;

            double det = (xm2*ym3-xm3*ym2)-(xm1*ym3-xm3*ym1)+(xm1*ym2-xm2*ym1);

            double a1 = (ym2-ym3)/det;
            double b1 = (xm3-xm2)/det;

            double a2 = (ym3-ym1)/det;
            double b2 = (xm1-xm3)/det;

            double a3 = (ym1-ym2)/det;
            double b3 = (xm2-xm1)/det;

            double *pRR, *pRC;
            switch (k)
            {
                case 0:
                    pRR = ro[iCell];
                    pRC = &ROc;
                    break;
                case 1:
                    pRR = ru[iCell];
                    pRC = &RUc;
                    break;
                case 2:
                    pRR = rv[iCell];
                    pRC = &RVc;
                    break;
                case 3:
                    pRR = re[iCell];
                    pRC = &REc;
                    break;
            }

            pRR[1] = (a1*d1n+a2*d2n+a3*d3n)*mesh->cells[iCell].HX;
            pRR[2] = (b1*d1n+b2*d2n+b3*d3n)*mesh->cells[iCell].HY;

            if (fabs(pRR[1]) < EPS) pRR[1] = 0.0;
            if (fabs(pRR[2]) < EPS) pRR[2] = 0.0;
        }

    }
}

void MethodGasGalerkin::calcNewValues() {
    double *aRO, *aRU, *aRV, *aRE;
    aRO = (double*)malloc(sizeof(double)*BASE_FUNC_COUNT);
    aRU = (double*)malloc(sizeof(double)*BASE_FUNC_COUNT);
    aRV = (double*)malloc(sizeof(double)*BASE_FUNC_COUNT);
    aRE = (double*)malloc(sizeof(double)*BASE_FUNC_COUNT);
    for (int iCell = 0; iCell < mesh->cCount; iCell++)
    {
        double ** A = invA[iCell];
        memset(aRO, 0, sizeof(double)*BASE_FUNC_COUNT);
        memset(aRU, 0, sizeof(double)*BASE_FUNC_COUNT);
        memset(aRV, 0, sizeof(double)*BASE_FUNC_COUNT);
        memset(aRE, 0, sizeof(double)*BASE_FUNC_COUNT);
        for (int j = 0; j < BASE_FUNC_COUNT; j++)
        {
            for (int k = 0; k < BASE_FUNC_COUNT; k++)
            {
                aRO[j] += A[j][k] * intRO[iCell][k];
                aRU[j] += A[j][k] * intRU[iCell][k];
                aRV[j] += A[j][k] * intRV[iCell][k];
                aRE[j] += A[j][k] * intRE[iCell][k];
            }
        }

        for (int j = 0; j < BASE_FUNC_COUNT; j++)
        {
            ro[iCell][j] += TAU * aRO[j];
            ru[iCell][j] += TAU * aRU[j];
            rv[iCell][j] += TAU * aRV[j];
            re[iCell][j] += TAU * aRE[j];

            if (fabs(ro[iCell][j]) <= EPS) ro[iCell][j] = 0.0;
            if (fabs(ru[iCell][j]) <= EPS) ru[iCell][j] = 0.0;
            if (fabs(rv[iCell][j]) <= EPS) rv[iCell][j] = 0.0;
            if (fabs(re[iCell][j]) <= EPS) re[iCell][j] = 0.0;

        }
    }

    free( aRO );
    free( aRU );
    free( aRV );
    free( aRE );
}

void MethodGasGalerkin::calcLimiter_II() {
    double xtt[9],ytt[9],t[3];
    double Rcur1, RUcur1, RVcur1, Ecur1;
    double Rcur2, RUcur2, RVcur2, Ecur2;
    double Rcur3, RUcur3, RVcur3, Ecur3;
    double sR, sRU, sRV, sE;
    double Pmin, Amin, Amax, Acur, fcur, fmax, ffcur, ffmax;
    int pointP;
    double gamma = GAM;



    for (int k = 0; k < mesh->cCount; k++)
    {
        double x1 = mesh->nodes[mesh->cells[k].nodesInd[0]].x;
        double x2 = mesh->nodes[mesh->cells[k].nodesInd[1]].x;
        double x3 = mesh->nodes[mesh->cells[k].nodesInd[2]].x;

        double y1 = mesh->nodes[mesh->cells[k].nodesInd[0]].y;
        double y2 = mesh->nodes[mesh->cells[k].nodesInd[1]].y;
        double y3 = mesh->nodes[mesh->cells[k].nodesInd[2]].y;

        double xc = (x1+x2+x3)/3.;
        double yc = (y1+y2+y3)/3.;

        Rcur1 = getRO(k, x1, y1, NEW); RUcur1 = getRU(k, x1, y1, NEW); RVcur1 = getRV(k, x1, y1, NEW); Ecur1 = getRE(k, x1, y1, NEW);
        Rcur2 = getRO(k, x2, y2, NEW); RUcur2 = getRU(k, x2, y2, NEW); RVcur2 = getRV(k, x2, y2, NEW); Ecur2 = getRE(k, x2, y2, NEW);
        Rcur3 = getRO(k, x3, y3, NEW); RUcur3 = getRU(k, x3, y3, NEW); RVcur3 = getRV(k, x3, y3, NEW); Ecur3 = getRE(k, x3, y3, NEW);

        if((Ecur1-(RUcur1*RUcur1+RVcur1*RVcur1)/(2.*Rcur1) > 0.) && (Ecur2-(RUcur2*RUcur2+RVcur2*RVcur2)/(2.*Rcur2) > 0.) && (Ecur3-(RUcur3*RUcur3+RVcur3*RVcur3)/(2.*Rcur3) > 0.)){
            goto lbl_Lim_1;

        } else{
            t[0] = -sqrt(3.0/5.0);
            t[1] = 0.;
            t[2] = -t[0];

            double feps = 1.e-13;

            for (int i = 0; i < 3; i++) {
                xtt[i] = (x2-x1)*(t[i]+1)/2.+x1;
                ytt[i] = (y2-y1)*(t[i]+1)/2.+y1;
            }
            for (int i = 3; i < 6; i++) {
                xtt[i] = (x3-x2)*(t[i-3]+1)/2.+x2;
                ytt[i] = (y3-y2)*(t[i-3]+1)/2.+y2;
            }

            for (int i = 6; i < 9; i++) {

                xtt[i] = (x1-x3)*(t[i-6]+1)/2.+x3;
                ytt[i] = (y1-y3)*(t[i-6]+1)/2.+y3;
            }

            double fRcur, Rmin;
            for (int j = 0; j < 9; j++)
            {
                fRcur = ro[k][0]+ro[k][1]*getF(1, k, xtt[j], ytt[j])+ro[k][2]*getF(2, k, xtt[j], ytt[j]);
                if (j == 0)
                {
                    Rmin = fRcur;
                }
                else
                {
                    if(fRcur <= Rmin) Rmin = fRcur;
                }
            }

            double alfR = (ro[k][0]==Rmin) ? 1.0 : MIN((ro[k][0]-feps)/(ro[k][0]-Rmin),1.);
            for (int i = 1; i < BASE_FUNC_COUNT; i++) ro[k][i] *= alfR;


            for (int j = 0; j < 9; j++)
            {
                double sR, sRU, sRV, sE;
                sR = getRO(k, xtt[j], ytt[j], NEW); sRU = getRU(k, xtt[j], ytt[j], NEW); sRV = getRV(k, xtt[j], ytt[j], NEW); sE = getRE(k, xtt[j], ytt[j], NEW);

                double sep = sE/sR-((sRU*sRU)/(sR*sR)+(sRV*sRV)/(sR*sR))/2.;

                double P = (gamma-1.)*sR*sep;

                if(j==0){
                    Pmin = P;
                    pointP = j;
                }else{
                    if(P<=Pmin){
                        Pmin = P;
                        pointP = j;
                    }
                }
            }


            if(Pmin>=feps){
                goto lbl_Lim_1;
            } else {
                Amin = 0.;
                Amax = 1.;
                Acur = (Amax - Amin)/2. + Amin;

                sR = 0.;
                sRU = 0.;
                sRV = 0.;
                sE = 0.;


                for (int i = 1; i < BASE_FUNC_COUNT; i++)
                {
                    sR  = sR  + ro[k][i]*getF(i,k,xtt[pointP],ytt[pointP]);
                    sRU = sRU + ru[k][i]*getF(i,k,xtt[pointP],ytt[pointP]);
                    sRV = sRV + rv[k][i]*getF(i,k,xtt[pointP],ytt[pointP]);
                    sE  = sE  + re[k][i]*getF(i,k,xtt[pointP],ytt[pointP]);
                }

                fcur = re[k][0]+Acur*sE-((ru[k][0]+Acur*sRU)*(ru[k][0]+Acur*sRU)+(rv[k][0]+Acur*sRV)*(rv[k][0]+Acur*sRV))/(2.*(ro[k][0]+Acur*sR));
                fmax = re[k][0]+Amax*sE-((ru[k][0]+Amax*sRU)*(ru[k][0]+Amax*sRU)+(rv[k][0]+Amax*sRV)*(rv[k][0]+Amax*sRV))/(2.*(ro[k][0]+Amax*sR));
                ffmax = (gamma - 1.)*fmax - feps;
                ffcur = (gamma - 1.)*fcur - feps;

                double tfl = 0.;

                while (fabs(ffmax)>feps) {
                    if(tfl>10) goto lbl_Lim_2;
                    if(ffmax*ffcur>0.)
                        Amax = Acur;
                    else
                        Amin = Acur;
                    Acur = (Amax - Amin)/2. + Amin;
                    for (int i = 1; i < BASE_FUNC_COUNT; i++)
                    {
                        sR  = sR  + ro[k][i]*getF(i,k,xtt[pointP],ytt[pointP]);
                        sRU = sRU + ru[k][i]*getF(i,k,xtt[pointP],ytt[pointP]);
                        sRV = sRV + rv[k][i]*getF(i,k,xtt[pointP],ytt[pointP]);
                        sE  = sE  + re[k][i]*getF(i,k,xtt[pointP],ytt[pointP]);
                    }

                    fcur = re[k][0]+Acur*sE-((ru[k][0]+Acur*sRU)*(ru[k][0]+Acur*sRU)+(rv[k][0]+Acur*sRV)*(rv[k][0]+Acur*sRV))/(2.*(ro[k][0]+Acur*sR));
                    fmax = re[k][0]+Amax*sE-((ru[k][0]+Amax*sRU)*(ru[k][0]+Amax*sRU)+(rv[k][0]+Amax*sRV)*(rv[k][0]+Amax*sRV))/(2.*(ro[k][0]+Amax*sR));
                    ffmax = (gamma - 1.)*fmax - feps;
                    ffcur = (gamma - 1.)*fcur - feps;


                    tfl = tfl + 1;
                }

                lbl_Lim_2:
                for (int  i = 1; i< BASE_FUNC_COUNT; i++)
                {
                    ro[k][i] *= Acur;
                    ru[k][i] *= Acur;
                    rv[k][i] *= Acur;
                    re[k][i] *= Acur;
                }
            }
        }
        lbl_Lim_1:;
    }
}

void MethodGasGalerkin::initLimiterParameters() {
    limAlfa = new double **[mesh->cCount];
    limNeigh = new int **[mesh->cCount];
    limLm = new Vector *[mesh->cCount];
    limLmN = new Vector *[mesh->cCount];
    limPm = new Point *[mesh->cCount];
    for (int i = 0; i < mesh->cCount; i++) {
        limAlfa[i] = new double *[3];
        limAlfa[i][0] = new double[2];
        limAlfa[i][1] = new double[2];
        limAlfa[i][2] = new double[2];

        limNeigh[i] = new int *[3];
        limNeigh[i][0] = new int[2];
        limNeigh[i][1] = new int[2];
        limNeigh[i][2] = new int[2];

        limLm[i] = new Vector[3];
        limLmN[i] = new Vector[3];
        limPm[i] = new Point[3];
    }

    deltaS1 = new double[4];
    deltaS2 = new double[4];
    deltaU1 = new double *[3];
    deltaU2 = new double *[3];
    for (int m = 0; m < 3; m++) {
        deltaU1[m] = new double[4];
        deltaU2[m] = new double[4];
    }

    matrL = new double *[4];
    matrR = new double *[4];
    for (int i = 0; i < 4; i++) {
        matrL[i] = new double[4];
        matrR[i] = new double[4];
    }

    for (int iCell = 0; iCell < mesh->cCount; iCell++) {
        int n0 = mesh->cells[iCell].neigh[0];
        int n1 = mesh->cells[iCell].neigh[1];
        int n2 = mesh->cells[iCell].neigh[2];
        for (int m = 0; m < 3; m++) {
            int iEdge = getEdgeByCells(iCell, mesh->cells[iCell].neigh[m]);
            limPm[iCell][m].x = mesh->edges[iEdge].c[0].x;
            limPm[iCell][m].y = mesh->edges[iEdge].c[0].y;
            limLm[iCell][m].x = mesh->edges[iEdge].c[0].x - mesh->cells[iCell].c.x;
            limLm[iCell][m].y = mesh->edges[iEdge].c[0].y - mesh->cells[iCell].c.y;
            double tmp = sqrt(limLm[iCell][m].x * limLm[iCell][m].x + limLm[iCell][m].y * limLm[iCell][m].y);
            limLmN[iCell][m].x = limLm[iCell][m].x / tmp;
            limLmN[iCell][m].y = limLm[iCell][m].y / tmp;
            choiseDirection(limNeigh[iCell][m][0], limNeigh[iCell][m][1], limAlfa[iCell][m][0], limAlfa[iCell][m][1], iCell, n0, n1, n2, limPm[iCell][m], m);
            if (limAlfa[iCell][m][0] < 0.0 || limAlfa[iCell][m][1] < 0.0) printf("ERROR!!!\n");
        }
    }
}

void MethodGasGalerkin::getMatrLR(double **L, double **R, double u, double v, double c) {
    double h = c*c/(GAM - 1.0)+0.5*(u*u+v*v);

    R[0][0] = 1.0;		R[0][1] = 1.0;		R[0][2] = 0.0;	R[0][3] = 1.0;
    R[1][0] = u-c;		R[1][1] = u;		R[1][2] = 0.0;	R[1][3] = u+c;
    R[2][0] = v;		R[2][1] = v;		R[2][2] = 1.0;	R[2][3] = v;
    R[3][0] = h-u*c;	R[3][1] = 0.5*v*v;	R[3][2] = v;	R[3][3] = h+u*c;

    inverseMatrix(R, L, 4);
}

int MethodGasGalerkin::getEdgeByCells(int c1, int c2) {
    for (int iEdge = 0; iEdge < mesh->eCount; iEdge++) {
        if ((mesh->edges[iEdge].c1 == c1 && mesh->edges[iEdge].c2 == c2) || (mesh->edges[iEdge].c1 == c2 && mesh->edges[iEdge].c2 == c1))
            return iEdge;
    }
}

void
MethodGasGalerkin::choiseDirection(int &nn1, int &nn2, double &a1, double &a2, int n0, int n1, int n2, int n3, Point pm,
                                   int mm) {
    int nn[3] = {n1, n2, n3};

    for (int m = 0; m < 3; m++)
    {
        nn1 = nn[m];
        nn2 = nn[(m+1)%3];
        if (nn1 >= 0 && nn2 >= 0)
        {
            double D  = (mesh->cells[nn1].c.x-mesh->cells[n0].c.x)*(mesh->cells[nn2].c.y-mesh->cells[n0].c.y)-(mesh->cells[nn1].c.y-mesh->cells[n0].c.y)*(mesh->cells[nn2].c.x-mesh->cells[n0].c.x);
            double D1 = (pm.x-mesh->cells[n0].c.x)*(mesh->cells[nn2].c.y-mesh->cells[n0].c.y)-(pm.y-mesh->cells[n0].c.y)*(mesh->cells[nn2].c.x-mesh->cells[n0].c.x);
            double D2 = (mesh->cells[nn1].c.x-mesh->cells[n0].c.x)*(pm.y-mesh->cells[n0].c.y)-(mesh->cells[nn1].c.y-mesh->cells[n0].c.y)*(pm.x-mesh->cells[n0].c.x);

            a1 = D1/D;
            a2 = D2/D;
            if (a1>=0.0 && a2>=0.0) return;
        }
    }
    if (n1 < 0 || n2 < 0 || n3 < 0)
    {
        a1 = 0.0;
        a2 = 0.0;
        nn1 = 0;
        nn2 = 0;
        return;
    } else {
        printf(" ERROR: choiseDirection()\n");
        printf("N0 = %d, N1 = %d, N2 = %d, N3 = %d, pm.x = %f, pm.y = %f, m = %d\n", n0, n1, n2, n3, pm.x, pm.y, mm);
        printf("p0.x = %f, p0.y = %f\n", mesh->cells[n0].c.x, mesh->cells[n0].c.y);
        printf("p1.x = %f, p1.y = %f\n", mesh->cells[n1].c.x, mesh->cells[n1].c.y);
        printf("p2.x = %f, p2.y = %f\n", mesh->cells[n2].c.x, mesh->cells[n2].c.y);
        printf("p3.x = %f, p3.y = %f\n", mesh->cells[n3].c.x, mesh->cells[n3].c.y);
        printf("\n\n");
        for (int m = 0; m < 3; m++)
        {
            nn1 = nn[m];
            nn2 = nn[(m+1)%3];
            printf("nn1 = %d, nn2 = %d\n", nn1 , nn2);
            if (nn1 >= 0 && nn2 >= 0)
            {
                double D  = (mesh->cells[nn1].c.x-mesh->cells[n0].c.x)*(mesh->cells[nn2].c.y-mesh->cells[n0].c.y)-(mesh->cells[nn1].c.y-mesh->cells[n0].c.y)*(mesh->cells[nn2].c.x-mesh->cells[n0].c.x);
                double D1 = (pm.x-mesh->cells[n0].c.x)*(mesh->cells[nn2].c.y-mesh->cells[n0].c.y)-(pm.y-mesh->cells[n0].c.y)*(mesh->cells[nn2].c.x-mesh->cells[n0].c.x);
                double D2 = (mesh->cells[nn1].c.x-mesh->cells[n0].c.x)*(pm.y-mesh->cells[n0].c.y)-(mesh->cells[nn1].c.y-mesh->cells[n0].c.y)*(pm.x-mesh->cells[n0].c.x);

                printf("D = %f, D1 = %f, D2 = %f\n");

                a1 = D1/D; if (fabs(a1)<EPS) a1 = 0.0;
                a2 = D2/D; if (fabs(a2)<EPS) a2 = 0.0;
                printf("a1 = %f, a2 = %f\n");
            }
        }
        exit(1);
    }
}

void MethodGasGalerkin::copyToOld() {
    for (int i = 0; i < mesh->cCount; i++)
    {
        memcpy(roOld[i], ro[i], BASE_FUNC_COUNT*sizeof(double));
        memcpy(ruOld[i], ru[i], BASE_FUNC_COUNT*sizeof(double));
        memcpy(rvOld[i], rv[i], BASE_FUNC_COUNT*sizeof(double));
        memcpy(reOld[i], re[i], BASE_FUNC_COUNT*sizeof(double));
    }
    memcpy(roOld, ro, mesh->cCount*sizeof(double));
    memcpy(ruOld, ru, mesh->cCount*sizeof(double));
    memcpy(rvOld, rv, mesh->cCount*sizeof(double));
    memcpy(reOld, re, mesh->cCount*sizeof(double));
}

void MethodGasGalerkin::flux(Param pl, Param pr, Vector n, double &fr, double &fu, double &fv, double &fe) {
    double rol, rul, rvl, rel;
    double ror, rur, rvr, rer;
    double frl, ful, fvl, fel;
    double frr, fur, fvr, fer;
    double alpha, unl, unr, q1,q2;

    unl = n.x*pl.u+n.y*pl.v;
    unr = n.x*pr.u+n.y*pr.v;

    rol = pl.r;
    rul = pl.r*pl.u;
    rvl = pl.r*pl.v;
    rel = pl.r*(pl.e+0.5*pl.U2());

    ror = pr.r;
    rur = pr.r*pr.u;
    rvr = pr.r*pr.v;
    rer = pr.r*(pr.e+0.5*pr.U2());

    frl = pl.r*unl;
    ful = frl*pl.u+pl.p*n.x;
    fvl = frl*pl.v+pl.p*n.y;
    fel = (rel+pl.p)*unl;

    frr = pr.r*unr;
    fur = frr*pr.u+pr.p*n.x;
    fvr = frr*pr.v+pr.p*n.y;
    fer = (rer+pr.p)*unr;

    q1 = sqrt(pl.p*GAM/pl.r)+fabs(pl.U2());
    q2 = sqrt(pr.p*GAM/pr.r)+fabs(pr.U2());
    alpha = (q1 > q2) ? q1 : q2;

    fr = 0.5*((frl+frr)-alpha*(ror-rol));
    fu = 0.5*((ful+fur)-alpha*(rur-rul));
    fv = 0.5*((fvl+fvr)-alpha*(rvr-rvl));
    fe = 0.5*((fel+fer)-alpha*(rer-rel));
}

void MethodGasGalerkin::bnd(Edge *e, Point pt, Param p1, Param &p2, Vector n) {
    switch (e->type) {
        case 1: // вытекание
            p2 = p1;
            break;
        case 2: // втекание
            p2.r = 1.;
            p2.u = 3.;
            p2.v = 0.0;
            p2.p = 1./1.4;
            p2.e = p2.p/p2.r/(GAM-1.0);
            p2.T = 0.0;
            p2.M = 0.0;
            break;
        case 3: // отражение
            p2 = p1;
            double Un = p1.u*e->n.x + p1.v*e->n.y;
            Vector V;
            V.x = e->n.x*Un*2.0;
            V.y = e->n.y*Un*2.0;
            p2.u = p1.u - V.x;
            p2.v = p1.v - V.y;
            break;
    }
}
