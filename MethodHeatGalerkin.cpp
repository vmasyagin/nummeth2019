#include "MethodHeatGalerkin.h"
#include <cstdio>
#include <cstring>
#include <cmath>
#include <stdlib.h>


void MethodHeatGalerkin::convertToParam(int i, Point pt, Param& p)
{
    p.r = 0.0;
    p.p = 0.0;
    p.u = 0.0;
    p.v = 0.0;
    p.e = 0.0;
    p.T = getT(i, pt);
    p.M = 0.0;
    p.qx = getQx(i, pt);
    p.qy = getQy(i, pt);
}

void MethodHeatGalerkin::calcGaussParam() {

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

double MethodHeatGalerkin::getF(int i, int iCell, Point pt)
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

double MethodHeatGalerkin::getDfDx(int i, int iCell, Point pt)
{
    if (BASE_FUNC_COUNT == 3)
    {
        switch(i) {
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
                return 2.0*(pt.x - mesh->cells[iCell].c.x) / (mesh->cells[iCell].HX*mesh->cells[iCell].HX);
            case 4:
                return 0.0;
            case 5:
                return (pt.y - mesh->cells[iCell].c.y) / (mesh->cells[iCell].HX*mesh->cells[iCell].HY);
        }
    }
    else {
        printf("ERROR: wrong basic functions count!\n");
        exit(1);
    }
}

double MethodHeatGalerkin::getDfDy(int i, int iCell, Point pt)
{
    if (BASE_FUNC_COUNT == 3)
    {
        switch(i) {
            case 0:
                return 0.0;
            case 1:
                return 0.0;
            case 2:
                return 1.0 / mesh->cells[iCell].HY;
        }
    }
    else if (BASE_FUNC_COUNT == 6)
    {
        switch(i) {
            case 0:
                return 0.0;
            case 1:
                return 0.0;
            case 2:
                return 1.0 / mesh->cells[iCell].HY;
            case 3:
                return 0.0;
            case 4:
                return 2.0*(pt.y - mesh->cells[iCell].c.y) / (mesh->cells[iCell].HY*mesh->cells[iCell].HY);
            case 5:
                return (pt.x - mesh->cells[iCell].c.x) / (mesh->cells[iCell].HX*mesh->cells[iCell].HY);
        }
    }
    else {
        printf("ERROR: wrong basic functions count!\n");
        exit(1);
    }
}

void MethodHeatGalerkin::init()
{
    mesh = new Mesh();
    mesh->initFromFiles((char*)"heat.1");

    T = new double*[mesh->cCount];
    qx = new double*[mesh->cCount];
    qy = new double*[mesh->cCount];

    intT = new double*[mesh->cCount];
    intQx = new double*[mesh->cCount];
    intQy = new double*[mesh->cCount];

    TMAX    = 1.0;
    TAU     = 1.e-5;
    BASE_FUNC_COUNT = 3;
    GP_CELL_COUNT = 3;
    GP_EDGE_COUNT = 2;
    C11 = 1.0;

    cellGP = new Point*[mesh->cCount];
    cellGW = new double*[mesh->cCount];
    cellJ = new double[mesh->cCount];
    edgeGP = new Point*[mesh->eCount];
    edgeGW = new double*[mesh->eCount];
    edgeJ = new double[mesh->eCount];

    for(int i = 0; i < mesh->cCount; i++){
        T[i] = new double[BASE_FUNC_COUNT];
        memset(T[i], 0.0, BASE_FUNC_COUNT* sizeof(double));

        qx[i] = new double[BASE_FUNC_COUNT];
        qy[i] = new double[BASE_FUNC_COUNT];

        intT[i] = new double[BASE_FUNC_COUNT];
        intQx[i] = new double[BASE_FUNC_COUNT];
        intQy[i] = new double[BASE_FUNC_COUNT];

        cellGP[i] = new Point[GP_CELL_COUNT];
        cellGW[i] = new double[GP_CELL_COUNT];
    }

    saveVTK(0);

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
}

static void inverseMatrix(double** a_src, double **am, int N)
{
    double **a = a_src;
    double detA = a[0][0] * a[1][1] * a[2][2] + a[1][0] * a[2][1] * a[0][2] + a[0][1] * a[1][2] * a[2][0]
                  - a[2][0] * a[1][1] * a[0][2] - a[1][0] * a[0][1] * a[2][2] - a[0][0] * a[2][1] * a[1][2];

    double m[3][3];
    m[0][0] = a[1][1] * a[2][2] - a[2][1] * a[1][2];
    m[0][1] = a[2][0] * a[1][2] - a[1][0] * a[2][2];
    m[0][2] = a[1][0] * a[2][1] - a[2][0] * a[1][1];
    m[1][0] = a[2][1] * a[0][2] - a[0][1] * a[2][2];
    m[1][1] = a[0][0] * a[2][2] - a[2][0] * a[0][2];
    m[1][2] = a[2][0] * a[0][1] - a[0][0] * a[2][1];
    m[2][0] = a[0][1] * a[1][2] - a[1][1] * a[0][2];
    m[2][1] = a[1][0] * a[0][2] - a[0][0] * a[1][2];
    m[2][2] = a[0][0] * a[1][1] - a[1][0] * a[0][1];

    am[0][0] = m[0][0] / detA;
    am[0][1] = m[1][0] / detA;
    am[0][2] = m[2][0] / detA;
    am[1][0] = m[0][1] / detA;
    am[1][1] = m[1][1] / detA;
    am[1][2] = m[2][1] / detA;
    am[2][0] = m[0][2] / detA;
    am[2][1] = m[1][2] / detA;
    am[2][2] = m[2][2] / detA;
}

void MethodHeatGalerkin::calcMassMatrix() {
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

void MethodHeatGalerkin::run()
{
    double t = 0.0;
    int step = 0;

    while (t < TMAX) {
        t += TAU;
        step++;


        for(int i = 0; i < mesh->cCount; i++){
            memset(intT[i], 0, sizeof(double)*BASE_FUNC_COUNT);
            memset(intQx[i], 0, sizeof(double)*BASE_FUNC_COUNT);
            memset(intQy[i], 0, sizeof(double)*BASE_FUNC_COUNT);
        }

        calcGradients();
        calcDiffusionVol();
        calcDiffusionSurf();

        calcNewValues();

        if (step % 10000 == 0)
        {
            saveVTK(step);
            printf("Calculation results for step %d are saved.\n", step);
        }
    }
}

MethodHeatGalerkin::~MethodHeatGalerkin() {

    for (int i = 0; i < mesh->cCount; i++) {
        delete[] T[i], qx[i], qy[i];
        delete[] intT[i], intQx[i], intQy[i];
        delete[] cellGP[i], cellGW[i];
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
    delete[] T, qx, qy;
    delete[] intT, intQx, intQy;
    delete[] cellGP, cellGW, cellJ;
    delete[] edgeGP, edgeGW, edgeJ;
    delete[] A, invA;
    delete mesh;
}

double MethodHeatGalerkin::getT(int iCell, Point pt) {
    double d = 0.0;
    for(int i = 0; i < BASE_FUNC_COUNT; i++){
        d += T[iCell][i]*getF(i, iCell, pt);
    }
    return d;
}

double MethodHeatGalerkin::getQx(int iCell, Point pt) {
    double d = 0.0;
    for(int i = 0; i < BASE_FUNC_COUNT; i++){
        d += qx[iCell][i]*getF(i, iCell, pt);
    }
    return d;
}

double MethodHeatGalerkin::getQy(int iCell, Point pt) {
    double d = 0.0;
    for(int i = 0; i < BASE_FUNC_COUNT; i++){
        d += qy[iCell][i]*getF(i, iCell, pt);
    }
    return d;
}

void MethodHeatGalerkin::bnd(Edge *e, Point pt, Param p1, Param &p2) {
    switch (e->type) {
        case 1: // тепловой поток = 0
            p2.T = 0.0;
            p2.qx = p1.qx;
            p2.qy = p1.qy;
            break;
        case 2: // постоянное значение
            p2.T = 1.0;
            p2.qx = p1.qx;
            p2.qy = p1.qy;
            break;
    }
}

void MethodHeatGalerkin::flux(Param pl, Param pr, Vector n, double &fT, double &fqx, double &fqy) {
    fT = 0.5*(pl.T + pr.T);
    fqx = 0.5*(pl.qx + pr.qx) + C11*(pr.T - pl.T)*n.x;
    fqy = 0.5*(pl.qy + pr.qy) + C11*(pr.T - pl.T)*n.y;
}

void MethodHeatGalerkin::calcGradients() {
    // поверхностные интегралы
    for (int i = 0; i < mesh->cCount; i++){
        for(int j = 0; j < BASE_FUNC_COUNT; j++) {
            double tmpIntQx = 0.0, tmpIntQy = 0.0;
            for (int z = 0; z < GP_CELL_COUNT; z++) {
                tmpIntQx += cellGW[i][z] * getT(i, cellGP[i][z]) * getDfDx(j, i, cellGP[i][z]);
                tmpIntQy += cellGW[i][z] * getT(i, cellGP[i][z]) * getDfDy(j, i, cellGP[i][z]);
            }
            tmpIntQx *= cellJ[i];
            tmpIntQy *= cellJ[i];
            intQx[i][j] -= tmpIntQx;
            intQy[i][j] -= tmpIntQy;
        }
    }

    // криволинейные интегралы
    for(int i = 0; i < mesh->eCount; i++) {
        Edge &e = mesh->edges[i];
        int c1 = e.c1;

        for (int j = 0; j < BASE_FUNC_COUNT; j++) {
            int c2 = e.c2;
            double tmpIntQx1 = 0.0, tmpIntQx2 = 0.0, tmpIntQy1 = 0.0, tmpIntQy2 = 0.0;
            for (int iGP = 0; iGP < GP_EDGE_COUNT; iGP++) {
                Point &pt = edgeGP[i][iGP];

                Param p1, p2;
                convertToParam(c1, pt, p1);
                if (e.c2 >= 0) {
                    convertToParam(e.c2, pt, p2);
                }
                else {
                    bnd(&e, pt, p1, p2);
                }
                double fT, fqx, fqy;
                flux(p1, p2, e.n, fT, fqx, fqy);

                tmpIntQx1 += e.n.x * fT * getF(j, c1, pt) * edgeGW[i][iGP];
                tmpIntQy1 += e.n.y * fT * getF(j, c1, pt) * edgeGW[i][iGP];
                if (e.c2 >= 0) {
                    tmpIntQx2 += e.n.x * fT * getF(j, c2, pt) * edgeGW[i][iGP];
                    tmpIntQy2 += e.n.y * fT * getF(j, c2, pt) * edgeGW[i][iGP];
                }
            }

            tmpIntQx1 *= edgeJ[i];
            tmpIntQy1 *= edgeJ[i];

            tmpIntQx2 *= edgeJ[i];
            tmpIntQy2 *= edgeJ[i];

            intQx[c1][j] += tmpIntQx1;
            intQy[c1][j] += tmpIntQy1;

            if (e.c2 >= 0) {
                intQx[c2][j] -= tmpIntQx2;
                intQy[c2][j] -= tmpIntQy2;
            }
        }

    }

    //вычисляем градиенты на новом шаге по времени
    for(int i = 0; i < mesh->cCount; i++) {
        for (int j = 0; j < BASE_FUNC_COUNT; j++) {
            double tmpQx = 0.0, tmpQy = 0.0;
            for (int z = 0; z < BASE_FUNC_COUNT; z++) {
                tmpQx += invA[i][j][z] * intQx[i][z];
                tmpQy += invA[i][j][z] * intQy[i][z];
            }
            qx[i][j] = tmpQx;
            qy[i][j] = tmpQy;
        }
    }
}

void MethodHeatGalerkin::calcDiffusionVol() {
    for (int i = 0; i < mesh->cCount; i++){
        for(int j = 0; j < BASE_FUNC_COUNT; j++) {
            double tmpIntT = 0.0;
            for (int z = 0; z < GP_CELL_COUNT; z++) {
                tmpIntT += cellGW[i][z] * getQx(i, cellGP[i][z]) * getDfDx(j, i, cellGP[i][z]);
                tmpIntT += cellGW[i][z] * getQy(i, cellGP[i][z]) * getDfDy(j, i, cellGP[i][z]);
            }
            tmpIntT *= cellJ[i];
            intT[i][j] -= tmpIntT;
        }
    }
}

void MethodHeatGalerkin::calcDiffusionSurf() {
    for(int i = 0; i < mesh->eCount; i++) {
        Edge &e = mesh->edges[i];
        int c1 = e.c1;

        for (int j = 0; j < BASE_FUNC_COUNT; j++) {
            int c2 = e.c2;
            double tmpIntT1 = 0.0, tmpIntT2 = 0.0;
            for (int iGP = 0; iGP < GP_EDGE_COUNT; iGP++) {
                Point &pt = edgeGP[i][iGP];

                Param p1, p2;
                convertToParam(c1, pt, p1);
                if (e.c2 >= 0) {
                    convertToParam(e.c2, pt, p2);
                }
                else {
                    bnd(&e, pt, p1, p2);
                }
                double fT, fqx, fqy;
                flux(p1, p2, e.n, fT, fqx, fqy);

                tmpIntT1 += e.n.x * fqx * getF(j, c1, pt) * edgeGW[i][iGP];
                tmpIntT1 += e.n.y * fqy * getF(j, c1, pt) * edgeGW[i][iGP];
                if (e.c2 >= 0) {
                    tmpIntT2 += e.n.x * fqx * getF(j, c2, pt) * edgeGW[i][iGP];
                    tmpIntT2 += e.n.y * fqy * getF(j, c2, pt) * edgeGW[i][iGP];
                }
            }

            tmpIntT1 *= edgeJ[i];

            tmpIntT2 *= edgeJ[i];

            intT[c1][j] += tmpIntT1;

            if (e.c2 >= 0) {
                intT[c2][j] -= tmpIntT2;
            }
        }

    }
}

void MethodHeatGalerkin::calcNewValues() {
    for(int i = 0; i < mesh->cCount; i++) {
        for (int j = 0; j < BASE_FUNC_COUNT; j++) {
            double tmpT = 0.0;
            for (int z = 0; z < BASE_FUNC_COUNT; z++) {
                tmpT += invA[i][j][z] * intT[i][z];
            }
            T[i][j] += TAU*tmpT;
        }
    }
}

