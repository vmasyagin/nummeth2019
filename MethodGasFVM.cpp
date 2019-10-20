#include "MethodGasFVM.h"
#include <cstdio>
#include <cstring>
#include <cmath>


void MethodGasFVM::convertToParam(int i, Point pt, Param& p)
{
    p.r = ro[i];
    p.u = ru[i]/p.r;
    p.v = rv[i]/p.r;
    p.e = re[i]/p.r-p.U2()/2.;
    p.p = p.r*p.e*(GAM-1.0);
    p.T = 0.0;
    p.M = p.U2()/sqrt(GAM*p.p/p.r);
}


void MethodGasFVM::init()
{
    mesh = new Mesh();
    mesh->initFromFiles((char*)"step.1");

    ro = new double[mesh->cCount];
    ru = new double[mesh->cCount];
    rv = new double[mesh->cCount];
    re = new double[mesh->cCount];

    initValues();
    saveVTK(0);

    intRO = new double[mesh->cCount];
    intRU = new double[mesh->cCount];
    intRV = new double[mesh->cCount];
    intRE = new double[mesh->cCount];

    TMAX    = 4.0;
    TAU     = 1.e-4;
}


void MethodGasFVM::initValues()
{
    for (int i = 0; i < mesh->cCount; i++) {
        ro[i] = 1.;
        ru[i] = 3.;
        rv[i] = 0.0;
        re[i] = (1./1.4)/(GAM-1.0)+0.5*(ru[i]*ru[i])/ro[i];
    }
}


void MethodGasFVM::run()
{
    double t = 0.0;
    int step = 0;
    while (t < TMAX) {
        t += TAU;
        step++;

        memset(intRO, 0, sizeof(double)*mesh->cCount);
        memset(intRU, 0, sizeof(double)*mesh->cCount);
        memset(intRV, 0, sizeof(double)*mesh->cCount);
        memset(intRE, 0, sizeof(double)*mesh->cCount);
        for (int ie = 0; ie < mesh->eCount; ie++) {
            Edge &e = mesh->edges[ie];
            int c1 = e.c1;
            Cell &cell1 = mesh->getCell(c1);

            Param p1, p2;
            convertToParam(c1, mesh->cells[c1].c, p1);
            if (e.c2 >= 0) {
                convertToParam(e.c2, mesh->cells[e.c2].c, p2);
            }
            else {
                bnd(&e, p1, p2);
            }

            double fr, fu, fv, fe;
            flux(p1, p2, e.n, fr, fu, fv, fe);
            fr *= e.l;
            fu *= e.l;
            fv *= e.l;
            fe *= e.l;
            intRO[c1] -= fr;
            intRU[c1] -= fu;
            intRV[c1] -= fv;
            intRE[c1] -= fe;
            if (e.c2 >= 0) {
                intRO[e.c2] += fr;
                intRU[e.c2] += fu;
                intRV[e.c2] += fv;
                intRE[e.c2] += fe;
            }
        }

        for (int i = 0; i < mesh->cCount; i++) {
            double CFL = TAU/mesh->cells[i].S;
            ro[i] += intRO[i]*CFL;
            ru[i] += intRU[i]*CFL;
            rv[i] += intRV[i]*CFL;
            re[i] += intRE[i]*CFL;
        }

        if (step % PRINT_STEP == 0)
        {
            saveVTK(step);
            printf("Calculation results for step %d are saved.\n", step);
        }
    }
}


void MethodGasFVM::flux(Param pl, Param pr, Vector n, double &fr, double &fu, double &fv, double &fe)
{
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


void MethodGasFVM::bnd(Edge *e, Param p1, Param &p2)
{
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


MethodGasFVM::~MethodGasFVM()
{
    delete mesh;
    delete[] ro, ru, rv, re;
    delete[] intRO, intRU, intRV, intRE;
}
