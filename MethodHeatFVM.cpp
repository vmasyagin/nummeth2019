#include "MethodHeatFVM.h"
#include <cstdio>
#include <cstring>
#include <cmath>


void MethodHeatFVM::convertToParam(int i, Point pt, Param& p)
{
    p.r = 0.0;
    p.p = 0.0;
    p.u = 0.0;
    p.v = 0.0;
    p.e = 0.0;
    p.T = T[i];
    p.M = 0.0;
}


void MethodHeatFVM::init()
{
    mesh = new Mesh();
    mesh->initFromFiles((char*)"heat.1");

    T = new double[mesh->cCount];
    memset(T, 0, sizeof(double)*mesh->cCount);

    saveVTK(0);

    intT = new double[mesh->cCount];

    TMAX    = 1.0;
    TAU     = 1.e-5;
}


void MethodHeatFVM::run()
{
}

MethodHeatFVM::~MethodHeatFVM()
{
    delete mesh;
    delete[] T, intT;
}
