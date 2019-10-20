#include "method.h"
#include "MethodHeatFVM.h"
#include "MethodHeatGalerkin.h"
#include "MethodGasFVM.h"
//#include "MethodGasGalerkin.h"
#include <cstdio>
#include <cstdlib>
#include "math.h"
#pragma warning(disable:4996)

Method* Method::create(int methodCode) {
    switch (methodCode) {
        case METHOD_CODE_HEAT_FVM:
            return new MethodHeatFVM();
        case METHOD_CODE_GAS_FVM:
            return new MethodGasFVM();
        case METHOD_CODE_HEAT_GALERKIN:
            return new MethodHeatGalerkin();
        //case METHOD_CODE_GAS_GALERKIN:
            //return new MethodGasGalerkin();
    }
}

void Method::saveVTK(int step) {
    char fName[50];

    sprintf(fName, "res_%010d.vtk", step);
    FILE * fp = fopen(fName, "w");
    if (!fp) {
        fprintf(stderr, "Can not open file '%s'\n", fName);
        exit(1);
    }
    fprintf(fp, "# vtk DataFile Version 2.0\n");
    fprintf(fp, "GASDIN data file\n");
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
    fprintf(fp, "POINTS %d float\n", mesh->nCount);
    for (int i = 0; i < mesh->nCount; i++)
    {
        fprintf(fp, "%f %f %f  ", mesh->nodes[i].x,  mesh->nodes[i].y, 0.0);
        if (i+1 % 8 == 0) fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
    fprintf(fp, "CELLS %d %d\n", mesh->cCount, 4*mesh->cCount);
    for (int i = 0; i < mesh->cCount; i++)
    {
        fprintf(fp, "3 %d %d %d\n", mesh->cells[i].nodesInd[0], mesh->cells[i].nodesInd[1], mesh->cells[i].nodesInd[2]);
    }
    fprintf(fp, "\n");

    fprintf(fp, "CELL_TYPES %d\n", mesh->cCount);
    for (int i = 0; i < mesh->cCount; i++) fprintf(fp, "5\n");
    fprintf(fp, "\n");

    fprintf(fp, "CELL_DATA %d\nSCALARS Density float 1\nLOOKUP_TABLE default\n", mesh->cCount);
    for (int i = 0; i < mesh->cCount; i++)
    {
        Param p;
        convertToParam(i, mesh->cells[i].c, p);
        fprintf(fp, "%25.16f ", p.r);
        if (i+1 % 8 == 0 || i+1 == mesh->cCount) fprintf(fp, "\n");
    }

    fprintf(fp, "SCALARS Pressure float 1\nLOOKUP_TABLE default\n");
    for (int i = 0; i < mesh->cCount; i++)
    {
        Param p;
        convertToParam(i, mesh->cells[i].c, p);
        fprintf(fp, "%25.16f ", p.p);
        if (i+1 % 8 == 0 || i+1 == mesh->cCount) fprintf(fp, "\n");
    }

    fprintf(fp, "SCALARS Temperature float 1\nLOOKUP_TABLE default\n");
    for (int i = 0; i < mesh->cCount; i++)
    {
        Param p;
        convertToParam(i, mesh->cells[i].c, p);
        fprintf(fp, "%25.16f ", p.T);
        if (i+1 % 8 == 0 || i+1 == mesh->cCount) fprintf(fp, "\n");
    }

    fprintf(fp, "VECTORS Velosity float\n");
    for (int i = 0; i < mesh->cCount; i++)
    {
        Param p;
        convertToParam(i, mesh->cells[i].c, p);
        fprintf(fp, "%25.16f %25.16f %25.16f ", p.u, p.v, 0.0);
        if (i+1 % 8 == 0 || i+1 == mesh->cCount) fprintf(fp, "\n");
    }

    fclose(fp);
}
