#pragma once
#include "mesh.h"

#define METHOD_CODE_HEAT_FVM 1
#define METHOD_CODE_GAS_FVM  2
#define METHOD_CODE_HEAT_GALERKIN  3
#define METHOD_CODE_GAS_GALERKIN  4

const double EPS = 1.0e-12;
const int PRINT_STEP = 1000;

enum DATA_LAYER
{
    OLD = 0,
    NEW = 1
};

struct Param
{
    double r;		//!< плотность
    double p;		//!< давление
    double u;		//!< первая  компонента  вектора  скорости
    double v;		//!< вторая  компонента  вектора  скорости
    double e;		//!< внутренняя   энергия
    double T;		//!< температура
    double M;		//!< число  Маха
    double qx;      //!< первая компонента градиента температуры
    double qy;      //!< вторая компонента градиента температуры

    inline double U2() { return u*u+v*v; };
};

class Method
{
protected:
    Mesh * mesh;
public:
    virtual void convertToParam(int, Point, Param&) = 0;
    void saveVTK(int step);
    virtual void init() = 0;
    virtual void run() = 0;

    static Method* create(int methodCode);

};

