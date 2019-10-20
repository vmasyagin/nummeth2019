#include <iostream>
#include "method.h"

using namespace std;

void printUsage(char * exeName)
{
    cout << "Usage:" << endl;
    cout << "  " << exeName << " <method_code>" << endl << endl;
    cout << "  <method_code>: " << METHOD_CODE_HEAT_FVM << " - heat eq. solver; " << endl <<
            "                 " << METHOD_CODE_GAS_FVM << " - gas dynamics eq. solver." << endl <<
            "                 " << METHOD_CODE_HEAT_GALERKIN << " - heat eq. Galerkin solver." << endl <<
            "                 " << METHOD_CODE_GAS_GALERKIN << " - gas dynamics eq. Galerkin solver." << endl;
}

int main(int argc, char** argv)
{
    if (argc < 2) {
        printUsage(argv[0]);
        return -1;
    }
    int mCode = atoi(argv[1]);
    Method * m = Method::create(mCode);
    if (!m) {
        return -2;
    }
    m->init();
    m->run();

    delete m;
}
