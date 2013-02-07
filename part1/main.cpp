#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <time.h>
#include <armadillo>

#include "lib.h"
//#include "varmc.h"
#include "variationalmc.h"

using namespace std;
using namespace arma;



int main() {

    int     N = 10;
    double  start = 2.0;
    double  end   = 2.9;
    double  E     = 0.0;
    double  minE  = 1e300;
    double  minA  = 0.0;
    double  minB  = 0.0;
    double  a;
    double  b;
    VariationalMC m;
    //VarMC m;

    //for (int i = 0; i < N; i++) {
        b = 1.0;

        for (int j = 0; j < N; j++) {
            a = start + 0.1 * j;
            E = m.runMetropolis(a,b);
            if (E < minE) {
                minE = E;
                minA = a;
                minB = b;
                cout << "minA=" << minA << endl;
                cout << "minB=" << minB << endl;
                cout << "minE=" << minE << endl;
            }
            cout << endl << endl;
      //  }
    }
    cout << "a " << minA << endl;
    cout << "b " << minB << endl;
    cout << "minE=" << minE << endl;
    return 0;
}


