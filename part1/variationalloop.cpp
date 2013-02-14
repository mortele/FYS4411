#include "variationalloop.h"

VariationalLoop::VariationalLoop() {
}

void VariationalLoop::initialize_helium() {
    this->N = 50;
    this->start = 2.0;
    this->end   = 2.9;
    this->E     = 0.0;
    this->minE  = 1e300;
    this->minA  = 0.0;
    this->minB  = 0.0;

}

void VariationalLoop::run() {
    b = 1.0;

    for (int j = 0; j < N; j++) {
        a = start + ((end-start)/N) * j;
        cout << "b = "<< b << endl << "a = " << a << endl;
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
    }
    cout << "a " << minA << endl;
    cout << "b " << minB << endl;
    cout << "minE=" << minE << endl;
}
