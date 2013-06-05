#include <iostream>
#include <fstream>
#include <iomanip>
#include "variationalloop.h"


VariationalLoop::VariationalLoop() {

}

void VariationalLoop::initialize_processes(int my_rank) {

    this->my_rank = my_rank;
}

void VariationalLoop::run() {
    double startClock, finishClock, energy, accepted;
    vec    energyVarGrad = zeros<vec>(4);
    vec    varGrad       = zeros<vec>(2);
    vec    alphaBeta     = zeros<vec>(2);

//    ofstream outFileAlpha;
//    char* fileName = "alphaHeMolecule_2";

//    cout << fileName << endl;
//    outFileAlpha.open(fileName, ios::out);

    alphaBeta(0) = 1.843;
    alphaBeta(1) = 0.3465;



    //    npart=2:  (1.84, 0.35), start: gaussian_deviate(&idum) * 10.0 / (sqrt(alph)), dt =  0.01, alphaBeta / 1
    //    npart=4:  (3.93, 0.1),  start: gaussian_deviate(&idum) * 10.0 / (sqrt(alph)), dt =  0.01, alphaBeta / 1
    //    npart=10: (9.0, 0.3),   start: gaussian_deviate(&idum) *  5.0 / (sqrt(alph)), dt = 0.001, alphaBeta / 50


    for (int i = 1; i < 2; i++) {
            startClock = clock();
//        outFileAlpha << alphaBeta(0) << " " << alphaBeta(1) << endl;

        // Run one metropolis loop.
        energyVarGrad = m.runMetropolis(alphaBeta(0), alphaBeta(1), my_rank);

        energy        = energyVarGrad(0);
        varGrad(0)    = energyVarGrad(1);
        varGrad(1)    = energyVarGrad(2);
        accepted      = energyVarGrad(3);

        if (accepted > 0.95) {
            // Compute new values of alpha / beta.
            alphaBeta     += (1/((double) i)) * varGrad * energy / 1;
        } else {
            cout << "\n\n\n\n\n\naccepted  = "<< accepted << " < 0.95 \n\n\n\n\n";
        }


        if (alphaBeta(0) < 0) {
            alphaBeta(0) = 0.01;
        } if (alphaBeta(1) < 0) {
            alphaBeta(1) = 0.01;
        }

        finishClock = clock();
        cout << "   * Time usage       = " << (finishClock - startClock) / 1000000.0 << " [s] "<< endl;
        cout << "   * Alpha            = " << alphaBeta(0) << endl;
        cout << "   * Beta             = " << alphaBeta(1) << endl;
        cout << "VarGrad=" << varGrad << endl;
        cout << endl << endl;
    }

}

