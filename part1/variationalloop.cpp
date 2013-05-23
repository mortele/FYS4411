#include <iostream>
#include <fstream>
#include <iomanip>
#include "variationalloop.h"


VariationalLoop::VariationalLoop() {
}

void VariationalLoop::initialize_helium() {
    this->N     = 40;
    this->start = 2.0;
    this->end   = 4.0;
    this->E     = 0.0;
    this->minE  = 1e300;
    this->minA  = 0.0;
    this->minB  = 0.0;
}

void VariationalLoop::run() {
    double startClock, finishClock, energy, accepted;
    vec    energyVarGrad = zeros<vec>(4);
    vec    varGrad       = zeros<vec>(2);
    vec    alphaBeta     = zeros<vec>(2);
    alphaBeta(0) = 3.7;
    alphaBeta(1) = 0.2;

    for (int i = 1; i < 20000; i++) {
        startClock = clock();

        // Run one metropolis loop.
        energyVarGrad = m.runMetropolis(alphaBeta(0), alphaBeta(1));
        energy        = energyVarGrad(0);
        varGrad(0)    = energyVarGrad(1);
        varGrad(1)    = energyVarGrad(2);
        accepted      = energyVarGrad(3);

        if (accepted > 0.95) {
            // Compute new values of alpha / beta.
            alphaBeta     += (1/((double) i + 100)) * varGrad * energy;
            //alphaBeta(0)   = 3.7;
            //alphaBeta(1) = alphaBeta(1) + 1/((double) i + 10) * varGrad(1) * energy;
        } else {
            cout << "\n\n\n\n\n\naccepted  = "<< accepted << " < 0.95 \n\n\n\n\n";
        }


        /*
        for (int k = 0; k < N; k++) {
            for (int j = 0; j < N; j++) {
                startClock = clock();
                b = start + ((end-start)/N) * k;
                a = start + ((end-start)/N) * j;
                b = 0.6;
                a = 3.6;

                //alphaBeta(1) -=  (1/((double) i)) * 0.01 * varGrad(1) * energy;


                if (alphaBeta(0) < 0) {
                    alphaBeta(0) = 0.01;
                } if (alphaBeta(1) < 0) {
                    alphaBeta(1) = 0.01;
                }
        */
        finishClock = clock();
        cout << "   * Time usage       = " << (finishClock - startClock) / 1000000.0 << " [s] "<< endl;
        cout << "   * Alpha            = " << alphaBeta(0) << endl;
        cout << "   * Beta             = " << alphaBeta(1) << endl;
        cout << "VarGrad=" << varGrad << endl;
        cout << "delta alphaBeta =" << (1/((double) i + 100)) * varGrad * energy << endl;
        cout << "   * i                = " << i << endl;
        cout << endl << endl;
        }

}




//  ofstream outFile("data_varying_a_and_b__N_100_points_from_0_9_to_2_0.dat");
//    for (int k = 0; k < N; k++) {
//        for (int j = 0; j < N; j++) {
//            startClock = clock();
//            b = start + ((end-start)/N) * k;
//            a = start + ((end-start)/N) * j;
//            b = 0.5;
//            a = 1.8;

//            cout << "/-------------------------------------------------------------\\" << endl;
//            printf("|   For parameters:    a = %8.4f,   and    b = %8.4f    |\n", a, b);
//            cout << "\\-------------------------------------------------------------/" << endl;

//            E = m.runMetropolis(a,b);
//            //acceptanceRatio = m.acceptanceRatio;

//            //outFile << setprecision(15) << a << " " << b << " " << E << " " << acceptanceRatio <<endl;

////            if (E < minE) {
////                minE = E;
////                minA = a;
////                minB = b;
//                //            cout << "minA=" << minA << endl;
//                //            cout << "minB=" << minB << endl;
//                //            cout << "minE=" << minE << endl;
////            }
//            finishClock = clock();
//            cout << "   * Time usage       = " << (finishClock - startClock) / 1000000.0 << " [s] "<< endl;
//            cout << endl << endl;
//        }
//    }

//    double Eexperimental = -2.9037;
//    cout << "/-------------------------------------------------------------\\" << endl;
//    printf("|                     Best parameters found                    |\n");
//    cout << "\\-------------------------------------------------------------/" << endl;
//    cout << "   * a                     = " << minA << endl;
//    cout << "   * b                     = " << minB << endl;
//    cout << "   * With corresponding E  = " << minE << endl;
//    cout << "   * Experimental value E' = " << Eexperimental << endl;
//    cout << "   * |E - E'|              = " << abs(minE - Eexperimental) << endl << endl << endl;

//    outFile.close();

