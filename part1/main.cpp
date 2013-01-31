#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <time.h>
#include <armadillo>

#include "lib.h"
#include "variationalmc.h"

using namespace std;
using namespace arma;


//double psi(double* coordinates) {
//    return 2.1;
//}

//double computeEnergy(double* c, double k) {
//    return c[0]*k;
//}


int main() {

    VariationalMC m;
    m.runMetropolis();

    return 0;






//    int     nDimensions     = 3;
//    int     nParticles      = 2;
//    int     nCycles         = 50;
//    int     N               = 5;
//    long    seed            = time(0);
//    double  newWaveFunction = 0.0;
//    double  oldWaveFunction = 0.0;
//    double  energy          = 0.0;
//    double  energySum       = 0.0;
//    double* suggestedMove   = new double[nDimensions * nParticles];
//    double* coordinates     = new double[nDimensions * nParticles];
//    double* newCoordinates  = new double[nDimensions * nParticles];

//    // Fill coordinates array with random values.
//    for (int i = 0; i < nDimensions * nParticles; i++) {
//        coordinates[i] = ran0(&seed);
//    }

//    // Metropolis.
//    for (int i = 0; i < N; i++) {

//        // Suggest change in state, i.e. change in variables x1,x2,..,xn.
//        for (int j = 0; j < nDimensions * nParticles; j++) {
//            suggestedMove[j] = ran0(&seed);
//        }

//        // Compute new wavefunction.
//        for (int j = 0; j < nDimensions * nParticles; j++) {
//            newCoordinates[j] = coordinates[j] + suggestedMove[j];
//        }
//        newWaveFunction = psi(newCoordinates);

//        // Metropolis step -- chose whether to accept step.
//        if (ran0(&seed) <= (newWaveFunction*newWaveFunction) / (oldWaveFunction*oldWaveFunction)) {

//            // Reject new move.
//            energySum += energy;
//        } else {

//            // Accept new move.
//            energy          = computeEnergy(newCoordinates, newWaveFunction);
//            oldWaveFunction = newWaveFunction;
//            coordinates     = newCoordinates;
//        }
//    }

//    VariationalMC m(1,2,3,4,5);
}


