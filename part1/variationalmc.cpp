#include <armadillo>
#include <time.h>

#include "variationalmc.h"
#include "lib.cpp"


using namespace std;
using namespace arma;

/* VMC constructor. */
VariationalMC::VariationalMC() :
    nParticles  (2),
    nDimensions (3),
    nCycles     (50),
    N       (200),
    idum    (10),
    charge  (2.0),
    h       (0.001),
    h2      (1000000),
    alph    (0.5*charge),
    beta    (20.0),
    Z       (20.0),
    stepSize(1.0) {
}

/* Runs the Metropolis algorithm nCycles times. */
void VariationalMC::runMetropolis() {

    mat coordinatesNew = zeros<mat>(nParticles, nDimensions);
    mat coordinatesOld = zeros<mat>(nParticles, nDimensions);

    double newWaveFunction = 0.0;
    double oldWaveFunction = 0.0;

    double energy          = 0.0;
    double energy2         = 0.0;
    double energySum       = 0.0;
    double energy2Sum      = 0.0;


    // Fill coordinates arrays with random values.
    for (int i = 0; i < nParticles; i++) {
        for (int j = 0; j < nDimensions; j++) {
            coordinatesNew(i,j) = coordinatesOld(i,j) = ran0(&idum);
        }
    }

    // Compute the wave function in this initial state.
    oldWaveFunction = computePsi(coordinatesNew);

    // Metropolis loop.
    for (int k = 0; k < nCycles; k++) {

        // Suggest new positions for all particles, i.e. new state.
        for (int i = 0; i < nParticles; i++) {
            for (int j = 0; j < nDimensions; j++) {
                coordinatesNew(i,j) += ran0(&idum) * stepSize;
            }
        }

        // Compute the wavefunction in this new state.
        newWaveFunction = computePsi(coordinatesNew);

        // Check if the suggested move is accepted.
        if (1==1) {
            coordinatesOld = coordinatesNew;

            // Energy changes from previous state.
            energy = computeEnergy(coordinatesNew);
        } else {
            coordinatesNew = coordinatesOld;

            // Energy remains unchanged.
        }

        // Add energy of this state to the energy sum.
        energySum  += energy;
        energy2Sum += energy * energy;
    }

    // Calculate the expected value of the energy, the energy squared, and the variance.
    energy  = energySum  / nCycles;
    energy2 = energy2Sum / nCycles;

    cout << "<E>  = " << energy << endl;
    cout << "<E²> = " << energy2 << endl;

}


/* Computes the wavefunction in a state defined by position matrix r. */
double VariationalMC::computePsi(const mat &r) {

    double r_ij, r_ii, r_jj;
    double returnVal = 0.0;

    for (int i = 0; i < nParticles; i++) {

        // Compute magnitude of r(i,:) position vector.
        r_ii = 0.0;
        for (int k = 0; k < nDimensions; k ++) {
            r_ii += r(i,k) * r(i,k);
        }
        r_ii = sqrt(r_ii);

            r_jj = r_ij = 0.0;
            for (int j = (i + 1); j < nParticles; j++) {

                // Compute magnitude of r(i,:) - r(j,:) position vector.
                r_ij = 0.0;
                for (int k = 0; k < nDimensions; k++) {
                    r_ij += (r(i,k) - r(j,k)) * (r(i,k) - r(j,k));
                }
                r_ij = sqrt(r_ij);

                // Compute magnitude of r(j,:) position vector.
                r_jj = 0.0;
                for (int k = 0; k < nDimensions; k ++) {
                    r_jj += r(j,k) * r(j,k);
                }
                r_jj = sqrt(r_jj);

                returnVal += (1 / r_ij) * exp(-alph * (r_ii + r_jj)) * (exp(r_ij) / (1 + beta * r_ij));
            }
        }

    return returnVal;

    // Pseudo-code outline of function:
    //
    //    for i=0..nParticles
    //        compute r_ii
    //
    //            for j = i+1..nParticles
    //                compute r_ij, and
    //                compute r_jj
    //
    //                compute contribution to wavefunction from interaction of particles i,j
    //
    //    return sum of all contributions to wavefunction
}




/* Computes the energy of a state defined by position matrix r. */
double VariationalMC::computeEnergy(const mat &r) {
    return r(1,0);
}

/* Computes a numerical approximation to the double derivative of psi. */
double VariationalMC::computeDoubleDerivative(double psiLow, double psi,double psiHigh) {
    return (psiLow - 2 * psi + psiHigh) / h2;

}

void VariationalMC::updateRmatrix(const mat &r, mat &R) {

    for (int i = 0; i < nParticles; i++) {

        // Compute magnitude of r(i,:) position vector.
        R(i,i) = 0.0; // = r_ii
        for (int k = 0; k < nDimensions; k ++) {
            R(i,i) += r(i,k) * r(i,k); // = r_ii
        }
        R(i,i) = sqrt(R(i,i));

        R(j,j) = R(i,j) = 0.0;
        for (int j = (i + 1); j < nParticles; j++) {

            // Compute magnitude of r(i,:) - r(j,:) position vector.
            r_ij = 0.0;
            for (int k = 0; k < nDimensions; k++) {
                r_ij += (r(i,k) - r(j,k)) * (r(i,k) - r(j,k));
            }
            r_ij = sqrt(r_ij);

            // Compute magnitude of r(j,:) position vector.
            r_jj = 0.0;
            for (int k = 0; k < nDimensions; k ++) {
                r_jj += r(j,k) * r(j,k);
            }
            r_jj = sqrt(r_jj);
        }
    }

    // Pseudo-code outline of function:
    //
    //    for i=0..nParticles
    //        compute R_ii
    //
    //            for j = i+1..nParticles
    //                compute R_ij
}


