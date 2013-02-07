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
    nCycles     (2000000),
    N       (2000000),
    idum    (10),
    charge  (2.0),
    h       (0.001),
    h2      (0.000001),
    alph    (0.5*charge),
    beta    (1.0),
    Z       (2.0),
    stepSize(0.5) {
}

/* Runs the Metropolis algorithm nCycles times. */
void VariationalMC::runMetropolis() {

    mat coordinatesNew = zeros<mat>(nParticles, nDimensions);
    mat coordinatesOld = zeros<mat>(nParticles, nDimensions);
    mat Rnew           = zeros<mat>(nParticles, nDimensions);   // Matrix of distances and magnitudes.
    mat Rold           = zeros<mat>(nParticles, nDimensions);

    double ecoeff          = 0.0;
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

    // Updates distances and magnitudes in this initial state.
    updateRmatrix(coordinatesNew, Rnew);
    updateRmatrix(coordinatesOld, Rold);

    // Compute the wave function in this initial state.
    oldWaveFunction = computePsi(Rnew);

    // Metropolis loop.
    for (int k = 0; k < nCycles; k++) {

        // Suggest new positions for all particles, i.e. new state.
        for (int i = 0; i < nParticles; i++) {
            for (int j = 0; j < nDimensions; j++) {
                coordinatesNew(i,j) += ran0(&idum) * stepSize;
            }
        }

        // Compute the wavefunction in this new state.
        updateRmatrix(coordinatesNew, Rnew);
        newWaveFunction = computePsi(Rnew);

        // Check if the suggested move is accepted.
        ecoeff = newWaveFunction / oldWaveFunction;
        if (ecoeff>ran0(&idum)) {   
            coordinatesOld = coordinatesNew;

            // Energy changes from previous state.
            energy = computeEnergy(Rnew, coordinatesNew, newWaveFunction);

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
double VariationalMC::computePsi(const mat &R) {

    double returnVal = 0.0;

    for (int i = 0; i < nParticles; i++) {
        for (int j = (i + 1); j < nParticles; j++) {
            returnVal += (1 / R(i,j)) * exp(-alph * (R(i,i) + R(j,j))) * (exp(R(i,j)) / (1 + beta * R(i,j)));
        }
    }
    return returnVal;
}




/* Computes the local energy of a state defined by position matrix r, and distance matrix R. EL = 1/psi * H * psi */
double VariationalMC::computeEnergy(mat &R, mat &r, double psi)
{
    double r12 = R(0,1);
    double r1 = R(0,0);
    double r2 = R(1,1);
    double E1 = (-Z*(1/r1 +1/r2) +1/r12);   //this is the commutative part of the hamiltonian (we can just multiply it with psi.)
    double E2 = 0;

    /*
                              for(int i=0; i< nDimensions*nParticles;i++){
                                  coordinates[i]-=h;
                                  psil = psi(coordinates);
                                  coordinates[i]+=2*h;
                                  psih = psi(coordinates);
                                  coordinates[i]-=h;
                                  E2-=derivative(psil, psi,psih);
                              }
                              */

    vec oldR(nParticles);
    double psil, psih;

    for(int i = 0; i<nParticles;i++){
        /*        
        r1 = R(i,0);

        for(int k=0; k<i; k++){
            oldR(k) = R(i,k); //R is the matrix of distances
        }


        for(int k = i + 1; k < nParticles; k++) { // nParticles
            oldR(k) = R(i,k);
        }
        */
        for(int j = 0; j<nDimensions;j++){
            r(i,j)-=h; //r is the array of coordinates

            updateForDerivative(R,r, i);
            psil = computePsi(R);

            r(i,j)+=2*h;
            updateForDerivative(R, r, i);
            psih = computePsi(R);
            r(i,j)-=h;
            E2-=computeDoubleDerivative(psil, psi,psih);
            updateForDerivative(R, r, i);   // set all values back to normal
        }
    }

    return E2/psi +E1;

}

/* Computes a numerical approximation to the double derivative of psi. */
double VariationalMC::computeDoubleDerivative(double psiLow, double psi,double psiHigh) {
    return (psiLow - 2 * psi + psiHigh) / h2;

}

void VariationalMC::updateRmatrix(const mat &r, mat &R) {

    for (int i = 0; i < nParticles; i++) {

        // Compute magnitude of r(i,:) position vector.
        R(i,i) = 0.0;
        for (int k = 0; k < nDimensions; k ++) {
            R(i,i) += r(i,k) * r(i,k);
        }
        R(i,i) = sqrt(R(i,i));

        for (int j = (i + 1); j < nParticles; j++) {
            R(j,j) = R(i,j) = 0.0;

            // Compute magnitude of r(i,:) - r(j,:) position vector.
            R(i,j) = 0.0;
            for (int k = 0; k < nDimensions; k++) {
                R(i,j) += (r(i,k) - r(j,k)) * (r(i,k) - r(j,k));
            }
            R(i,j) = sqrt(R(i,j));
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


/*Updates the distance matrix when we have just changed one coordinate of particle "i" (like we do when computing the derivative)*/
void VariationalMC::updateForDerivative(mat &R, const mat &r, int i){
    vec dx(nDimensions);
    dx.zeros();
    double dxx;

    for(int k=0; k<i; k++) {
        for(int l =0;l<nDimensions;l++){
            dxx = r(i,l) - r(k,l); // [l+i*nDimensions]-r[l+k*nDimensions];   //this may need to be changed to r[i,l] - r[k,l]
                                                          // (likewise for the next loops)
            dx(l) = dxx*dxx;
        }

        R(i,k) = sqrt(sum(dx)); //R is the matrix of distances
        dx.zeros();
    }


    for(int k=i+1;k<nParticles;k++){
        for(int l =0;l<nDimensions;l++){

            dxx = r(i,l) - r(k,l);;
            dx(l) = dxx * dxx;

        }

        R(i,k) = sqrt(sum(dx)); //R is the matrix of distances
        dx.zeros();
    }

    for(int l =0;l<nDimensions;l++){
        dxx = r(i,l); //r[l+i*nDimensions]*r[l+i*nDimensions];
        dx(l) = dxx*dxx;
    }

    R(i,i) = sqrt(sum(dx));
    dx.zeros();
}

