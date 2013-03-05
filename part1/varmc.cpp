#include <armadillo>
#include <time.h>

#include "varmc.h"
#include "lib.h"


using namespace std;
using namespace arma;

/* VMC constructor. */
VarMC::VarMC() :
    nParticles  (2),
    nDimensions (3),
    nCycles     (100000000),
    N       (nCycles / 10),
    idum    (17),
    charge  (2.0),
    h       (0.001),
    h2      (h * h),
    alph    (1.0),
    alph2   (alph * alph),
    beta    (1.0),
    Z       (2.0),
    stepSize(0.07),
    dx      (zeros(nDimensions)),
    acceptanceRatio(0.0) {
}

/* Runs the Metropolis algorithm nCycles times. */
double VarMC::runMetropolis(double alpha, double beta) {
    this->alph = alpha;
    this->alph2 = alph * alph;
    this->beta = beta;

    mat coordinatesNew = zeros<mat>(nParticles, nDimensions);
    mat coordinatesOld = zeros<mat>(nParticles, nDimensions);
    mat Rnew           = zeros<mat>(nParticles, nParticles);   // Matrix of distances and magnitudes.
    mat Rold           = zeros<mat>(nParticles, nParticles);

    double ecoeff          = 0.0;
    double newWaveFunction = 0.0;
    double oldWaveFunction = 0.0;

    double energy          = 0.0;
    double energy2         = 0.0;

    double energySum       = 0.0;
    double energy2Sum      = 0.0;

    int    accepted        = 0;

    double randI;
    int    iRand;

//    mat test  = zeros<mat>(2,3);
//    mat test2 = zeros<mat>(2,3);

//    test(0,0) = 1;
//    test(0,1) = 2;
//    test(0,2) = 3;
//    test(1,0) = 4;
//    test(1,1) = 5;
//    test(1,2) = 6;

//    updateForDerivative(test2, test, 1);

//    cout << test2 << endl;

    // Fill coordinates arrays with random values.
    for (int i = 0; i < nParticles; i++) {
        for (int j = 0; j < nDimensions; j++) {
            coordinatesNew(i,j) = coordinatesOld(i,j) = (ran0(&idum)-0.5) / (0.5*alph);
        }
    }

    // Updates distances and magnitudes in this initial state.
    updateRmatrix(coordinatesNew, Rnew);
    updateRmatrix(coordinatesOld, Rold);

    // Compute the wave function in this initial state.
    oldWaveFunction = computePsi(Rnew);

    // The Metropolis loop itself.
    for (int k = 0; k < nCycles; k++) {

        // Suggest new positions for one particle, i.e. new state.
        randI = ran0(&idum) * nParticles;
        iRand = floor(randI);


        for (int j = 0; j < nDimensions; j++) {
            coordinatesNew(iRand,j) += (ran0(&idum)-0.5) * stepSize;
        }

        // Compute the wavefunction in this new state.
        updateForDerivative(Rnew, coordinatesNew,iRand);
        newWaveFunction = computePsi(Rnew);

        // Check if the suggested move is accepted.
        ecoeff = newWaveFunction * newWaveFunction / (oldWaveFunction * oldWaveFunction);
        if (ecoeff > ran0(&idum)) {
            accepted++;
            coordinatesOld  = coordinatesNew;
            oldWaveFunction = newWaveFunction;

            // c = a
            // a = b
            // b = c

            // Energy changes from previous state.
            energy = computeEnergy(Rnew, coordinatesNew, newWaveFunction);
        } else {
            coordinatesNew  = coordinatesOld;

            // Energy remains unchanged.
        }

        // Add energy of this state to the energy sum.
        energySum  += energy;
        energy2Sum += energy * energy;

        if (k == N) {
            energySum  = 0.0;
            energy2Sum = 0.0;
        }
    }

    // Calculate the expected value of the energy, the energy squared, and the variance.
    energy  = energySum  / (nCycles - N);
    energy2 = energy2Sum / (nCycles - N);

    cout << "   * <E>              = " << energy << endl;
    cout << "   * <E²>             = " << energy2 << endl;
    cout << "   * Variance         = " << (energy2 - energy*energy)/(nCycles-N) << endl;
    cout << "   * Std. dev.        = " << sqrt((energy2 - energy*energy)/(nCycles-N)) << endl;
    cout << "   * Acceptance ratio = " << ((double) accepted) / nCycles << endl;
    this->acceptanceRatio = ((double) accepted) / nCycles;

    return energy;
}


/* Computes the wavefunction in a state defined by position matrix r. */
double VarMC::computePsi(const mat &R) {

    double returnVal = 0.0;
    returnVal = exp(-alph * (R(0,0) + R(1,1)) + R(0,1) / (2 * (1  + beta * R(0,1)))); // ) * exp(

    //    for (int i = 0; i < nParticles; i++) {
    //        for (int j = (i + 1); j < nParticles; j++) {
    //            returnVal += (1 / R(i,j)) * exp(-alph * (R(i,i) + R(j,j))) * (exp(R(i,j)) / (1 + beta * R(i,j)));
    //        }
    //    }
    return returnVal;
}




/* Computes the local energy of a state defined by position matrix r, and distance matrix R.
 * EL = 1/psi * H * psi */
double VarMC::computeEnergy(mat &R, mat &r, double psi) {
//    double b1 = beta * R(0,1);
//    double b2 = 1 + b1;
//    double b3 = 1/(2 * b2 * b2);
//    double prikk = r(0,0) * r(1,0) +  r(0,1) * r(1,1) + r(0,2) * r(1,2);

//    double E_L1 = (alph - Z) * (1 / R(0,0)  + 1 / R(1,1)) + 1 / R(0,1) - alph2;
//    double E_L2 = E_L1 + b3 * ( (alph * (R(0,0) + R(1,1))) / (R(0,1))  * (1 - (prikk / (R(0,0) * R(1,1)))) - b3 - 2 / R(0,1) + ((2*beta) / b2));
//    return E_L2;

    double b1 = beta * R(0,1);
    double b2 = 1 + b1;
    double b3 = 1/(2 * b2 * b2);
    double prikk = r(0,0) * r(1,0) +  r(0,1) * r(1,1) + r(0,2) * r(1,2);

    double E_L1 = (alph - Z) * (1 / R(0,0)  + 1 / R(1,1)) + 1 / R(0,1) - alph2;
    double E_L2 = E_L1 + b3 * ( (alph * (R(0,0) + R(1,1))) / (R(0,1))  * (1 - (prikk / (R(0,0) * R(1,1)))) - b3 - 2 / R(0,1) + ((2*beta) / b2));
    return E_L2;
}

/* Computes a numerical approximation to the double derivative of psi. */
double VarMC::computeDoubleDerivative(double psiLow, double psi,double psiHigh) {
    return (psiLow - 2 * psi + psiHigh) / h2;
}

void VarMC::updateRmatrix(const mat &r, mat &R) {

    for (int i = 0; i < nParticles; i++) {

        // Compute magnitude of r(i,:) position vector.
        R(i,i) = 0.0;
        for (int k = 0; k < nDimensions; k ++) {
            R(i,i) += r(i,k) * r(i,k);
        }
        R(i,i) = sqrt(R(i,i));

        for (int j = (i + 1); j < nParticles; j++) {
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


/* Updates the distance matrix when we have just changed one coordinate of particle "i"
 * (like we do when computing the derivative)*/
void VarMC::updateForDerivative(mat &R, const mat &r, int i){

    double dxx, sum;
    for(int k=0; k<i; k++) {
        sum = 0;
        for(int l =0;l<nDimensions;l++){
            dxx = r(i,l) - r(k,l); // [l+i*nDimensions]-r[l+k*nDimensions];   //this may need to be changed to r[i,l] - r[k,l]
                                                         // (likewise for the next loops)
            sum += dxx*dxx;
        }

        R(k,i) = sqrt(sum); //R is the matrix of distances

    }


    for(int k=i+1;k<nParticles;k++){
        sum = 0;
        for(int l =0;l<nDimensions;l++){

            dxx = r(i,l) - r(k,l);;
            sum += dxx*dxx;

        }

        R(i,k) = sqrt(sum); //R is the matrix of distances

    }
    sum = 0;
    for(int l =0;l<nDimensions;l++){
        dxx = r(i,l); //r[l+i*nDimensions]*r[l+i*nDimensions];
        sum += dxx*dxx;
    }

    R(i,i) = sqrt(sum);
}

