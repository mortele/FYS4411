#include <armadillo>
#include <time.h>
#include <iomanip>
#include <cmath>

#include "variationalmc.h"
#include "wavefunctions.h"
#include "lib.cpp"


using namespace std;
using namespace arma;

/* VMC constructor. */
VariationalMC::VariationalMC() :
    nParticles  (8),
    nDimensions (3),
    nCycles     (50000),
    N       (2 * nCycles / 10),
    idum    (time(0)),
    h       (0.00001),
    h2      (h * h),
    alph    (1.0),
    alph2   (alph * alph),
    beta    (1.0),
    MolecDist (4.6),
//    Z       (nParticles),
    Z       (nParticles/2), //for molecules
    stepSize(0.1),
    D       (0.5),
    dt      (0.01), // 0.0007
    dx      (zeros(nDimensions)),
    spins   (zeros(nParticles,nParticles)) {
}


/* Runs the Metropolis algorithm nCycles times. */
vec VariationalMC::runMetropolis(double alpha, double beta) {
    this->alph  = alpha;
    this->alph2 = alph * alph;
    this->beta  = beta;

    mat correlationsOld = zeros<mat>(nParticles, nParticles);
    mat correlationsNew = zeros<mat>(nParticles, nParticles);
    mat coordinatesNew  = zeros<mat>(nParticles, nDimensions);
    mat coordinatesOld  = zeros<mat>(nParticles, nDimensions);
//    mat Rnew            = zeros<mat>(nParticles, nParticles);   // Matrix of distances and magnitudes.
//    mat Rold            = zeros<mat>(nParticles, nParticles);
    mat Rnew            = zeros<mat>(nParticles+1, nParticles);   // Matrix of distances and magnitudes. For molecules
    mat Rold            = zeros<mat>(nParticles+1, nParticles);
    mat slaterOldUp     = zeros<mat>(nParticles/2, nParticles/2);
    mat slaterOldDown   = zeros<mat>(nParticles/2, nParticles/2);
    mat slaterNewUp     = zeros<mat>(nParticles/2, nParticles/2);
    mat slaterNewDown   = zeros<mat>(nParticles/2, nParticles/2);
    mat slaterGradient  = zeros<mat>(nParticles, nDimensions);
    mat slaterGradientOld   = zeros<mat>(nParticles, nDimensions);
    mat jastrowGradient     = zeros<mat>(nParticles, nParticles);
    mat jastrowGradientOld  = zeros<mat>(nParticles, nParticles);
    mat jastrowLaplacian    = zeros<mat>(nParticles, nParticles);
    mat jastrowLaplacianOld = zeros<mat>(nParticles, nParticles);
    mat quantumForceOld     = zeros<mat>(nParticles, nDimensions);
    mat quantumForceNew     = zeros<mat>(nParticles, nDimensions);

    vec variationalGrad     = zeros<vec>(2);
    vec variationalGradE    = zeros<vec>(2);
    vec variationalGradSum  = zeros<vec>(2);
    vec variationalGradESum = zeros<vec>(2);

    double ecoeff          = 0.0;
    double Rsd             = 0.0;
    double R               = 0.0;
    double newWaveFunction = 0.0;
    double oldWaveFunction = 0.0;

    double correlationOld  = 0.0;
    double correlationNew  = 0.0;

    double energy          = 0.0;
    double energy2         = 0.0;

    double energySum       = 0.0;
    double energy2Sum      = 0.0;

    double energyUp        = 0.0;
    double energyDown      = 0.0;
    double energyPot       = 0.0;
    double energycrossterm = 0.0;
    double energyJas       = 0.0;

    double greensFunction  = 0.0;

    int    accepted        = 0;
    int    t               = 0;

    double randI;
    int    iRand;

    // Fills the spin matrix with relevant values for Beryllium.
    fillSpinMatrix(this->spins);

    // Fill coordinates matrix with random values.
    for (int i = 0; i < nParticles; i++) {
        for (int j = 0; j < nDimensions; j++) {
            coordinatesNew(i,j) = gaussian_deviate(&idum) * 10.0/ (alph);
            coordinatesOld(i,j) = coordinatesNew(i,j);
        }
    }

//    // Updates distances and magnitudes in this initial state.
//    updateRmatrix(coordinatesNew, Rnew);
//    updateRmatrix(coordinatesOld, Rold);


    for(int i=0; i<nParticles;i++) {
        updateForDerivative(Rnew, coordinatesNew,i);
        updateForDerivative(Rold, coordinatesNew,i);
    }

    // Evaluate the slater determinand, and find its inverse.
    evaluateSlater(slaterOldUp, Rnew, coordinatesNew, 0);  // Up --> k=0.
    evaluateSlater(slaterOldDown, Rnew, coordinatesNew, 1); // Down --> k=1.

    // Invert the slater derterminands.
    slaterOldUp   = slaterOldUp.i();
    slaterOldDown = slaterOldDown.i();
    slaterNewUp   = slaterOldUp;
    slaterNewDown = slaterOldDown;


    // Compute the correlations between particles.
    fillCorrelationsMatrix(correlationsOld, spins, Rnew);
    fillCorrelationsMatrix(correlationsNew, spins, Rnew);



    energyPot = computePotentialEnergyClosedForm(Rnew);
    energy = energyUp + energyDown + energyPot;

    updateVariationalGradient(variationalGrad,
                              variationalGradE,
                              Rnew,
                              slaterNewUp,
                              slaterNewDown,
                              correlationsNew,
                              energy,
                              variationalGradSum,
                              variationalGradESum);



    // Compute the correlation factor in the initial state.
    correlationOld  = computeCorrelation(Rnew);



    for (int i=0; i<nParticles; i++) {
        if (i>=nParticles/2) {
            computeSlaterGradient(Rnew, coordinatesNew,slaterOldDown, slaterGradient,1, i);
        }
        else {
            computeSlaterGradient(Rnew, coordinatesNew,  slaterOldUp, slaterGradient,1, i);
        }

        computeJastrowGradient(Rnew, jastrowGradient, i);
        computeJastrowLaplacian(Rnew, jastrowLaplacian, i);
    }
    jastrowGradientOld = jastrowGradient;
    jastrowLaplacianOld = jastrowLaplacian;
    slaterGradientOld = slaterGradient;
    computeQuantumForce(quantumForceNew, Rnew,coordinatesNew, jastrowGradient, slaterGradient, energycrossterm);
    computeQuantumForce(quantumForceOld, Rnew, coordinatesNew, jastrowGradient, slaterGradient, energycrossterm);
//    cout << Rnew << endl;

    // Compute the intial energy.
    energyUp   = computeKineticEnergyClosedForm(Rnew,coordinatesNew,slaterOldUp, 0);
    energyDown = computeKineticEnergyClosedForm(Rnew,coordinatesNew,slaterOldDown, 1);
    energyJas  = computeJastrowEnergy(Rnew, jastrowLaplacian, jastrowGradient);
    energyPot = computePotentialEnergyClosedForm(Rnew);
    energy = energyUp + energyDown + energyPot + energyJas + energycrossterm;;
    //cout << energy << endl;

    // Metropolis loop.
    for (int k = 0; k < nCycles; k++) {

        // Suggest new positions for a single particle, i.e. new state.
        randI = ran0(&idum) * nParticles;
        iRand = floor(randI);

        for (int j = 0; j < nDimensions; j++) {
            // Brute force way:
//            coordinatesNew(iRand,j) += (ran0(&idum)-0.5) * stepSize;

            // Importance sampled way:
            coordinatesNew(iRand, j) += gaussian_deviate(&idum) * sqrt(dt) +
                                        quantumForceOld(iRand,j) * dt * D; // * 2 * D;
        }


        // Compute the wavefunction in this new state.
        updateForDerivative(Rnew, coordinatesNew,iRand); // updates R.


        // Compute the quantum force in this new state.
//        quantumForceNew = computeQuantumForce(Rnew, coordinatesNew, newWaveFunction);

//        // Compute the inside of the exponential term of the difference between Greens functions.
//        greensFunction = 0.0;
//        for (int i = 0; i < nParticles; i++) {
//            for (int j = 0; j < nDimensions; j++) {
//                greensFunction += 0.5 * (quantumForceOld(nDimensions*i+j) +
//                                  quantumForceNew(nDimensions*i+j)) *
//                                  (D * dt * 0.5 * (quantumForceOld(nDimensions*i+j) -
//                                  quantumForceNew(nDimensions*i+j)) -
//                                  coordinatesNew(i,j) + coordinatesOld(i,j));
//            }
//        }
//        // Compute the fraction GreensF(new) / GreensF(old).
//        greensFunction = exp(greensFunction);


        updateCorrelationsMatrix(correlationsNew, Rnew, spins, iRand);
        double Rc = computeRc(correlationsOld, correlationsNew, iRand);


        //cout << Rc << endl;

        if (iRand >= (nParticles/2)) {
            // Down part.
            Rsd = computeSlaterRatio(slaterOldDown, Rnew,coordinatesNew, iRand, iRand-(nParticles/2));
            computeSlaterGradient(Rnew, coordinatesNew,slaterOldDown, slaterGradient,Rsd, iRand);
        } else {
            // Up part.
            Rsd = computeSlaterRatio(slaterOldUp, Rnew,coordinatesNew, iRand, iRand);
            computeSlaterGradient(Rnew, coordinatesNew,slaterOldUp, slaterGradient,Rsd, iRand);
        }
        double screwyou = 0;
        computeJastrowGradient(Rnew, jastrowGradient, iRand);
        computeJastrowLaplacian(Rnew, jastrowLaplacian, iRand);
        computeQuantumForce(quantumForceNew, Rnew, coordinatesNew, jastrowGradient, slaterGradient, energycrossterm);
        computeQuantumForce(quantumForceOld, Rold, coordinatesOld, jastrowGradient, slaterGradientOld, screwyou);

        // Compute the inside of the exponential term of the difference between Greens functions.
        greensFunction = 0.0;

        for (int i = 0;i < nParticles; i++) {
            for (int j = 0; j < nDimensions; j++) {
                greensFunction += 0.5* (quantumForceOld(i,j) +
                                        quantumForceNew(i,j)) *
                        (D * dt * 0.5* (quantumForceOld(i,j) -
                                        quantumForceNew(i,j)) -
                         coordinatesNew(i,j) + coordinatesOld(i,j));
                //cout << greensFunction << endl;
            }
        }
        //cout << "Rold : \n \n \n"<< Rold << endl;
        // Compute the fraction GreensF(new) / GreensF(old).
        greensFunction = exp(greensFunction);
//        cout <<greensFunction << endl;
        R= Rsd*Rc;
        // Check if the suggested move is accepted
        ecoeff =  R * R * greensFunction;
//        if (greensFunction<1e-5) {
//            cout << iRand << endl;
//            cout << greensFunction << endl;
//            cout << slaterGradient << endl;
//            cout << quantumForceNew << endl;
//            cout << coordinatesNew << endl;
//            cout << coordinatesOld << endl;
//            cout << slaterGradientOld << endl;
//            cout << quantumForceOld << endl;
//        }

        if (ecoeff > ran0(&idum)) {
            // Accept new step, calculate new energy.
            accepted++;
            coordinatesOld.row(iRand) = coordinatesNew.row(iRand);
            quantumForceOld = quantumForceNew;           
            correlationsOld  = correlationsNew;
            slaterGradientOld = slaterGradient;
            jastrowGradientOld = jastrowGradient;
            jastrowLaplacianOld = jastrowLaplacian;


            // Check which slater determinand we need to change -- up or down.
            if (iRand >= (nParticles/2)) {
                // Down part.
                updateSlaterInverse(slaterNewDown, slaterOldDown, Rnew, Rold, coordinatesNew, iRand -nParticles/2, 1, Rsd);
                slaterOldDown = slaterNewDown;
            } else {
                // Up part.
                updateSlaterInverse(slaterNewUp, slaterOldUp, Rnew, Rold, coordinatesNew, iRand, 0, Rsd);
                slaterOldUp = slaterNewUp;
//                mat slatertest = zeros<mat>(nParticles/2, nParticles/2);
//                evaluateSlater(slatertest, Rnew,0);
//                cout << " identity : "<<det(slaterNewUp*slatertest) << endl;
//                if (abs(det(slaterNewUp*slatertest)-1)>1e-12) {
//                    cout << det(slaterNewUp*slatertest) << endl;
//                    exit(1);
//                }
            }
            Rold = Rnew;
            double energy1;
            // Compute the energy.
            energyUp   = computeKineticEnergyClosedForm(Rnew,coordinatesNew,slaterNewUp, 0);
            energyDown = computeKineticEnergyClosedForm(Rnew,coordinatesNew,slaterNewDown, 1);
            energyPot  = computePotentialEnergyClosedForm(Rnew);
            energyJas  = computeJastrowEnergy(Rnew, jastrowLaplacian, jastrowGradient);

            energy    =  energyUp + energyDown + energyPot + energyJas + energycrossterm; //
            //energy     = computeEnergy(Rnew, coordinatesNew, R);
            //cout << energyUp << " " << energyDown << " " << energyPot << endl;

            //energy = energyUp + energyDown + energyPot + energyJas + energycrossterm;
            //energy = computeEnergy(Rnew, coordinatesNew, R);

            updateVariationalGradient(variationalGrad,
                                      variationalGradE,
                                      Rnew,
                                      slaterNewUp,
                                      slaterNewDown,
                                      correlationsNew,
                                      energy,
                                      variationalGradSum,
                                      variationalGradESum);

        } else {


            // Reject suggested step, energy remains as before.
            coordinatesNew.row(iRand) = coordinatesOld.row(iRand);
            correlationsNew = correlationsOld;
            slaterGradient = slaterGradientOld;
            jastrowGradient = jastrowGradientOld;
            jastrowLaplacian = jastrowLaplacianOld;
            quantumForceNew = quantumForceOld;
            Rnew = Rold;


            updateVariationalGradientSum(variationalGradSum,variationalGradESum,variationalGrad,variationalGradE);

            //cout << "not accepted bitch!!" << endl;
            // Check which slater determinand we need to change -- up or down.
//            if (iRand >= (nParticles/2)) {
//                // Down part.
//               slaterNewDown = slaterOldDown;
//            } else {
//                // Up part.
//                slaterNewUp = slaterOldUp;
//            }
        }


        // Add energy of this state to the energy sum.
        energySum  += energy;
        energy2Sum += energy * energy;
        //cout << energyUp + energyDown << endl;

        // Throw away the first N samples.
        if (k == N) {
            energySum  = 0.0;
            energy2Sum = 0.0;
            accepted   = 0;
            variationalGradSum  = zeros(2);
            variationalGradESum = zeros(2);
        }
    }

    // Calculate the expected value of the energy, the energy squared, and the variance.
    energy  = energySum  / ((double) (nCycles - N-1));
    energy2 = energy2Sum / ((double) (nCycles - N-1));


    variationalGradSum  /= ((double) (nCycles - N-1));
    variationalGradESum /= ((double) (nCycles - N-1));

    variationalGrad = 2*(variationalGradESum - energy * variationalGradSum);

    cout << "<E>  = " << setprecision(15) << energy << endl;
    cout << "<E²> = " << setprecision(15) << energy2 << endl;
    cout << "Variance  = " << (energy2 - energy*energy)/((double) nCycles) << endl;
    cout << "Std. dev. = " << sqrt((energy2 - energy*energy)/sqrt(nCycles)) << endl;
    cout << "Accepted steps / total steps = " << ((double) accepted) / (nCycles - N-1) << endl;
    cout << "Total steps, nCycles = " << nCycles << endl;

    vec returnVec = zeros<vec>(3);
    returnVec(0) = energy;
    returnVec(1) = variationalGrad(0);
    returnVec(2) = variationalGrad(1);
    return returnVec;
}


/* Old version of compute psi. */
double VariationalMC::computePsi2(const mat &R) {
    double returnVal = exp(-alph * ( R(0,0) + R(1,1) ) + R(0,1) / (2 * (1  + beta * R(0,1))));
    return returnVal;
}


/* New version of compute psi. Slater determinands, bitches! */
/*
double VariationalMC::computePsi(const mat &R) {
    // Compute the correlation part of the wave function.
    double correlation = 1.0; //computeCorrelation(R);

    // Compute the Slater determinand part of the wave function.
    double detup   = psi_s1(R(0,0)) * psi_s2(R(1,1)) - psi_s1(R(1,1)) * psi_s2(R(0,0));
    double detdown = psi_s1(R(2,2)) * psi_s2(R(3,3)) - psi_s1(R(3,3)) * psi_s2(R(2,2));

    return detup * detdown * correlation;
}
*/

/* Computes the local energy of a state defined by position matrix r, and distance matrix R.
 * EL = 1/psi * H * psi */
double VariationalMC::computeEnergy(mat &R, mat &r, double psi) {

    double b1 = beta * R(0,1);
    double b2 = 1 + b1;
    double b3 = 1/(2 * b2 * b2);
    double prikk = r(0,0) * r(1,0) +  r(0,1) * r(1,1) + r(0,2) * r(1,2);
    double E_L1 = (alph - Z) * (1 / R(0,0)  + 1 / R(1,1)) + 1 / R(0,1) - alph2; // (alph - Z) +  + 1 / R(0,1)
    double E_L2 = E_L1  + b3 * ((alph * (R(0,0) + R(1,1))) / (R(0,1))  * (1  - (prikk / (R(0,0) * R(1,1)))) - b3 - 2 / R(0,1) + ((2*beta) / b2)); //
    //double E_L2 =  b3 * ((alph * (R(0,0) + R(1,1))) / (R(0,1))  * (1 - prikk / (R(0,0) * R(1,1)))); //
    return E_L2;
}

//We need to add energycrossterm and the correlation part
double VariationalMC::computeKineticEnergyClosedForm( mat& R, mat& r, const mat& slater, int k) {
    double returnVal = 0.0;
    int o = nParticles/2*k;
    for (int i = 0; i < nParticles/2; i++) {
        for (int j = 0; j < nParticles/2; j++) {
            returnVal += psiDoubleDerivative(R,r,i+o, j) * slater(j,i);
        }
    }
    return returnVal / (-2.0);
}

//double VariationalMC::computePotentialEnergyClosedForm(const mat& R) {
//    double returnVal = 0.0;
//    double E1 = 0.0;
//    // Compute the commutative part of H.
//    for(int i = 0; i < nParticles; i++) {
//        returnVal -= Z/R(i,i);

//        for(int j = i+1; j<nParticles; j++) {
//            E1 += 1/R(i,j);
//        }
//    }

//    return returnVal; // + E1;
//}

double VariationalMC::computePotentialEnergyClosedForm(const mat& R) {
    double returnVal = 0.0;
    double E1 = 0.0;
    returnVal = Z*Z/MolecDist;
    // Compute the commutative part of H.
    for(int i = 0; i < nParticles; i++) {
        returnVal -= Z/R(i,i);
        returnVal -= Z/R(i+1,i);

        for(int j = i+1; j<nParticles; j++) {
            E1 += 1/R(i,j);
        }
    }

    return returnVal + E1;
}

/* Computes the local energy, by numerical differentiation, of a state defined by position matrix r, and distance matrix R.
 * EL = 1/psi * H * psi */
/*
double VariationalMC::computeEnergyNumerical(mat &R, mat &r, double psi) {
    double psil, psih;
    double E1  = 0;
    double E2  = 0;
    double psi2, psi3;

    // Compute the commutative part of H.
    for(int i = 0;i<nParticles; i++) {
        E1 -= Z/R(i,i);
        //for(int j = i+1; j<nParticles; j++) {
        //    E1 += 1/R(i,j);
        //}
    }

    for(int i = 0; i<nParticles;i++){
        for(int j = 0; j<nDimensions;j++){
            psi3 = computePsi(R);
            r(i,j) += h;
            updateForDerivative(R, r, i);
            psih = computePsi(R);

            r(i,j) -= 2 * h;
            updateForDerivative(R, r, i);
            psil = computePsi(R);

            // Reset R back to oringial positions.
            r(i,j) += h;
            updateForDerivative(R, r, i);

            psi2 = computePsi(R);
            E2 -= computeDoubleDerivative(psil, psi2, psih);
        }
    }
    if (fabs(E2 / (2 * psi) + E1)>80) {
        cout << "hei";
        //return -20;
    }
    return E2 / (2 * psi) + E1;
}
*/

/* Computes a numerical approximation to the double derivative of psi. */
double VariationalMC::computeDoubleDerivative(double psiLow, double psi,double psiHigh) {
    return (psiLow - 2 * psi + psiHigh) / h2;
}

double VariationalMC::computeFirstDerivative(double psiLow, double psiHigh) {
    return (psiHigh - psiLow) / (2 * h);
}


/* Takes two uniformely distributed random numbers rand1, and rand2, and transforms
 * them into two gaussian distributed random numbers. */
double VariationalMC::gaussian_deviate(long * idum) {
    static int iset = 0;
    static double gset;
    double fac, rsq, v1, v2;

    if ( idum < 0) iset =0;
    if (iset == 0) {
        do {
            v1 = 2.*ran0(idum) -1.0;
            v2 = 2.*ran0(idum) -1.0;
            rsq = v1*v1+v2*v2;
        } while (rsq >= 1.0 || rsq == 0.);
        fac = sqrt(-2.*log(rsq)/rsq);
        gset = v1*fac;
        iset = 1;
        return v2*fac;
    } else {
        iset =0;
        return gset;
    }
} // end function for gaussian deviates


void VariationalMC::updateRmatrix(const mat &r, mat &R) {

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


///* Updates the distance matrix when we have just changed one coordinate of particle "i"
// * (like we do when computing the derivative)*/
//void VariationalMC::updateForDerivative(mat &R, const mat &r, int i){

//    double dxx, sum;
//    for(int k=0; k<i; k++) {
//        sum = 0;
//        for(int l =0;l<nDimensions;l++){
//            dxx = r(i,l) - r(k,l);
//            sum += dxx*dxx;
//        }
//        R(k,i) = sqrt(sum); //R is the matrix of distances
//    }

//    for(int k=i+1;k<nParticles;k++){
//        sum = 0;
//        for(int l =0;l<nDimensions;l++){
//            dxx = r(i,l) - r(k,l);;
//            sum += dxx*dxx;
//        }
//        R(i,k) = sqrt(sum); //R is the matrix of distances
//    }

//    sum = 0;
//    for(int l =0;l<nDimensions;l++){
//        dxx = r(i,l); //r[l+i*nDimensions]*r[l+i*nDimensions];
//        sum += dxx*dxx;
//    }
//    R(i,i) = sqrt(sum);
//}

/* Updates the distance matrix when we have just changed one coordinate of particle "i"
 * (like we do when computing the derivative)*/

//version for calculating molecules
void VariationalMC::updateForDerivative(mat &R, const mat &r, int i){
    double dxx, sum;
    double dist2 = MolecDist/2.0;
    for(int k=0; k<i; k++) {
        sum = 0;
        for(int l =0;l<nDimensions;l++){
            dxx = r(i,l) - r(k,l);
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
        if (l==0) {
            sum +=(dxx-dist2)*(dxx-dist2);
        } else {
            sum += dxx*dxx;
        }
    }
    R(i,i) = sqrt(sum);
    sum = 0;
    for(int l =0;l<nDimensions;l++){
        dxx = r(i,l); //r[l+i*nDimensions]*r[l+i*nDimensions];
        if (l==0) {
            sum +=(dxx+dist2)*(dxx+dist2);
        } else {
            sum += dxx*dxx;
        }
    }
    R(i+1,i) = sqrt(sum);
}




/* Computes the quantum force. Uses numerical differentiation to find the gradient of Psi. */
//vec VariationalMC::computeQuantumForce(mat &R, mat &r, double psi) {

//    vec gradient(nDimensions*nParticles);
//    double psiHigh, psiLow;

//    for (int k = 0; k < nParticles; k++) {
//        for (int o = 0; o < nDimensions; o++) {

//            // Compute psi(R - dR).
//            r(k,o) -= h;
//            updateForDerivative(R, r, k);
//            psiLow = computePsi(R);

//            // Compute psi(R + dR).
//            r(k,o) += 2 * h;
//            updateForDerivative(R, r, k);
//            psiHigh = computePsi(R);

//            // Return to original configuration, R.
//            r(k,o) -= h;
//            updateForDerivative(R, r, k);

//            gradient(nDimensions*k+o) = computeFirstDerivative(psiLow, psiHigh);
//        }
//    }
//    return gradient;
//}
void VariationalMC::computeQuantumForce(mat & QuantumForce, mat &R, mat &r, mat & jastrowGradient, mat & slaterGradient, double & energycrossterm) {
    energycrossterm = 0.0;
    for (int k = 0; k<nParticles; k++) {
        for (int j =0 ; j<nDimensions; j++) {//we call the function jastrowgradient far too many times here, this should be done more effectively!!!!!
            double sum = 0;
            for (int i = 0; i<k; i++) {
                sum += (r(k,j)-r(i,j))/R(i,k)*jastrowGradient(i,k);
            }
            for (int i = k+1; i < nParticles; i++) {
                sum -= (r(i,j)-r(k,j))/R(k,i)*jastrowGradient(k,i);
            }

            QuantumForce(k,j) = 2*(slaterGradient(k,j) + sum);
//            QuantumForce(k,j) = 2*slaterGradient(k,j);
            energycrossterm  -=  sum*sum/2.0 + (slaterGradient(k,j) * sum); //;
        }
    }
}



/*
The hydrogenlike wavefunctions and their derivatives:
*/



double VariationalMC::psi_s1(double distance){
    return exp(-alph*distance);
}

double VariationalMC::psi_s2(double distance){
    return (1-alph*distance/2.)*exp(-alph*distance/2);
}

double VariationalMC::psi_p2z(double distance, mat &r1, int i){
    double r = distance;
//    double z = r1(i,2);
    return r1(i,2)*exp(-alph*r/2);
}

//double VariationalMC::psi_p2z(double distance, mat &r1, int i){
//    double r = distance;
//    double z = r1(i,2);
//    return r*z*exp(-alph*r/2)/12;
//}

double VariationalMC::psi_p2y(double distance, mat &r1, int i){
    double r = distance;
//    double y = r1(i,1);
    return r1(i,1)*exp(-alph*r/2);
}

double VariationalMC::psi_p2x(double distance, mat &r1, int i){
    double r = distance;
//    double x = r1(i,0);
    return r1(i,0)*exp(-alph*r/2);
}

double VariationalMC::psi_s1_derivative(double distance, double coord) { //distance is r, coord is x,y or z
    return -alph*coord*exp(-alph*distance)/distance;
}

double VariationalMC::psi_s2_derivative(double distance,double coord) { //distance is r, coord is x,y or z
    return 1/4.*alph*(alph*distance - 4)*coord *exp(-alph*distance/2)/distance;
}

double VariationalMC::psi_p2z_derivative(double distance, mat &r1, int i, int j) { //distance is r, i is which particle, j is which coordinate
    double r = distance;
    double x = r1(i,0);
    double y = r1(i,1);
    double z = r1(i,2);
    if (j == 0) {
        return -alph*x*z*exp(-alph*r/2)/(2*r);
    }
    else if (j == 1) {
        return -alph*y*z*exp(-alph*r/2)/(2*r);
    }
    else if (j == 2){
        return -(alph*pow(z, 2) - 2*r)*exp(-alph*r/2)/(2*r);
    }
}

//double VariationalMC::psi_p2z_derivative(double distance, mat &r1, int i, int j) { //distance is r, coord is x,y or z
//    double r = distance;
//    double x = r1(i,0);
//    double y = r1(i,1);
//    double z = r1(i,2);
//    if (j == 0) {
//        return -x*z*(alph*r - 2)*exp(-alph*r/2)/(24*r);
//    }
//    else if (j == 1) {
//        return -y*z*(alph*r - 2)*exp(-alph*r/2)/(24*r);
//    }
//    else if (j == 2){
//        return -(alph*r*pow(z, 2) - 2*pow(r, 2) - 2*pow(z, 2))*exp(-alph*r/2)/(24*r);
//    }
//}




double VariationalMC::psi_p2y_derivative(double distance, mat &r1, int i, int j) { //distance is r, i is  which particle, j is which coordinate
    double r = distance;
    double x = r1(i,0);
    double y = r1(i,1);
    double z = r1(i,2);
    if (j == 0) {
        return -alph*x*y*exp(-alph*r/2)/(2*r);
    }
    else if (j == 1) {
        return -(alph*pow(y, 2) - 2*r)*exp(-alph*r/2)/(2*r);
    }
    else if (j == 2){
        return -alph*y*z*exp(-alph*r/2)/(2*r);
    }
}


double VariationalMC::psi_p2x_derivative(double distance, mat &r1, int i, int j) { //distance is r, i is  which particle, j is which coordinate
    double r = distance;
    double x = r1(i,0);
    double y = r1(i,1);
    double z = r1(i,2);
    if (j == 0) {
        return -(alph*pow(x, 2) - 2*r)*exp(-alph*r/2)/(2*r);
    }
    else if (j == 1) {
        return -alph*x*y*exp(-alph*r/2)/(2*r);
    }
    else if (j == 2){
        return -alph*x*z*exp(-alph*r/2)/(2*r);
    }
}


double VariationalMC::psi_s1_doubleDerivative(double distance) {
    return (alph2 - 2 * alph / distance) * exp(-alph * distance);
}

double VariationalMC::psi_s2_doubleDerivative(double distance) {
    return (5.0/4.0 * alph2 - 2 * alph / distance - alph2*alph * distance / 8.) * exp(-alph * distance / 2.);
}

double VariationalMC::psi_p2z_doubleDerivative(double distance, mat &r1, int i){ //we can save alot here since this is the same expression for each of the coordinates
    double r = distance;
    double z = r1(i,2);
    return alph*z*(alph*r - 8)*exp(-alph*r/2)/(4*r);
}
//double VariationalMC::psi_p2z_doubleDerivative(double distance, mat &r1, int i){ //we can save alot here since this is the same expression for each of the coordinates
//    double r = distance;
//    double z = r1(i,2);
//    return z*(alph*alph*r*r - 12*alph*r + 16)*exp(-alph*r/2)/(48*r);
//}

double VariationalMC::psi_p2y_doubleDerivative(double distance, mat &r1, int i){
    double r = distance;
    double y = r1(i,1);
    return alph*y*(alph*r - 8)*exp(-alph*r/2)/(4*r);
}

double VariationalMC::psi_p2x_doubleDerivative(double distance, mat &r1, int i){
    double r = distance;
    double x = r1(i,0);
    return alph*x*(alph*r - 8)*exp(-alph*r/2)/(4*r);
}




void VariationalMC::computeSlaterGradient(mat &Rnew,
                                          mat &r,
                                          mat& slater,
                                          mat& gradient,
                                          double R,
                                          int particle) {

    int k = particle %(nParticles / 2); //the remainder
    for (int o = 0; o < nDimensions; o++) {
        double sum = 0;
        for (int j = 0; j < (nParticles / 2); j++) {
            sum += psiDerivative(Rnew,r,particle, o, j)*slater(j,k);
        }
        gradient(particle,o) = 1.0/R * sum;
//        if (abs(gradient(particle,o))> 20 ) {
//        cout << gradient(particle,o) << endl;
//        cout << slater << endl;
//        cout << Rnew << endl;
//        cout << r << endl;
//        }
    }
}


double VariationalMC::psiDerivative(mat &R, mat &r, int i, int j, int k) {
    double distance = R(i,i);
    double distance2 = R(i+1,i);
    double R2 = MolecDist/2;
    if (k == 0) {
        if (j==0)
            return psi_s1_derivative(distance,r(i,j)-R2) + psi_s1_derivative(distance2,r(i,j)+R2);
        else
            return psi_s1_derivative(distance,r(i,j)) + psi_s1_derivative(distance2,r(i,j));

    } else if (k == 1) {
        if (j==0)
            return psi_s1_derivative(distance,r(i,j)-R2) - psi_s1_derivative(distance2,r(i,j)+R2);
        else
            return psi_s1_derivative(distance,r(i,j)) - psi_s1_derivative(distance2,r(i,j));
    } else if (k == 2) {
        if (j==0)
            return psi_s2_derivative(distance,r(i,j)-R2) + psi_s2_derivative(distance2,r(i,j)+R2);
        else
            return psi_s2_derivative(distance,r(i,j)) + psi_s2_derivative(distance2,r(i,j));
    } else if (k == 3) {
        if (j==0)
            return psi_s2_derivative(distance,r(i,j)-R2) - psi_s2_derivative(distance2,r(i,j)+R2);
        else
            return psi_s2_derivative(distance,r(i,j)) - psi_s2_derivative(distance2,r(i,j));
    } else if (k == 4) {
        return psi_p2z_derivative(distance,r,i,j);
    } else {
        return 0;
        cout << "Error in psiDerivative: k = " << k  << endl;
    }
}

double VariationalMC::psiDoubleDerivative(mat &R, mat &r, int i, int j) {
    double distance = R(i,i);
    double distance2 = R(i+1,i);
    if (j == 0) {
        return psi_s1_doubleDerivative(distance) + psi_s1_doubleDerivative(distance2);
    } else if (j == 1) {
        return psi_s1_doubleDerivative(distance) - psi_s1_doubleDerivative(distance2);
    } else if (j == 2) {
        return psi_s2_doubleDerivative(distance) + psi_s2_doubleDerivative(distance2);
    } else if (j == 3) {
        return psi_s2_doubleDerivative(distance) - psi_s2_doubleDerivative(distance2);
    } else if (j == 4) {
        return psi_p2z_doubleDerivative(distance,r,i);
    } else {
        return 0;
        cout << "Error in psiDoubleDerivative: j = " <<  j  << endl;
    }
}

//double VariationalMC::psiDerivative(mat &R, mat &r, int i, int j, int k) {
//    double distance = R(i,i);
//    if (k == 0) {

//        return psi_s1_derivative(distance,r(i,j));

//    } else if (k == 1) {
//        return psi_s2_derivative(distance,r(i,j));
//    } else if (k == 2) {
//        return psi_p2x_derivative(distance,r,i,j);
//    } else if (k == 3) {
//        return psi_p2y_derivative(distance,r,i,j);
//    } else if (k == 4) {
//        return psi_p2z_derivative(distance,r,i,j);
//    } else {
//        return 0;
//        cout << "Error in psiDerivative: k = " << k  << endl;
//    }
//}

//double VariationalMC::psiDoubleDerivative(mat &R, mat &r, int i, int j) {
//    double distance = R(i,i);

//    if (j == 0) {
//        return psi_s1_doubleDerivative(distance);
//    } else if (j == 1) {
//        return psi_s2_doubleDerivative(distance);
//    } else if (j == 2) {
//        return psi_p2x_doubleDerivative(distance,r,i);
//    } else if (j == 3) {
//        return psi_p2y_doubleDerivative(distance,r,i);
//    } else if (j == 4) {
//        return psi_p2z_doubleDerivative(distance,r,i);
//    } else {
//        return 0;
//        cout << "Error in psiDoubleDerivative: j = " <<  j  << endl;
//    }
//}


double VariationalMC::computeSlaterRatio(const mat& slaterInverseOld, mat &R, mat &r, int i, int p) {
    double Rsd = 0;
    for (int j = 0; j < (nParticles / 2); j++) {
        Rsd += slaterPsi(R,r,i, j) * slaterInverseOld(j,p);
    }
    return Rsd;
}

void VariationalMC::computeJastrowGradient(const mat& R, mat& jastrowGradient, int particle) {
    for (int i = 0; i<particle; i++) {
        double factor = 1+beta*R(i,particle);
        jastrowGradient(i,particle) = spins(i,particle)/(factor*factor);
    }
    for (int i = particle+1; i < nParticles; i++) {
        double factor = 1+beta*R(particle,i);
        jastrowGradient(particle,i) = spins(particle,i)/(factor*factor);
    }
}

void VariationalMC::computeJastrowLaplacian(const mat& R, mat& jastrowLaplacian, int particle) {
    for (int i = 0; i<particle; i++) {
        double factor = 1 + beta*R(i,particle);
        jastrowLaplacian(i,particle) = -2*spins(i,particle)*beta/(factor*factor*factor);
    }
    for (int i = particle+1; i < nParticles; i++) {
        double factor = 1 + beta*R(particle,i);
        jastrowLaplacian(particle,i) = -2*spins(particle,i)*beta/(factor*factor*factor);
    }
}

double VariationalMC::computeJastrowEnergy(const mat& R, mat& jastrowLaplacian, mat& jastrowGradient) {
    double sum = 0.0;
    for (int k = 0; k<nParticles; k++) {
        for (int i = 0; i<k; i++) {
            sum +=(nDimensions -1)/R(i,k)*jastrowGradient(i,k) +jastrowLaplacian(i,k);
        }
        for (int i = k+1; i < nParticles; i++) {
            sum +=(nDimensions -1)/R(k,i)*jastrowGradient(k,i) +jastrowLaplacian(k,i);
        }
    }
    return sum/(-2.0);
}

///* Choses which wave function to calculate from the int j, then calls the corresponding psi_j. */
//double VariationalMC::slaterPsi(mat &R, mat &r, int i, int j) {
//    double distance = R(i,i);
//    if (j == 0) {
//        return psi_s1(distance);
//    } else if (j == 1) {
//        return psi_s2(distance);
//    } else if (j == 2) {
//        return psi_p2x(distance,r,i);
//    } else if (j == 3) {
//        return psi_p2y(distance,r,i);
//    } else if (j == 4) {
//        return psi_p2z(distance,r,i);
//    } else {
//        return 0;
//        cout << "Error in slaterPsi" << j<< endl;
//    }
//}

/* Choses which wave function to calculate from the int j, then calls the corresponding psi_j. */ //FOR MOLECULES
double VariationalMC::slaterPsi(mat &R, mat &r, int i, int j) {
    double distance = R(i,i);
    double distance2 = R(i+1,i);
    if (j == 0) {
        return psi_s1(distance) + psi_s1(distance2);
    } else if (j == 1) {
        return psi_s1(distance) - psi_s1(distance2);
    } else if (j == 2) {
        return psi_s2(distance) + psi_s2(distance2);
    } else if (j == 3) {
        return psi_s2(distance) - psi_s2(distance2);
    } else if (j == 4) {
        return psi_p2z(distance,r,i);
    } else {
        return 0;
        cout << "Error in slaterPsi" << j<< endl;
    }
}




/* Computes the correlation factor of the wave function. */
double VariationalMC::computeCorrelation(const mat& R) {
    // Compute the correlation part of the wave function.
    double correlation = 0.0;
    for (int i = 0; i < nParticles; i++) {
        for (int j = (i+1); j < nParticles; j++) {
            correlation += spins(i,j) * R(i,j) / (1 + beta * R(i,j));
        }
    }
    correlation = exp(correlation);
    return correlation;
}


/* Updates the inverse of one part of the slater determinand (up, or down part). */
void VariationalMC::updateSlaterInverse(mat&        slaterNew,
                                        const mat&  slaterOld,
                                        mat&  Rnew,
                                        const mat&  Rold,
                                        mat&        r,
                                        int         i,
                                        int         p,
                                        double      R) {
    int o = i+nParticles/2*p;
    for (int k = 0;k<nParticles/2; k++) {
        for (int j = 0; j<nParticles/2; j++) {
            if(j != i) {
                double sum = 0;
                for( int l = 0; l< nParticles/2; l++) {
                    sum += slaterOld(l,j) * slaterPsi(Rnew,r,o, l);
                }
                slaterNew(k,j) = slaterOld(k,j) - slaterOld(k,i) * sum / R;
            }
            else {
                slaterNew(k,j) = slaterOld(k,i)/ R;
            }
        }
    }
}

// Evaluates the slater determinand. k=0 means the spin up determinand, k=1 is spin down;
void VariationalMC::evaluateSlater(mat& slater, mat& R,mat &r, int k) {
    // Offset is N/2 for the spin down determinand, 0 for spin up.
    int offset = (k*nParticles/2);

    for (int i = 0; i < (nParticles/2); i++) {
        for (int j = 0; j < (nParticles/2); j++) {
            slater(i,j) = slaterPsi(R,r,i+offset, j);
        }
    }
}



/* Fills the matrix of spin values, which is used to compute the correlation
part of the wavefunction. */
void VariationalMC::fillSpinMatrix(mat& spin) {
    // Set which electrons have spins up, and which have spins down.
    vec electronSpins(nParticles);
    electronSpins.zeros();
    for (int i = 0; i < (nParticles / 2); i++) {
        electronSpins(i) = 1;
    }

    // Set the electron-electron spin interaction matrix.
    for (int i = 0; i < nParticles; i++) {
        for (int j = i+1; j < nParticles; j++) {
            if (electronSpins(i) != electronSpins(j)) {
                spin(i,j) = 0.5;
            } else {
                spin(i,j) = 0.25;
            }
        }
    }
}



/* Compute all the correlation factors. */
void VariationalMC::fillCorrelationsMatrix( mat& correlationsMat,
                                            const mat& spins,
                                            const mat& R) {

  for (int i = 0; i < nParticles; i++) {
    for (int j = i+1; j < nParticles; j++) {
      correlationsMat(i,j) = spins(i,j) * R(i,j) / (1 + beta * R(i,j));
    }
  }
}


/* Computes the ratio of the new correlations to the old, Rc. */
double VariationalMC::computeRc( const mat& correlationsOld,
                                 const mat& correlationsNew,
                                 int particle) {

  int    k  = particle;
  double Rc = 0;

  for (int i = 0; i < k; i++) {
    Rc += correlationsNew(i,k) - correlationsOld(i,k);
  }
  for (int i = k+1;i<nParticles; i++) {
    Rc += correlationsNew(k,i) - correlationsOld(k,i);
  }
  return exp(Rc);
}



/* Updates the new correlations matrix. */
void VariationalMC::updateCorrelationsMatrix( mat& correlationsMat,
                                              const mat& R,
                                              const mat& spins,
                                              int  particle) {

  int k = particle;

  for (int i = 0; i < k; i++) {
    correlationsMat(i,k) = spins(i,k) * R(i,k) / (1 + beta * R(i,k));
  }
  for (int i = (k+1); i < nParticles; i++) {
    correlationsMat(k,i) = spins(k,i) * R(k,i) / (1 + beta * R(k,i));
  }
}


/*
 * Computes the derivative of the logarithm of the jastrow part of the total wave function
 * with regards to beta.
 */
double VariationalMC::computeJastrowBetaDerivative(const mat& R,
                                                   const mat& correlationsMat) {
  double sum = 0;
  for (int i = 0; i < nParticles; i++) {
    for (int j = i+1; j < nParticles; j++) {
        sum -= correlationsMat(i,j) * R(i,j) / (1 + beta * R(i,j));
    }
  }
  return sum;
}


/*
 * Computes the derivative of the the slater determinands with regards to alpha.
 */
double VariationalMC::computeSlaterAlphaDerivative(const mat& R,
                                                   const mat& slaterInvUp,
                                                   const mat& slaterInvDown) {

  double sum = 0;
  for (int i = 0; i < nParticles/2; i++) {
    for (int j = 0; j < nParticles/2; j++) {
      sum += slaterAlphaDerivative(i, R(j,j)) * slaterInvUp(i,j);
    }
  }
  for (int i = 0; i < nParticles/2; i++) {
    for (int j = 0; j < nParticles/2; j++) {
        sum += slaterAlphaDerivative(i, R(j+nParticles/2,j+nParticles/2)) * slaterInvDown(i,j);
    }
  }
  return sum;
}


/*
 * Choses which wave function to calculate the alpha derivative of from the int
 * j, then calls the corresponding psi_alphaDerivative.
 */
double VariationalMC::slaterAlphaDerivative(int j, double r) {
    if (j == 0) {
        return psi_s1_alphaDerivative(r);
    } else if (j == 1) {
        return psi_s2_alphaDerivative(r);
    } else {
        return 0;
        cout << "Error in slaterAlphaDerivative" << j << endl;
    }
}


double VariationalMC::psi_s1_alphaDerivative(double r) {
    return -r * exp(-alph * r);
}


double VariationalMC::psi_s2_alphaDerivative(double r) {
    return (r * r * alph / 4.0 - r) * exp(- alph * r / 2.0);
}


/*
 * Updates the gradient of the trial wave function with regards to the variational
 * parameters.
 */
void VariationalMC::updateVariationalGradient(mat&        varGradient,
                                              mat&        varGradientE,
                                              const mat&  R,
                                              const mat&  slaterInvUp,
                                              const mat&  slaterInvDown,
                                              const mat&  correlationsMat,
                                              double      E,
                                              mat&        varGradientSum,
                                              mat&        varGradientESum) {
  varGradient(0)   = computeSlaterAlphaDerivative(R, slaterInvUp, slaterInvDown);
  varGradient(1)   = computeJastrowBetaDerivative(R, correlationsMat);
  varGradientE(0)  = varGradient(0) * E;
  varGradientE(1)  = varGradient(1) * E;

  // Updates the sums of the gradient and the gradient * E.
  updateVariationalGradientSum(varGradientSum, varGradientESum, varGradient, varGradient);
}


/*
 * Updates the sum of the gradients of the trial wave function with regards
 * to the variational parameters, alpha and beta.
 */
void VariationalMC::updateVariationalGradientSum(mat& varGradientSum,
                                                 mat& varGradientESum,
                                                 const mat& varGradient,
                                                 const mat& varGradientE) {
  varGradientSum  += varGradient;
  varGradientESum += varGradientE;
}
