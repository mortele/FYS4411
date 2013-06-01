

//#include <armadillo>
//#include <time.h>
//#include <iomanip>
//#include <cmath>

//#include "variationalmc.h"
//#include "wavefunctions.h"
//#include "lib.cpp"


//using namespace std;
//using namespace arma;

///* VMC constructor. */
//VariationalMC::VariationalMC():
//    nParticles  (2),
//    nDimensions (3),
//    nCycles     (500000),
//    N       (2 * nCycles / 10),
//    idum    (time(0)),
//    h       (0.00001),
//    h2      (h * h),
//    alph    (1.0),
//    alph2   (alph * alph),
//    beta    (1.0),
//    Z       (nParticles),
//    stepSize(0.01),
//    D       (0.5),
//    dt      (0.01), // 0.0007
//    dx      (zeros(nDimensions)),
//    spins   (zeros(nParticles,nParticles)) {
//}


///* Runs the Metropolis algorithm nCycles times. */
///*
//int VariationalMC::runMetropolis(double alpha, double beta) {
//    this->alph  = alpha;
//    this->alph2 = alph * alph;
//    this->beta  = beta;

//    mat correlationsOld = zeros<mat>(nParticles, nParticles);
//    mat correlationsNew = zeros<mat>(nParticles, nParticles);
//    mat coordinatesNew  = zeros<mat>(nParticles, nDimensions);
//    mat coordinatesOld  = zeros<mat>(nParticles, nDimensions);
//    mat Rnew            = zeros<mat>(nParticles+1, nParticles);   // Matrix of distances and magnitudes.
//    mat Rold            = zeros<mat>(nParticles+1, nParticles);
//    mat quantumForceOld     = zeros<mat>(nParticles, nDimensions);
//    mat quantumForceNew     = zeros<mat>(nParticles, nDimensions);

//    double ecoeff          = 0.0;
//    double Rsd             = 0.0;
//    double R               = 0.0;
//    double newWaveFunction = 0.0;
//    double oldWaveFunction = 0.0;

//    double correlationOld  = 0.0;
//    double correlationNew  = 0.0;

//    double energy          = 0.0;
//    double energy2         = 0.0;

//    double energySum       = 0.0;
//    double energy2Sum      = 0.0;

//    double greensFunction  = 0.0;

//    int    accepted        = 0;
//    int    t               = 0;

//    double randI;
//    int    iRand;


//    // Fill coordinates matrix with random values.
//    for (int i = 0; i < nParticles; i++) {
//        for (int j = 0; j < nDimensions; j++) {
//            coordinatesNew(i,j) = gaussian_deviate(&idum) * 2.0/ (alph);
//            coordinatesOld(i,j) = coordinatesNew(i,j);
//        }
//    }

//    // Updates distances and magnitudes in this initial state.
//    updateRmatrix(coordinatesNew, Rnew);
//    updateRmatrix(coordinatesOld, Rold);



//    energyPot = computePotentialEnergyClosedForm(Rnew);
//    energy = energyUp + energyDown + energyPot;

//    // Compute the correlation factor in the initial state.
//    correlationOld  = computeCorrelation(Rnew);

//    for (int i=0; i<nParticles; i++) {
//        if (i>=nParticles/2) {
//            computeSlaterGradient(Rnew, coordinatesNew,slaterOldDown, slaterGradient,1, i);
//        }
//        else {
//            computeSlaterGradient(Rnew, coordinatesNew,  slaterOldUp, slaterGradient,1, i);
//        }

//        computeJastrowGradient(Rnew, jastrowGradient, i);
//        computeJastrowLaplacian(Rnew, jastrowLaplacian, i);
//    }
//    jastrowGradientOld = jastrowGradient;
//    jastrowLaplacianOld = jastrowLaplacian;
//    slaterGradientOld = slaterGradient;
//    computeQuantumForce(quantumForceNew, Rnew,coordinatesNew, jastrowGradient, slaterGradient, energycrossterm);
//    computeQuantumForce(quantumForceOld, Rnew, coordinatesNew, jastrowGradient, slaterGradient, energycrossterm);



//    // Compute the intial energy.
//    energyKin   = computeKineticEnergyClosedForm(Rnew,coordinatesNew,slaterOldUp, 0);

//    energyPot = computePotentialEnergyClosedForm(Rnew);
//    energy = energyUp + energyDown + energyPot;// + energyJas + energycrossterm;;
//    //cout << energy << endl;

//    // Metropolis loop.
//    for (int k = 0; k < nCycles; k++) {

//        // Suggest new positions for a single particle, i.e. new state.
//        randI = ran0(&idum) * nParticles;
//        iRand = floor(randI);

//        for (int j = 0; j < nDimensions; j++) {
//            // Brute force way:
//            //coordinatesNew(iRand,j) += (ran0(&idum)-0.5) * stepSize;

//            // Importance sampled way:
//            coordinatesNew(iRand, j) += gaussian_deviate(&idum) * sqrt(dt) +
//                                        quantumForceOld(iRand,j) * dt * D; // * 2 * D;
//        }


//        // Compute the wavefunction in this new state.
//        updateForDerivative(Rnew, coordinatesNew,iRand); // updates R.


//        updateCorrelationsMatrix(correlationsNew, Rnew, spins, iRand);
//        double Rc = computeRc(correlationsOld, correlationsNew, iRand);


//        if (iRand >= (nParticles/2)) {
//            // Down part.
//            Rsd = computeSlaterRatio(slaterOldDown, Rnew(iRand,iRand),coordinatesNew, iRand, iRand-(nParticles/2));
//            computeSlaterGradient(Rnew, coordinatesNew,slaterOldDown, slaterGradient,Rsd, iRand);
//        } else {
//            // Up part.
//            Rsd = computeSlaterRatio(slaterOldUp, Rnew(iRand,iRand),coordinatesNew, iRand, iRand);
//            computeSlaterGradient(Rnew, coordinatesNew,slaterOldUp, slaterGradient,Rsd, iRand);
//        }
//        double screwyou = 0;
//        computeJastrowGradient(Rnew, jastrowGradient, iRand);
//        computeJastrowLaplacian(Rnew, jastrowLaplacian, iRand);
//        computeQuantumForce(quantumForceNew, Rnew, coordinatesNew, jastrowGradient, slaterGradient, energycrossterm);
//        computeQuantumForce(quantumForceOld, Rold, coordinatesOld, jastrowGradient, slaterGradientOld, screwyou);

//        // Compute the inside of the exponential term of the difference between Greens functions.
//        greensFunction = 0.0;

//        for (int i = 0;i < nParticles; i++) {
//            for (int j = 0; j < nDimensions; j++) {
//                greensFunction += 0.5* (quantumForceOld(i,j) +
//                                        quantumForceNew(i,j)) *
//                        (D * dt * 0.5* (quantumForceOld(i,j) -
//                                        quantumForceNew(i,j)) -
//                         coordinatesNew(i,j) + coordinatesOld(i,j));
//                //cout << greensFunction << endl;
//            }
//        }
//        //cout << Rnew << endl;
//        // Compute the fraction GreensF(new) / GreensF(old).
//        greensFunction = exp(greensFunction);
//        //cout << greensFunction << endl;
//        R= Rsd;//*Rc;
//        // Check if the suggested move is accepted
//        ecoeff =  R * R * greensFunction;
////        if (greensFunction<1e-5) {
////            cout << iRand << endl;
////            cout << greensFunction << endl;
////            cout << slaterGradient << endl;
////            cout << quantumForceNew << endl;
////            cout << coordinatesNew << endl;
////            cout << coordinatesOld << endl;
////            cout << slaterGradientOld << endl;
////            cout << quantumForceOld << endl;
////        }

//        if (ecoeff > ran0(&idum)) {
//            // Accept new step, calculate new energy.
//            accepted++;
//            coordinatesOld.row(iRand) = coordinatesNew.row(iRand);
//            quantumForceOld = quantumForceNew;
//            correlationsOld  = correlationsNew;
//            slaterGradientOld = slaterGradient;
//            jastrowGradientOld = jastrowGradient;
//            jastrowLaplacianOld = jastrowLaplacian;


//            // Check which slater determinand we need to change -- up or down.
//            if (iRand >= (nParticles/2)) {
//                // Down part.
//                updateSlaterInverse(slaterNewDown, slaterOldDown, Rnew, Rold, coordinatesNew, iRand -nParticles/2, 1, Rsd);
//                slaterOldDown = slaterNewDown;
//            } else {
//                // Up part.
//                updateSlaterInverse(slaterNewUp, slaterOldUp, Rnew, Rold, coordinatesNew, iRand, 0, Rsd);
//                slaterOldUp = slaterNewUp;
//            }
//            Rold = Rnew;
//            double energy1;
//            // Compute the energy.
//            energyUp   = computeKineticEnergyClosedForm(Rnew,coordinatesNew,slaterNewUp, 0);
//            energyDown = computeKineticEnergyClosedForm(Rnew,coordinatesNew,slaterNewDown, 1);
//            energyPot  = computePotentialEnergyClosedForm(Rnew);
//            energyJas  = computeJastrowEnergy(Rnew, jastrowLaplacian, jastrowGradient);

//            energy    =  energyUp + energyDown + energyPot; // + energyJas + energycrossterm; //


//            updateVariationalGradient(variationalGrad,
//                                      variationalGradE,
//                                      Rnew,
//                                      slaterNewUp,
//                                      slaterNewDown,
//                                      correlationsNew,
//                                      energy,
//                                      variationalGradSum,
//                                      variationalGradESum);

//        } else {


//            // Reject suggested step, energy remains as before.
//            coordinatesNew.row(iRand) = coordinatesOld.row(iRand);
//            correlationsNew = correlationsOld;
//            slaterGradient = slaterGradientOld;
//            jastrowGradient = jastrowGradientOld;
//            jastrowLaplacian = jastrowLaplacianOld;
//            quantumForceNew = quantumForceOld;
//            Rnew = Rold;


//            updateVariationalGradientSum(variationalGradSum,variationalGradESum,variationalGrad,variationalGradE);
//        }



//        // Add energy of this state to the energy sum.
//        energySum  += energy;
//        energy2Sum += energy * energy;
//        //cout << energyUp + energyDown << endl;

//        // Throw away the first N samples.
//        if (k == N) {
//            energySum  = 0.0;
//            energy2Sum = 0.0;
//            accepted   = 0;
//            variationalGradSum  = zeros(2);
//            variationalGradESum = zeros(2);
//        }
//    }

//    // Calculate the expected value of the energy, the energy squared, and the variance.
//    energy  = energySum  / ((double) (nCycles - N-1));
//    energy2 = energy2Sum / ((double) (nCycles - N-1));


//    cout << "<E>  = " << setprecision(15) << energy << endl;
//    cout << "<E²> = " << setprecision(15) << energy2 << endl;
//    cout << "Variance  = " << (energy2 - energy*energy)/((double) nCycles) << endl;
//    cout << "Std. dev. = " << sqrt((energy2 - energy*energy)/sqrt(nCycles)) << endl;
//    cout << "Accepted steps / total steps = " << ((double) accepted) / (nCycles - N-1) << endl;
//    cout << "Total steps, nCycles = " << nCycles << endl;

//    return 0;
//}


///* Old version of compute psi. */

//double VariationalMC::computePsiH2(const mat &R) {
//    R12 = R(0,1);
//    R1P1 = R(0,0);
//    R1P2 = R(1,0);
//    R2P1 = R(1,1);
//    R2P2 = R(2,1);
//    double returnVal = (exp(-alph * R1P1 ) + exp(-alph*R1P2))*(exp(-alph * R2P1 ) + exp(-alph*R2P2))*(exp( R12/ (2 * (1  + beta * R12))));
//    return returnVal;
//}


////We need to add energycrossterm and the correlation part
//double VariationalMC::computeKineticEnergyClosedForm(const mat& R, mat& r, const mat& slater, int k) {
//    double returnVal = 0.0;
//    int o = nParticles/2*k;
//    for (int i = 0; i < nParticles/2; i++) {
//        for (int j = 0; j < nParticles/2; j++) {
//            returnVal += psiDoubleDerivative(R(i+o,i+o),r,i, j) * slater(j,i);
//        }
//    }
//    return returnVal / (-2.0);
//}

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

//    return returnVal;// + E1;
//}





///* Takes two uniformely distributed random numbers rand1, and rand2, and transforms
// * them into two gaussian distributed random numbers. */
//double VariationalMC::gaussian_deviate(long * idum) {
//    static int iset = 0;
//    static double gset;
//    double fac, rsq, v1, v2;

//    if ( idum < 0) iset =0;
//    if (iset == 0) {
//        do {
//            v1 = 2.*ran0(idum) -1.0;
//            v2 = 2.*ran0(idum) -1.0;
//            rsq = v1*v1+v2*v2;
//        } while (rsq >= 1.0 || rsq == 0.);
//        fac = sqrt(-2.*log(rsq)/rsq);
//        gset = v1*fac;
//        iset = 1;
//        return v2*fac;
//    } else {
//        iset =0;
//        return gset;
//    }
//} // end function for gaussian deviates


//void VariationalMC::updateRmatrix(const mat &r, mat &R) {

//    for (int i = 0; i < nParticles; i++) {

//        // Compute magnitude of r(i,:) position vector.
//        R(i,i) = 0.0;
//        for (int k = 0; k < nDimensions; k ++) {
//            R(i,i) += r(i,k) * r(i,k);
//        }
//        R(i,i) = sqrt(R(i,i));

//        for (int j = (i + 1); j < nParticles; j++) {

//            // Compute magnitude of r(i,:) - r(j,:) position vector.
//            R(i,j) = 0.0;
//            for (int k = 0; k < nDimensions; k++) {
//                R(i,j) += (r(i,k) - r(j,k)) * (r(i,k) - r(j,k));
//            }
//            R(i,j) = sqrt(R(i,j));
//        }
//    }


//    // Pseudo-code outline of function:
//    //
//    //    for i=0..nParticles
//    //        compute R_ii
//    //
//    //            for j = i+1..nParticles
//    //                compute R_ij
//}


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
//        if (l==0) {
//            sum +=(dxx-dist2)*(dxx-dist2);
//        } else {
//            sum += dxx*dxx;
//        }
//    }
//    R(i,i) = sqrt(sum);
//    sum = 0;
//    for(int l =0;l<nDimensions;l++){
//        dxx = r(i,l); //r[l+i*nDimensions]*r[l+i*nDimensions];
//        if (l==0) {
//            sum +=(dxx+dist2)*(dxx+dist2);
//        } else {
//            sum += dxx*dxx;
//        }
//    }
//    R(i+1,i) = sqrt(sum);
//}


///* Computes the quantum force. Uses numerical differentiation to find the gradient of Psi. */
////vec VariationalMC::computeQuantumForce(mat &R, mat &r, double psi) {

////    vec gradient(nDimensions*nParticles);
////    double psiHigh, psiLow;

////    for (int k = 0; k < nParticles; k++) {
////        for (int o = 0; o < nDimensions; o++) {

////            // Compute psi(R - dR).
////            r(k,o) -= h;
////            updateForDerivative(R, r, k);
////            psiLow = computePsi(R);

////            // Compute psi(R + dR).
////            r(k,o) += 2 * h;
////            updateForDerivative(R, r, k);
////            psiHigh = computePsi(R);

////            // Return to original configuration, R.
////            r(k,o) -= h;
////            updateForDerivative(R, r, k);

////            gradient(nDimensions*k+o) = computeFirstDerivative(psiLow, psiHigh);
////        }
////    }
////    return gradient;
////}
//void VariationalMC::computeQuantumForce(mat & QuantumForce, mat &R, mat &r, mat & jastrowGradient, mat & slaterGradient, double & energycrossterm) {
//    energycrossterm = 0.0;
//    for (int k = 0; k<nParticles; k++) {
//        for (int j =0 ; j<nDimensions; j++) {//we call the function jastrowgradient far too many times here, this should be done more effectively!!!!!
//            double sum = 0;
//            for (int i = 0; i<k; i++) {
//                sum += (r(k,j)-r(i,j))/R(i,k)*jastrowGradient(i,k);
//            }
//            for (int i = k+1; i < nParticles; i++) {
//                sum -= (r(i,j)-r(k,j))/R(k,i)*jastrowGradient(k,i);
//            }
//            //QuantumForce(k,j) = 2*(slaterGradient(k,j) + sum);
//            QuantumForce(k,j) = 2*slaterGradient(k,j);
//            energycrossterm  -=  sum*sum/2.0 + (slaterGradient(k,j) * sum); //;
//        }
//    }
//}



///*
//The hydrogenlike wavefunctions and their derivatives:
//*/



//double VariationalMC::psi_s1(double distance){
//    return exp(-alph*distance);
//}

//double VariationalMC::psi_s2(double distance){
//    return (1-alph*distance/2.)*exp(-alph*distance/2);
//}

//double VariationalMC::psi_p2z(double distance, mat &r1, int i){
//    double r = distance;
//    double z = r1(i,2);
//    return z*exp(-alph*r/2);
//}

////double VariationalMC::psi_p2z(double distance, mat &r1, int i){
////    double r = distance;
////    double z = r1(i,2);
////    return r*z*exp(-alph*r/2)/12;
////}

//double VariationalMC::psi_p2y(double distance, mat &r1, int i){
//    double r = distance;
//    double y = r1(i,1);
//    return y*exp(-alph*r/2);
//}

//double VariationalMC::psi_p2x(double distance, mat &r1, int i){
//    double r = distance;
//    double x = r1(i,0);
//    return x*exp(-alph*r/2);
//}

//double VariationalMC::psi_s1_derivative(double distance, double coord) { //distance is r, coord is x,y or z
//    return -alph*coord*exp(-alph*distance)/distance;
//}

//double VariationalMC::psi_s2_derivative(double distance,double coord) { //distance is r, coord is x,y or z
//    return 1/4.*alph*(alph*distance - 4)*coord *exp(-alph*distance/2)/distance;
//}

//double VariationalMC::psi_p2z_derivative(double distance, mat &r1, int i, int j) { //distance is r, coord is x,y or z
//    double r = distance;
//    double x = r1(i,0);
//    double y = r1(i,1);
//    double z = r1(i,2);
//    if (j == 0) {
//        return -alph*x*z*exp(-alph*r/2)/(2*r);
//    }
//    else if (j == 1) {
//        return -alph*y*z*exp(-alph*r/2)/(2*r);
//    }
//    else if (j == 2){
//        return -(alph*pow(z, 2) - 2*r)*exp(-alph*r/2)/(2*r);
//    }
//}

////double VariationalMC::psi_p2z_derivative(double distance, mat &r1, int i, int j) { //distance is r, coord is x,y or z
////    double r = distance;
////    double x = r1(i,0);
////    double y = r1(i,1);
////    double z = r1(i,2);
////    if (j == 0) {
////        return -x*z*(alph*r - 2)*exp(-alph*r/2)/(24*r);
////    }
////    else if (j == 1) {
////        return -y*z*(alph*r - 2)*exp(-alph*r/2)/(24*r);
////    }
////    else if (j == 2){
////        return -(alph*r*pow(z, 2) - 2*pow(r, 2) - 2*pow(z, 2))*exp(-alph*r/2)/(24*r);
////    }
////}




//double VariationalMC::psi_p2y_derivative(double distance, mat &r1, int i, int j) { //distance is r, coord is x,y or z
//    double r = distance;
//    double x = r1(i,0);
//    double y = r1(i,1);
//    double z = r1(i,2);
//    if (j == 0) {
//        return -alph*x*y*exp(-alph*r/2)/(2*r);
//    }
//    else if (j == 1) {
//        return -(alph*pow(y, 2) - 2*r)*exp(-alph*r/2)/(2*r);
//    }
//    else if (j == 2){
//        return -alph*y*z*exp(-alph*r/2)/(2*r);
//    }
//}

//double VariationalMC::psi_p2x_derivative(double distance, mat &r1, int i, int j) { //distance is r, coord is x,y or z
//    double r = distance;
//    double x = r1(i,0);
//    double y = r1(i,1);
//    double z = r1(i,2);
//    if (j == 0) {
//        return -(alph*pow(x, 2) - 2*r)*exp(-alph*r/2)/(2*r);
//    }

//    else if (j == 1) {
//        return -alph*x*y*exp(-alph*r/2)/(2*r);
//    }
//    else if (j == 2){
//        return -alph*x*z*exp(-alph*r/2)/(2*r);
//    }
//}


//double VariationalMC::psi_s1_doubleDerivative(double distance) {
//    return (alph2 - 2 * alph / distance) * exp(-alph * distance);
//}

//double VariationalMC::psi_s2_doubleDerivative(double distance) {
//    return (5.0/4.0 * alph2 - 2 * alph / distance - alph2*alph * distance / 8.) * exp(-alph * distance / 2.);
//}

//double VariationalMC::psi_p2z_doubleDerivative(double distance, mat &r1, int i){ //we can save alot here since this is the same expression for each of the coordinates
//    double r = distance;
//    double z = r1(i,2);
//    return alph*z*(alph*r - 8)*exp(-alph*r/2)/(4*r);
//}
////double VariationalMC::psi_p2z_doubleDerivative(double distance, mat &r1, int i){ //we can save alot here since this is the same expression for each of the coordinates
////    double r = distance;
////    double z = r1(i,2);
////    return z*(alph*alph*r*r - 12*alph*r + 16)*exp(-alph*r/2)/(48*r);
////}

//double VariationalMC::psi_p2y_doubleDerivative(double distance, mat &r1, int i){
//    double r = distance;
//    double y = r1(i,1);
//    return alph*y*(alph*r - 8)*exp(-alph*r/2)/(4*r);
//}

//double VariationalMC::psi_p2x_doubleDerivative(double distance, mat &r1, int i){
//    double r = distance;
//    double x = r1(i,0);
//    return alph*x*(alph*r - 8)*exp(-alph*r/2)/(4*r);
//}


//double VariationalMC::psiDerivative(double distance, mat &r, int i, int j, int k) {
//    if (k == 0) {
//        return psi_s1_derivative(distance,r(i,j));
//    } else if (k == 1) {
//        return psi_s2_derivative(distance,r(i,j));
//    } else if (k == 2) {
//        return psi_p2z_derivative(distance,r,i,j);
//    } else if (k == 3) {
//        return psi_p2y_derivative(distance,r,i,j);
//    } else if (k == 4) {
//        return psi_p2x_derivative(distance,r,i,j);
//    } else {
//        return 0;
//        cout << "Error in psiDerivative: k = " << k  << endl;
//    }
//}

//double VariationalMC::psiDoubleDerivative(double distance, mat &r, int i, int j) {
//    if (j == 0) {
//        return psi_s1_doubleDerivative(distance);
//    } else if (j == 1) {
//        return psi_s2_doubleDerivative(distance);
//    } else if (j == 2) {
//        return psi_p2z_doubleDerivative(distance,r,i);
//    } else if (j == 3) {
//        return psi_p2y_doubleDerivative(distance,r,i);
//    } else if (j == 4) {
//        return psi_p2x_doubleDerivative(distance,r,i);
//    } else {
//        return 0;
//        cout << "Error in psiDoubleDerivative: j = " << j  << endl;
//    }
//}



///* Choses which wave function to calculate from the int j, then calls the corresponding psi_j. */
//double VariationalMC::slaterPsi(double distance, mat &r, int i, int j) {
//    if (j == 0) {
//        return psi_s1(distance1) + psi_s1(distance2);
//    } else if (j == 1) {
//        return psi_s2(distance);
//    } else if (j == 2) {
//        return psi_p2z(distance,r,i);
//    } else if (j == 3) {
//        return psi_p2y(distance,r,i);
//    } else if (j == 4) {
//        return psi_p2x(distance,r,i);
//    } else {
//        return 0;
//        cout << "Error in slaterPsi" << j<< endl;
//    }
//}
//*/
