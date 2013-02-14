#ifndef VARIATIONALMC_H
#define VARIATIONALMC_H

#include <armadillo>

using namespace std;
using namespace arma;


class VariationalMC {
    private:
        int     nParticles;
        int     nDimensions;
        int     nCycles;
        int     N;
        long    idum;
        double  charge;
        double  h;
        double  h2;
        double  alph;
        double  beta;
        double  Z;
        double  stepSize;
        double  D;

        double computePsi(const mat&);
        double computeEnergy(mat&, mat&, double);
        double computeDoubleDerivative(double, double, double);
        double computeFirstDerivative (double, double);
        void   updateRmatrix(const mat&, mat&);
        void   updateForDerivative(mat&, const mat&, int );
        vec    computeQuantumForce(mat&, mat&, double);

    public:
        VariationalMC();
        double runMetropolis(double, double);
};

#endif // VARIATIONALMC_H
