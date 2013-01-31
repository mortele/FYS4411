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

        double computePsi(const mat &);
        double computeEnergy(mat &, mat &, double);
        double computeDoubleDerivative(double, double, double);
        void   updateRmatrix(const mat &, mat &);
        void   updateForDerivative(mat &, const mat &, int );

    public:
        VariationalMC();
        void runMetropolis();
};

#endif // VARIATIONALMC_H
