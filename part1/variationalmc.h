#ifndef VARIATIONALMC_H
#define VARIATIONALMC_H

class VariationalMC {
    private:
        double h;
        double h2;
        double alph;
        double beta;
        double Z;

    public:
        VariationalMC(double, double, double, double, double);

};

#endif // VARIATIONALMC_H
