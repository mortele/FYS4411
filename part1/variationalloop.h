#ifndef VARIATIONALLOOP_H
#define VARIATIONALLOOP_H
#include "variationalmc.h"


class VariationalLoop {
    private:
        VariationalMC m;    // Importance sampled with closed form energy expression.
        //VarMC   m;          // Brute force with closed form energy expression.


        int     N;
        double  start;
        double  end;
        double  E;
        double  minE;
        double  minA;
        double  minB;
        double  a;
        double  b;
        int     my_rank;


    public:
        VariationalLoop();
        void initialize_helium(int);
        void run();
};

#endif // VARIATIONALLOOP_H
