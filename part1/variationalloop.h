#ifndef VARIATIONALLOOP_H
#define VARIATIONALLOOP_H
#include "variationalmc.h"

class VariationalLoop {
    private:
        VariationalMC m;
        int     N;
        double  start;
        double  end;
        double  E;
        double  minE;
        double  minA;
        double  minB;
        double  a;
        double  b;


    public:
        VariationalLoop();
        void initialize_helium();
        void run();
};

#endif // VARIATIONALLOOP_H
