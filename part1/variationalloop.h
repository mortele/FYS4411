#ifndef VARIATIONALLOOP_H
#define VARIATIONALLOOP_H
#include "variationalmc.h"


class VariationalLoop {
    private:
        VariationalMC m;    // Importance sampled with closed form energy expression.

        double  E;

        int     my_rank;


    public:
        VariationalLoop();
        void initialize_processes(int);
        void run();
};

#endif // VARIATIONALLOOP_H
