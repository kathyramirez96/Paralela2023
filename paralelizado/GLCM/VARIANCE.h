#ifndef VARIANCE_H
#define VARIANCE_H

#include <iostream>
#include <vector>

using namespace std;

class VARIANCE;

class VARIANCE
{
public:
    VARIANCE();
    double f7_svar(double **P, int Ng, double S) const;
    double f10_dvar (double **P, int Ng) const;
};

#endif
