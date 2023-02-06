#ifndef AVEREGEP_H
#define AVEREGEP_H

#include <iostream>
#include <vector>

using namespace std;

class AVEREGEP;

class AVEREGEP
{
public:
    AVEREGEP();
    double f6_savg (double **P, int Ng, int hilos) const;
};

#endif
