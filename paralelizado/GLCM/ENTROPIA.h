#ifndef ENTROPIA_H
#define ENTROPIA_H

#include <iostream>
#include <vector>

using namespace std;

class ENTROPIA;

class ENTROPIA
{
public:
    ENTROPIA();
    double f8_sentropy (double **P, int Ng) const;
    double f9_entropy (double **P, int Ng) const;
    double f11_dentropy (double **P, int Ng) const;
};

#endif
