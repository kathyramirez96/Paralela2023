#ifndef SQUAREVARIANCE_H
#define SQUAREVARIANCE_H

#include <iostream>
#include <vector>

using namespace std;

class SQUAREVARIANCE;

class SQUAREVARIANCE
{
public:
    SQUAREVARIANCE();
    double f4_var (double **P, int Ng) const;
};

#endif
