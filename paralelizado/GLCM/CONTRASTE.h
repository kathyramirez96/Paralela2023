#ifndef CONTRASTE_H
#define CONTRASTE_H

#include <iostream>
#include <vector>

using namespace std;

class CONTRASTE;

class CONTRASTE
{
public:
    CONTRASTE();
    double f2_contrast (double **P, int Ng) const;
};

#endif
