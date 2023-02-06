#ifndef IMAGENES_H
#define IMAGENES_H

#include <iostream>
#include <vector>

using namespace std;

class IMAGENES;

class IMAGENES
{
public:
    IMAGENES();
    void mensaje(string data) const;
    double ** ESCALAGRISES(vector<vector<int>>  imagen, int min, int max) const;
    int ObtenertoneCount(vector<vector<int>>  imagen, int min, int max) const;
    void generarExcel(double m_asm,
                      double m_contrast,
                      double m_corr,
                      double m_var,
                      double m_idm,
                      double m_savg,
                      double m_svar,
                      double m_sentropy,
                      double m_entropy,
                      double m_dvar,
                      double m_dentropy,
                      double m_icorr1,
                      double m_icorr2,
                      double m_maxcorr) const;
};

#endif // IMAGENES_H
