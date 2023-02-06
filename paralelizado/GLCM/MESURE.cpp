#include "MESURE.h"
#include <iostream>
#include <fstream>
#include <QCoreApplication>
#include <opencv2/opencv.hpp>
#include <opencv2/core/mat.hpp>
#include <opencv2/highgui.hpp>
#include <iostream>
#include "xlsxwriter.h"


using namespace cv;
using namespace std;

//DEFINICIONES DE PARAMETRICAS
#define PGM_MAXMAXVAL 255
#define EPSILON 0.000000001
#define RADIX 2.0
#define SIGN(x,y) ((y)<0 ? -fabs(x) : fabs(x))
#define SWAP(a,b) {y=(a);(a)=(b);(b)=y;}

MESURE::MESURE(){

}

//Initialize functions that have been used for measuring co-occurence matrixes for 0,45,90,135 degree angle
double** CoOcMat_Angle_0   (int distance, u_int8_t **grays, int rows, int cols, int* tone_LUT, int tone_count);
double** CoOcMat_Angle_45  (int distance, u_int8_t **grays, int rows, int cols, int* tone_LUT, int tone_count);
double** CoOcMat_Angle_90  (int distance, u_int8_t **grays, int rows, int cols, int* tone_LUT, int tone_count);
double** CoOcMat_Angle_135 (int distance, u_int8_t **grays, int rows, int cols, int* tone_LUT, int tone_count);

//INICIALIZAR FUNCIONES
double f1_asm (double **P, int Ng);

//RECURSOS DE LOS METODOS
double *allocate_vector (int nl, int nh);
double **allocate_matrix (int nrl, int nrh, int ncl, int nch);
void free_matrix(double **matrix,int nrh);


//6 INFORMATION MESURE CORRELATION 1 Y 2
double MESURE::f12_icorr (double **P, int Ng) const {
    cout << "f12 mesure funciona";
    int i, j;
    double *px, *py;
    double hx = 0, hy = 0, hxy = 0, hxy1 = 0, hxy2 = 0;
    px = allocate_vector (0, Ng);
    py = allocate_vector (0, Ng);
    for (i = 0; i < Ng; ++i) {
        for (j = 0; j < Ng; ++j) {
            px[i] += P[i][j];
            py[j] += P[i][j];
        }
    }
    for (i = 0; i < Ng; ++i)
        for (j = 0; j < Ng; ++j) {
            hxy1 -= P[i][j] * log10 (px[i] * py[j] + EPSILON)/log10(2.0);
            hxy2 -= px[i] * py[j] * log10 (px[i] * py[j] + EPSILON)/log10(2.0);
            hxy -= P[i][j] * log10 (P[i][j] + EPSILON)/log10(2.0);
        }
    for (i = 0; i < Ng; ++i) {
        hx -= px[i] * log10 (px[i] + EPSILON)/log10(2.0);
        hy -= py[i] * log10 (py[i] + EPSILON)/log10(2.0);
    }

    free(px);
    free(py);
        if ((hx > hy ? hx : hy)==0) return(1);
        else
    return ((hxy - hxy1) / (hx > hy ? hx : hy));
}
double MESURE::f13_icorr (double **P, int Ng) const {
    cout << "f13 mesure funciona";
    int i, j;
    double *px, *py;
    double hx = 0, hy = 0, hxy = 0, hxy1 = 0, hxy2 = 0;
    px = allocate_vector (0, Ng);
    py = allocate_vector (0, Ng);
    for (i = 0; i < Ng; ++i) {
        for (j = 0; j < Ng; ++j) {
          px[i] += P[i][j];
          py[j] += P[i][j];
        }
    }
    for (i = 0; i < Ng; ++i)
        for (j = 0; j < Ng; ++j) {
            hxy1 -= P[i][j] * log10 (px[i] * py[j] + EPSILON)/log10(2.0);
            hxy2 -= px[i] * py[j] * log10 (px[i] * py[j] + EPSILON)/log10(2.0);
            hxy -= P[i][j] * log10 (P[i][j] + EPSILON)/log10(2.0);
        }
    for (i = 0; i < Ng; ++i) {
        hx -= px[i] * log10 (px[i] + EPSILON)/log10(2.0);
        hy -= py[i] * log10 (py[i] + EPSILON)/log10(2.0);
    }
    free(px);
    free(py);
    return (sqrt (fabs (1 - exp (-2.0 * (hxy2 - hxy)))));
}

