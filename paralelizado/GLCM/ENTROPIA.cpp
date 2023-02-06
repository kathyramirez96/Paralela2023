#include "ENTROPIA.h"
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

ENTROPIA::ENTROPIA(){

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

//SUMA ENTROPIA
double ENTROPIA::f8_sentropy (double **P, int Ng) const {
    cout << "Entroipia funcion";
    int i, j;
    double sentropy = 0;
    double *Pxpy = allocate_vector (0, 2*Ng);
    for (i = 0; i <= 2 * Ng; ++i)
        Pxpy[i] = 0;
    for (i = 0; i < Ng; ++i)
        for (j = 0; j < Ng; ++j)
          Pxpy[i + j + 2] += P[i][j];
    for (i = 2; i <= 2 * Ng; ++i)
        sentropy -= Pxpy[i] * log10 (Pxpy[i] + EPSILON)/log10(2.0) ;
    free (Pxpy);
    return sentropy;
}

//9 ENTROPIA
double ENTROPIA::f9_entropy (double **P, int Ng) const {
    int i, j;
    double entropy = 0;
    for (i = 0; i < Ng; ++i)
        for (j = 0; j < Ng; ++j)
            entropy += P[i][j] * log10 (P[i][j] + EPSILON)/log10(2.0) ;
    return -entropy;
}


//11 DIFERENCIA ENTROPIA
double ENTROPIA::f11_dentropy (double **P, int Ng) const {
    int i, j;
    double sum = 0;
    double *Pxpy = allocate_vector (0, 2*Ng);

    for (i = 0; i <= 2 * Ng; ++i)
        Pxpy[i] = 0;

    for (i = 0; i < Ng; ++i)
        for (j = 0; j < Ng; ++j)
            Pxpy[abs (i - j)] += P[i][j];

    for (i = 0; i < Ng; ++i)
        /*    sum += Pxpy[i] * log10 (Pxpy[i] + EPSILON); */
        sum += Pxpy[i] * log10 (Pxpy[i] + EPSILON)/log10(2.0) ;

    free (Pxpy);
    return -sum;
}


