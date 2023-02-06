#include "CONTRASTE.h"
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

CONTRASTE::CONTRASTE(){

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



//7 CONTRASTE
double CONTRASTE::f2_contrast (double **P, int Ng) const {
    cout << "contraste funciona";
    int i, j, n;
    double sum = 0, bigsum = 0;
    for (n = 0; n < Ng; ++n) {
        for (i = 0; i < Ng; ++i)
            for (j = 0; j < Ng; ++j) {
                if ((i - j) == n || (j - i) == n)
                    sum += P[i][j];
                }
        bigsum += n * n * sum;
        sum = 0;
    }
    return bigsum;
}

