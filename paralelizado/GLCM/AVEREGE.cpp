#include "AVEREGE.h"
#include <iostream>
#include <fstream>
#include <QCoreApplication>
#include <opencv2/opencv.hpp>
#include <opencv2/core/mat.hpp>
#include <opencv2/highgui.hpp>
#include <iostream>
#include "xlsxwriter.h"
#include <omp.h>

using namespace cv;
using namespace std;

//DEFINICIONES DE PARAMETRICAS
#define PGM_MAXMAXVAL 255
#define EPSILON 0.000000001
#define RADIX 2.0
#define SIGN(x,y) ((y)<0 ? -fabs(x) : fabs(x))
#define SWAP(a,b) {y=(a);(a)=(b);(b)=y;}

AVEREGE::AVEREGE(){

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



//4 SUM AVEREGE
double AVEREGE::f6_savg (double **P, int Ng, int hilos) const {
        cout << "AVEREGE funcionando";
        int i, j;
        double savg = 0;
        double *Pxpy = allocate_vector (0, 2*Ng);
#pragma omp parallel num_threads(hilos)
    {

        for (i = 0; i <= 2 * Ng; ++i)
            Pxpy[i] = 0;

        for (i = 0; i < Ng; ++i)
            for (j = 0; j < Ng; ++j)
              Pxpy[i + j] += P[i][j];

        for (i = 0; i <= (2 * Ng - 2); ++i)
            savg += i * Pxpy[i];
    }
        free (Pxpy);
        return savg;

}
