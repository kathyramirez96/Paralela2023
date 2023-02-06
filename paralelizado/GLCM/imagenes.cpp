#include "imagenes.h"
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

IMAGENES::IMAGENES(){

}

//Initialize functions that have been used for measuring co-occurence matrixes for 0,45,90,135 degree angle
double** CoOcMat_Angle_0   (int distance, u_int8_t **grays, int rows, int cols, int* tone_LUT, int tone_count);
double** CoOcMat_Angle_45  (int distance, u_int8_t **grays, int rows, int cols, int* tone_LUT, int tone_count);
double** CoOcMat_Angle_90  (int distance, u_int8_t **grays, int rows, int cols, int* tone_LUT, int tone_count);
double** CoOcMat_Angle_135 (int distance, u_int8_t **grays, int rows, int cols, int* tone_LUT, int tone_count);

//-----------------------------CALCULOS-----------------------------
// LOCACIONES DE MATRIZ CON RANGO
double **allocate_matrix (int nrl, int nrh, int ncl, int nch)
{
    int i;
    double **m;
    m = (double **) malloc ((unsigned) (nrh - nrl + 1) * sizeof (double *));
    if (!m) fprintf (stderr, "memoria local no encontrada "), exit (1);
    m -= ncl;
    for (i = nrl; i <= nrh; i++) {
        m[i] = (double *) malloc ((unsigned) (nch - ncl + 1) * sizeof (double));
        if (!m[i]) fprintf (stderr, "memory allocation failure (allocate_matrix 2) "), exit (2);
        m[i] -= ncl;
    }
    return m;
}


//VECTOR DE LOACION
double *allocate_vector (int nl, int nh) {
    double *v;

    v = (double *) calloc (1, (unsigned) (nh - nl + 1) * sizeof (double));
    if (!v) fprintf (stderr, "memory allocation failure (allocate_vector) "), exit (1);

    return v - nl;
}


//MATRIZ DE CONCURRENCIA CON ANGULO 0
double** CoOcMat_Angle_0 (int distance, u_int8_t **grays, int rows, int cols, int* tone_LUT, int tone_count)
{
    int d = distance;
    int x, y;
    int row, col, itone, jtone;
    double count=0.0;
    double** matrix = allocate_matrix (0, tone_count, 0, tone_count);

    for (itone = 0; itone < tone_count; ++itone)
        for (jtone = 0; jtone < tone_count; ++jtone)
            matrix[itone][jtone] = 0.0;

    for (row = 0; row < rows; ++row)
        for (col = 0; col < cols; ++col) {
            if (col + d < cols) {
                x = tone_LUT[grays[row][col]];
                y = tone_LUT[grays[row][col + d]];
                matrix[x][y]++;
                matrix[y][x]++;
                count += 2.0 ;
            }
        }

    for (itone = 0; itone < tone_count; ++itone){
          for (jtone = 0; jtone < tone_count; ++jtone){
            if (count==0.0)   /* protect from error */
               matrix[itone][jtone]=0.0;
               else matrix[itone][jtone] /= count;
          }
    }

    return matrix;
}

//INICIALIZAR FUNCIONES
double f1_asm (double **P, int Ng);

//RECURSOS DE LOS METODOS
double *allocate_vector (int nl, int nh);
double **allocate_matrix (int nrl, int nrh, int ncl, int nch);
void free_matrix(double **matrix,int nrh);


// CREAR ARRAY 2D CON LAS DIMENSIONES DEL DEGRADADO A GRISES
unsigned char** generarEscala(vector<vector<int>>  imagen, int rows, int cols){
    unsigned char **pGray;
       pGray = new unsigned char *[rows];
       for(int i = 0; i < rows; i++){
           pGray[i] = new unsigned char[cols];
       }
       for (int y = 0; y < rows; y++)
       {
           for (int x = 0; x < cols; x++)
           {
              if(imagen.data()[x][y] > 125){
                     pGray[y][x] = (u_int8_t)255;
              }
              else{
                    pGray[y][x] = (u_int8_t)0;
              }
           }
       }
    return pGray;
}



//0 ESCALA DE GRISES
double ** IMAGENES::ESCALAGRISES(vector<vector<int>>  imagen, int min, int max) const{
    int row, col, rows, cols;
    cols = max; //weigth - col - max
    rows = min; //heigth - fil - min
    // CREAR ARRAY 2D CON LAS DIMENSIONES DEL DEGRADADO A GRISES
    unsigned char **pGray = generarEscala(imagen, rows,cols);
    int toneLUT[PGM_MAXMAXVAL + 1];		// toneLUT is an array that can hold 256 values
    int toneCount = 0;
    int iTone;

    //RELLENAR CON -1
    for(row = PGM_MAXMAXVAL; row >= 0; --row)
          toneLUT[row] = -1;
    for(row = rows - 1; row >= 0; --row){
          for(col = 0; col < cols; ++col){
              std::vector<int, std::allocator<int>> const *aux = imagen.data();
              toneLUT[(u_int8_t)aux[row][imagen.size() * row + col * 1]] = (u_int8_t)aux[row][imagen.size() * row + col * 1];
          }
    }
    for (row = PGM_MAXMAXVAL, toneCount = 0; row >= 0; --row){
          if (toneLUT[row] != -1){
              toneCount++;
          }
      }
    for (row = 0, iTone = 0; row <= PGM_MAXMAXVAL; row++){
          if (toneLUT[row] != -1)
            toneLUT[row] = iTone++;
      }
    double **pMatriz;
    int distancia = 1;
    pMatriz = CoOcMat_Angle_0(distancia, pGray, rows, col, toneLUT, toneCount);
    return pMatriz;
}

//0 CONTADOR MATRIZ
int IMAGENES::ObtenertoneCount(vector<vector<int>>  imagen, int min, int max) const{
    int row, col, rows, cols;
    cols = max; //weigth - col - max
    rows = min; //heigth - fil - min
    int toneLUT[PGM_MAXMAXVAL + 1];		// toneLUT is an array that can hold 256 values
    int toneCount = 0;
    int iTone;
    //RELLENAR CON -1
    for(row = PGM_MAXMAXVAL; row >= 0; --row)
          toneLUT[row] = -1;
    for(row = rows - 1; row >= 0; --row){
          for(col = 0; col < cols; ++col){
              std::vector<int, std::allocator<int>> const *aux = imagen.data();
              toneLUT[(u_int8_t)aux[row][imagen.size() * row + col * 1]] = (u_int8_t)aux[row][imagen.size() * row + col * 1];
          }
    }
    for (row = PGM_MAXMAXVAL, toneCount = 0; row >= 0; --row){
          if (toneLUT[row] != -1){
              toneCount++;
          }
      }
    for (row = 0, iTone = 0; row <= PGM_MAXMAXVAL; row++){
          if (toneLUT[row] != -1)
            toneLUT[row] = iTone++;
      }
    return toneCount;
}


void IMAGENES::generarExcel(double m_asm,
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
                            double m_maxcorr) const
{
    lxw_workbook  *workbook  = workbook_new("/home/user/PARALELA/RESULTADOS/resultados.xlsx");
    lxw_worksheet *worksheet = workbook_add_worksheet(workbook, NULL);
    worksheet_write_string(worksheet, 0, 0, "SECOND ANGULAR MOMENT", NULL);
    worksheet_write_number(worksheet, 0, 1, m_asm, NULL);
    worksheet_write_string(worksheet, 1, 0, "CORRELACION", NULL);
    worksheet_write_number(worksheet, 1, 1, m_corr, NULL);
    worksheet_write_string(worksheet, 2, 0, "DIFERENCIA ENTROPIA", NULL);
    worksheet_write_number(worksheet, 2, 1, m_dentropy, NULL);
    worksheet_write_string(worksheet, 3, 0, "DIFERENCIA VARIANZA", NULL);
    worksheet_write_number(worksheet, 3, 1, m_dvar, NULL);
    worksheet_write_string(worksheet, 4, 0, "ENTROPIA", NULL);
    worksheet_write_number(worksheet, 4, 1, m_entropy, NULL);
    worksheet_write_string(worksheet, 5, 0, "INVERSE DIFERENCE MOEMNT", NULL);
    worksheet_write_number(worksheet, 5, 1, m_idm, NULL);
    worksheet_write_string(worksheet, 6, 0, "CONTRASTE", NULL);
    worksheet_write_number(worksheet, 6, 1, m_contrast, NULL);
    worksheet_write_string(worksheet, 7, 0, "INFORMATION MESURE CORRELATION 1", NULL);
    worksheet_write_number(worksheet, 7, 1, m_icorr1, NULL);
    worksheet_write_string(worksheet, 8, 0, "INFORMATION MESURE CORRELATION 2", NULL);
    worksheet_write_number(worksheet, 8, 1, m_icorr2, NULL);
    worksheet_write_string(worksheet, 9, 0, "MAXIMO CORRELATION COEFICIENT", NULL);
    worksheet_write_number(worksheet, 9, 1, m_maxcorr, NULL);
    worksheet_write_string(worksheet, 10, 0, "SUM AVEREGE", NULL);
    worksheet_write_number(worksheet, 10, 1, m_savg, NULL);
    worksheet_write_string(worksheet, 11, 0, "SUMA ENTROPIA", NULL);
    worksheet_write_number(worksheet, 11, 1, m_sentropy, NULL);
    worksheet_write_string(worksheet, 12, 0, "SUM OF SQUARE VARIANCE", NULL);
    worksheet_write_number(worksheet, 12, 1, m_var, NULL);
    worksheet_write_string(worksheet, 13, 0, "SUM VARIANCE", NULL);
    worksheet_write_number(worksheet, 13, 1, m_svar, NULL);
    workbook_close(workbook);
    cout << "DATOS ALMACENADIOS";
}

void IMAGENES::mensaje(string data) const
{
    cout << "\n" << data <<"\n";
}



