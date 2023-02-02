#include <QCoreApplication>
#include <opencv4/opencv2/opencv.hpp>
#include <opencv4/opencv2/core/mat.hpp>
#include <opencv2/opencv.hpp>
#include <opencv4/opencv2/highgui.hpp>
#include <iostream>
#include "xlsxwriter.h"
#include <stdio.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <stdlib.h>


using namespace cv;
using namespace std;


//DEFINICIONES DE PARAMETRICAS
#define PGM_MAXMAXVAL 255
#define EPSILON 0.000000001
#define RADIX 2.0
#define SIGN(x,y) ((y)<0 ? -fabs(x) : fabs(x))
#define SWAP(a,b) {y=(a);(a)=(b);(b)=y;}


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
double** CoOcMat_Angle_0 (int distance, u_int8_t **grays,
                         int rows, int cols, int* tone_LUT, int tone_count)
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


//MAX CORREALCION HEIS

int hessenberg (double **a, int n, double wr[], double wi[])
{
  int nn, m, l, k, j, its, i, mmin;
  double z, y, x, w, v, u, t, s, r=0.0, q=0.0, p=0.0, anorm;

  anorm = fabs (a[1][1]);
  for (i = 2; i <= n; i++)
    for (j = (i - 1); j <= n; j++)
      anorm += fabs (a[i][j]);
  nn = n;
  t = 0.0;
  while (nn >= 1)
  {
    its = 0;
    do
    {
      for (l = nn; l >= 2; l--)
      {
    s = fabs (a[l - 1][l - 1]) + fabs (a[l][l]);
    if (s == 0.0)
      s = anorm;
    if ((double) (fabs (a[l][l - 1]) + s) == s)
      break;
      }
      x = a[nn][nn];
      if (l == nn)
      {
    wr[nn] = x + t;
    wi[nn--] = 0.0;
      }
      else
      {
    y = a[nn - 1][nn - 1];
    w = a[nn][nn - 1] * a[nn - 1][nn];
    if (l == (nn - 1))
    {
      p = 0.5 * (y - x);
      q = p * p + w;
      z = sqrt (fabs (q));
      x += t;
      if (q >= 0.0)
      {
        z = p + SIGN (z, p);
        wr[nn - 1] = wr[nn] = x + z;
        if (z)
          wr[nn] = x - w / z;
        wi[nn - 1] = wi[nn] = 0.0;
      }
      else
      {
        wr[nn - 1] = wr[nn] = x + p;
        wi[nn - 1] = -(wi[nn] = z);
      }
      nn -= 2;
    }
    else
    {
      if (its == 30)
        {
         return 0;
         }
      if (its == 10 || its == 20)
      {
        t += x;
        for (i = 1; i <= nn; i++)
          a[i][i] -= x;
        s = fabs (a[nn][nn - 1]) + fabs (a[nn - 1][nn - 2]);
        y = x = 0.75 * s;
        w = -0.4375 * s * s;
      }
      ++its;
      for (m = (nn - 2); m >= l; m--)
      {
        z = a[m][m];
        r = x - z;
        s = y - z;
        p = (r * s - w) / a[m + 1][m] + a[m][m + 1];
        q = a[m + 1][m + 1] - z - r - s;
        r = a[m + 2][m + 1];
        s = fabs (p) + fabs (q) + fabs (r);
        p /= s;
        q /= s;
        r /= s;
        if (m == l)
          break;
        u = fabs (a[m][m - 1]) * (fabs (q) + fabs (r));
        v = fabs (p) * (fabs (a[m - 1][m - 1]) +
                fabs (z) + fabs (a[m + 1][m + 1]));
        if ((double) (u + v) == v)
          break;
      }
      for (i = m + 2; i <= nn; i++)
      {
        a[i][i - 2] = 0.0;
        if (i != (m + 2))
          a[i][i - 3] = 0.0;
      }
      for (k = m; k <= nn - 1; k++)
      {
        if (k != m)
        {
          p = a[k][k - 1];
          q = a[k + 1][k - 1];
          r = 0.0;
          if (k != (nn - 1))
        r = a[k + 2][k - 1];
          if ( (x = fabs (p) + fabs (q) + fabs (r)) )
          {
        p /= x;
        q /= x;
        r /= x;
          }
        }
        if ( (s = SIGN (sqrt (p * p + q * q + r * r), p)) )
        {
          if (k == m)
          {
        if (l != m)
          a[k][k - 1] = -a[k][k - 1];
          }
          else
        a[k][k - 1] = -s * x;
          p += s;
          x = p / s;
          y = q / s;
          z = r / s;
          q /= p;
          r /= p;
          for (j = k; j <= nn; j++)
          {
        p = a[k][j] + q * a[k + 1][j];
        if (k != (nn - 1))
        {
          p += r * a[k + 2][j];
          a[k + 2][j] -= p * z;
        }
        a[k + 1][j] -= p * y;
        a[k][j] -= p * x;
          }
          mmin = nn < k + 3 ? nn : k + 3;
          for (i = l; i <= mmin; i++)
          {
        p = x * a[i][k] + y * a[i][k + 1];
        if (k != (nn - 1))
        {
          p += z * a[i][k + 2];
          a[i][k + 2] -= p * r;
        }
        a[i][k + 1] -= p * q;
        a[i][k] -= p;
          }
        }
      }
    }
      }
    } while (l < nn - 1);
  }
return 1;
}

//MAX CORRELACION OPERACION
void mkbalanced (double **a, int n)
{
  int last, j, i;
  double s, r, g, f, c, sqrdx;

  sqrdx = RADIX * RADIX;
  last = 0;
  while (last == 0)
  {
    last = 1;
    for (i = 1; i <= n; i++)
    {
      r = c = 0.0;
      for (j = 1; j <= n; j++)
    if (j != i)
    {
      c += fabs (a[j][i]);
      r += fabs (a[i][j]);
    }
      if (c && r)
      {
    g = r / RADIX;
    f = 1.0;
    s = c + r;
    while (c < g)
    {
      f *= RADIX;
      c *= sqrdx;
    }
    g = r * RADIX;
    while (c > g)
    {
      f /= RADIX;
      c /= sqrdx;
    }
    if ((c + r) / f < 0.95 * s)
    {
      last = 0;
      g = 1.0 / f;
      for (j = 1; j <= n; j++)
        a[i][j] *= g;
      for (j = 1; j <= n; j++)
        a[j][i] *= f;
    }
      }
    }
  }
}

//MAX CORRELACION REDUCCION

void reduction (double **a, int n)
{
  int m, j, i;
  double y, x;

  for (m = 2; m < n; m++)
  {
    x = 0.0;
    i = m;
    for (j = m; j <= n; j++)
    {
      if (fabs (a[j][m - 1]) > fabs (x))
      {
    x = a[j][m - 1];
    i = j;
      }
    }
    if (i != m)
    {
      for (j = m - 1; j <= n; j++)
    SWAP (a[i][j], a[m][j])
    for (j = 1; j <= n; j++)
      SWAP (a[j][i], a[j][m])
      a[j][i] = a[j][i];
    }
    if (x)
    {
      for (i = m + 1; i <= n; i++)
      {
    if ( (y = a[i][m - 1]) )
    {
      y /= x;
      a[i][m - 1] = y;
      for (j = m; j <= n; j++)
        a[i][j] -= y * a[m][j];
      for (j = 1; j <= n; j++)
        a[j][m] += y * a[j][i];
    }
      }
    }
  }
}



//SECOND ANGULAR MOMENT
double f1_asm (double **P, int Ng) {
    int i, j;
    double sum = 0;

    for (i = 0; i < Ng; ++i)
        for (j = 0; j < Ng; ++j)
            sum += P[i][j] * P[i][j];

    return sum;
}

//CORRELACION
double f3_corr (double **P, int Ng) {
    int i, j;
    double sum_sqrx = 0, sum_sqry = 0, tmp, *px;
    double meanx =0 , meany = 0 , stddevx, stddevy;
    px = allocate_vector (0, Ng);
    for (i = 0; i < Ng; ++i)
        px[i] = 0;
    for (i = 0; i < Ng; ++i)
        for (j = 0; j < Ng; ++j)
            px[i] += P[i][j];
    for (i = 0; i < Ng; ++i) {
        meanx += px[i]*i;
        sum_sqrx += px[i]*i*i;
    }
    meany = meanx;
    sum_sqry = sum_sqrx;
    stddevx = sqrt (sum_sqrx - (meanx * meanx));
    stddevy = stddevx;
    for (tmp = 0, i = 0; i < Ng; ++i)
        for (j = 0; j < Ng; ++j)
              tmp += i*j*P[i][j];

    free(px);
        if (stddevx * stddevy==0) return(1);
        else return (tmp - meanx * meany) / (stddevx * stddevy);
}


//DIFERENCIA ENTROPIA
double f11_dentropy (double **P, int Ng) {
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




// DIFERENCIA VARIANZA
double f10_dvar (double **P, int Ng) {
    int i, j;
    double sum = 0, sum_sqr = 0, var = 0;
    double *Pxpy = allocate_vector (0, 2*Ng);

    for (i = 0; i <= 2 * Ng; ++i)
        Pxpy[i] = 0;

    for (i = 0; i < Ng; ++i)
        for (j = 0; j < Ng; ++j)
            Pxpy[abs (i - j)] += P[i][j];

    for (i = 0; i < Ng; ++i) {
        sum += i * Pxpy[i] ;
        sum_sqr += i * i * Pxpy[i] ;
    }
    var = sum_sqr - sum*sum ;
    free (Pxpy);
    return var;
}


//ENTROPIA
double f9_entropy (double **P, int Ng) {
    int i, j;
    double entropy = 0;
    for (i = 0; i < Ng; ++i)
        for (j = 0; j < Ng; ++j)
            entropy += P[i][j] * log10 (P[i][j] + EPSILON)/log10(2.0) ;
    return -entropy;
}

//INVERSE DIFERENCE MOEMNT
double f5_idm (double **P, int Ng) {
    int i, j;
    double idm = 0;

    for (i = 0; i < Ng; ++i)
        for (j = 0; j < Ng; ++j)
            idm += P[i][j] / (1 + (i - j) * (i - j));

    return idm;
}



//CONTRASTE
double f2_contrast (double **P, int Ng) {
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



//INFORMATION MESURE CORRELATION 1 Y 2
double f12_icorr (double **P, int Ng) {
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


double f13_icorr (double **P, int Ng) {
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



//MAXIMO CORRELATION COEFICIENT
double f14_maxcorr (double **P, int Ng) {
    int i, j, k;
    double *px, *py, **Q;
    double *x, *iy, tmp;
    double f=0.0;
    px = allocate_vector (0, Ng);
    py = allocate_vector (0, Ng);
    Q = allocate_matrix (1, Ng + 1, 1, Ng + 1);
    x = allocate_vector (1, Ng);
    iy = allocate_vector (1, Ng);

    for (i = 0; i < Ng; ++i) {
        for (j = 0; j < Ng; ++j) {
            px[i] += P[i][j];
            py[j] += P[i][j];
        }
    }

    for (i = 0; i < Ng; ++i) {
        for (j = 0; j < Ng; ++j) {
            Q[i + 1][j + 1] = 0;
            for (k = 0; k < Ng; ++k)
                          if (px[i] && py[k])
                Q[i + 1][j + 1] += P[i][k] * P[j][k] / px[i] / py[k];
        }
    }
    mkbalanced (Q, Ng);
    reduction (Q, Ng);
    if (!hessenberg (Q, Ng, x, iy)) {
        for (i=1; i<=Ng+1; i++) free(Q[i]+1);
        free(Q+1);
        free((char *)px);
        free((char *)py);
        free((x+1));
        free((iy+1));
        return 0.0;
    }
    for (i = 2, tmp = x[1]; i <= Ng; ++i)
        tmp = (tmp > x[i]) ? tmp : x[i];

    if (x[Ng - 1]>=0)
      f = sqrt(x[Ng - 1]);

    for (i=1; i<=Ng+1; i++) free(Q[i]+1);
    free(Q+1);
    free((char *)px);
    free((char *)py);
    free((x+1));
    free((iy+1));

    return f;
}

//SUM AVEREGE
double f6_savg (double **P, int Ng) {
    int i, j;
    double savg = 0;
    double *Pxpy = allocate_vector (0, 2*Ng);

    for (i = 0; i <= 2 * Ng; ++i)
        Pxpy[i] = 0;

    for (i = 0; i < Ng; ++i)
        for (j = 0; j < Ng; ++j)
          Pxpy[i + j] += P[i][j];

    for (i = 0; i <= (2 * Ng - 2); ++i)
        savg += i * Pxpy[i];

    free (Pxpy);
    return savg;
}

//SUMA ENTROPIA
double f8_sentropy (double **P, int Ng) {
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


//SUM OF SQUARE VARIANCE
double f4_var (double **P, int Ng) {
    int i, j;
    double mean = 0, var = 0;
    for (i = 0; i < Ng; ++i)
        for (j = 0; j < Ng; ++j)
            mean += i * P[i][j];

    for (i = 0; i < Ng; ++i)
        for (j = 0; j < Ng; ++j)
          var += (i - mean) * (i - mean) * P[i][j];

    return var;
}

//SUM VARIANCE
double f7_svar (double **P, int Ng, double S) {
    int i, j;
    double var = 0;
    double *Pxpy = allocate_vector (0, 2*Ng);

    for (i = 0; i <= 2 * Ng; ++i)
        Pxpy[i] = 0;

    for (i = 0; i < Ng; ++i)
        for (j = 0; j < Ng; ++j)
          Pxpy[i + j] += P[i][j];

    for (i = 0; i <= (2 * Ng - 2); ++i)
        var += (i - S) * (i - S) * Pxpy[i];

    free (Pxpy);
    return var;
}



std::string to_stringC(double x)
{
    std::ostringstream ss;
    ss << x;
    return ss.str();
}



int main()
{
    void *mem1, *mem2;
    struct rusage uso;
    struct rlimit limite;
    getrusage(RUSAGE_SELF, &uso);
    printf("Uso de RAM: %ld KB\n", (long)uso.ru_maxrss);

    std::string image_path = samples::findFile("/home/user/PARALELA/pro1/sky.png");
    Mat img = imread(image_path, IMREAD_GRAYSCALE);
    if(img.empty())
    {
        std::cout << "Could not read the image: " << image_path << std::endl;
        return 1;
    }

    // IMPRIMIR VALORES IMAGEN
    unsigned char *p;
    uchar intensity;
    p = img.data;
    int row, col, rows, cols;
    cols = img.cols; //weigth - col
    rows = img.rows; //heigth - fil

    for(int y = 0; y < img.rows; y++){
            for(int x = 0; x < img.cols; x++){
                intensity = img.data[img.step * y + x * 1];
            }
            cout << endl;
        }

    // CREAR ARRAY 2D CON LAS DIMENSIONES DEL DEGRADADO A GRISES
    unsigned char **pGray;
        pGray = new unsigned char *[rows];

        for(int i = 0; i < rows; i++){
            pGray[i] = new unsigned char[cols];
        }

        for (int y = 0; y < rows; y++)
        {
            for (int x = 0; x < cols; x++)
            {
                pGray[y][x] = img.data[img.step * y + x * 1];
            }
        }
    // IMPRIMIR LA ESCALA DE GRISES
        cout << endl << "GRISES: " << endl << endl;
        for(int y = 0; y < img.cols; y++){
            for(int x = 0; x < img.rows; x++){
                intensity = pGray[y][x];
                cout << (unsigned int)intensity << "   ";
            }
            cout << endl;
        }


        int toneLUT[PGM_MAXMAXVAL + 1];		// toneLUT is an array that can hold 256 values
        int toneCount = 0;
        int iTone;

        //RELLENAR CON -1
        for(row = PGM_MAXMAXVAL; row >= 0; --row)
                toneLUT[row] = -1;

        for(row = rows - 1; row >= 0; --row){
                for(col = 0; col < cols; ++col){
                    toneLUT[(u_int8_t)img.data[img.step * row + col * 1]] = (u_int8_t)img.data[img.step * row + col * 1];
                }
         }
        for (row = PGM_MAXMAXVAL, toneCount = 0; row >= 0; --row){
                if (toneLUT[row] != -1)
                    toneCount++;
                else
                    ;
            }
        for (row = 0, iTone = 0; row <= PGM_MAXMAXVAL; row++){
                if (toneLUT[row] != -1)
                  toneLUT[row] = iTone++;
            }

        double **pMatriz;
        int distancia = 1;
        pMatriz = CoOcMat_Angle_0(distancia, pGray, rows, col, toneLUT, toneCount);
        double m_asm, m_contrast, m_corr, m_var, m_idm, m_savg, m_svar, m_sentropy, m_entropy, m_dvar, m_dentropy, m_icorr1, m_icorr2, m_maxcorr;
        //SECOND ANGULAR MOMENT
        m_asm = f1_asm (pMatriz , toneCount);
        //CORRELACION
        m_corr = f3_corr(pMatriz, toneCount);
        //DIFERENCIA ENTROPIA
        m_dentropy = f11_dentropy(pMatriz, toneCount);
        // DIFERENCIA VARIANZA
        m_dvar = f10_dvar(pMatriz, toneCount);
        //ENTROPIA
        m_entropy = f9_entropy(pMatriz, toneCount);
        //INVERSE DIFERENCE MOEMNT
        m_idm = f5_idm(pMatriz, toneCount);
        //CONTRASTE
        m_contrast = f2_contrast(pMatriz , toneCount);
        //INFORMATION MESURE CORRELATION 1 Y 2
        m_icorr1 = f12_icorr(pMatriz, toneCount);
        m_icorr2 = f13_icorr(pMatriz, toneCount);
        //MAXIMO CORRELATION COEFICIENT
        m_maxcorr = f14_maxcorr(pMatriz, toneCount);
        //SUM AVEREGE
        m_savg = f6_savg(pMatriz, toneCount);
        //SUMA ENTROPIA
        m_sentropy = f8_sentropy(pMatriz, toneCount);
        //SUM OF SQUARE VARIANCE
        m_var = f4_var(pMatriz, toneCount);
        //SUM VARIANCE
        m_svar = f7_svar(pMatriz, toneCount, m_sentropy);


            imshow("BIenvenido", img);
            waitKey(0);


            lxw_workbook  *workbook  = workbook_new("/home/user/pro1/resultados.xlsx");
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
        return 0;
}




