#include "CORRELATIONCOEFICIENT.h"
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

CORRELATIONCOEFICIENT::CORRELATIONCOEFICIENT(){

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

//12 CORRELACION
double CORRELATIONCOEFICIENT::f3_corr (double **P, int Ng) const {
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
    stddevx = sqrt (sum_sqry - (meanx * meanx));
    stddevy = stddevx;
    for (tmp = 0, i = 0; i < Ng; ++i)
        for (j = 0; j < Ng; ++j)
              tmp += i*j*P[i][j];

    free(px);
        if (stddevx * stddevy==0) return(1);
        else return (tmp - meanx * meany) / (stddevx * stddevy);
}


//5 MAXIMO CORRELATIONCOEFICIENT
double CORRELATIONCOEFICIENT::f14_maxcorr (double **P, int Ng) const {
    cout << "\ncorrelacion coeficiente funciona \n";
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

