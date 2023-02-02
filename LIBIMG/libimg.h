#ifndef LIBIMG_H
#define LIBIMG_H

#include <QCoreApplication>
#include <opencv4/opencv2/opencv.hpp>
#include <opencv4/opencv2/core/mat.hpp>
#include <opencv2/opencv.hpp>
#include <opencv4/opencv2/highgui.hpp>
#include <iostream>
#include "xlsxwriter.h"
#include <omp.h>

class LIBIMG
{
public:
    LIBIMG();
    void ExecHilo1() const;
    void ExecHilo4() const;
    void ExecHilo8() const;
    void ExecHilo16() const;
    void ExecHilo32() const;
};

#endif // LIBIMG_H
