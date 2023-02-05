#include <iostream>
#include <vector>
#include <string>

#include "dicom/dicomutils.h"
#include "dicom/DicomReader.h"
#include "fcm/fcm.h"

using namespace std;

int main()
{

    DicomReader dicomObj("/home/will/Projects/dicom-classifier/data/20586908_6c613a14b80a8591_MG_R_CC_ANON.dcm");
    int rows = dicomObj.getHeight();
    int cols = dicomObj.getWidth();

    int **arr = dicomObj.getImageArray(12);

    double **data = DicomUtils::asFCMPointsData(arr, rows, cols);
    FCM *fcm;
    fcm = new FCM(2.0, 0.00001); // Fuzzines and epsilon

    int clusters = 3;

    fcm->init(data, clusters, rows, cols);

    fcm->eval();
    double **centers;
    centers = fcm->getCenters();

    return 0;
}

