#include <iostream>
#include <vector>
#include <string>

#include "dicom_read/dicomutils.h"
#include "dicom_read/DicomReader.h"
#include "fcm/fcm.h"

using namespace std;

int main()
{

    DicomReader dicomObj("/home/user/PARALELA/dataset/MasaMicro1.dcm");
    int rows = dicomObj.getHeight();
    int cols = dicomObj.getWidth();

    //int **arr = dicomObj.getImageArray(12);

    //double **data = DicomUtils::asFCMPointsData(arr, rows, cols);
    //FCM *fcm;
    //fcm = new FCM(2.0, 0.00001); // Fuzzines and epsilon

    //int clusters = 3;

    //fcm->init(data, clusters, rows, cols);

    //fcm->eval();
    //double **centers;
    //centers = fcm->getCenters();
    cout << "FUNCIONO!" << rows << " and " << cols;

    return 0;
}

